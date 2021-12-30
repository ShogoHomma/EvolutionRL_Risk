#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 23:31:12 2018

@author: shogohomma
"""

# import libraries
import csv, datetime, sys, math, os, copy, zipfile, random, itertools #, seaborn
import numpy as np
import pandas as pd

#%%

# os.getcwd()

"""
Set Parameters
"""

#-------- Run Condition --------#

RunCond                 = 'Parallel' # ['Parallel', 'Local']

if RunCond == 'Parallel':
    
    command_arg             = sys.argv
    script_name             = command_arg[0]
    
else:
    
    pass


## Replication数は全てterminalで指定する
replication_i           = int(command_arg[1])
Population              = 30000
Trial                   = 500
Offline_trial           = 100
round_thres             = 4

options                 = 2
choice                  = 'Normal' # either ['Normal', 'Binom']

# Normal
payoff_norm_mean1       = float(command_arg[2])
payoff_norm_sigma1      = float(command_arg[3])
payoff_norm_mean2       = float(command_arg[4])
payoff_norm_sigma2      = float(command_arg[5])

# Binomial
payoff_bi1_x            = 80 # Risky Option
payoff_bi1_y            = 0
payoff_bi1_prob         = 0.5 # prob of x

payoff_bi2_x            = 40 # Sure Option
payoff_bi2_y            = 0
payoff_bi2_prob         = 1 # prob of x


#-- Q-learning --#
ap                      = 'Uniform' # either ['Random', 'Uniform', real number[0,1]]
an                      = 'Uniform' # either ['Random', 'Uniform', real number[0,1]]
beta                    = 'Uniform' # either ['Random', 'Uniform', real number[0, limit_beta]]
limit_beta              = 0.5



#-- parameter set --#

parameters_dict = {
                    "Run_condition": RunCond, "Population": Population, "Trial": Trial, "Offline_trial": Offline_trial, 
                    "round_threshold": round_thres,
                    "number_options": options, "type_op": choice, 
                    "payoff_norm_mean1": payoff_norm_mean1, "payoff_norm_sigma1": payoff_norm_sigma1, "payoff_norm_mean2": payoff_norm_mean2, "payoff_norm_sigma2": payoff_norm_sigma2,
                    "payoff_bi1_x": payoff_bi1_x, "payoff_bi1_y": payoff_bi1_y, "payoff_bi1_prob": payoff_bi1_prob,
                    "payoff_bi2_x": payoff_bi2_x, "payoff_bi2_y": payoff_bi2_y, "payoff_bi2_prob": payoff_bi2_prob,
                    "alpha_positive": ap, "alpha_negative": an, "inverse_temp": beta, "limit_inverse_temp": limit_beta
                    }

"""
    Define Classes & Methods 
"""


class SimulationInfo:
    
    """
    This class contains information about setting of simulation
    """
    
    def __init__(self, parameter):
        
        """
        引数 parameterには、pamameters_dictを与える
        """
        
        # Simulation Setting
        self.Population = parameter["Population"]
        self.Trial = parameter["Trial"]
        self.Offline_trial = parameter["Offline_trial"]
                
        # Choices & Generating rewards
        self.number_options = parameter["number_options"]
        self.type_option = parameter["type_op"]
        
        self.round_threshold = parameter["round_threshold"]

        
        ## Normal Distribution
        self.payoff_norm_mean1 = parameter["payoff_norm_mean1"] # Risky
        self.payoff_norm_sigma1 = parameter["payoff_norm_sigma1"]
        self.payoff_norm_mean2 = parameter["payoff_norm_mean2"] # Sure
        self.payoff_norm_sigma2 = parameter["payoff_norm_sigma2"]
        
        ## Binomial Distribution
        self.payoff_bi1_x = parameter["payoff_bi1_x"] # Risky
        self.payoff_bi1_y = parameter["payoff_bi1_y"]
        self.payoff_bi1_prob = parameter["payoff_bi1_prob"]
        self.payoff_bi2_x = parameter["payoff_bi2_x"] # Sure
        self.payoff_bi2_y = parameter["payoff_bi2_y"]
        self.payoff_bi2_prob = parameter["payoff_bi2_prob"]
        
        # Learning Parameters
        self.alpha_positive = parameter["alpha_positive"]
        self.alpha_negative = parameter["alpha_negative"]
        self.inverse_temp = parameter["inverse_temp"]
        self.limit_inverse_temp = parameter["limit_inverse_temp"]
        
        
class SaveData:
    
    def __init__(self, sim_info):
        
        population = sim_info.Population
        trial = sim_info.Trial
        offline_trial = sim_info.Offline_trial
        
        self.path = "./Analysis/Results/"
        self.dir_name = None  # define later in this attribute
        self.file_name_trial = None # define later in this 
        self.file_name_offline = None # define late in this attribute
        self.file_name_nat= None # define later in this attribute
        self.header_trial = ['agent', 'trial', 'ap', 'an', 'bt', 'Q1', 'Q0', 'p', 'choice', 'get_payoff', 'delta']
        self.header_offline = ['agent', 'trial','ap', 'an', 'bt', 'Q1', 'Q0', 'p', 'offline_choice', 'offline_get_payoff']
        self.header_nat = ['agent', 'average_payoff', 'average_pvalue', 'choice_rate1', 'rate0to0', 'rate0to1', 'rate1to0', 'rate1to1', 'ap', 'an', 'bt']
        self.output_trial = np.empty((population*trial, len(self.header_trial)), dtype = 'float32') # empty matrix for output to main csv
        self.output_offline = np.empty((population*offline_trial, len(self.header_offline)), dtype = 'float32') # empty matrix for output to offline csv
        self.output_nat = np.empty((population, len(self.header_nat))) # empty matrix for output to natural selection csv
        self.time_stamp = datetime.datetime.now()
        
        ### create directory name ###
        
        pop = 'N' + f'{sim_info.Population:0=4d}'
        trial = 'T' + f'{sim_info.Trial:0=4d}'
        ap = 'ap' + f'{sim_info.alpha_positive}' 
        an = 'an' + f'{sim_info.alpha_negative}'
        beta = 'bt' + f'{sim_info.inverse_temp}'
        
        if sim_info.type_option == 'Normal':
            choice = 'C' + f'{sim_info.type_option}' + '_1mean' + f'{sim_info.payoff_norm_mean1}' + '_1sigma' + f'{sim_info.payoff_norm_sigma1}' + '_2mean' + f'{sim_info.payoff_norm_mean2}' + '_2sigma' + f'{sim_info.payoff_norm_sigma2}'
        elif sim_info.type_option == 'Binom':
            choice = 'C' + f'{sim_info.type_option}' + '_1x' + f'{sim_info.payoff_bi1_x}' + '_1y' + f'{sim_info.payoff_bi1_y}' + '_2x'+ f'{sim_info.payoff_bi2_x}'+ '_2y'+ f'{sim_info.payoff_bi2_y}'
        else:
            print("Error in name of choice in directory name. Check class saveData")
            sys.exit()
        
        time_stamp_chr = f'{self.time_stamp:%Y%m%d%H%M%S}'
        
        dir_name_list = [pop, trial, ap, an, beta, choice, time_stamp_chr]
        self.dir_name = self.path +  "_".join(dir_name_list)
        self.path_text = self.dir_name + '/_Parameters_memo.txt'
        
        # if there is no focul dir, create the dir
        
        if os.path.exists(self.dir_name):
            pass
        else:
            print(f"make new dir: {self.dir_name}")
            os.mkdir(self.dir_name)
            
        
    def writeParameterInfo(self, sim_info, parameters_dict):
        
        para_dict = parameters_dict
        
        with open(self.path_text, mode='w') as f:
    
            for key, value in para_dict.items():
    
                f.write(f'{key} : {value}\n')
            
            # finally write the date info
            f.write(f'\n\n{self.time_stamp:%Y%m%d%H%M%S}\n\n')
        
    def createFile(self, sim_info, replication_i):
        
        ### create file name ###
        
        #replication_i = 0
        rep = 'r' + f'{replication_i:0=3d}'
        
        self.file_name_trial = f"./{self.dir_name}/{rep}.csv"
        self.file_name_offline = f"./{self.dir_name}/{rep}_offline.csv"
        self.file_name_nat = f"./{self.dir_name}/{rep}_nat.csv"
        
        # stop the program if there is the focal file.
        """
        if os.path.exists(self.file_name_trial):
            print("same name file exists! stop process.")
            sys.exit()
        else:
            #print(f"make new file: {self.file_name_trial}")
            pd.DataFrame(columns = self.header_trial).to_csv(self.file_name_trial, index = False)
        """
        
        if os.path.exists(self.file_name_nat):
            print("same name file exists! stop process.")
            sys.exit()
        else:
            #print(f"make new file: {self.file_name_nat}")
            pd.DataFrame(columns = self.header_nat).to_csv(self.file_name_nat, index = False)
    
    
    def recordData_trial(self, agent, trial_i, pop_i, sim_info):
        
        output = [agent[pop_i].name, trial_i, agent[pop_i].ap, agent[pop_i].an, agent[pop_i].beta, agent[pop_i].Q[trial_i, 1], 
                    agent[pop_i].Q[trial_i, 0], agent[pop_i].p_history[trial_i], agent[pop_i].choice_history[trial_i], agent[pop_i].payoff_history[trial_i], agent[pop_i].delta[trial_i]]
        self.output_trial[pop_i*sim_info.Trial + trial_i, :] = output
        
        
    def writeCSV_trial(self, sim_info):
        
        f = open(self.file_name_trial, "a")
        writecsv = csv.writer(f, lineterminator="\n")
        
        for data_i in self.output_trial:    
            writecsv.writerow(data_i)
                
        f.close()
    
    
    def recordData_offline(self, agent, off_trial_i, off_pop_i):
        
        output = [agent[off_pop_i].name, off_trial_i, agent[off_pop_i].ap, agent[off_pop_i].an, agent[off_pop_i].beta, agent[off_pop_i].Q[Trial-1, 1], 
                    agent[off_pop_i].Q[Trial-1, 0], agent[off_pop_i].p_history[Trial-1], agent[off_pop_i].offline_choice[off_trial_i], agent[off_pop_i].offline_perform[off_trial_i]]
        self.output_offline[off_pop_i*sim_info.Offline_trial + off_trial_i, :] = output
    
    
    def writeCSV_offline(self, sim_info):
        
        f = open(self.file_name_offline, "a")
        writecsv = csv.writer(f, lineterminator="\n")
                    
        for data_i in self.output_offline:
            writecsv.writerow(data_i)
        
        f.close()
    
    
    def recordData_NaturalSelection(self, agent, sim_info):
        
        population = sim_info.Population
        trial = sim_info.Trial
        
        for i in range(population):
            
            rate_choice1 = np.mean(agent[i].choice_history)
            mean_payoff = np.mean(agent[i].payoff_history)
            mean_p = np.mean(agent[i].p_history)
            
            sum00 = 0
            sum01 = 0
            sum10 = 0
            sum11 = 0
            for j in range(1, trial):
                
                if (agent[i].choice_history[j-1] == 0) & (agent[i].choice_history[j] == 0):
                    
                    sum00 += 1
                
                elif (agent[i].choice_history[j-1] == 0) & (agent[i].choice_history[j] == 1):
                    
                    sum01 += 1
                
                elif (agent[i].choice_history[j-1] == 1) & (agent[i].choice_history[j] == 0):
                
                    sum10 += 1
                    
                else: # in the case of 11
                    sum11 += 1
            
            rate00 = sum00/trial
            rate01 = sum01/trial
            rate10 = sum10/trial
            rate11 = sum11/trial
            
            output = [i, mean_payoff, mean_p, rate_choice1, rate00, rate01, rate10 , rate11, agent[i].ap, agent[i].an, agent[i].beta]
            self.output_nat[i, :] = output

    
    def writeCSV_NaturalSelection(self, sim_info):
        
        f = open(self.file_name_nat, "a")
        writecsv = csv.writer(f, lineterminator="\n")
        
        for data_i in self.output_nat:
            writecsv.writerow(data_i)
        f.close()
    
    """
    def zipFile(self, sim_info):
        
        for i in range(len(self.file_name_trial)):
            
            if os.path.exists(self.file_name_trial[i]):
                
                rep = 'r' + f'{replication_i:0=3d}'
                
                output_zip_name = f"./{self.dir_name}/{rep}_{i}.zip"
                archive_name = f"./archive_dir/{rep}_{i}.csv"
                
                output_file = zipfile.ZipFile(output_zip_name, 'w', zipfile.ZIP_DEFLATED)
                output_file.write(filename=self.file_name_trial[i], arcname=archive_name)
                output_file.close()
                
            else:
                print("the file does not exist! stop process.")
                sys.exit()
            
            if os.path.exists(output_zip_name):
                
                os.remove(self.file_name_trial[i]) # delete original csv
            
            else:
                
                print("the zip file does not exist! stop process.")
                sys.exit()
    """
        
        
class Agent:
    
    """
     # This class contains information of Agent and methods worked to agents #
        setGenes: set parameters at 1st generation
        selectChoice_Norm: select choices according to softmax function and get payoffs generated from Normal
        selectChoice_Binom: select choices according to softmax function and get binary payoffs
        qlearning: update Q values
    """
    
    def __init__(self, name, trial, offline_trial, options):
        
        self.name = name
        self.Q = np.zeros((trial + 1, options))
        self.delta = np.empty(trial)
        self.choice_history = np.empty(trial, dtype = "int")
        self.payoff_history = np.empty(trial)
        self.offline_choice = np.empty(offline_trial, dtype = "int")
        self.offline_perform = np.empty(offline_trial)
        self.p_history = np.empty(trial)
        self.ap = 0
        self.an = 0
        self.beta = 0
        self.fitness = 0
     
    def initialize(self, name, trial, offline_trial, options):
        
        self.name = name
        self.Q = np.zeros((trial + 1, options))
        self.delta = np.empty(trial)
        self.choice_history = np.empty(trial, dtype = "int")
        self.payoff_history = np.empty(trial)
        self.offline_choice = np.empty(offline_trial, dtype = "int")
        self.offline_perform = np.empty(offline_trial)
        self.p_history = np.empty(trial)
        self.fitness = 0
    
    def setGenes(self, ap, an, beta):
        
        self.ap = ap
        self.an = an
        self.beta = beta
    
    def selectChoice_Norm(self, trial_i, sim_info, payoffgenerate):
    
        """
            Calculate Prob of choice 1 by softmax function
            payoffgenerate = ["positive", "negative", "both"]
        """
        
        type_op = sim_info.type_option
        round_thres = sim_info.round_threshold
        payoff_norm_mean1 = sim_info.payoff_norm_mean1
        payoff_norm_sigma1 = sim_info.payoff_norm_sigma1
        payoff_norm_mean2 = sim_info.payoff_norm_mean2
        payoff_norm_sigma2 = sim_info.payoff_norm_sigma2
        
        if type_op == 'Normal':
            #print('type option is Normal.')
            pass
        else:
            print('Error in type option. Check function selectChoice_Norm.')
            sys.exit()
        

        self.p_history[trial_i] = 1 /(1+ math.exp(- self.beta * (self.Q[trial_i, 1] - self.Q[trial_i, 0])))
        self.p_history = self.p_history.round(round_thres)
        p_trial = self.p_history[trial_i]
        
        self.choice_history[trial_i] = np.random.choice([1,0], p = [p_trial, 1.0 - p_trial])
        
        choice_trial = self.choice_history[trial_i]
        
        
        """
            Acquire Payoffs
        """
        
        if choice_trial == 1: # high-risk choice
            
            if payoffgenerate == "positive":
                
                payoff_trial = -1 # set any minus value
                
                while payoff_trial < 0: # loop until the value gets > 0
                    
                    payoff_trial = np.random.normal(payoff_norm_mean1, payoff_norm_sigma1)
                
                
            elif payoffgenerate == "negative":
                
                payoff_trial = 1 # set any plus value
                
                while payoff_trial > 0: # loop until the value gets < 0
                    
                    payoff_trial = np.random.normal(payoff_norm_mean1, payoff_norm_sigma1)
            
            
            elif payoffgenerate == "both":
                
                payoff_trial = np.random.normal(payoff_norm_mean1, payoff_norm_sigma1)
                
            else:
               # print error
               print("Error in payoffgenerate. Check function selectChoice_Norm.")
               sys.exit() 
            
            
        elif choice_trial == 0: # low-risk choice
            
            if payoffgenerate == "positive":
                
                payoff_trial = -1 # set any minus value
                
                while payoff_trial < 0: # loop until the value gets > 0
                    
                    payoff_trial = np.random.normal(payoff_norm_mean2, payoff_norm_sigma2)
                   
            elif payoffgenerate == "negative":
                
                payoff_trial = 1 # set any plus value
                
                while payoff_trial > 0: # loop until the value gets < 0
                    
                    payoff_trial = np.random.normal(payoff_norm_mean2, payoff_norm_sigma2)
                
            elif payoffgenerate == "both":
                
                payoff_trial = np.random.normal(payoff_norm_mean2, payoff_norm_sigma2)
                
            else:
               # print error
               print("Error in payoffgenerate. Check function selectChoice_Norm.")
               sys.exit() 
            
        else:
            # print error
            print("Error in the type of Choices. Check function SelectChoice_Norm.")
            sys.exit()
    
        self.payoff_history[trial_i] = payoff_trial
    
    
    def selectChoice_Binom(self, trial_i, sim_info):
        
        type_op = sim_info.type_option
        round_thres = sim_info.round_threshold
        payoff_bi1_x = sim_info.payoff_bi1_x
        payoff_bi1_y = sim_info.payoff_bi1_y
        payoff_bi1_prob = sim_info.payoff_bi1_prob
        
        payoff_bi2_x = sim_info.payoff_bi2_x
        payoff_bi2_x = sim_info.payoff_bi2_x
        payoff_bi2_prob = sim_info.payoff_bi2_prob
        
        if type_op == 'Binom':
            # print('type option is Binomial.')
            pass
        else:
            print('Error in type option. Check function selectChoice_Binom')
            sys.exit()
        
        self.p_history[trial_i] = 1 /(1+ math.exp(- self.beta * (self.Q[trial_i, 1] - self.Q[trial_i, 0])))
        self.p_history = self.p_history.round(round_thres)
        p_trial = self.p_history[trial_i]
        
        self.choice_history[trial_i] = np.random.choice([1,0], p = [p_trial, 1.0 - p_trial])
        
        choice_trial = self.choice_history[trial_i]
        
        if choice_trial == 1: # high-risk choice
            
            payoff_trial = np.random.choice([payoff_bi1_x, payoff_bi1_y], p = [payoff_bi1_prob, 1 - payoff_bi1_prob])
            
        elif choice_trial == 0: # low-risk choice
            
            payoff_trial = np.random.choice([payoff_bi2_x, payoff_bi2_y], p = [payoff_bi2_prob, 1 - payoff_bi2_prob])
            
        else:
            # print error
            print("Error in the kind of Choices. Check function selectChoice_Binom.")
            sys.exit()
            
        self.payoff_history[trial_i] = payoff_trial
        self.payoff_history = self.payoff_history.round(round_thres)
        
        
    def qlearning(self, trial_i):
        
        self.delta[trial_i] = self.payoff_history[trial_i] - self.Q[trial_i, self.choice_history[trial_i]]
        
        if self.delta[trial_i] >= 0: # if prediction error >= 0
        
            self.Q[trial_i+1, self.choice_history[trial_i]] = self.Q[trial_i, self.choice_history[trial_i]] + self.ap * self.delta[trial_i]
            self.Q[trial_i+1, 1-self.choice_history[trial_i]] = self.Q[trial_i, 1-self.choice_history[trial_i]]
    
        else: # if prediction error < 0
        
            self.Q[trial_i+1, self.choice_history[trial_i]] = self.Q[trial_i, self.choice_history[trial_i]] + self.an * self.delta[trial_i]
            self.Q[trial_i+1, 1-self.choice_history[trial_i]] = self.Q[trial_i, 1-self.choice_history[trial_i]]
        
        self.Q = self.Q.round(round_thres)
        
    
    def offlinePerform(self, off_trial_i, sim_info, payoffgenerate):
                    
        type_option = sim_info.type_option
        round_thres = sim_info.round_threshold
        
        # Normal
        payoff_norm_mean1 = sim_info.payoff_norm_mean1
        payoff_norm_sigma1 = sim_info.payoff_norm_sigma1
        payoff_norm_mean2 = sim_info.payoff_norm_mean2
        payoff_norm_sigma2 = sim_info.payoff_norm_sigma2
        
        # Binomial
        payoff_bi1_x = sim_info.payoff_bi1_x
        payoff_bi1_y = sim_info.payoff_bi1_y
        payoff_bi1_prob = sim_info.payoff_bi1_prob
        
        payoff_bi2_x = sim_info.payoff_bi2_x
        payoff_bi2_x = sim_info.payoff_bi2_x
        payoff_bi2_prob = sim_info.payoff_bi2_prob

        
        p_offline = self.p_history[Trial-1]
        
        if type_option == "Binom":
            
            choice_offline_i = np.random.choice([1,0], p = [p_offline, 1.0 - p_offline])
            
            if choice_offline_i == 1: # high-risk choice
                
                payoff_offline_i = np.random.choice([payoff_bi1_x, payoff_bi1_y], p = [payoff_bi1_prob, 1 - payoff_bi1_prob])
                
            elif choice_offline_i == 0: # low-risk choice
                
                payoff_offline_i = np.random.choice([payoff_bi2_x, payoff_bi2_y], p = [payoff_bi2_prob, 1 - payoff_bi2_prob])
            
        
        elif type_option == "Normal":
                                        
            choice_offline_i = np.random.choice([1,0], p = [p_offline, 1.0 - p_offline])
            
            if choice_offline_i == 1: # high-risk choice
                
                if payoffgenerate == "positive":
                    
                    payoff_offline_i = -1 # set any minus value
                    
                    while payoff_offline_i < 0: # loop until the value gets > 0
                        
                        payoff_offline_i = np.random.normal(payoff_norm_mean1, payoff_norm_sigma1)
                
                
                elif payoffgenerate == "negative":
                    
                    payoff_offline_i = 1 # set any plus value
                    
                    while payoff_offline_i > 0: # loop until the value gets < 0
                        
                        payoff_offline_i = np.random.normal(payoff_norm_mean1, payoff_norm_sigma1)
                    
                elif payoffgenerate == "both":
                    
                    payoff_offline_i = np.random.normal(payoff_norm_mean1, payoff_norm_sigma1)
                    
                else:
                   # print error
                   print("Error in payoffgenerate. Check function offlinePerform.")
                   sys.exit() 
                   
            elif choice_offline_i == 0: # low-risk choice
                
                if payoffgenerate == "positive":
                    
                    payoff_offline_i = -1 # set any minus value
                    
                    while payoff_offline_i < 0: # loop until the value gets > 0
                        
                        payoff_offline_i = np.random.normal(payoff_norm_mean2, payoff_norm_sigma2)
                
                
                elif payoffgenerate == "negative":
                    
                    payoff_offline_i = 1 # set any plus value
                    
                    while payoff_offline_i > 0: # loop until the value gets < 0
                        
                        payoff_offline_i = np.random.normal(payoff_norm_mean2, payoff_norm_sigma2)
                    
                elif payoffgenerate == "both":
                    
                    payoff_offline_i = np.random.normal(payoff_norm_mean2, payoff_norm_sigma2)
                    
                else:
                   # print error
                   print("Error in payoffgenerate. Check function offlinePerform.")
                   sys.exit() 
                
            else:
                # print error
                print("Error in the type of Choices. Check function offlinePerform.")
                sys.exit()
                
        else:
            
            # print error
            print("Error in the kind of Options. Check function offlinePerform.")
            sys.exit()
        
        self.offline_choice[off_trial_i] = choice_offline_i
        self.offline_perform[off_trial_i] = payoff_offline_i
        self.offline_perform = self.offline_perform.round(round_thres)
        
    
"""
 # Define Functions #
    Make agent instance, 
    Initialize Agent at start of generation, 
    Assign Genes to Agents, 
    roulette Wheel Selection,
    cause mutation
"""


def createAgent(sim_info):
    
    """
    create agent instance
    """
    
    population = sim_info.Population
    trial = sim_info.Trial
    offline_trial = sim_info.Offline_trial
    options = sim_info.number_options
    
    agent = np.empty(population, dtype = 'object')
    agent = [Agent(name = i, trial = trial, offline_trial = offline_trial, options = options) for i in range(population)]
    return agent
    
    
def assignGenes(agent, sim_info):
    
    """
    Assign Genes to each agent. Values are randomly generated.
    """
        
    ap_num = 100
    an_num = 100
    bt_num = 3
    
    ap_list = np.linspace(start=0.01, stop=1.0, num=ap_num, dtype='float')
    an_list = np.linspace(start=0.01, stop=1.0, num=an_num, dtype='float')
    bt_list = np.linspace(start=0.1, stop=0.4, num=bt_num, dtype='float')
    
    product_gene = list(itertools.product(ap_list, an_list, bt_list))
    if len(product_gene) == Population:
            
        for i,j in enumerate(product_gene):
        
            agent[i].setGenes(j[0], j[1], j[2])
        
    else:
        print('Error in length of Uniform Genes. Check function assignGenes')
        sys.exit()

   

#----------------------------------------------------------------------------# 
###-------------------------- Run simulations -----------------------------###



###-------------- Run on Local or parallel -------------------###


print('start:', datetime.datetime.now())

sim_info = SimulationInfo(parameters_dict) # make instance of SimulationInfo()
save_data = SaveData(sim_info) # make instance of SaveData()
save_data.writeParameterInfo(sim_info, parameters_dict)

save_data.createFile(sim_info, replication_i) # create csv files
agent = createAgent(sim_info) # make instance of createAgent()
assignGenes(agent, sim_info) 

population = sim_info.Population
trial = sim_info.Trial
#offline_trial = sim_info.Offline_trial

        
for pop_i in range(population):
    
    for trial_i in range(trial):
        
        agent[pop_i].selectChoice_Norm(trial_i, sim_info, "both") # Selection & get Payoffs # [Binom, Norm]
        agent[pop_i].qlearning(trial_i) # Q Learning
        save_data.recordData_trial(agent, trial_i, pop_i, sim_info) # write Results

"""
for off_pop_i in range(population):
    
    for off_trial_i in range(offline_trial):
        
        agent[off_pop_i].offlinePerform(off_trial_i, sim_info, 'both') # compute offline Performance
        save_data.recordData_offline(agent, off_trial_i, off_pop_i, sim_info) 
"""


save_data.recordData_NaturalSelection(agent, sim_info) # write csv of each fitness

save_data.writeCSV_NaturalSelection(sim_info) 
#save_data.writeCSV_trial(sim_info)
#save_data.writeCSV_offline(sim_info)

            
#save_data.zipFile(sim_info) # Zip trial_file
    
print('finish:', datetime.datetime.now())

    
###---------------------------------------------------###
