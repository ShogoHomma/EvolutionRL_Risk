#include <vector>
#include <string>
#include <random>
#include <sstream>
#include <iostream>
#include <fstream> 
#include <iomanip>
#include <algorithm>
#include <exception>
#include <time.h>
#include <ctime>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

// parameters for running setting
const char Run_condition[10] = "parallel";
const int Replication = 10;

//---------------------- //
// The values of these parameters are put into by command line args

string script_name;
int replication_i;
int PopSize;
int Generation;
int Trial;
int Task_N;

// parameters for task setting
// a vector which contains a sequence of mean or sd of tasks (e.g., "20,10,5")
// the number of elements of a vector shoule be equal to Task_N
string pay1; // payoff_norm_mean1
string pay2; // payoff_norm_sigma1
string pay3; // payoff_norm_mean2
string pay4; // payoff_norm_sigma2

vector<double> payoff_norm_mean1_list; // risky option
vector<double> payoff_norm_sigma1_list; 
vector<double> payoff_norm_mean2_list; // sure option
vector<double> payoff_norm_sigma2_list;

//---------------------- //


// parameters for Natural Selection & mutation

const double base_fitness = 20.0;
const double mut_ap = 0.01;
const double mut_an = 0.01;
const double mut_bt = 0.005;
const double bt_upper = 0.5;

const string init_ap = "R"; // {"R" : Random, "F" : Fixed}
const double fixed_ap = 0.5; // if init_ap = "F", fixed_ap is assigned to all agents. 
const string init_an = "R";
const double fixed_an = 0.5;
const string init_bt = "R";
const double fixed_bt = 0.2;

// This function was cited from https://www.sejuku.net/blog/49378
// This function is used for spliting the character lists by the specific character (given by arg del). 
// In this simulation, it is used for spliting the payoff_structure by character ","
// e.g., "20,10,5" is converted to "20", "10", "5". and they are later converted to float 20, 10, 5 by atof() functon.
vector<string> split(string str, char del) {
    
    int first = 0;
    int last = str.find_first_of(del); //find_first_of() は先頭から文字列がどこにあるのかを検索
    //cout << "first 1: " << first << endl;
    //cout << "last 1: " << last << endl;

    vector<string> result;
    
    while (first < str.size()) { // strの長さになるまでループ
        
        string subStr(str, first, last - first); // strをfirstからlastまでで切り取り、subStrとする
        //cout << "subStr: " << subStr << endl;

        result.push_back(subStr);
        
        first = last + 1;
        last = str.find_first_of(del, first); // 2個目の引数は、どこから検索を開始するか（デフォは0）
        // cout below for checking how this function works, never mind (you can delete them). 
        //cout << "first 2: " << first << endl; // 1回目 9, 2回目19(最後) 
        //cout << "last 2: " << last << endl; // 2回目 -1, 2回目 -1

        if (last == string::npos) {
            // nposはfindで見つからなかった時の返り値で、値としては-1
            // つまり、lastが何もなかった場合（=終わりの場合）, 最後の場所(size)がlastに入る
            last = str.size();
        }
    }
 
    return result;
}

// classes
class Agent{

    int ID;
    vector<vector<double> > Qvalue_list; // two dimension, store by row major
    vector<vector<int> > Choice_list; // two dimension, store by row major
    vector<vector<double> > Payoff_list; // two dimension, store by row major
    vector<vector<double> > Delta_list; // two dimension, store by row major
    vector<vector<long double> > Pvalue_list; // two dimension, store by row major
    vector<double> mean_payoff; 
    vector<double> mean_p; 
    double fitness;
    int offspring;
    vector<double> rateChoice1; // rate of choosing option 1 (risky option)
    vector<double> rate00; // transition rate from option 0 to 0 (stay 0)
    vector<double> rate01; // transition rate from option 0 to 1
    vector<double> rate10; // transition rate from option 1 to 0
    vector<double> rate11; // transition rate from option 1 to 1 (stay 1)

    // genes 
    double ap;
    double an;
    double bt;

    public:

        // constructor
        Agent(); 

        // initialize functions
        void initialize();

        // defining funcitions for displaying / replacing the above values
        void set_ID(int n) {ID = n;}
        int return_ID() {return ID;}
        double return_Choice(int task_i, int trial_i) {return Choice_list[task_i][trial_i];}
        double return_Payoff(int task_i, int trial_i) {return Payoff_list[task_i][trial_i];}
        double return_Delta(int task_i, int trial_i) {return Delta_list[task_i][trial_i];}
        double return_P(int task_i, int trial_i) {return Pvalue_list[task_i][trial_i];}
        double return_Q1(int trial_i) {return Qvalue_list[trial_i][1];}
        double return_Q0(int trial_i) {return Qvalue_list[trial_i][0];}

        void set_Genes(double alp, double aln, double beta) {ap = alp; an = aln; bt = beta; }
        void show_Genes() {cout << "ap: " << ap << "; an: " << an << "; bt: " << bt << endl;}
        double return_ap() {return ap;}
        double return_an() {return an;}
        double return_bt() {return bt;}

        void calculate_MeanPayoff(int task_i);
        double return_MeanPayoff(int task_i) {return mean_payoff[task_i];}
        void calculate_rateChoice1(int task_i);
        double return_rateChoice1(int task_i) {return rateChoice1[task_i];}
        void calculate_MeanP(int task_i);
        double return_MeanP(int task_i) {return mean_p[task_i];}
        void calculate_rateTransition(int task_i);
        double return_rate00(int task_i) {return rate00[task_i];}
        double return_rate01(int task_i) {return rate01[task_i];}
        double return_rate10(int task_i) {return rate10[task_i];}
        double return_rate11(int task_i) {return rate11[task_i];}
        void set_Fitness(double n) {fitness = n;}
        double return_Fitness() {return fitness;}
        void set_Offspring(int n) {offspring = n;}
        int return_Offspring() {return offspring;}
        
        // defining select option & Q-learning functions
        void select_Option(int task_i, int trial_i, vector<double> &payoff_structure);
        void qLearning(int task_i, int trial_i);

};

// constructor of Agent
Agent::Agent()
{
    // initialize each member lists
    Qvalue_list.assign(Trial+1, vector<double>(2, 0.0));
    Choice_list.assign(Task_N, vector<int>(Trial, 0));
    Payoff_list.assign(Task_N, vector<double>(Trial,0.0));
    Delta_list.assign(Task_N, vector<double>(Trial, 0.0));
    Pvalue_list.assign(Task_N, vector<long double>(Trial, 0.0));

    // initialize other params
    fitness = 0.0;
    offspring = 0;
    mean_payoff.assign(Task_N, 0.0);
    mean_p.assign(Task_N, 0.0);
    rateChoice1.assign(Task_N, 0.0);
    rate00.assign(Task_N, 0.0);
    rate01.assign(Task_N, 0.0);
    rate10.assign(Task_N, 0.0);
    rate11.assign(Task_N, 0.0);

}

// functions to initialize // maybe a redundant function
void Agent::initialize() {

    // initialize each member lists
    Qvalue_list.assign(Trial+1, vector<double>(2, 0.0));
    Choice_list.assign(Task_N, vector<int>(Trial, 0));
    Payoff_list.assign(Task_N, vector<double>(Trial,0.0));
    Delta_list.assign(Task_N, vector<double>(Trial, 0.0));
    Pvalue_list.assign(Task_N, vector<long double>(Trial, 0.0));

    // initialize other params
    fitness = 0.0;
    offspring = 0;
    mean_payoff.assign(Task_N, 0.0);
    mean_p.assign(Task_N, 0.0);
    rateChoice1.assign(Task_N, 0.0);
    rate00.assign(Task_N, 0.0);
    rate01.assign(Task_N, 0.0);
    rate10.assign(Task_N, 0.0);
    rate11.assign(Task_N, 0.0);

}

// This function for calculating agent' mean payoff within a task
// The calculation for a mean payoff of all tasks (payoff_mean_all) is conducted later 
void Agent::calculate_MeanPayoff(int task_i) {

    double sum = 0.0;
    for (int i = 0; i < Trial; ++i)
        sum += Payoff_list[task_i][i];

    mean_payoff[task_i] = sum/Trial;
    //cout << "mean : " << mean_payoff << endl;

}

void Agent::calculate_rateChoice1(int task_i) {

    int sum = 0;
    for (int i = 0; i < Trial; ++i)
        sum += Choice_list[task_i][i];

    rateChoice1[task_i] = (double)sum/Trial;

}

void Agent::calculate_MeanP(int task_i) {

    double sum = 0.0;
    for (int i = 0; i < Trial; ++i)
        sum += Pvalue_list[task_i][i];
    
    mean_p[task_i] = (double)sum/Trial;

}

// caluculate frequency of choice transition
void Agent::calculate_rateTransition(int task_i) {

    int sum00 = 0; int sum01 = 0; int sum10 = 0; int sum11 = 0;
    // start from trial 1 (not from trial 0)
    for (int i = 1; i < Trial; ++i) {

        // choice transition 0 to 0 (stay 0)
        if (Choice_list[task_i][i-1] == 0 && Choice_list[task_i][i] == 0) {
            sum00 += 1;
        }
        // choice transition 0 to 1
        else if (Choice_list[task_i][i-1] == 0 && Choice_list[task_i][i] == 1) {
            sum01 += 1;
        }
        // choice transition 1 to 0
        else if (Choice_list[task_i][i-1] == 1 && Choice_list[task_i][i] == 0) {
            sum10 += 1;
        }
        // choice transition 1 to 1 (stay 1)
        else if (Choice_list[task_i][i-1] == 1 && Choice_list[task_i][i] == 1) {
            sum11 += 1;
        }
    
    }

    rate00[task_i] = (double)sum00/(Trial-1);
    rate01[task_i] = (double)sum01/(Trial-1);
    rate10[task_i] = (double)sum10/(Trial-1);
    rate11[task_i] = (double)sum11/(Trial-1);    

}

// select_Option function (calculate p value from Qvalues and select a option, get payoffs)
// payoff_sturucture should be given by pointer as an argument
void Agent::select_Option(int task_i, int trial_i, vector<double> &payoff_structure) {

    // set random seed every trial
    random_device rd;
    mt19937_64 get_rand_mt(rd());

    //cout << "\ntrial : " << trial_i << "\n" << endl;

    double payoff_norm_mean1 = payoff_structure[0]; // risky option
    double payoff_norm_sigma1 = payoff_structure[1];
    double payoff_norm_mean2 = payoff_structure[2]; // sure option
    double payoff_norm_sigma2 = payoff_structure[3];

    // calculate p value by Softmax function
    // Qvalue_list[trial_i][1]: prob of choosing option 1 (risky option)
    // Qvalue_list[trial_i][0]: prob of choosing option 0 (sure option)
    Pvalue_list[task_i][trial_i] = 1/(1 + exp(- bt * (Qvalue_list[trial_i][1]-(Qvalue_list[trial_i][0]))));
    //cout << "p value is : " << Pvalue_list[trial_i] << endl;

    // create a binomial distirbution 
    binomial_distribution<int> bnm_dist(1, Pvalue_list[task_i][trial_i]);

    // sample either choice 1 or 0 from binomial distribution
    int choice = bnm_dist(get_rand_mt);
    Choice_list[task_i][trial_i] = choice;
    //cout << "Choice 1 or 0: " << Choice_list[trial_i] << endl;

    // get payoffs depending on the choice
    // if the agent chose risky option
    if (Choice_list[task_i][trial_i] == 1) {
        
        // create normal distribution to generate payoffs
        normal_distribution<> nrml_dist_risk(payoff_norm_mean1, payoff_norm_sigma1);
        double payoff = nrml_dist_risk(get_rand_mt);

        Payoff_list[task_i][trial_i] = payoff;
        //cout << "Payoff risky option: " << Payoff_list[trial_i] << endl;
        
    } 
    // if the agent chose sure option
    else if (Choice_list[task_i][trial_i] == 0) {

        // create normal distribution to generate payoffs
        normal_distribution<> nrml_dist_sure(payoff_norm_mean2, payoff_norm_sigma2);
        double payoff = nrml_dist_sure(get_rand_mt);

        Payoff_list[task_i][trial_i] = payoff;
        //cout << "Payoff sure option: " << Payoff_list[trial_i] << endl;

    }
    // otherwise, display error
    else {

        terminate();
    
    }

}

// update Q values by calculating reward prediction error
void Agent::qLearning(int task_i, int trial_i) {

    double delta = Payoff_list[task_i][trial_i] - Qvalue_list[trial_i][Choice_list[task_i][trial_i]];
    Delta_list[task_i][trial_i] = delta;
    //cout << "delta : " << delta << endl;

    // if RPE is positive
    if (delta >= 0) {

        Qvalue_list[trial_i + 1][Choice_list[task_i][trial_i]] = Qvalue_list[trial_i][Choice_list[task_i][trial_i]] + ap * delta; // update for chosen option 
        Qvalue_list[trial_i + 1][1 - Choice_list[task_i][trial_i]] = Qvalue_list[trial_i][1 - Choice_list[task_i][trial_i]]; // update for non-chosen option 

    }
    // if RPE is negative
    else if (delta < 0) {

        Qvalue_list[trial_i + 1][Choice_list[task_i][trial_i]] = Qvalue_list[trial_i][Choice_list[task_i][trial_i]] + an * delta; // update for chosen option 
        Qvalue_list[trial_i + 1][1 - Choice_list[task_i][trial_i]] = Qvalue_list[trial_i][1 - Choice_list[task_i][trial_i]]; // update for non-chosen option 

    }
    // otherwise, display error
    else {

        terminate();

    }

    // display Q values
    //cout << "Q value 0 : " << Qvalue_list[trial_i + 1][0] << endl;
    //cout << "Q value 1 : " << Qvalue_list[trial_i + 1][1] << endl;

}

class saveData {

    char time_stamp_str[100]; // time to start
    time_t start_time;
    string path_name;
    string folder_name;
    string file_name_params;
    string file_name_trial; // to record first & last generation
    vector<string> file_name_trial_list; // to record all generation
    string file_name_gene;
    string file_name_nat;

    string system_order_st; // for file copy order of system()
    
    // pointers to output lists (define in constructor)
    vector<vector<double> > output_trial; // 2 dimension list for by-trial data: 14 columns
    vector<vector<double> > output_nat; // 2 dimension list for by-agent data: 16 columns (full ver.) or 7 columns (short ver.)
    vector<double> output_gene; // 1 dimension list for by-generation data: 8 columns
    
    public:

        int show_gene_time[10]; // to record time
        
        saveData(); // constructor

        // functions to create directory and csv files to save
        void create_dir();
        void create_file();
        
        // functions to record data to lists
        void recordData_trial(Agent agent, int replication_i, int gene_i, int pop_i, int task_i, int trial_i);
        void recordData_nat(Agent agent, int replication_i, int gene_i, int pop_i, int task_i);
        void recordData_gene(Agent agent, int replication_i, int gene_i, int pop_i);
        
        // functions to write output to files
        void writeParameterInfo();
        void writeTaskInfo();
        void writeElapsedTime(int replication_i, int gene_i, int check_i);
        void writeCSV_trial(int gene_i);
        void writeCSV_nat();
        void writeCSV_gene();

};

// constructor of saveData
saveData::saveData()
{
    // initialize start time stamp
    time_t now = time(NULL); struct tm tm; 
    strftime(time_stamp_str, sizeof(time_stamp_str), "%Y%m%d%H%M%S", localtime_r(&now, &tm));
    cout << "\n--- start time is ... " << time_stamp_str << endl;
    
    start_time = now;
    //cout << "start time : " << start_time << endl;

    // set generation number to record clocks
    // if Generation = 100, ten_perc_gene = 10, then, show_gene_time = [9, 19, ..., 99](record data every 10 generation)
    int ten_perc_gene = Generation/10;
    for (int i = 0; i < 10; ++i) {

        if (i == 0) {
            show_gene_time[i] = ten_perc_gene - 1;
        }
        else {
            show_gene_time[i] = show_gene_time[i - 1] + ten_perc_gene;
        }
    }

    // initialize output
    output_trial.assign(PopSize*Task_N*Trial, vector<double>(14, 0.0)); // 14 columns
    //output_nat.assign(PopSize, vector<double>(16, 0.0)); // full output ver. 16 columns
    output_nat.assign(PopSize*Task_N, vector<double>(7, 0.0)); // simple output ver. 7 columns
    output_gene.assign(8, 0.0); // 8 columns

}

// function to create directory to make files
void saveData::create_dir() {

    path_name = "./Analysis/Results";
    stringstream ss;
    
    // info about generation, popSize, trial
    ss << "G" << setw(5) << setfill('0') << Generation;
    folder_name += ss.str(); ss.str(""); ss.clear(stringstream::goodbit);
        // ss.str(""); ss.clear(stringstream::goodbit); -> clear baffer & clear the state of ss
    ss << "_N" << setw(5) << setfill('0') << PopSize;
    folder_name += ss.str(); ss.str(""); ss.clear(stringstream::goodbit);
    ss << "_T" << setw(4) << setfill('0') << Trial;
    folder_name += ss.str(); ss.str(""); ss.clear(stringstream::goodbit);

    // info about how to assign initial genes (Random or Fixed), and type of choices (Norm or Binom)
    folder_name += "_ap"; folder_name += init_ap; folder_name += "_an"; folder_name += init_an; 
    folder_name += "_bt"; folder_name += init_bt;
    folder_name += "_CN"; 
    ss << "_taskN" << setw(2) << setfill('0') << Task_N;
    folder_name += ss.str(); ss.str(""); ss.clear(stringstream::goodbit);
    
    // info about environmental change
    folder_name += "_V";
    

    /* ------------------------------------------------------------- */
    // create path to save files
     path_name += "/"; path_name += folder_name; path_name += '_'; path_name += time_stamp_str; path_name += "/"; 

    // create the directory
    if (mkdir("./Analysis", 0777) == 0) 
        cout << "---- create the directory ./Analysis/" << endl; 
    if (mkdir("./Analysis/Results", 0777) == 0) 
        cout << "---- create the directory ./Analysis/Results ." << endl; 
    if (mkdir(path_name.c_str(), 0777) == 0) 
        cout << "---- create the directory ..." << endl << path_name.c_str() << endl; 
    
    // Copy this .cpp file to the directory
    system_order_st = "cp EvolQL_multiTask_parallel.cpp ";
    system_order_st += path_name;

    system(system_order_st.c_str()); // run copy

}

// function to create csv files to save data
void saveData::create_file() {

    file_name_params += path_name; file_name_params += "_Parameters_memo.txt"; 

    stringstream rep;
    rep << "r" << setw(3) << setfill('0') << replication_i; 
    file_name_gene += path_name; file_name_gene += rep.str(); file_name_gene += "_gene.csv";
    file_name_nat += path_name; file_name_nat += rep.str(); file_name_nat += "_nat.csv";

    // ------------ set each file for output ---------- //
    ofstream ofs0(file_name_params);
    if (!ofs0) {
        cout << "cannot open file!. Check file_name_params ..." << endl << file_name_params << endl;
    }
    else {
        cout << "succeed! created " << file_name_params << endl;
    }
    
    /*
    // at the moment, we do not create file_trials (because they are large files & unnecessary for analysis)
    // create by-trial csv files
    // for file_trial, we devide into 10 files (because they are large files)
    for (int file_trial_i = 0; file_trial_i < 10; ++file_trial_i) {

        string file_num = to_string(file_trial_i);
        file_name_trial += path_name; file_name_trial += rep.str();
        file_name_trial += "_trial_"; file_name_trial += file_num; file_name_trial += ".csv";
        ofstream ofs1(file_name_trial);
        if (!ofs1) {
            cout << "cannot open file!. Check file_name_trial ..." << endl << file_name_trial << endl;
        }
        else {
            cout << "succeed! created " << file_name_trial << endl;
            // 14 columns, write the header names
            ofs1 << "replication,generation,agent,task,trial,ap,an,bt,Q1,Q0,p,choice,get_payoff,delta" << endl;
        }

        file_name_trial_list.push_back(file_name_trial);
        file_name_trial.clear(); // reset 

    }
    */

    // create by-trial csv files for the first & last generation
    file_name_trial += path_name; file_name_trial += rep.str(); file_name_trial += "_trial.csv";
    ofstream ofs1(file_name_trial);
    if (!ofs1) {
        cout << "cannot open file!. Check file_name_trial ..." << endl << file_name_trial << endl;
    }
    else {
        cout << "succeed! created " << file_name_trial << endl;
        // 14 columns, write the header names
        ofs1 << "replication,generation,agent,task,trial,ap,an,bt,Q1,Q0,p,choice,get_payoff,delta" << endl;
    }
    
    // create by-agent csv file (_nat.csv)
    ofstream ofs2(file_name_nat);
    if (!ofs2) {
        cout << "cannot open file!. Check file_name_nat ..." << endl << file_name_nat << endl;
    }
    else {
        cout << "succeed! created " << file_name_nat << endl;
        
        // full output ver. (16 columns), write the header names
        //ofs2 << "replication,generation,agent,task,ancestor,mean_payoff,mean_p,fitness,rate_choice1,rate0to0,rate0to1,rate1to0,rate1to1,ap,an,bt" << endl;
        
        // short output ver. (7 columns), write the header names
        // cf. BasicEvolQL_parallel.cpp: This short ver mainly omit 3 columns: replication, generation, and agent
        // because they make file size too large and redundant info.
        ofs2 << "task,mean_payoff,mean_p,rate_choice1,ap,an,bt" << endl; 
    }

    // create by-generation csv file (_gene.csv)
    ofstream ofs3(file_name_gene);
    if (!ofs3) {
        cout << "cannot open file!. Check file_name_gene ..." << endl << file_name_gene << endl;
    }
    else {
        cout << "succeed! created " << file_name_gene << endl;
        // 8 columns, write the header names
        ofs3 << "replication,generation,mean_ap,mean_an,mean_bt,sd_ap,sd_an,sd_bt" << endl;
    }

}

// function to write Sim Parameters info to text file 
void saveData::writeParameterInfo() {

    ofstream fout_params(file_name_params, ios::app);
    fout_params << "Script Name : " << script_name << endl << "Run_condition : " << Run_condition << endl << "Number of Replication : " << Replication << endl << \
     "Population : " << PopSize << endl << "Generation : " << Generation << endl << "Trial : " << Trial << endl << "Task_N: " << Task_N << endl << \
     "mutation_ap : " << mut_ap << endl << "mutation_an : " << mut_an << endl << "mutation_bt : " << mut_bt << endl << \
     "volatility_prob : " << "None" << endl << "volatility_env : " << "None" << endl << "type_change_option : " << "None" << endl << "volatility_sigma : " << "None" << endl << \
     "type_op : " << "Normal" << endl << "payoff_norm_mean1 : " << pay1 << endl << "payoff_norm_sigma1 : " << pay2 << endl << \
     "payoff_norm_mean2 : " << pay3 << endl << "payoff_norm_sigma2 : " << pay4 << endl << \
     "init_alpha_positive : " << init_ap << endl << "init_alpha_negative : " << init_an << endl << "init_inverse_temp : " << init_bt << endl << \
     "limit_inverse_temp : " << bt_upper << endl << "base_fitness : " << base_fitness << \
    endl << endl <<  "start : " << time_stamp_str << endl;

}

// This function writes task informations to file_name_params (_Parameters_memo.txt)
void saveData::writeTaskInfo() {

    ofstream fout_params(file_name_params, ios::app);

    for (int i = 0; i < Task_N; ++i) {

        fout_params << "Task " << i << ": N(" << payoff_norm_mean1_list[i] << ", " << payoff_norm_sigma1_list[i] << "), N(" << \
        payoff_norm_mean2_list[i] << ", " << payoff_norm_sigma2_list[i] << ")" << endl;

    }

}

// function to record elapsed time in 10% of generation
void saveData::writeElapsedTime(int replication_i, int gene_i, int check_i) {

    int gene_bool;
    char now_time_str[200];
    char elapsed_time_str[200];

    // if gene_i is 10% of Generation, record clock  
    time_t now = time(NULL); struct tm tm; 
    strftime(now_time_str, sizeof(now_time_str), "%Y-%m-%d %H:%M:%S", localtime_r(&now, &tm));
    cout << "show gene : " << show_gene_time[check_i] << ", time: " << now_time_str << endl;
    time_t elapsed = difftime(now, start_time);
    //strftime(elapsed_time_str, sizeof(elapsed_time_str), "%H:%M:%S", localtime_r(&elapsed, &tm));
    
    ofstream fout_params(file_name_params, ios::app);
    fout_params << "rep_i : " << replication_i << ", gene_i : " << gene_i << \
    ", end time : " << now_time_str << ", elapsed time : " << elapsed << " sec" << endl;
    
    start_time = now;

}

// function to record data to output_trial
void saveData::recordData_trial(Agent agent, int replication_i, int gene_i, int pop_i, int task_i, int trial_i) {

    double output[14] = {(double)replication_i, (double)gene_i, (double)agent.return_ID(), (double)task_i, (double)trial_i, \
                        agent.return_ap(), agent.return_an(), agent.return_bt(), \
                        agent.return_Q1(trial_i), agent.return_Q0(trial_i), (double)agent.return_P(task_i, trial_i), \
                        (double)agent.return_Choice(task_i,trial_i), agent.return_Payoff(task_i,trial_i), agent.return_Delta(task_i,trial_i)};
    
    for (int data_i = 0; data_i < 14; ++data_i) {
        output_trial[pop_i*Trial*Task_N + task_i*Trial + trial_i][data_i] = output[data_i];
    }

}

// function to record data to output_nat
void saveData::recordData_nat(Agent agent, int replication_i, int gene_i, int pop_i, int task_i) {

    /*
    // full output version (16 columns)
    double output[16] = {(double)replication_i, (double)gene_i, (double)agent.return_ID(), (double)task_i, (double)agent.return_Offspring(), agent.return_MeanPayoff(task_i), \
                        (double)agent.return_MeanP(task_i), agent.return_Fitness(), agent.return_rateChoice1(task_i), \
                        agent.return_rate00(task_i), agent.return_rate01(task_i), agent.return_rate10(task_i), agent.return_rate11(task_i), \
                        agent.return_ap(), agent.return_an(), agent.return_bt()};

    for (int data_i = 0; data_i < 16; ++data_i) {
        output_nat[pop_i][data_i] = output[data_i];
    }
    */

    // short output version (7 colums)
    double output[7] = {(double)task_i, agent.return_MeanPayoff(task_i),\
                        (double)agent.return_MeanP(task_i), agent.return_rateChoice1(task_i), agent.return_ap(), agent.return_an(), agent.return_bt()};
    
    for (int data_i = 0; data_i < 7; ++data_i) {
        output_nat[pop_i*Task_N + task_i][data_i] = output[data_i];
    }

}


// function to record data to output_gene
void saveData::recordData_gene(Agent agent, int replication_i, int gene_i, int pop_i) {

    // lists for calculating mean & SD of parameters
    static vector<double> ap_list;
    static vector<double> an_list;
    static vector<double> bt_list;
    if (pop_i == 0) { // only if the first agent, the lists are initialized (maybe unnecessary process)

        ap_list.assign(PopSize, 0.0);
        an_list.assign(PopSize, 0.0);
        bt_list.assign(PopSize, 0.0);

    }
    
    ap_list[pop_i] = agent.return_ap(); 
    an_list[pop_i] = agent.return_an(); 
    bt_list[pop_i] = agent.return_bt();
    //cout << "Name : " << pop_i << " ap : " << agent.return_ap() << " an : " << agent.return_an() << " bt : " << agent.return_bt() << endl;

    if (pop_i == PopSize - 1) {
        
        // display
        /*
        cout << endl;
        cout << "ap_list: ";
        for (int j = 0; j < PopSize; j++) 
            cout << ap_list[j] << ", ";
        cout << endl;
        cout << "an_list: ";
        for (int j = 0; j < PopSize; j++) 
            cout << an_list[j] << ", ";
        cout << endl;
        cout << "bt_list: ";
        for (int j = 0; j < PopSize; j++) 
            cout << bt_list[j] << ", ";
        cout << endl;
         */
        
        // calculate mean of parameters
        double mean_ap, mean_an, mean_bt;
        double sum_ap = 0.0; double sum_an = 0.0; double sum_bt = 0.0;

        for (int i = 0; i < PopSize; ++i) {

            //cout << "Name : " << i << " ap : " << ap_list[i] << " an : " << an_list[i] << " bt : " << bt_list[i] << endl;
            sum_ap += ap_list[i]; 
            sum_an += an_list[i]; 
            sum_bt += bt_list[i];

        }
    
        mean_ap = sum_ap/PopSize; mean_an = sum_an/PopSize; mean_bt = sum_bt/PopSize;
        //cout << "mean_ap: " << mean_ap << ", mean_an: " << mean_an << ", mean_bt: " << mean_bt << endl;

        // calculate SD of parameters
        double sd_ap, sd_an, sd_bt;
        double var_ap = 0.0; double var_an = 0.0; double var_bt = 0.0;

        for (int i = 0; i < PopSize; ++i) {
            
            var_ap += (ap_list[i] - mean_ap) * (ap_list[i] - mean_ap);
            var_an += (an_list[i] - mean_an) * (an_list[i] - mean_an);
            var_bt += (bt_list[i] - mean_bt) * (bt_list[i] - mean_bt);

        }
        
        sd_ap = sqrt(var_ap/PopSize); sd_an = sqrt(var_an/PopSize); sd_bt = sqrt(var_bt/PopSize);
        //cout << "sd_ap: " << sd_ap << ", sd_an: " << sd_an << ", sd_bt: " << sd_bt << endl;

        double output[8] = {(double)replication_i, (double)gene_i, mean_ap, mean_an, mean_bt, sd_ap, sd_an, sd_bt};
        
        for (int data_i = 0; data_i < 8; ++data_i) {
            output_gene[data_i] = output[data_i];
        }
        
    }
    
}

// --- functions to write csv (they are run every generation)

// function to data to trial csv
void saveData::writeCSV_trial(int gene_i) {

    /*
    // write trial Data (all generation)

    int file_index = gene_i/(Generation/10); // file_index : 0 ~ 9
    //cout << "file_index : " << file_index << endl;
    
    ofstream fout_trial(file_name_trial_list[file_index], ios::app);
    for (size_t data_i = 0; data_i < (PopSize*Trial*Task_N); ++data_i) {

        for (size_t data_j = 0; data_j < 14; ++data_j) {
            
            if (data_j != 13) {
                fout_trial << output_trial[data_i][data_j] << ", ";
            }
            else { // if last data
                fout_trial << output_trial[data_i][data_j] << endl;
            }

        }
        
    }
    */

    // write trial Data (first & last generation)
    ofstream fout_trial(file_name_trial, ios::app);
    for (size_t data_i = 0; data_i < (PopSize*Trial*Task_N); ++data_i) {

        for (size_t data_j = 0; data_j < 14; ++data_j) {
            
            if (data_j != 13) {
                fout_trial << output_trial[data_i][data_j] << ", ";
            }
            else { // if last data
                fout_trial << output_trial[data_i][data_j] << endl;
            }

        }
        
    }

}

// function to data to nat csv
void saveData::writeCSV_nat() {

    ofstream fout_nat(file_name_nat, ios::app);
    /*
    // full output ver.
    for (size_t data_i = 0; data_i < PopSize; ++data_i) {

        for (size_t data_j = 0; data_j < 16; ++data_j) {

            if (data_j != 15) {
                fout_nat << output_nat[data_i][data_j] << ", ";
            } 
            else { // if last data
                fout_nat << output_nat[data_i][data_j] << endl;
            }

        }
        
    }
    */
    
    // short output ver.
    for (size_t data_i = 0; data_i < (PopSize*Task_N); ++data_i) {

        for (size_t data_j = 0; data_j < 7; ++data_j) {

            if (data_j != 6) {
                fout_nat << output_nat[data_i][data_j] << ", ";
            } 
            else { // if last data
                fout_nat << output_nat[data_i][data_j] << endl;
            }

        }
        
    }

}

// function to data to gene csv
void saveData::writeCSV_gene() {

    ofstream fout_gene(file_name_gene, ios::app);
    for (size_t data_i = 0; data_i < 8; ++data_i) {
        
        if (data_i != 7) {
            fout_gene << output_gene[data_i] << ", ";
        }
        else { // if last data
            fout_gene << output_gene[data_i] << endl;
        }

    }
    
}


// ----------------- MAIN PROCESS ------------------- //
int main(int argc, char *argv[]){

    //------ get parameters set from command -----//

    script_name += argv[0];
    replication_i = atoi(argv[1]);
    PopSize = atoi(argv[2]);
    Generation = atoi(argv[3]);
    Trial = atoi(argv[4]);
    
    // parameters for choice setting
    char del = ',';
    pay1 = argv[5]; // strings of payoff_norm_mean1 (eg, "10,20,40")
    pay2 = argv[6]; // strings of payoff_norm_sigma1
    pay3 = argv[7]; // strings of payoff_norm_mean2
    pay4 = argv[8]; // strings of payoff_norm_sigma2

    // ------------------------------------------//
    
    for (const auto str_i : split(pay1, del)){
        payoff_norm_mean1_list.push_back(atof(str_i.c_str())); // add payoff_norm_mean1 (double vector)
    }
    
    for (const auto str_i : split(pay2, del)){
        payoff_norm_sigma1_list.push_back(atof(str_i.c_str())); // add payoff_norm_sigma1 (double vector)
    }
    
    for (const auto str_i : split(pay3, del)){
        payoff_norm_mean2_list.push_back(atof(str_i.c_str())); // add payoff_norm_mean2 (double vector)
    }
    
    for (const auto str_i : split(pay4, del)){
        payoff_norm_sigma2_list.push_back(atof(str_i.c_str())); // add payoff_norm_sigma2 (double vector)
    }

    Task_N = payoff_norm_mean1_list.size();

    // check equal to Task_N
    if ((payoff_norm_mean1_list.size() == Task_N) && (payoff_norm_sigma1_list.size() == Task_N) && (payoff_norm_mean2_list.size() == Task_N) && (payoff_norm_sigma2_list.size() == Task_N)) {
        // cout << "proper task setting." << endl;
    } else {
        cout << "Stop. task setting is improper." << endl;
        terminate(); // terminate programs
    }

    //initialzing random seed for shuffle using <random>
    //mt19937_64 get_rand_mt(time(0));
    random_device rd;
    mt19937_64 get_rand_mt_1(rd());
    
    saveData savedat;
    savedat.create_dir();
    savedat.create_file();
    savedat.writeParameterInfo();
    savedat.writeTaskInfo();
    
    // Agent objects
    Agent nAgent[PopSize];

    // assign initial Genes
    uniform_real_distribution<double> ap_dist(0.0, 1.0);
    uniform_real_distribution<double> an_dist(0.0, 1.0);
    uniform_real_distribution<double> bt_dist(0.0, bt_upper);

    for (int pop_i = 0; pop_i < PopSize; ++pop_i){

        double ap_i;
        double an_i;
        double bt_i;

        if (init_ap == "R") {
            ap_i = ap_dist(get_rand_mt_1);
        }
        else { // "Fixed"
            ap_i = fixed_ap;
        }

        if (init_an == "R") {
            an_i = an_dist(get_rand_mt_1);
        }
        else { // "Fixed"
            an_i = fixed_an;
        }
        
        if (init_bt == "R") {
            bt_i = bt_dist(get_rand_mt_1);
        }
        else { // "Fixed"
            bt_i = fixed_bt;
        }
        
        nAgent[pop_i].set_Genes(ap_i, an_i, bt_i);
        //cout << pop_i << ": ";
        //nAgent[pop_i].show_Genes();

    }
    
    // ---------------------------------------------------------------------------------- //
    // loop by gene_i
    for (int gene_i = 0; gene_i < Generation; ++gene_i) {

        //cout << "\n ---- Generation : " << gene_i << " ----- " << endl;

        // assign ID & initialize params
        for (int pop_i = 0; pop_i < PopSize; ++pop_i){
            nAgent[pop_i].set_ID(pop_i);
            nAgent[pop_i].initialize();
        }

        // ------------------------------------------ //
        // loop by pop_i: select_option + qlearning
        for (int pop_i = 0; pop_i < PopSize; ++pop_i) {

            //cout << "\n ---- Agent : " << pop_i << " ---- " << endl;

            // loop by task_i:  
            for (int task_i = 0; task_i < Task_N; ++task_i) {

                //cout << "\n ---- Task : " << task_i << " ---- " << endl;
                
                vector<double> payoff_structure;
                payoff_structure.push_back(payoff_norm_mean1_list[task_i]);
                payoff_structure.push_back(payoff_norm_sigma1_list[task_i]);
                payoff_structure.push_back(payoff_norm_mean2_list[task_i]);
                payoff_structure.push_back(payoff_norm_sigma2_list[task_i]);

                for (int trial_i = 0; trial_i < Trial; ++trial_i) {
                    
                    nAgent[pop_i].select_Option(task_i, trial_i, payoff_structure); // select options
                    nAgent[pop_i].qLearning(task_i, trial_i); // qlearning

                    // record by-trial data to a vector (output_trial)
                    if ((gene_i == 0)||(gene_i == Generation-1)) {
                        savedat.recordData_trial(nAgent[pop_i], replication_i, gene_i, pop_i, task_i, trial_i); // * WARNING * This line should be usually commented out because we dont need by-trial csv
                    }
                    
                }

            }

        }
        // ------------------------------------------ //
        // the end of loop pop_i

        // record by-generation data to a vector (output_gene)
        for (int pop_i = 0; pop_i < PopSize; ++pop_i) 
            savedat.recordData_gene(nAgent[pop_i], replication_i, gene_i, pop_i);


        //------------ Calculating Fitness ------------//
        //cout << "\nFitness\n" << endl;

        vector<double> payoff_base_vec;
        vector<double> vecFitness;
        for (int pop_i = 0; pop_i < PopSize; ++pop_i) {

            double payoff_all = 0.0;
            double payoff_mean_all = 0.0;
            //cout << "agent " << pop_i << ": " ;
            for (int task_i = 0; task_i < Task_N; ++task_i) {
                
                nAgent[pop_i].calculate_MeanPayoff(task_i);
                payoff_all += nAgent[pop_i].return_MeanPayoff(task_i);
                
                /*
                if (task_i != Task_N -1) {
                    cout << "task_i : " << task_i << " : " << nAgent[pop_i].return_MeanPayoff(task_i) << ", ";
                } else {
                    cout << "task_i : " << task_i << " : " << nAgent[pop_i].return_MeanPayoff(task_i) << endl;
                }
                */

            }

            // calculate agent's mean payoff of all tasks
            payoff_mean_all = payoff_all/Task_N;
            // always add base fitness 20 to mean_payoff of all agents
            payoff_base_vec.push_back(payoff_mean_all + base_fitness);

        }

        // get minimum payoff among all agents
        double min_payoff = *min_element(payoff_base_vec.begin(), payoff_base_vec.end());
        
        // if minimum payoff is negative, raise the bottom by fabs(min_payoff)
        if (min_payoff <= 0) {
            
            //cout << "\nMinimum payoff was negative.\n" << endl;
            for (int pop_i = 0; pop_i < PopSize; ++pop_i) {

                nAgent[pop_i].set_Fitness(payoff_base_vec[pop_i] + fabs(min_payoff) + 0.1);
                vecFitness.push_back(nAgent[pop_i].return_Fitness());
                //cout << "payoff_base_vec agent " << pop_i << " : " << payoff_base_vec[pop_i] << endl;
                //cout << "fitness agent " << pop_i << " : " << vecFitness[pop_i] << endl;

            }

        }
        // if minimum payoff is positive, continue to calculate without any adjustment
        else if (min_payoff > 0) {

            //cout << "\nMinimum payoff was positive.\n" << endl;
            for (int pop_i = 0; pop_i < PopSize; ++pop_i) {

                nAgent[pop_i].set_Fitness(payoff_base_vec[pop_i]);
                vecFitness.push_back(nAgent[pop_i].return_Fitness());
                //cout << "payoff_base_vec agent " << pop_i << " : " << payoff_base_vec[pop_i] << endl;
                //cout << "fitness agent " << pop_i << " : " << vecFitness[pop_i] << endl;

            }

        }

        //------------ Natural Selection ------------//
        // create a multinomial distribution (discrete_distribution), 
        // and determine the number of offspring an agent produces, proportional to its fitness

        random_device rd;
        mt19937_64 get_rand_mt_2(rd());
        discrete_distribution<> prob(vecFitness.begin(), vecFitness.end());

        // repeatedly sampling agents' ID from multinomial distribution
        // vecOffspring will contains IDs (0 ~ PopSize-1, allowing duplication) who are successful in producing offspring
        vector<int> vecOffspring;

        for (size_t n = 0; n < PopSize; ++n) {
            size_t result = prob(get_rand_mt_2);
            vecOffspring.push_back(result);
        }

        /*
        // displaying ID number of agents who succeeded in producing offspring
        for (int i = 0; i < PopSize; i++)
            cout << vecOffspring[i] << ", ";
        cout << endl;
        */

        // calculate the number of offsprings of each agents
        int bool_list[PopSize];
        int count_list[PopSize];
        // for agent i
        for (int pop_i = 0; pop_i < PopSize; ++pop_i) {
            
            int count = 0;
            //cout << "-- agent : " << pop_i << endl;
            for (int j = 0; j < PopSize; ++j) {
                
                // check whether agent i is found in vecOffspring or not
                bool_list[j] = vecOffspring[j] == pop_i; // return true (1) or not (0)
                // cout << bool_list[j];

                if (j == PopSize - 1) {
                    
                    for (int k = 0; k < PopSize; ++k)
                        count += bool_list[k];

                }
            
            }
            nAgent[pop_i].set_Offspring(count);
            //cout << " offspring : " << nAgent[pop_i].return_Offspring() << endl;

        }

        // record by-agent data to a vector (output_nat)
        for (int pop_i = 0; pop_i < PopSize; ++pop_i) {
            
            for (int task_i = 0; task_i < Task_N; ++task_i) {

                nAgent[pop_i].calculate_MeanP(task_i);
                nAgent[pop_i].calculate_rateChoice1(task_i);
                nAgent[pop_i].calculate_rateTransition(task_i);
                savedat.recordData_nat(nAgent[pop_i], replication_i, gene_i, pop_i, task_i);

            }
            
        }
        
        //--------------- Copy Genes ---------------//
        // new_ap[PopSize] .. temp list for storing gene ap
        //cout << "\nCopy genes.\n" << endl;

        double new_ap[PopSize];
        double new_an[PopSize];
        double new_bt[PopSize];
        for (int pop_i = 0; pop_i < PopSize; ++pop_i) {

            new_ap[pop_i] = nAgent[vecOffspring[pop_i]].return_ap();
            new_an[pop_i] = nAgent[vecOffspring[pop_i]].return_an();
            new_bt[pop_i] = nAgent[vecOffspring[pop_i]].return_bt();
            //cout << "pop:" << vecOffspring[pop_i] << "; ap: " << new_ap[pop_i] << "; an: " << new_an[pop_i] << "; bt: " << new_bt[pop_i] << endl;

        }

        // set random seed
        mt19937_64 get_rand_mt_3(rd());

        //--------------- Mutation ---------------//
        //cout << "\nMutation\n" << endl;

        for (int pop_i = 0; pop_i < PopSize; ++pop_i) {
        
            double next_ap_i;
            double next_an_i;
            double next_bt_i;
            normal_distribution<double> nrml_dist_ap(new_ap[pop_i], mut_ap);
            normal_distribution<double> nrml_dist_an(new_an[pop_i], mut_an);
            normal_distribution<double> nrml_dist_bt(new_bt[pop_i], mut_bt);

            // sample repeatedly as long as the value within parameter range
            if (init_ap == "R") {

                do {
                    next_ap_i = nrml_dist_ap(get_rand_mt_3);
                    //cout << "next ap : " << next_ap_i << endl;
                } while((next_ap_i < 0)||(next_ap_i > 1));

            } 
            else { // any fixed value

                next_ap_i = fixed_ap;

            }
            
            if (init_an == "R") {

                do {
                    next_an_i = nrml_dist_an(get_rand_mt_3);
                    //cout << "next an : " << next_an_i << endl;
                } while((next_an_i < 0)||(next_an_i > 1));
            
            }
            else { // any fixed value

                next_an_i = fixed_an;

            }

            if (init_bt == "R") {

                do {
                    next_bt_i = nrml_dist_bt(get_rand_mt_3);
                    //cout << "next bt : " << next_bt_i << endl;
                } while((next_bt_i < 0)||(next_bt_i > bt_upper));
            
            }
            else { // any fixed value

                next_bt_i = fixed_bt;

            }
            
            // set new genes to a new agent (offspring)
            nAgent[pop_i].set_Genes(next_ap_i, next_an_i, next_bt_i);
        }

        /*
        cout << "\n Display new agents genes\n" << endl;
        for (int pop_i = 0; pop_i < PopSize; pop_i++) {
            nAgent[pop_i].show_Genes();
        }
        */

        // ------ write output data to CSV files ----- //
        if ((gene_i == 0) || (gene_i == Generation-1)) {
            cout << "write CSV trial at gene : " << gene_i << endl;
            savedat.writeCSV_trial(gene_i); // * WARNING * This line should be usually commented out because we dont need by-trial csv
        }
        
        savedat.writeCSV_gene();
        savedat.writeCSV_nat();

        // ----------------------------------- //

        int gene_bool; // record finishing time every 10% of all generations
        for (int check_i = 0; check_i < 10; ++check_i) {
            
            // check whether the current generation matches with 10% of generation
            gene_bool = gene_i == savedat.show_gene_time[check_i]; // gene_bool = true (1) or not (0)
            if (gene_bool == 1) { // if it is 10% of generations, let's record the time
                
                savedat.writeElapsedTime(replication_i, gene_i, check_i);

            }
        }

    // -------------------------------------------------------------------------------- //
    // the end of loop gene_i 
    }

    return 0;

}
