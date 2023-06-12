# Risk Preference and Evolution of Reinforcement Learning

This repository contains the scripts to run simulations and analyze the data in

Homma, S. & Takezawa, M. (in prep) Risk preference as an outcome of evolutionarily adaptive learning mechanisms: an evolutionary simulation under diverse risky environments.

  
## Simulation Code

  
All of codes are assumed to be run in GNU parallel.
  
  
When you run a code, you should specify parameters by `argv` (command line).
  
  
For example, to run BasicEvolQL_parallel.cpp in GNU parallel, 


```shell
clang++ -O2 BasicEvolQL_parallel.cpp -o BasicEvolQL
```
  
 ```shell
parallel ./BasicEvolQL ::: {1..10} ::: 10000 ::: 5000 ::: 500 ::: 40 ::: 30 25 20 15 10 ::: 20 ::: 5
```

  
 The simulation codes generate a save directory in a path Analysis/Results/. The save directory is usually specified by a suffix '_YYYYMMDDHHMMSS', which is used for reading data in the analysis script.
  
  
 - BasicEvolQL_parallel.cpp : A script to run a single task simulation of RL model with the asymmetric learning rates
   
   
 - BasicQL.py : A script to run the alpha poisitve and alpha negative effect on behavior (for making a heatmap)
   
   
 - generate_RandMultiTask.cpp : A script to generate a task set for multiple-task simulations. You should run this script before running the multiple-task simulations below.
   
   
 - EvolQL_RandMultiTask_parallel.cpp : A script to run a multiple-task simulation of RL model with the asymmetric learning rates
    
   
 - EvolQL_RandMultiTask_parallel_singleLR.cpp : A script to run a multiple-task simulation of RL model with the single learning rate
  

 - EvolQL_RandMultiTask_parallel_Perseverance.cpp : A script to run a multiple-task simulation of RL model with perseverance parameter



## Analysis Script


- Figures_for_Paper.Rmd : An analysis script to get figures in the main text.
  
  
- Figures_for_Supplementary.Rmd : An analysis script to get figures in the supplementary information.


- _functions_common and _functions_project : directories which contains .R functions. These codes are loaded in the .Rmd file.

