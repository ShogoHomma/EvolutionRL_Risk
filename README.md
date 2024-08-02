# Risk Preference and Evolution of Reinforcement Learning

This repository contains the scripts to run simulations and analyze the data in

Homma S, Takezawa M (2024) Risk preference as an outcome of evolutionarily adaptive learning mechanisms: An evolutionary simulation under diverse risky environments. PLoS ONE 19(8): e0307991. https://doi.org/10.1371/journal.pone.0307991


## Data

Data is available at OSF https://osf.io/skd9r/. (Because the data size was large, the upload to this repository was rejected.)

When analyzing the data with the analysis scripts, place the zipped data in Data directory, and unzip them.


## Simulation Code
  
  
When you run a code, you should specify parameters as command line arguments.

For example, to run `BasicEvolQL_parallel.cpp`, the parameter settings below,

- Population Size : 10000
- Generation : 5000
- Trial : 500
- Task : Normal(40, 20), Normal(20, 5)

```
clang++ -O2 BasicEvolQL_parallel.cpp -o BasicEvolQL
```

```
./BasicEvolQL 1 10000 5000 500 40 20 20 5
```

The first argument `1` is a run (replication) number, which is given as a vector when running in GNU parallel.

When you run a multiple-task simulation, you should run `generate_RandMultiTask.cpp`, and then `EvolQL_RandMultiTask_*.cpp`. 

If you run in the parameter settings:

- Population Size : 1000
- Generation : 3000
- Trial : 500
- Number of tasks : 4
- Number of risk-seeking tasks : 2
- The normal distribution which generates random tasks : Normal(0, 20) (usually fixed)
- The range of SD : [5, 30] (usually fixed)

``` 
./generate_RandMultiTask 1000 3000 500 4 2 100 0 20 30 5
```

```
./EvolQL_RandMultiTask_parallel 1 1000 3000 500 4 2 100 "path_to_save_directory"
```

The last argument should be the path to save directory generated by `generate_RandMultiTask.cpp`. You should copy the first line in `_path_info.txt` in the generated directory.

  
 The simulation codes generate a save directory in a path Analysis/Results/. The save directory is usually specified by a suffix '_YYYYMMDDHHMMSS', which is used for reading data in the analysis script.
  
  
 - `BasicEvolQL_parallel.cpp` : A script to run a single task simulation of RL model with the asymmetric learning rates
   
   
 - `BasicQL.py` : A script to run a simulation of the alpha poisitve and alpha negative effect on behavior (for making a heatmap)
   
   
 - `generate_RandMultiTask.cpp` : A script to generate a task set for multiple-task simulations. You should run this script before running the multiple-task simulations below.
   
   
 - `EvolQL_RandMultiTask_parallel.cpp` : A script to run a multiple-task simulation of RL model with the asymmetric learning rates
    
   
 - `EvolQL_RandMultiTask_parallel_singleLR.cpp` : A script to run a multiple-task simulation of RL model with the single learning rate
  
  
 - `EvolQL_RandMultiTask_parallel_Perseverance.cpp` : A script to run a multiple-task simulation of RL model with perseverance parameter



## Analysis Script


- `Figures_for_Paper.Rmd` : An analysis script to get figures in the main text.
  
  
- `Figures_for_Supplementary.Rmd` : An analysis script to get figures in the supporting information.


- _functions_common and _functions_project : directories which contain .R functions. These codes are loaded and used in the .Rmd file.

