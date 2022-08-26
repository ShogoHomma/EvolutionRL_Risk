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

const int SimN = 100;

//---------------------- //
// The values of these parameters are put into by command line args

string script_name;
//int replication_i;
int PopSize;
int Generation;
int Trial;
int Task_N;
int Riskseek_N;
int Taskgroup_N;

double generate_norm_mean;
double generate_norm_sigma;
double norm_sd1_upper;
double norm_sd1_lower; // set more than norm2_sd

const double norm_sd2 = 5;

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


//---------------------- //

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

// ------------- define TaskEnvironment Class ------------- //

class TaskEnvironment {

    vector<double> taskgroup_norm_mean1_list; // risky option (option 1)
    vector<double> taskgroup_norm_sd1_list;
    vector<double> taskgroup_norm_mean2_list; // non-risky option (option 0)
    vector<double> taskgroup_norm_sd2_list;
    vector<double> effectsize_list;

    public:

        TaskEnvironment();
        
        // functions to generate random tasks (task group)
        void generate_random_tasks();
        double calculate_effectsize(double m1, double sd1, double m2, double sd2);
        double return_m1(int task_i) {return taskgroup_norm_mean1_list[task_i]; }
        double return_sd1(int task_i) {return taskgroup_norm_sd1_list[task_i]; }
        double return_m2(int task_i) {return taskgroup_norm_mean2_list[task_i]; }
        double return_sd2(int task_i) {return taskgroup_norm_sd2_list[task_i]; }
        double return_eff_size(int task_i) {return effectsize_list[task_i]; }

};

// constructor
TaskEnvironment::TaskEnvironment() {

    taskgroup_norm_mean1_list.assign(Taskgroup_N, 0.0); // risky option (option 1)
    taskgroup_norm_sd1_list.assign(Taskgroup_N, 0.0);
    taskgroup_norm_mean2_list.assign(Taskgroup_N, 0.0); // non-risky option (option 0)
    taskgroup_norm_sd2_list.assign(Taskgroup_N, 0.0);
    effectsize_list.assign(Taskgroup_N, 0.0);

}

double TaskEnvironment::calculate_effectsize(double m1, double sd1, double m2, double sd2) {

    double target_eff_size = (m1 - m2)/sqrt((sd1 * sd1 + sd2 * sd2)/2);

    return target_eff_size;

}

void TaskEnvironment::generate_random_tasks() {
    
    random_device rd;
    mt19937_64 get_rand_mt_rt(rd());

    for (int task_i = 0; task_i < (Taskgroup_N/2); ++task_i) {
        
        uniform_real_distribution<double> unfm_dist_sd1(norm_sd1_lower, norm_sd1_upper);
        normal_distribution<double> nrml_dist_m1(generate_norm_mean, generate_norm_sigma);
        normal_distribution<double> nrml_dist_m2(generate_norm_mean, generate_norm_sigma);
        
        double sampled_sd1;
        double sampled_m1;
        double sampled_m2;
        double target_eff_size;

        do {

            sampled_sd1 = unfm_dist_sd1(get_rand_mt_rt);
            sampled_m1 = nrml_dist_m1(get_rand_mt_rt);
            sampled_m2 = nrml_dist_m2(get_rand_mt_rt);

            target_eff_size = calculate_effectsize(sampled_m1, sampled_sd1, sampled_m2, norm_sd2);

        } while((sampled_m1 < sampled_m2)||(target_eff_size > 5)); // m1 > m2, effect_size < 5となるまでサンプリング

        cout << "set task_i: " << task_i << ", N(" << sampled_m1 << ", " << sampled_sd1 << "), N(" << sampled_m2 << ", " << norm_sd2 << "), " << "eff_size: " << target_eff_size << endl;
        
        taskgroup_norm_mean1_list[task_i] = sampled_m1;
        taskgroup_norm_sd1_list[task_i] = sampled_sd1;
        taskgroup_norm_mean2_list[task_i] = sampled_m2;
        taskgroup_norm_sd2_list[task_i] = norm_sd2;
        effectsize_list[task_i] = target_eff_size;
        
    }

    // the rest of half
    for (int task_i = 0; task_i < (Taskgroup_N/2); ++ task_i) {
        
        double target_eff_size;

        taskgroup_norm_mean1_list[task_i + (Taskgroup_N/2)] = (-1) * taskgroup_norm_mean1_list[task_i];
        taskgroup_norm_sd1_list[task_i + (Taskgroup_N/2)] = taskgroup_norm_sd1_list[task_i];
        taskgroup_norm_mean2_list[task_i + (Taskgroup_N/2)] = (-1) * taskgroup_norm_mean2_list[task_i];
        taskgroup_norm_sd2_list[task_i + (Taskgroup_N/2)] = norm_sd2;

        target_eff_size = calculate_effectsize(taskgroup_norm_mean1_list[task_i + (Taskgroup_N/2)], taskgroup_norm_sd1_list[task_i + (Taskgroup_N/2)], taskgroup_norm_mean2_list[task_i + (Taskgroup_N/2)], taskgroup_norm_sd2_list[task_i + (Taskgroup_N/2)]);

        effectsize_list[task_i + (Taskgroup_N/2)] = target_eff_size;

    }
    
}

// ------------- define saveData Class ------------- //

class saveData {

    char time_stamp_str[100];
    time_t start_time;
    string path_name;
    string folder_name;
    string file_name_pathinfo;
    string file_name_taskgroup;

    string system_order_st_cp; // for file copy order of system()
    
    public:

        saveData(); // constructor

        // functions to create directory and csv files to save
        void create_dir();
        void create_file();
        
        // functions to write output to files
        void writePathInfo();
        void writeTaskGroup(TaskEnvironment taskenv);

};

// constructor of saveData
saveData::saveData()
{
    // initialize start time stamp
    time_t now = time(NULL); struct tm tm; 
    strftime(time_stamp_str, sizeof(time_stamp_str), "%Y%m%d%H%M%S", localtime_r(&now, &tm));
    cout << "\n--- generating tasks time is ... " << time_stamp_str << endl;
    
    //start_time = now;

}

// function to create directory to make files
void saveData::create_dir() {

    path_name = "./Analysis/Results";
    stringstream ss;
    
    // info about generation, popSize, trial
    ss << "Sim" << setw(4) << setfill('0') << SimN;
    folder_name += ss.str(); ss.str(""); ss.clear(stringstream::goodbit);

    ss << "_G" << setw(5) << setfill('0') << Generation;
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
    ss << "_RiskN" << setw(2) << setfill('0') << Riskseek_N;
    folder_name += ss.str(); ss.str(""); ss.clear(stringstream::goodbit);
    ss << "_Taskg" << setw(4) << setfill('0') << Taskgroup_N;
    folder_name += ss.str(); ss.str(""); ss.clear(stringstream::goodbit);


    /* ------------------------------------------------------------- */
    // create path to save files
    path_name += "/"; path_name += folder_name; path_name += '_'; path_name += time_stamp_str; path_name += "/"; 

    // create the directory
    if (mkdir("./Analysis", 0777) == 0) 
        cout << "---- create the directory ./Analysis/" << endl; 
    if (mkdir("./Analysis/Results", 0777) == 0) 
        cout << "---- create the directory ./Analysis/Results." << endl;
    if (mkdir(path_name.c_str(), 0777) == 0) 
        cout << "---- create the directory ..." << endl << path_name.c_str() << endl; 
    
    // Copy this .cpp file to the directory
    system_order_st_cp = "cp generate_RandMultiTask.cpp ";
    system_order_st_cp += path_name;
    
    system(system_order_st_cp.c_str()); // run copy

}

// function to create csv files to save data
void saveData::create_file() {

    // ------------ set file pathinfo ---------- //

    file_name_pathinfo += path_name; file_name_pathinfo += "_pathinfo_"; file_name_pathinfo += time_stamp_str; file_name_pathinfo += ".txt";
 
    ofstream ofs_pathinfo(file_name_pathinfo);
    if (!ofs_pathinfo) {
        cout << "cannot open file!. Check file_name: " << endl << file_name_pathinfo << endl;
    }
    else {
        cout << "succeed! created " << file_name_pathinfo << endl;
    }

    // ------------ set file taskgroup ---------- //

    file_name_taskgroup += path_name; file_name_taskgroup += "_taskgroup.csv";

    ofstream ofs_taskgroup(file_name_taskgroup);
    if (!ofs_taskgroup) {
        cout << "cannot open file!. Check file_name: " << endl << file_name_taskgroup << endl;
    }
    else {
        cout << "succeed! created " << file_name_taskgroup << endl;
        // 6 columns, write the header names
        ofs_taskgroup << "task_i,m1,sd1,m2,sd2,eff_size" << endl;
    }

}

// This function writes task informations to file_name_params (_Parameters_memo.txt)
void saveData::writePathInfo() {

    ofstream fout_pathinfo(file_name_pathinfo, ios::app);

    fout_pathinfo << path_name << endl;
    fout_pathinfo << "start time: " << time_stamp_str << endl;

    fout_pathinfo << "Script Name : " << script_name << endl << \
    "Task Group_N: " << Taskgroup_N << endl << \
    "generate Norm mean : " << generate_norm_mean << endl << "generate Norm sigma : " << generate_norm_sigma << endl << \
    "Norm sd upper : " << norm_sd1_upper << endl << "Norm sd lower : " << norm_sd1_lower << endl;

}

// function to record data to output_trial
void saveData::writeTaskGroup(TaskEnvironment taskenv) {

    ofstream fout_taskgroup(file_name_taskgroup, ios::app);

    for (int task_i = 0; task_i < Taskgroup_N; ++task_i) {

        fout_taskgroup << task_i << ", " << taskenv.return_m1(task_i) << ", "  << taskenv.return_sd1(task_i) << ", " << \
            taskenv.return_m2(task_i) << ", " << taskenv.return_sd2(task_i) << ", " << taskenv.return_eff_size(task_i) << endl;

    }

}


// ----------------- MAIN PROCESS ------------------- //

int main(int argc, char *argv[]) {

    //------ get parameters set from command -----//

    script_name += argv[0];
    
    PopSize = atoi(argv[1]);
    Generation = atoi(argv[2]);
    Trial = atoi(argv[3]);
    Task_N = atoi(argv[4]);
    Riskseek_N = atoi(argv[5]);
    Taskgroup_N = atoi(argv[6]);
    
    generate_norm_mean = atof(argv[7]);
    generate_norm_sigma = atof(argv[8]);
    norm_sd1_upper = atof(argv[9]);
    norm_sd1_lower = atof(argv[10]);

    /*
    Taskgroup_N = 100;
    int Riskseek_N;
    PopSize = 10000;
    Generation = 3000;
    Trial = 500;

    generate_norm_mean = 0;
    generate_norm_sigma = 20;
    norm_sd1_upper = 30;
    norm_sd1_lower = 5;
    */

    saveData savedat;
    TaskEnvironment taskenv;

    savedat.create_dir();
    savedat.create_file();
    taskenv.generate_random_tasks();
    
    savedat.writePathInfo();
    savedat.writeTaskGroup(taskenv);
    
    return 0;

}
