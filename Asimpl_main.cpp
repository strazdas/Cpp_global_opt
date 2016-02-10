#include <iostream>
#include <getopt.h>
// getopt example: http://stackoverflow.com/questions/8793020/using-getopt-long-c-how-do-i-code-up-a-long-short-option-to-both-require-a
#include "Asimpl.h"
#include <stdio.h>
#include <fstream>

#define no_argument 0
#define required_argument 1
#define optional_argument 2


using namespace std;

void print_stats(int calls[], int subregions[], int n) {
    sort(calls, calls+n);
    sort(subregions, subregions+n);
    int calls_sum = 0;
    for (int i=0; i < n; i++) {
        calls_sum += calls[i];
    };
    cout << "Calls50: " << calls[49] << " Calls100: " << calls[99] << " Average: " << calls_sum/100. <<
        " Subregions50: " << subregions[49] << " Subregions100: " << subregions[99] << endl; 
};


int main(int argc, char* argv[]) {
    // Parse parameters
    const struct option longopts[] = {
        {"func_cls", required_argument, 0, 'c'},
        {"func_id", required_argument, 0, 'f'},
        {"task_id", required_argument, 0, 't'},
        {"callback", required_argument, 0, 'b'},
        {"max_duration", optional_argument, 0, 'd'},
        {"max_calls", optional_argument, 0, 'i'},
    };
    int cls;
    int fid;
    int task_id;
    char* callback;
    int max_calls = 40000;
    int max_duration = 3600;

    int opt_id;
    int iarg = 0;
    while(iarg != -1) {
        iarg = getopt_long(argc, argv, "cftbdi", longopts, &opt_id);
        switch (iarg) {
            case 'c':
                cls = strtoul(optarg, 0, 0);
                break;
            case 'f':
                fid = strtoul(optarg, 0, 0);
                break;
            case 't':
                task_id = strtoul(optarg, 0, 0);
                break;
            case 'b':
                callback = strdup(optarg);
                break;
            case 'd':
                max_duration = strtoul(optarg, 0, 0);
                break;
            case 'i':
                max_calls = strtoul(optarg, 0, 0);
                break;
        };
    };

    // Minimize 
    GKLSFunction* func;
    Asimpl* alg;

    alg = new Asimpl(max_calls, max_duration);
    func = new GKLSFunction(cls, fid);

    alg->minimize(func);

    // Save results
    string cmd;
    stringstream cmd_ss; 
    cmd_ss << callback
           << " --calls=" << func->_calls
           << " --subregions=" << alg->_partition.size()
           << " --duration=" << alg->_duration  
           << " --task_id=" << task_id  
           << " --status=" << alg->_status  
           << " -exe=" << argv[0];  
    cmd = cmd_ss.str();
    popen(cmd.c_str(), "r");

    // Free memory
    delete alg;
    delete func;
    return 0;
};
