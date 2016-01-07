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
    char* callback = {'\0'};
    int max_duration = -1;
    int max_calls = -1;

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

    if ((max_duration >= 0) && (max_calls >= 0)) {
        alg = new Asimpl(max_calls=max_calls, max_duration=max_duration);
    } else if (max_duration >= 0) {
        alg = new Asimpl(max_duration=max_duration);
    } else if (max_calls >= 0) {
        alg = new Asimpl(max_calls=max_calls);
    } else {
        alg = new Asimpl();
    };
    func = new GKLSFunction(cls, fid);

    alg->minimize(func);

    // Save results
    cout << "Calls " << func->_calls << endl;
    if (callback != '\0') {
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
    };

    // Free memory
    delete alg;
    delete func;
    return 0;
};
