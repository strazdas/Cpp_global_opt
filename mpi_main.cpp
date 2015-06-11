#include <iostream>
#include <sstream>
#include <fstream>
#include "mpi.h"
#include "Disimplv.h"
#include <sys/stat.h>

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

int main (int argc, char *argv[]) {
    int cls = atoi(argv[1]);
    int fid_from = 1;
    int fid_till = 100;
    int pid;   // procesor id
    int pool;  // total processors in this pool
    ofstream result_file; 

    if (cls < 5) {
        fid_till = 100;
    } else {
        MPI::Init();
        pool = MPI::COMM_WORLD.Get_size();
        pid = MPI::COMM_WORLD.Get_rank();
        fid_from = 1 + pid * 10;
        fid_till = (pid + 1) * 10; 
    };
        
    GKLSFunction* func;
    Disimplv* alg;
    for (int fid=fid_from; fid <= fid_till; fid++) {
        // Initialise function and algorithm
        alg = new Disimplv(1.0, 1000000);
        func = new GKLSFunction(cls, fid);

        // Minimize function
        alg->minimize(func);

        // Save results to file
        stringstream filename; 
        mkdir(("/home/albertas/C++_global_opt/results/" + alg->_name).c_str(), 0777);  // Tries to create directory 
        filename << "/home/albertas/C++_global_opt/results/" << alg->_name << "/" << cls << "_" << fid;
        result_file.open(filename.str().c_str(), ios::app);
        result_file << (*alg) << endl;
        result_file.close();

        delete alg;
        delete func;
    };

    if (cls >= 5) { MPI::Finalize(); };
    return 0;
};
