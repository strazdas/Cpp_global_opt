#include <iostream>
#include "Disimplv.h"

using namespace std;


int main() {
    GKLSFunction* func;
    Disimplv* alg;
    for (int cls=1; cls <= 1; cls++) {
        int calls[100];
        for (int fid=1; fid <= 100; fid++) {
            fid = 3;
            alg = new Disimplv(1.0, 30000);
            func = new GKLSFunction(cls, fid);    // Class1, id80
            alg->minimize(func);
            cout << fid << ". ";
            func->print();
            calls[fid-1] = func->_calls;
            delete alg;
            delete func;
            break;
        };
        int calls_sum = 0;
        for (int i=0; i< 100; i++) {
            calls_sum += calls[i];
        };
        // cout << "Class " << cls << " average: " << calls_sum/100. << endl; 
        break;
    };
    GKLS_free();
    GKLS_domain_free();
    return 0;
};
