#include <iostream>
#include "Disimplv.h"

using namespace std;


int main() {
    // Class1, id80
    GKLSFunction* func = new GKLSFunction(1, 80);
    // Branin* func = new Branin();
    Disimplv* alg = new Disimplv(0.01, 100000);
    alg->minimize(func);
    // alg->print();
    func->print(); 

    delete func;
    delete alg;
    return 0;
}
