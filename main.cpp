#include <iostream>
#include "Disimplv.h"

using namespace std;


int main() {
    Branin* func = new Branin();
    Disimplv* alg = new Disimplv(0.01, 100000);
    alg->minimize(func);
    // alg->print();
    func->print(); 

    delete func;
    delete alg;
    return 0;
}
