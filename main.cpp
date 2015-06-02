#include <iostream>
#include "Disimplv.h"
#include <math.h> 

using namespace std;


int main() {
    Branin* func = new Branin();
    Disimplv* alg = new Disimplv(1.0, 100000);
    alg->minimize(func);
    alg->print();
    func->print(); 

    delete func;
    delete alg;
    return 0;
}
