#include <iostream>
#include "Disimplv.h"
#include <math.h> 

using namespace std;


int main() {
    cout << "Start!" << endl;
    // // Initialize algorithm
    // Disimplv alg();
    // alg._max_calls = 10000;
    // alg._error = 0.01;
    // alg._mirror_division = false;

    // Initialize function
    // Function* func = new Function();
    // Point* point = func->get_point(0, 1);
    // point->add_value(1.2);
    // func->_points.print();
    // point->add_value(1.5);
    // Point* point2 = func->get_point(0, 1);
    // point->print();
    // point2->print();
    // Point* test = func->_ub;
    // delete func;

    cout << "Branin function test" << endl;
    Branin* func2 = new Branin();
    // Point* p = func2->get(-M_PI, 12.275);
    // Point* p = func2->get(M_PI, 2.275);
    // Point* p = func2->get(9.42478, 2.475);
    // Point* p = func2->get(-5, 0);
    Point* p = func2->get(0., 0.);
    p->print();

    // func._name = "branin";
    // func._D = 2;
    // func._lb = Point(-5, 0);
    // func._ub = Point(10, 15);
    // func._x_min = Point(-pi, 12.275);
    // func._f_min = 0.397887;

    // Triangulate feasable region 
    return 0;
}
