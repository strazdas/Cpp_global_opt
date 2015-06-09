#include <iostream>
#include "Disimplv.h"

using namespace std;


int main() {
    PointTree* points = new PointTree();
    points->add(1.);
    cout << "1:  "; points->print();
    points->add(3.);
    cout << "3:  "; points->print();
    points->add(7.);
    cout << "7:  "; points->print();
    points->add(2.);
    cout << "2:  "; points->print();
    points->add(15.);
    cout << "15: "; points->print();
    points->add(16.);
    cout << "16: "; points->print();
    points->add(11.);
    cout << "11: "; points->print();
    points->add(18.);
    cout << "18: "; points->print();
    points->add(13.);
    cout << "13: "; points->print();
    points->add(14.);
    cout << "14: "; points->print();
    points->add(12.);
    cout << "12: "; points->print();
    points->add(17.);
    cout << "17: "; points->print();
    points->add(5.);
    cout << "5:  "; points->print();
    points->add(4.);
    cout << "4:  "; points->print();
    points->add(6.);
    cout << "6:  "; points->print();
    points->add(8.);
    cout << "8:  "; points->print();
    points->add(9.);
    cout << "9:  "; points->print();
    points->add(19.);
    points->add(21.);
    points->add(22.);
    points->add(23.);
    points->add(24.);
    points->add(25.);
    points->add(26.);
    points->add(27.);
    points->add(28.);
    cout << "19: "; points->print();
    points->print();
    return 0;
};
