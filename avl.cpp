/* Testing PointTree class implementation */
#include <iostream>
#include "Disimplv.h"

using namespace std;


int main() {
    PointTree* points = new PointTree();
    double x0[3] = {1., 21., 31.};   Point* p0 = new Point(x0, 3);  points->add(p0);
    double x1[3] = {2., 21., 31.};   Point* p1 = new Point(x1, 3);  points->add(p1);
    double x2[3] = {3., 21., 31.};   Point* p2 = new Point(x2, 3);  points->add(p2);
    double x3[3] = {3., 22., 31.};   Point* p3 = new Point(x3, 3);  points->add(p3);
    double x4[3] = {3., 23., 31.};   Point* p4 = new Point(x4, 3);  points->add(p4);
    double x5[3] = {4., 23., 31.};   Point* p5 = new Point(x5, 3);  Point* rp5 = points->add(p5);
    double x6[3] = {3., 23., 31.};   Point* p6 = new Point(x6, 3);  Point* rp6 = points->add(p6); p6->add_value(-9.);
    cout << "p3" << p3 << "   p5 " << p5 << "   rp5 " << rp5 << endl;
    cout << "p3" << p3 << "   p6 " << p6 << "   rp6 " << rp6 << endl;
    points->print();
    delete points;
    return 0;
};
