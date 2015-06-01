#ifndef DISIMPLV_H
#define DISIMPLV_H_ 
#include <vector>
#include <string>
#include <iostream>
#include <math.h> 

using namespace std;


class Point {
    Point(const Point& other){}
    Point& operator=(const Point& other){};
public:
    Point(){};
    Point(double c1){};
    Point(double c1, double c2){
        _X.push_back(c1);
        _X.push_back(c2);
    };

    vector<double> _X;  // Coordinates in normalised [0,1]^n space  
    vector<double> _values;
         
    void add_value(double value) {
        _values.push_back(value);
    };

    void print(){
        for (int i=0; i < _X.size(); i++){
            cout << _X[i] << "  ";
        };
        for (int i=0; i < _values.size(); i++){
            if (i == 0) { cout << "->  "; };
            cout << _values[i] << "  ";
        };
        cout << endl;
    };

    virtual ~Point(){
        vector<double>::iterator cit = _X.begin();
        while (cit != _X.end()) {
            cit = _X.erase(cit);    
        };
        vector<double>::iterator vit = _values.begin();
        while (vit != _values.end()) {
            vit = _values.erase(vit);    
        };
    };
};


class Points {  // Binary balancing tree or simply linked list
    Points(const Points& other){}
    Points& operator=(const Points& other){}
public:                
    Points(){};
    vector<Point*> _points;
   
    void add_point(Point* point){
        _points.push_back(point);
    };

    Point* get_point(double c1, double c2){
        for (int i=0; i < _points.size(); i++){
            if (_points[i]->_X[0] == c1 && _points[i]->_X[1] == c2){
                return _points[i];
            };
        };
        return 0;
    };

    void print(){
        for (int i=0; i < _points.size(); i++) {
            _points[i]->print();
        };
    };

    virtual ~Points(){
        for (int i=0; i < _points.size(); i++) {
            delete _points[i];
        }
        _points.clear();
    };
};


//// typedef vector<Point*> Points;
// 
// class Simplex {
//     Simplex(const Simplex& other){}
//     Simplex& operator=(const Simplex& other){}
// public:
//     Simplex(){}
//     Points* _vertexes;
//     double _tolerance;
//     double _hash;
//     double _parent_hash;
//     Points* _min_AB;
// 
//     string toString(); 
//     virtual ~Simplex();
// }
// 
class Function {
    Function(const Function& other){};
    Function& operator=(const Function& other){};
public:
    Function(){
        _calls = 0;
    };
    string _name;
    int _D;
    Point* _lb;
    Point* _ub;
    Point* _glob_x;  // Point where global function minimum is (should be list)
    double _glob_f;  // Predefined global function minimum

    int _calls;
    double _f_min;
    Point* _x_min;
    Points _points;

    Point* get(double c1, double c2){
        Point* cached_point = _points.get_point(c1, c2);
        if (cached_point) {
            return cached_point;
        } else {
            Point* new_point = new Point(c1, c2);
            new_point->add_value(value(new_point));
            _calls += 1;
            _points.add_point(new_point);
            return new_point;
        };
    };

    double transform(Point* point, int i) {  // Transforms single point coordinate from [0,1] to [l,u]
        cout << "Transformation ub, lb: " << _lb->_X[i]<< ", "<< _ub->_X[i] << endl;
        cout << "Transformation value: " << point->_X[i] * (_ub->_X[i]-_lb->_X[i]) + _lb->_X[i] << endl;
        return point->_X[i] * (_ub->_X[i]-_lb->_X[i]) + _lb->_X[i];  
    };

    virtual double value(Point* point) = 0;

    virtual ~Function(){};
};

class Branin : public Function {
public:
    Branin(){
        _calls = 0;
        _name = "Branin";
        _D = 2;
        _lb = new Point(-5, 0);
        _ub = new Point(10, 15);
        _glob_x = new Point(M_PI, 12.275);  // Point where global function minimum is (should be list)
        double _glob_f = 0.397887;  // Predefined global function minimum
    };

    double value(Point* point) {
        double x1 = transform(point, 0); 
        double x2 = transform(point, 1); 
        double part1 = pow((x2 - 5/(4*pow(M_PI, 2))*pow(x1,2) + 5/M_PI*x1 -6), 2);
        double part2 = 10*(1 - 1/(8*M_PI))*cos(x1) + 10;
        return part1 + part2;
    };
};



class Disimplv {
    Disimplv(const Disimplv& other){}
    Disimplv& operator=(const Disimplv& other){}
public:
    Disimplv(){}
//     Disimplv(Function f, Point lb, Point ub, double error, )
//     Disimplv(int N, int max_calls, , Function f, double L, Points& D){
//         _N = N; _max_calls = max_calls; _error = error; _f = f; _L = L; _D = D;
//     }
//     int _N;         // Variable space dimension  
//     int _max_calls; // Max function calls
//     float _error;     // Error
//     Function _f;    // Function
//     double _L;      // Lipschitz constant
//     Points* _D;
// 
//     void toString();
//     void optimize();
    virtual ~Disimplv(); // Destruktorius
};

#endif
