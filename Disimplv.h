#ifndef DISIMPLV_H
#define DISIMPLV_H_ 
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <math.h> 
#include <limits>

using namespace std;




class Point {
    Point(const Point& other){}
    Point& operator=(const Point& other){};
public:
    Point(){};
    Point(int *c, int argc){
        for (int i=0; i<argc ; i++){
            _X.push_back(double(c[i]));
        };
    };
    Point(double *c, int argc){
        for (int i=0; i<argc ; i++){
            _X.push_back(c[i]);
        };
    };
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

    int size(){
        return _X.size();
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
   
    void add(Point* point){
        _points.push_back(point);
    };

    Point* get(double *c, int argc){
        for (int i=0; i < _points.size(); i++){
            bool matches = true;
            for (int j = 0; j < argc; j++){
                if (_points[i]->_X[j] != c[j]) {
                    matches = false;
                    break;
                };
            };
            if (matches) {
                return _points[i];
            };
        };
        return 0;
    };

    Point* get(Point* p){
        for (int i=0; i < _points.size(); i++){
            bool matches = true;
            for (int j = 0; j < p->size(); j++){
                if (_points[i]->_X[j] != p->_X[j]) {
                    matches = false;
                    break;
                };
            };
            if (matches) {
                return _points[i];
            };
        };
        return 0;
    };

    // Point* get_point(double c1, double c2){
    //     for (int i=0; i < _points.size(); i++){
    //         if (_points[i]->_X[0] == c1 && _points[i]->_X[1] == c2){
    //             return _points[i];
    //         };
    //     };
    //     return 0;
    // };

    int size() {
        return _points.size();
    };

    void print(){
        cout << "Points: " << endl;
        for (int i=0; i < _points.size(); i++) {
            _points[i]->print();
        };
    };

    virtual ~Points(){
        for (int i=0; i < _points.size(); i++) {
            delete _points[i];
        };
        _points.clear();
    };
};


class Simplex {
    Simplex(const Simplex& other){}
    Simplex& operator=(const Simplex& other){}
public:
    Simplex(){
        _is_in_partition = true;
    };
    vector<Point*> _verts;
    bool _is_in_partition;
//     double _tolerance;
//     double _hash;
//     double _parent_hash;
//     Points* _min_AB;

    void add_vertex(Point* vertex){
        _verts.push_back(vertex);
    };

    int size() {
        return _verts.size();
    };

    void print(){
        cout << "Simplex: " << endl;
        for (int i=0; i < _verts.size(); i++){
            _verts[i]->print();
        }
    };
    static void print(vector<Simplex*> simplexes, string label="Printing simplexes:"){
        cout << label << endl;
        for (int i=0; i < simplexes.size(); i++){
            simplexes[i]->print();
        };
    };

    virtual ~Simplex(){
        _verts.clear();
    };  
};
// typedef vector<Simplex*> Simplexes;


class Function {
    Function(const Function& other){};
    Function& operator=(const Function& other){};
public:
    Function(){
        _calls = 0;
        _f_min = numeric_limits<int>::max();
    };
    string _name;
    int _D;
    Point* _lb;
    Point* _ub;
    Point* _glob_x;  // Point where global function minimum is (should be list)
    double _glob_f;  // Predefined global function minimum
    double _L;       // Lipschitz constant

    int _calls;
    double _f_min;
    Point* _x_min;
    Points _points;


    Point* get(double *c, int argc){
        Point* cached_point = _points.get(c, argc);
        if (cached_point) {
            return cached_point;
        } else {
            Point* p = new Point(c, argc);
            double val = value(p);
            if (_f_min > val) {
                _f_min = val;
                _x_min = p;
            };
            p->add_value(val);
            _calls += 1;
            _points.add(p);
            return p;
        };
    };

    Point* get(Point* p){
        Point* cached_point = _points.get(p);
        if (cached_point) {
            return cached_point;
        } else {
            double val = value(p);
            if (_f_min > val) {
                _f_min = val;
                _x_min = p;
            };
            p->add_value(val);
            _calls += 1;
            _points.add(p);
            return p;
        };
    };

    void print(){
        cout << "\nCalls: " << _calls << "   f_min: " << _f_min << endl;
    };

    double transform(Point* point, int i) {     // Transforms single point coordinate from [0,1] to [l,u]
        return point->_X[i] * (_ub->_X[i]-_lb->_X[i]) + _lb->_X[i];  
    };

    double pe(){
        if (_glob_f != 0) {
            return (_f_min - _glob_f) / fabs(_glob_f) * 100.; 
        };
        return _f_min * 100.; 
    };

    virtual double value(Point* point) = 0;

    virtual ~Function(){};
};

class Branin : public Function {
    Branin(const Branin& other){};
    Branin& operator=(const Branin& other){};
public:
    Branin(): Function(){
        _name = "Branin";
        _D = 2;
        _lb = new Point(-5, 0);
        _ub = new Point(10, 15);
        _glob_x = new Point(M_PI, 12.275);  // Point where global function minimum is (should be list)
        _glob_f = 0.397887;          // Predefined global function minimum
        _L = 109.94813585;
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
    Disimplv(const Disimplv& other){};
    Disimplv& operator=(const Disimplv& other){};
public:
    // Disimplv(){};
    Disimplv(double min_pe, int max_calls){
        _min_pe = min_pe;
        _max_calls = max_calls;
    };
    double _min_pe;
    int _max_calls;
    vector<Simplex*> _partition;
    vector<Simplex*> _all_simplexes;
    Function* _func;
//     Disimplv(Function f, Point lb, Point ub, double error, )
//     Disimplv(int N, int max_calls, , Function f, double L, Points& D){
//         _N = N; _max_calls = max_calls; _error = error; _f = f; _L = L; _D = D;
//     }
//     double _L;      // Lipschitz constant
//     Points* _D;
// 
    void partition_feasable_region(){
        int n = _func->_D;
        int number_of_simpleces = 1;
        for (int i = 1; i <= n; i++) {
            number_of_simpleces *= i;
        };

        int teta[n];
        for (int i=0; i < n; i++){
            teta[i] = i;
        };

        do {
            int triangle[n+1][n];

            for (int k = 0; k < n; k++) {
                triangle[0][k] = 0;
            };

            for (int vertex=0; vertex < n; vertex++) {
                for (int j = 0; j < n + 1; j++) {
                    triangle[vertex + 1][j] = triangle[vertex][j];
                };
                triangle[vertex + 1][teta[vertex]] = 1;
            }

            Simplex* simpl = new Simplex();
            for (int i=0; i < n + 1; i++){
                Point* point = _func->get(new Point(triangle[i], n)); 
                simpl->add_vertex(point);
            };
            _partition.push_back(simpl);
            _all_simplexes.push_back(simpl);

        } while (next_permutation(teta, teta+n));
    };

    void print(){
        cout << endl << "Disimpl-v" << endl;
        cout << " Simplexes in partition: " << _partition.size() << endl;
        for (int i=0; i < _partition.size(); i++){
            _partition[i]->print();
        };
    };

    double l2norm(Point* p1, Point* p2) {
        double squared_sum = 0;
        for (int i=0; i < p1->size(); i++){
            squared_sum += pow(p1->_X[i] - p2->_X[i], 2);
        };
        return sqrt(squared_sum);
    };

    vector<Simplex*> select_simplexes_to_divide(string strategy="all"){
        vector<Simplex*> selected_simplexes;
        // if (strategy == "min_vert") {
        //     selected_simplexes;
        // };
        if (strategy == "all") {
            selected_simplexes = _partition;
        };
        for (int i=0; i < selected_simplexes.size(); i++){
            selected_simplexes[i]->_is_in_partition = false;
        };
        return selected_simplexes;
    };

    vector<Simplex*> divide_simplex(Simplex* simplex, string strategy="longest_half") {
        vector<Simplex*> divided_simplexes;
        if (strategy== "longest_half") {
            double max_edge_length = 0;
            int v1_idx = 0; 
            int v2_idx = 0; 

            for (int a=0; a < simplex->size(); a++) {
                for (int b=0; b < simplex->size(); b++){
                    if (b > a) {
                        double edge_length = l2norm(simplex->_verts[a], simplex->_verts[b]); 
                        if (edge_length > max_edge_length) {
                            max_edge_length = edge_length;
                            v1_idx = a;
                            v2_idx = b;
                        }; }; }; }; 

            // Find middle point
            int n = _func->_D;
            double c[n];
            for (int i=0; i < n; i++) {
                c[i] = (simplex->_verts[v1_idx]->_X[i] + simplex->_verts[v2_idx]->_X[i]) / 2.;
            };
            Point* middle_point = _func->get(c, n);

            // Construct two new simplexes using this middle point.
            Simplex* left_simplex = new Simplex();
            Simplex* right_simplex = new Simplex();
            for (int i=0; i < simplex->size(); i++){
                // Point* point = _func->get(new Point(triangle[i], n)); 
                if (i != v1_idx){
                    right_simplex->add_vertex(simplex->_verts[i]);
                } else {
                    right_simplex->add_vertex(middle_point);
                };
                if (i != v2_idx) {
                    left_simplex->add_vertex(simplex->_verts[i]);
                } else {
                    left_simplex->add_vertex(middle_point);
                };
            };

            divided_simplexes.push_back(left_simplex);
            divided_simplexes.push_back(right_simplex);
            return divided_simplexes;
        };
    };

    static bool should_remove_from_partition(Simplex* s) {
        return !s->_is_in_partition;
    };

    void minimize(Function* func){
        _func = func;
        partition_feasable_region();

        int iterations = 0;
        while (_func->_calls <= _max_calls && _func->pe() > _min_pe){
            // Selects simplexes to divide
            vector<Simplex*> simplexes_to_divide = select_simplexes_to_divide();

            // Divides selected simplexes
            vector<Simplex*> new_simplexes;
            for (int i=0; i < simplexes_to_divide.size(); i++){
                vector<Simplex*> divided_simplexes = divide_simplex(simplexes_to_divide[i]);
                for (int j=0; j < divided_simplexes.size(); j++) {
                    new_simplexes.push_back(divided_simplexes[j]);
                };
            };

            // Remove partitioned simplexes from _partition
            _partition.erase(remove_if(_partition.begin(), _partition.end(), should_remove_from_partition), _partition.end());

            // Add new simplexes to _partition and _all_simplexes
            for (int i=0; i < new_simplexes.size(); i++) {
                _partition.push_back(new_simplexes[i]);
                _all_simplexes.push_back(new_simplexes[i]);
            };

            iterations += 1;
            // Simplex::print(_partition, "Partition_:");
            // Simplex::print(_all_simplexes, "\nAll simplexes:");
            break;
        };

        // Draw partitioning: output simplex coordinates to file and draw it with Python
    };


    virtual ~Disimplv(){
        for (int i=0; i < _all_simplexes.size(); i++) {
            delete _all_simplexes[i];
        };
        _all_simplexes.clear();
        _partition.clear();
    };
};

#endif
