#ifndef FUNCTIONS_H
#define FUNCTIONS_H 
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h> 
#include <limits>

#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "gkls.h"
#include "rnd_gen.h"

using namespace std;

class Simplex;

class Point {
    Point(const Point& other){}
    Point& operator=(const Point& other){};
public:
    Point(int D){
        _D = D;
        _X = (double*) malloc((D)*sizeof(double));
    };
    Point(int *c, int argc){
        _D = argc;
        _X = (double*) malloc((argc)*sizeof(double));
        for (int i=0; i<argc ; i++){
            _X[i] = double(c[i]);
        };
    };
    Point(double *c, int argc){
        _D = argc;
        _X = (double*) malloc((argc)*sizeof(double));
        for (int i=0; i<argc ; i++){
            _X[i] = c[i];
        };
    };
    Point(double c, int argc){
        _D = argc;
        _X = (double*) malloc((argc)*sizeof(double));
        for (int i=0; i<argc ; i++){
            _X[i] = c;
        };
    };
    Point(double c1, double c2){
        _D = 2;
        _X = (double*) malloc(2*sizeof(double));
        _X[0] = c1;
        _X[1] = c2;
    };

    int _D;
    double* _X;  // Coordinates in normalised [0,1]^n space  
    vector<double> _values;
    vector<Simplex*> _simplexes;  // Simplexes, which have this point as vertex
         
    void add_value(double value) {
        _values.push_back(value);
    };

    Point* copy() {
        Point* point_copy = new Point(_X, _D);
        for (int i=0; i < _values.size(); i++) {
            point_copy->add_value(_values[i]);
        };
        return point_copy;
    };

    int size(){
        return _D;
    };

    void _neighbours_estimates_should_be_updated();
    // {
    //     for (int sid=0; sid < _simplexes.size(); sid++) {
    //         _simplexes[sid]->_should_estimates_be_updated = true;
    //     };
    // };

    void print(){
        // cout.precision(17);
        cout << "       ";
        for (int i=0; i < size(); i++){
            // cout << fixed << _X[i] << "  \t";
            cout << _X[i] << "  \t";
        };
        for (int i=0; i < _values.size(); i++){
            if (i == 0) { cout << "->\t"; };
            // cout << fixed << _values[i] << "  ";
            cout << _values[i] << "  ";
        };
        cout << endl;
    };

    static bool compare_by_value(Point* p1, Point* p2) {
        return p1->_values[0] < p2->_values[0];
    };

    friend ostream& operator<<(ostream& o, const Point& p){
        for (int i=0; i < p._D; i++) {
            o << p._X[i];
            if (i != p._D - 1) { o << ","; };
        };
        return o;
    };

    virtual ~Point(){
        free(_X);
        vector<double>::iterator vit = _values.begin();
        while (vit != _values.end()) {
            vit = _values.erase(vit);    
        };
    };
};


class PointTree;

class PointTreeNode {  // Binary balancing tree or simply linked list
    PointTreeNode(const PointTreeNode& other){}
    PointTreeNode& operator=(const PointTreeNode& other){}
public:                
    PointTreeNode(double value=numeric_limits<double>::max()){
        _height = 1;
        _value = value;
        _parent = 0;
        _left = 0;
        _right = 0;
        _subtree = 0;
        _point = 0;
    };
    int _height;
    double _value;
    PointTreeNode* _parent;
    PointTreeNode* _left;
    PointTreeNode* _right;
    PointTree* _subtree;     // Next dimension head
    Point* _point;           // Only last dimension node will have _point != 0;

    void print(){
        if (_left != 0) { cout << "l"; _left->print(); };
        cout << _value << "("<< _height << ")";
        if (_right != 0) { cout << "r"; _right->print(); };
    };

    virtual ~PointTreeNode(); // {
    //     if (_left != 0) { delete _left; };
    //     if (_right != 0) { delete _right; };
    //     if (_point != 0) { delete _point; };
    //     if (_subtree != 0) { delete _subtree; };
    // };
};

class PointTree { // Head of the tree
    PointTree(const PointTree& other){}
    PointTree& operator=(const PointTree& other){}
public:
    PointTree(){
         _tree_root = 0;
         _dim = 1;
    };
    PointTree(int dim){
         _tree_root = 0;
         _dim = dim;
    };
    PointTreeNode* _tree_root;
    int _dim;

    void update_height(PointTreeNode* node) {
        int lh = 0;
        int rh = 0;
        if (node->_left != 0) { lh = node->_left->_height; }; 
        if (node->_right != 0) { rh = node->_right->_height; };
        if (lh > rh) {
            node->_height = lh + 1;
        } else {
            node->_height = rh + 1;
        };
        // Also update all ancestors heights
        if (node->_parent != 0) {
            update_height(node->_parent);
        };
    };

    void left_right_rebalance(PointTreeNode* node) {
        PointTreeNode* diatteched_node;
        // node left right  <-  node left right left 
        diatteched_node = node->_left->_right;
        node->_left->_right = node->_left->_right->_left;
        if (node->_left->_right != 0) { node->_left->_right->_parent = node->_left; };
        // Diatteched left = node->_left
        diatteched_node->_left = node->_left;
        node->_left->_parent = diatteched_node;
        // node left  <-  node left right
        node->_left = diatteched_node;
        diatteched_node->_parent = node;
        // Update heights
        update_height(node);
        update_height(diatteched_node);
        update_height(diatteched_node->_left);
    };
    void left_left_rebalance(PointTreeNode* node) {
        PointTreeNode* diatteched;
        diatteched = node->_left;
        node->_left = node->_left->_right;
        if (node->_left != 0) { node->_left->_parent = node; };  
        diatteched->_parent = node->_parent;
        if (node->_parent != 0) {
            if (node->_parent->_left == node) {
                node->_parent->_left = diatteched;
            } else {
                node->_parent->_right = diatteched;
            };
        } else {
            _tree_root = diatteched;
        };
        diatteched->_right = node;
        node->_parent = diatteched;
        // Update heights
        update_height(node);
        update_height(diatteched);
    };
    void right_left_rebalance(PointTreeNode* node) {
        PointTreeNode* diatteched_node;
        // node left right  <-  node left right left 
        diatteched_node = node->_right->_left;
        node->_right->_left = node->_right->_left->_right;
        if (node->_right->_left != 0) { node->_right->_left->_parent = node->_right; };
        // Diatteched left = node->_left
        diatteched_node->_right = node->_right;
        node->_right->_parent = diatteched_node;
        // node left  <-  node left right
        node->_right = diatteched_node;
        diatteched_node->_parent = node;
        // Update heights
        update_height(node);
        update_height(diatteched_node);
        update_height(diatteched_node->_right);
    };
    void right_right_rebalance(PointTreeNode* node) {
        PointTreeNode* diatteched;
        diatteched = node->_right;
        node->_right = node->_right->_left;
        if (node->_right != 0) { node->_right->_parent = node; };
        diatteched->_parent = node->_parent;                       
        if (node->_parent != 0) {
            if (node->_parent->_left == node) {
                node->_parent->_left = diatteched;
            } else {
                node->_parent->_right = diatteched;
            };
        } else {
            _tree_root = diatteched;
        };
        diatteched->_left = node;
        node->_parent = diatteched;
        // Update heights
        update_height(node);
        update_height(diatteched);
    };

    void check_if_balanced(PointTreeNode* node) {
        int lh = 0;
        int rh = 0;
        int llh = 0;
        int lrh = 0;
        int rlh = 0;
        int rrh = 0;
        if (node->_left != 0) {
            lh = node->_left->_height;
            if (node->_left->_left != 0) { llh = node->_left->_left->_height; };
            if (node->_left->_right != 0) { lrh = node->_left->_right->_height; };
        };
        if (node->_right != 0) {
            rh = node->_right->_height;
            if (node->_right->_left != 0) { rlh = node->_right->_left->_height; };
            if (node->_right->_right != 0) { rrh = node->_right->_right->_height; };
        };
        if (abs(rh - lh) > 1) {
            // Not balanced, so rebalance
            if (rh > lh) {
                if (rrh > rlh) {
                    right_right_rebalance(node);
                } else {
                    right_left_rebalance(node);
                    right_right_rebalance(node);
                };
            };
            if (rh < lh) {
                if (llh > lrh) {
                    left_left_rebalance(node);
                } else {
                    left_right_rebalance(node);
                    left_left_rebalance(node);
                };
            };
        };
        if (node->_parent != 0) {
            check_if_balanced(node->_parent);
        };
    };

    Point* process_next_dimension(PointTreeNode* node, Point* point){
        // Creates next dimension tree if needed and adds point to it
        // Its last dimension
        if (point->_D == _dim) {  // Don't need next dimension
            if (node->_point == 0) {    // Save or return the point
                node->_point = point;
                return 0;
            } else {
                return node->_point;
            };
        };
        // Its not last dimension
        if (node->_subtree == 0) {  // Create subtree if it doesn't already exist
            node->_subtree = new PointTree(_dim + 1);
        };
        Point* found_point = node->_subtree->add(point);  // Get or insert point to the subtree
        if (found_point != 0) {  // We got point so return it 
            return found_point;
        } else {
            return 0;  // We inserted point
        };
    };

    Point* add(Point* point){
        // Get same point or insert given (if inserted returns 0)
        PointTreeNode* node = _tree_root;
        double value = point->_X[_dim -1];
        if (_tree_root == 0) {  // Create first tree node
            _tree_root = new PointTreeNode(value);
            node = _tree_root;
            process_next_dimension(node, point);
        } else {
            while (true) {  // Walk through tree
                if (value > node->_value) {
                    if (node->_right == 0) {
                        node->_right = new PointTreeNode(value);
                        node->_right->_parent = node;
                        update_height(node->_right);
                        process_next_dimension(node->_right, point);
                        check_if_balanced(node->_right);
                        return 0;
                    };
                    node = node->_right;
                } else if (value < node->_value) {
                    if (node->_left == 0) {
                        node->_left = new PointTreeNode(value);
                        node->_left->_parent = node;
                        update_height(node->_left);
                        process_next_dimension(node->_left, point);
                        check_if_balanced(node->_left);
                        return 0;
                    };
                    node = node->_left;
                } else {
                    // Node value matches given point value, move to next dimension.
                    return process_next_dimension(node, point);
                };
            };
        };
    };

    void print(){
        _tree_root->print();
        cout << endl;
    };

    virtual ~PointTree(){
        delete _tree_root;
    };
};

PointTreeNode::~PointTreeNode() {
    if (_left != 0) { delete _left; };
    if (_right != 0) { delete _right; };
    if (_point != 0) { delete _point; };
    if (_subtree != 0) { delete _subtree; };
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


class Function {
    Function(const Function& other){};
    Function& operator=(const Function& other){};
public:
    Function(){
        _calls = 0;
        _f_min = numeric_limits<int>::max();
        _points = new PointTree(); // new Points();
    };
    string _name;
    int _D;
    Point* _lb;
    Point* _ub;
    Point* _glob_x;  // Point where global function minimum is (should be list)
    double _delta;   // Accuracy for stoping criteria based on distance from glob_x
    double _glob_f;  // Predefined global function minimum
    double _L;       // Global Lipschitz constant

    int _calls;
    double _f_min;  // Best known function value
    Point* _x_min;  // Point where best known function value is
    double _distance_to_glob_x;  // Infinity distance to _glob_x from nearest known point
    Point* _x_nearest_to_glob_x;  // Nearest to _glob_x known point (infinity distance) 
    PointTree* _points;  // Points* _points;

    double get_distance_to_glob_x(Point* p) {
        double max_distance = 0;
        double dist = 0;
        for (int i=0; i<_D; i++) {
            dist = fabs(p->_X[i] - _glob_x->_X[i]);
            if (dist > max_distance) {
                max_distance = dist;
            };
        };
        return max_distance;
    };

    void update_meta(Point* p) {
        double val = value(p);
        if (_f_min > val) {
            _f_min = val;
            _x_min = p;
        };
        double distance_to_glob_x = get_distance_to_glob_x(p);
        if (distance_to_glob_x < _distance_to_glob_x) {
            _distance_to_glob_x = distance_to_glob_x;
            _x_nearest_to_glob_x = p;
        };

        p->add_value(val);
        _calls += 1;
    };

    Point* get(double *c, int argc){
        Point* p = new Point(c, argc);
        Point* cached_point = _points->add(p);
        if (cached_point) {
            delete p;
            return cached_point;
        } else {
            update_meta(p);
            return p;
        };
    };

    Point* get(Point* p){  // Get Point with its function value
        Point* cached_point = _points->add(p);
        if (cached_point) {
            return cached_point;
        } else {
            update_meta(p);
            return p;
        };
    };

    void print(){
        cout << _name << "  calls: " << _calls << "   f_min: " << _f_min << endl;
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

    bool is_accurate_enougth(){
        for (int i=0; i<_D; i++) {
            if (_delta * (_ub->_X[i] - _lb->_X[i]) < fabs(_x_nearest_to_glob_x->_X[i] - _glob_x->_X[i])) { // Infinity norm
                return false;
            };
        };
        return true; 
    };

    virtual double value(Point* point) = 0;

    virtual ~Function(){
        delete _points;
    };
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
        double part1 = pow((x2 - 5./(4*pow(M_PI, 2))*pow(x1,2) + 5./M_PI*x1 -6), 2);
        double part2 = 10.*(1. - 1./(8*M_PI))*cos(x1) + 10.;
        return part1 + part2;
    };
};

class RastriginShrinked: public Function {
    RastriginShrinked(const RastriginShrinked& other){};
    RastriginShrinked& operator=(const RastriginShrinked& other){};
public:
    RastriginShrinked(): Function(){
        _name = "RastriginShrinked";
        _D = 2;
        _lb = new Point(-0.5, -0.5);
        _ub = new Point(1.25, 1.25);
        _glob_x = new Point(0., 0.);
        _glob_f = 0;          
        _L = 15.6204993518;
    };

    double value(Point* point) {
        double x1 = transform(point, 0); 
        double x2 = transform(point, 1); 
        return 2*10 + 4*pow(x1,2) + 4*pow(x2,2) - 10*(cos(2*M_PI*x1) + cos(2*M_PI*x2));
    };

    virtual ~RastriginShrinked(){
        delete _lb;
        delete _ub;
        delete _glob_x;
    };
};


class GKLSFunction: public Function {
    GKLSFunction(const GKLSFunction& other){};
    GKLSFunction& operator=(const GKLSFunction& other){};
public:
    GKLSFunction(int cls, int function_id): Function() {
        int _GKLS_class_D[] = {2, 2, 3, 3, 4, 4, 5, 5};
        double _GKLS_class_global_dists[] = {0.9, 0.9, 0.66, 0.9, 0.66, 0.9, 0.66, 0.66};
        double _GKLS_class_global_radiuses[] = {0.2, 0.1, 0.2, 0.2, 0.2, 0.2, 0.3, 0.2};
        double _GKLS_class_detlas[] = {1e-4, 1e-4, 1e-6, 1e-6, 1e-6, 1e-6, 1e-7, 1e-7};

        stringstream function_name; 
        function_name << cls << "_" << function_id;  
        _name =  function_name.str();
        _cls = cls;
        _fid = function_id;

        cls -= 1;
        _D = _GKLS_class_D[cls];
        _global_dist = _GKLS_class_global_dists[cls];
        _global_radius = _GKLS_class_global_radiuses[cls]; 
        _delta = pow(_GKLS_class_detlas[cls], 1./_D);   // _delta = _GKLS_class_detlas[cls];

        _lb = new Point(-1., _D);
        _ub = new Point(1., _D);
        _x_nearest_to_glob_x = _lb;
        _distance_to_glob_x = numeric_limits<double>::max();

        _glob_f = -1.;          // Predefined global function minimum
        // _L = ;

        // assert(GKLS_set_default()== GKLS_OK);     // Standartiniai nustatymai

        GKLS_dim = _D;
        GKLS_global_dist = _global_dist;
        GKLS_global_radius = _global_radius;
        GKLS_num_minima = 10;
        GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;
        assert(GKLS_domain_alloc() == GKLS_OK);
        // for (unsigned int i = 0; i < GKLS_dim; i++) {
        //     GKLS_domain_left[i] = -1;
        //     GKLS_domain_right[i] = 1;
        // };
        assert(GKLS_parameters_check() == GKLS_OK);
        assert(GKLS_arg_generate(function_id) == GKLS_OK);
        int n = GKLS_glob.num_global_minima;
        assert(n == 1);
        int glob_idx = 1;
        assert(GKLS_minima.f[glob_idx] == GKLS_global_value);
        _glob_x = new Point(_D);  // Point where global function minimum is (should be list)
        for (int i=0; i < _D; i++) {
            _glob_x->_X[i] = (GKLS_minima.local_min[glob_idx][i] - _lb->_X[i]) / (_ub->_X[i]-_lb->_X[i]);
        };
    };

    int _fid;
    int _cls;
    double _global_dist;
    double _global_radius;


    double transform(Point* point, int i) {     // Transforms single point coordinate from [0,1] to [l,u]
        return point->_X[i] * (_ub->_X[i]-_lb->_X[i]) + _lb->_X[i];  
    };

    double value(Point* point) {
        assert(GKLS_arg_generate(_fid) == GKLS_OK);   // Needed for multicriteria problems. Slows down function evaluations.
        double transformed_point[_D];
        for (int i=0; i<_D; i++){
            transformed_point[i] = transform(point,i);
        };
        double value = GKLS_D_func(transformed_point);
        return value;
    };

    virtual ~GKLSFunction(){
        delete _lb;
        delete _ub;
        delete _glob_x;
        GKLS_free();
        GKLS_domain_free();
    };
};

#endif
