#ifndef DISIMPLV_H
#define DISIMPLV_H 
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <math.h> 
#include <string>
#include <sys/time.h>
#include <vector>
#include "Eigen/Dense"
#include "functions.h"

using namespace std;


/* Utility functions */
double l2norm(Point* p1, Point* p2) {
    double squared_sum = 0;
    for (int i=0; i < p1->size(); i++){
        squared_sum += pow(p1->_X[i] - p2->_X[i], 2);
    };
    return sqrt(squared_sum);
};

double Determinant(double **a, int n) {
   /* Taken from http://paulbourke.net/miscellaneous/determinant/ */
    int i, j, j1, j2;
    double det = 0;
    double **m = NULL;

    if (n < 1) { /* Error */ cout << "Determinant cannot be calculated for empty matrix" << endl;
    } else if (n == 1) { /* Shouldn't get used */
        det = a[0][0];
    } else if (n == 2) {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {
        det = 0;
        for (j1=0;j1<n;j1++) {
            m = (double**) malloc((n-1)*sizeof(double *));
            for (i=0;i<n-1;i++)
                m[i] = (double*) malloc((n-1)*sizeof(double));
            for (i=1; i<n; i++) {
                j2 = 0;
                for (j=0; j<n; j++) {
                    if (j == j1) continue;
                    m[i-1][j2] = a[i][j];
                    j2++;
                }
            }
            det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant(m,n-1);
            for (i=0;i<n-1;i++) free(m[i]);
            free(m);
        }
    }
    return(det);
};

typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp() {
    struct timeval now;
    gettimeofday(&now, NULL);
    return now.tv_usec + (timestamp_t) now.tv_sec * 1000000;
};

/* Class definitions */
enum LowerBoundStrategy { MinVert, LongestEdgeLB, LowestEdgeLB, All };
const char* LBS[] = { "Min vert", "Longest edge LB", "Lowest edge LB", "All", 0 };
enum LStrategy { Self, ParentSelf, Neighbours };
const char* LS[] = { "Self", "Parent+Self", "Neighbours", 0 };
enum DivisionStrategy { LongestHalf };
const char* DS[] = { "Longest Half", 0 };
enum SimplexGradientStrategy { FFMinVert, FFMaxVert, FFAllVertMean };
const char* SGS[] = { "FF min vert", "FF max vert", "FF all vert mean", 0 };   // "All vert max" should match "Max vert"

class SimplexTree;
class SimplexTreeNode;


class Simplex {
    Simplex(const Simplex& other){}
    Simplex& operator=(const Simplex& other){}
public:
    Simplex(LowerBoundStrategy lower_bound_strategy,
            LStrategy L_strategy,
            double parent_L_part,
            SimplexGradientStrategy simplex_gradient_strategy
        ){  // Should accept Lower_Bound_strategy and L_strategy parameter.
        _lower_bound_strategy = lower_bound_strategy;
        _L_strategy = L_strategy;
        _parent_L_part = parent_L_part;
        _simplex_gradient_strategy = simplex_gradient_strategy;
        _is_in_partition = true;
        _parent = 0;
        _diameter = 0;
        _le_v1 = 0;
        _le_v2 = 0;
        _min_vert = 0;
        _max_vert = 0;
        _max_vert_value = numeric_limits<double>::min();
        _min_vert_value = numeric_limits<double>::max();
        _should_be_divided = false;
    };
    LowerBoundStrategy _lower_bound_strategy;
    LStrategy _L_strategy;
    SimplexGradientStrategy _simplex_gradient_strategy;

    vector<Point*> _verts;
    bool _is_in_partition;
    bool _should_be_divided;  // Should be divided in next iteration
    Simplex* _parent;

    Point* _le_v1;      // Longest edge vertex1
    Point* _le_v2; 
    double _diameter;   // Longest edge length

    double _L;          // Cumulative estimate of Lipschitz constant 
    double _grad_norm;  // Lipschitz constant estimate calculated by Simplex Gradient Euclidean norm.
    double _parent_L_part;

    Point* _min_vert;   // Pointer to vertex with lowest function value 
    double _min_vert_value;  // _min_vert function value 
    Point* _max_vert;
    double _max_vert_value;
    double _metric__vert_min_value;     // _f_min - glob_f / _diameter

    Point* _longest_edge_lb; 
    double _longest_edge_lb_value; 
    double _metric__longest_edge_lb;     // _f_min - glob_f / _diameter
    
    void init_parameters(Function* func){   // Called when all verts have been added

        // Note: claculating metrics needed by algorithm would reduce calculations
        double edge_length;  // Temporary variable
        for (int a=0; a < _verts.size(); a++) {
            // Finds _diameter
            for (int b=0; b < _verts.size(); b++){
                if (b > a) {
                    edge_length = l2norm(_verts[a], _verts[b]); 
                    if (edge_length > _diameter) {
                        _diameter = edge_length;
                        _le_v1 = _verts[a];
                        _le_v2 = _verts[b];
                    };
                };
            };
        }; 

        // Sort vertexes and set _min_vert, _max_vert
        sort(_verts.begin(), _verts.end(), Point::compare_by_value);
        _min_vert = _verts[0];
        _min_vert_value = _min_vert->_values[0];
        _max_vert = _verts[_verts.size()-1];
        _max_vert_value = _max_vert->_values[0];

        // Find adaptive Lipschitz constant
        _grad_norm = find_simplex_gradient_norm(_simplex_gradient_strategy);      // Check in the article if global Lipschitz constant is defined
        if (_L_strategy == Self) { _L = _grad_norm; };
        if (_L_strategy == ParentSelf) {
            if (_parent != 0) {
                _L = _parent_L_part * _parent->_grad_norm + (1 - _parent_L_part) * _grad_norm;
            } else {
                _L = _grad_norm;
            };
        };
        if (_L_strategy == Neighbours) {
            // _L and other parameters will be updated later after all simplexes were initialized
            return;
        };

        // Find longest edge lower bound
        _longest_edge_lb_value = find_edge_lb_value(_le_v1, _le_v2, _L);

        // epsilon = 1e-4
        double E;
        if (1e-4 * fabs(func->_glob_f) > 1e-8) {
            E = 1e-4 * fabs(func->_glob_f);
        } else {
            E = 1e-8;
        };
        _metric__vert_min_value = (_min_vert_value - (func->_glob_f + E)) / _diameter;

        // Note: Bus labiau dalinami prasti, bet dideli simplexai
        _metric__longest_edge_lb = (_longest_edge_lb_value - (func->_glob_f + E)) / _diameter;

        // Note: gali būti, kad slope apibrėžimas pas mane netinkamas atmetant
        // simpleksus su epsilon (potencialiai optimalių simpleksų parinkimo metu).



        // All needed parameters should be calculated here.
        // These are: 
        // self.D = len(verts) - 1
        // self.verts = verts
        // self.values = values
        // self.L = L
        // self.min_value = min(values)
        // self.diameter, self.le_vert1, self.le_vert2 = self.find_longest_edge(verts)
        // self.min_lb = self.find_lower_bound_minimum()
        // self.tolerance = abs(self.min_lb - self.min_value)
        // self.min_value_id = values.index(self.min_value)
        // self.max_value_id = values.index(max(values))
    };

    static void extend_region_with_vertex_neighbours(Point* vertex, SimplexTree* region, int depth);

    static void update_estimates(vector<Simplex*> simpls, Function* func);

    double find_simplex_gradient_norm(SimplexGradientStrategy simplex_gradient_strategy){  // Estimates L constant at minimum Simplex vertex point.
        int D = _verts.size() - 1;
        double L_estimate = 0;
        Eigen::VectorXd f_diff(D);
        Eigen::MatrixXd x_diff(D, D);
        Eigen::MatrixXd x_diff_inv_T(D, D);
        Eigen::VectorXd grad(D);


        if (simplex_gradient_strategy == FFMinVert) {  // Gradient at min vertex
            for (int i=1; i < D+1; i++) { 
                f_diff(i - 1) = _verts[i]->_values[0] - _min_vert->_values[0];
            }; 

            for (int i=1; i < D+1; i++) {
                for (int j=0; j < D; j++) {
                    x_diff(i-1, j) = _verts[i]->_X[j] - _min_vert->_X[j];
                };
            };
        };
        if (simplex_gradient_strategy == FFMaxVert) {  // Gradient at max vertex
            for (int i=0; i < D; i++) { 
                f_diff(i) = _verts[i]->_values[0] - _max_vert->_values[0];
            }; 

            for (int i=0; i < D; i++) {
                for (int j=0; j < D; j++) {
                    x_diff(i, j) = _verts[i]->_X[j] - _max_vert->_X[j];
                };
            };
        };
        if (simplex_gradient_strategy == FFMaxVert) {
            throw "FFMaxVert gradient strategy not implemented yet";
            // FFAllVertMean
        };

        // Find gradient at lowest point
        x_diff_inv_T = x_diff.inverse().transpose();
        for (int i=0; i < D; i++) {
            grad[i] = x_diff_inv_T.row(i).dot(f_diff);
        };
        // cout << grad << endl;

        // Find norm of gradient at _min_vert 
        for (int i=0; i < grad.size(); i++){
            L_estimate += pow(grad[i], 2);
        };
        L_estimate = sqrt(L_estimate);
        return L_estimate;
    };


    double find_edge_lb_value(Point* _le_v1, Point* _le_v2, double _L) {   // Needs testing
        double dist = l2norm(_le_v1, _le_v2);
        double x1[2] = {0., _le_v1->_values[0]};             // (0, simplex[0][-1]['obj'][0])
        double x2[2] = {dist, _le_v1->_values[0] - _L*dist}; // (dist, simplex[0][-1]['obj'][0] - L*dist)
        double x3[2] = {dist, _le_v2->_values[0]};           // (dist, simplex[1][-1]['obj'][0])
        double x4[2] = {0, _le_v1->_values[0] - _L*dist};    // (0, simplex[1][-1]['obj'][0] - L*dist)

        // 2D line intersection based on  http://mathworld.wolfram.com/Line-LineIntersection.html
        double av[2] = {x2[0] - x1[0], x2[1] - x1[1]};
        double bv[2] = {x4[0] - x3[0], x4[1] - x3[1]};
        double cv[2] = {x3[0] - x1[0], x3[1] - x1[1]};

        // cross_product(v1, v2) = (v1.X * v2.Y) - (v1.Y * v2.X)
        double s = ((cv[0]*bv[1] - cv[1]*bv[0]) * (av[0]*bv[1] - av[1]*bv[0])/ (pow(av[0]*bv[1] - av[1]*bv[0], 2) ));

        double intersection[2] = {x1[0] + (x2[0]-x1[0])*s, x1[1] + (x2[1]-x1[1])*s};
        // X = a(simplex[0][:-1]) + s[0]/float(dist) * (a(simplex[1][:-1]) - a(simplex[0][:-1]))
        // return [list(X) + [s[1]]]
        return intersection[1];
    };

    static bool wont_be_divided(Simplex* s) {
        return !s->_should_be_divided;
    };
    static bool not_in_partition(Simplex* s) {
        return !s->_is_in_partition;
    };

    static double compare_diameter(Simplex* s1, Simplex* s2) {
        return s1->_diameter < s2->_diameter; 
    };

    // static double compare_metric(Simplex* s1, Simplex* s2) {
    //     return s1->_metric__vert_min_value < s2->_metric__vert_min_value;
    // };

    void add_vertex(Point* vertex){
        _verts.push_back(vertex);
        vertex->_simplexes.push_back(this);
    };

    int size() {
        return _verts.size();
    };

    void print(){
        cout << " Simplex  (" << _diameter << ", " << _min_vert_value << "):" << endl;
        for (int i=0; i < _verts.size(); i++){
            _verts[i]->print();
        };
    };

    static void print(vector<Simplex*> simplexes, string label="Printing simplexes:"){
        cout << label << endl;
        for (int i=0; i < simplexes.size(); i++){
            simplexes[i]->print();
        };
    };

    static void log_partition(vector<Simplex*> simplexes,
                              vector<Simplex*> selected,
                              string label="Partition:",
                              int iteration=0){
       ofstream log_file; 
       log_file.open("log/partition.txt", ios::app);
       log_file << label << iteration << ":" << endl;
       for (int i=0; i < simplexes.size(); i++) {
           for (int j=0; j < simplexes[i]->_verts.size(); j++) {
               for (int k=0; k < simplexes[i]->_verts[j]->size(); k++){
                    log_file << simplexes[i]->_verts[j]->_X[k] << " ";
               };
               log_file << " (" << simplexes[i]->_verts[j]->_values[0]<<"); ";
           };
           if (simplexes[i]->_L_strategy == Self) {
                log_file << " ("<< simplexes[i]->_diameter << "," << simplexes[i]->_min_vert_value << ")" << endl;
           };
           if (simplexes[i]->_L_strategy == Neighbours) {
                log_file << " ("<< simplexes[i]->_diameter << "," << simplexes[i]->_longest_edge_lb_value << ")" << endl;
           };
       };
       log_file << "Selected:" << endl;
       for (int i=0; i < selected.size(); i++) {
           for (int j=0; j < selected[i]->_verts.size(); j++) {
               for (int k=0; k < selected[i]->_verts[j]->size(); k++){
                    log_file << selected[i]->_verts[j]->_X[k] << " ";
               };
               log_file << " (" << selected[i]->_verts[j]->_values[0]<<"); ";
           };
           if (selected[i]->_L_strategy == Self) {
               log_file << " ("<< selected[i]->_diameter << "," << selected[i]->_min_vert_value << ")" << endl;
           };
           if (selected[i]->_L_strategy == Neighbours) {
               log_file << " ("<< selected[i]->_diameter << "," << selected[i]->_longest_edge_lb_value << ")" << endl;
           };
       };

       log_file.close();
    };

    virtual ~Simplex(){
        _verts.clear();
    };  
};


class SimplexTreeNode {
    SimplexTreeNode(const SimplexTreeNode& other){}
    SimplexTreeNode& operator=(const SimplexTreeNode& other){}
public:                
    SimplexTreeNode(Simplex* value){
        _height = 1;
        _value = value;
        _parent = 0;
        _left = 0;
        _right = 0;
    };
    int _height;
    Simplex* _value;
    SimplexTreeNode* _parent;
    SimplexTreeNode* _left;
    SimplexTreeNode* _right;

    void print(){
        if (_left != 0) { cout << "l"; _left->print(); };
        cout << _value << "("<< _height << ")";
        if (_right != 0) { cout << "r"; _right->print(); };
    };

    virtual ~SimplexTreeNode();
};

class SimplexTree {  // Binary balancing tree 
    SimplexTree(const SimplexTree& other){}
    SimplexTree& operator=(const SimplexTree& other){}
public:
    SimplexTree(){
         _max_grad_norm = numeric_limits<double>::min();
         _tree_root = 0;
    };
    double _max_grad_norm;
    SimplexTreeNode* _tree_root;

    void update_height(SimplexTreeNode* node) {
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
    void left_right_rebalance(SimplexTreeNode* node) {
        SimplexTreeNode* diatteched_node;
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
    void left_left_rebalance(SimplexTreeNode* node) {
        SimplexTreeNode* diatteched;
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
    void right_left_rebalance(SimplexTreeNode* node) {
        SimplexTreeNode* diatteched_node;
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
    void right_right_rebalance(SimplexTreeNode* node) {
        SimplexTreeNode* diatteched;
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

    void check_if_balanced(SimplexTreeNode* node) {  // Rebalances tree if its not balanced
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

    Simplex* add(Simplex* value) {
        // Get same point or insert given (if inserted returns 0)
        SimplexTreeNode* node = _tree_root;
        double grad_norm = value->_grad_norm;
        if (grad_norm > _max_grad_norm) {
            _max_grad_norm = grad_norm;
        };
        if (_tree_root == 0) {  // Create first tree node
            _tree_root = new SimplexTreeNode(value);
        } else {
            while (true) {  // Walk through tree
                if (value > node->_value) {
                    if (node->_right == 0) {
                        node->_right = new SimplexTreeNode(value);
                        node->_right->_parent = node;
                        update_height(node->_right);
                        check_if_balanced(node->_right);
                        return 0;
                    };
                    node = node->_right;
                } else if (value < node->_value) {
                    if (node->_left == 0) {
                        node->_left = new SimplexTreeNode(value);
                        node->_left->_parent = node;
                        update_height(node->_left);
                        check_if_balanced(node->_left);
                        return 0;
                    };
                    node = node->_left;
                } else {
                    // Node value matches given simplex value
                    return value;
                };
            };
        };
    };

    void print() {
        _tree_root->print();
        cout << endl;
    };

    virtual ~SimplexTree();

};

SimplexTreeNode::~SimplexTreeNode() {
    if (_left != 0) { delete _left; };
    if (_right != 0) { delete _right; };
};

SimplexTree::~SimplexTree() {
    delete _tree_root;
};


void Simplex::extend_region_with_vertex_neighbours(Point* vertex, SimplexTree* region, int depth) {
    // Recursively adds vertex neighbours to region
    for (int sid=0; sid < vertex->_simplexes.size(); sid++) {
        Simplex* simpl = vertex->_simplexes[sid];
        Simplex* result = region->add(simpl);
        if (depth != 0) {
            if (result == 0 && simpl->_is_in_partition) {
                for (int vid=0; vid < simpl->_verts.size(); vid++) {
                    if (simpl->_verts[vid] != vertex) {
                        extend_region_with_vertex_neighbours(simpl->_verts[vid], region, depth-1);
                    };
                };
            };
        };    
    };
};

void Simplex::update_estimates(vector<Simplex*> simpls, Function* func) {   // Neighbours strategy - updates estimates
    SimplexTree* region;
    int depth = 1;  // What about higher depth?
    for (int sid=0; sid < simpls.size(); sid++) {
        region = new SimplexTree();
        for (int vid=0; vid < simpls[sid]->_verts.size(); vid++) {
            extend_region_with_vertex_neighbours(simpls[sid]->_verts[vid], region, depth);
        };
        simpls[sid]->_L = region->_max_grad_norm;
        simpls[sid]->_longest_edge_lb_value  = simpls[sid]->find_edge_lb_value(simpls[sid]->_le_v1, simpls[sid]->_le_v2, simpls[sid]->_L); 
        
        double E;
        if (1e-4 * fabs(func->_glob_f) > 1e-8) {
            E = 1e-4 * fabs(func->_glob_f);
        } else {
            E = 1e-8;
        };
        simpls[sid]->_metric__longest_edge_lb = (simpls[sid]->_longest_edge_lb_value - (func->_glob_f + E)) / simpls[sid]->_diameter;
        delete region;
    };
//         simpls[sid]->_lb = simpls[sid].find_lb();
//     };
//     // for (int i=0; i < _verts.size(); i++) {
//     //     cout << _verts[i]->_simplexes.size() << endl;
//     // };
//     
};

class Algorithm {
    Algorithm(const Algorithm& other){};
    Algorithm& operator=(const Algorithm& other){};
public:
    Algorithm(){};
    string _name;
    int _max_calls;
    string _stop_criteria;
    double _min_pe;
    double _duration;   // Duration in seconds
    double _epsilon;    // How small simplexes should still be partitioned 
    vector<Simplex*> _partition;
    vector<Simplex*> _all_simplexes;
    Function* _func;

    LowerBoundStrategy _lower_bound_strategy;     // Simplex lower bound estimation strategy 
    LStrategy _L_strategy;  // Lipschitz constant estimation strategy
    DivisionStrategy _division_strategy;      // Simplex edge division strategy
    SimplexGradientStrategy _simplex_gradient_strategy;
    double _parent_L_part;

    /* Global optimization strategies */
    void partition_feasable_region_combinatoricly(){
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

            Simplex* simpl = new Simplex(_lower_bound_strategy, _L_strategy, _parent_L_part, _simplex_gradient_strategy);
            for (int i=0; i < n + 1; i++){
                Point* tmp_point = new Point(triangle[i], n);
                Point* point = _func->get(tmp_point); 
                if (tmp_point != point) {
                    delete tmp_point;
                };
                simpl->add_vertex(point);
            };
            simpl->init_parameters(_func);
            _partition.push_back(simpl);
            _all_simplexes.push_back(simpl);

        } while (next_permutation(teta, teta+n));
    };

    vector<Simplex*> divide_simplex(Simplex* simplex, string strategy="longest_half") {
        vector<Simplex*> divided_simplexes;
        if (strategy== "longest_half") {
            // Find middle point
            int n = _func->_D;
            double c[n];
            for (int i=0; i < n; i++) {
                c[i] = (simplex->_le_v1->_X[i] + simplex->_le_v2->_X[i]) / 2.;
            };
            Point* middle_point = _func->get(c, n);
            // Construct two new simplexes using this middle point.
            Simplex* left_simplex = new Simplex(_lower_bound_strategy, _L_strategy, _parent_L_part, _simplex_gradient_strategy);
            Simplex* right_simplex = new Simplex(_lower_bound_strategy, _L_strategy, _parent_L_part, _simplex_gradient_strategy);

            for (int i=0; i < simplex->size(); i++){
                // Point* point = _func->get(new Point(triangle[i], n)); 
                if (simplex->_verts[i] != simplex->_le_v1){
                    right_simplex->add_vertex(simplex->_verts[i]);
                } else {
                    right_simplex->add_vertex(middle_point);
                };
                if (simplex->_verts[i] != simplex->_le_v2) {
                    left_simplex->add_vertex(simplex->_verts[i]);
                } else {
                    left_simplex->add_vertex(middle_point);
                };
            };
            left_simplex->_parent = simplex;
            right_simplex->_parent = simplex;
            left_simplex->init_parameters(_func);
            right_simplex->init_parameters(_func);
            simplex->_is_in_partition = false;

            divided_simplexes.push_back(left_simplex);
            divided_simplexes.push_back(right_simplex);
            return divided_simplexes;
        };
    };

    friend ostream& operator<<(ostream& o, const Algorithm& a){
        o << "alg: '" << a._name << "', func: '" << a._func->_name << 
           "', calls: " << a._func->_calls << ", subregions: " << a._partition.size() <<
           ", duration: " << a._duration << ", stop_criteria: '" << a._stop_criteria << 
           "', f_min: " << a._func->_f_min << ", x_min: " << (*a._func->_x_min);
        return o;
    };

    /* Test functions */
    void test_unique_simplexes(){
        //// Test to check if simplexes are unique in _partition.
        for (int i=0; i < _partition.size(); i++){
            for (int j=0; j < _partition.size(); j++){
                if (i > j) {
                    // Compare all vertexes
                    bool same = true;
                    for (int k=0; k < _partition[i]->_verts.size(); k++){
                        bool found_same_vert = false;
                        for (int l=0; l < _partition[j]->_verts.size(); l++){
                            if (_partition[i]->_verts[k] == _partition[j]->_verts[l]){
                                found_same_vert = true;
                            };
                        };
                        if (!found_same_vert) {
                            same = false;
                        };
                    };
                    if (same){
                        cout << "NOT UNIQUE PARTITION  " << i << "==" << j << endl;
                        _partition[i]->print();
                        _partition[j]->print();
                    };
                };
            };
        };
    };

    virtual ~Algorithm(){
        for (int i=0; i < _all_simplexes.size(); i++) {
            delete _all_simplexes[i];
        };
        _all_simplexes.clear();
        _partition.clear();
    };
};


class Disimplv : public Algorithm {
    Disimplv(const Disimplv& other){};
    Disimplv& operator=(const Disimplv& other){};
public:
    Disimplv(LowerBoundStrategy lower_bound_strategy=MinVert,
             LStrategy L_estimation_strategy=Self,
             DivisionStrategy division_strategy=LongestHalf,
             SimplexGradientStrategy simplex_gradient_strategy=FFMinVert,
             double parent_L_part=0.1,
             double epsilon=0.0001,
             int max_calls=1000000){
        _lower_bound_strategy = lower_bound_strategy;
        _L_strategy = L_estimation_strategy;
        _simplex_gradient_strategy = simplex_gradient_strategy;
        _division_strategy = division_strategy;
        _stop_criteria = "x_dist_Serg";
        _max_calls = max_calls;
        _epsilon = epsilon;
        _parent_L_part = parent_L_part;
        // _min_pe = min_pe;

        // Construct algorithm name
        if (_lower_bound_strategy == MinVert) { _name = "Disimpl-v"; };
        if (_lower_bound_strategy == LongestEdgeLB) { _name = "Disimpl-2v-le"; }; 
        stringstream alg_name; 
        alg_name << _name << "_e" << epsilon;  
        _name =  alg_name.str();

        // Clean partition log file
        ofstream log_file; 
        log_file.open("log/partition.txt");
        log_file.close();
    };

    int nextv(int v, int m) {
        if (v == m) {
            return 0;
        };
        return v + 1;
    };
    int predv(int v, int m) {
        if (v == 0) {
            return m;
        };
        return v - 1;
    };

    vector<Simplex*> convex_hull(vector<Simplex*> simplexes) {
        int m = simplexes.size() - 1;
        if (m <= 1) { return simplexes; };
        int START = 0;
        int v = START;
        int w = m;
        bool flag = false;
        bool leftturn = false;
        int a, b, c;
        double det_val;
        while ((nextv(v, m) != START) or (flag == false)) {
            if (nextv(v, m) == w) {
                flag = true;
            }
            a = v;
            b = nextv(v, m);
            c = nextv(nextv(v, m), m);   // d = x = _diameter;  f = y = _min_vert_value;

            double* matrix[3];
            double line1[3] = {simplexes[a]->_diameter, simplexes[a]->_min_vert_value, 1.};
            double line2[3] = {simplexes[b]->_diameter, simplexes[b]->_min_vert_value, 1.};
            double line3[3] = {simplexes[c]->_diameter, simplexes[c]->_min_vert_value, 1.};
            matrix[0] = line1;
            matrix[1] = line2;
            matrix[2] = line3;
            det_val = Determinant(matrix, 3);

            if (det_val >= 0){
                leftturn = 1;
            } else {
                leftturn = 0;
            };
            if (leftturn) {
                v = nextv(v, m);
            } else {
                simplexes.erase(simplexes.begin() + nextv(v, m));
                m -= 1;
                w -= 1;
                v = predv(v, m);
            };
        };
        return simplexes;
    };

    vector<Simplex*> select_simplexes_by_min_vert() {
        vector<Simplex*> selected_simplexes;

        // Sort simplexes by their diameter
        vector<Simplex*> sorted_partition = _partition;   // Note: Could sort globally, resorting would take less time

        // Simplex::print(sorted_partition, "Selecting for division from: ");

        sort(sorted_partition.begin(), sorted_partition.end(), Simplex::compare_diameter);
        double f_min = _func->_f_min;

        // Find simplex with  minimum metric  and  unique diameters
        Simplex* min_metric_simplex = sorted_partition[0]; // Initial value
        vector<double> diameters;
        vector<Simplex*> best_for_size;

        bool unique_diameter;
        bool found_with_same_size;
        for (int i=0; i < sorted_partition.size(); i++) {
            if (sorted_partition[i]->_metric__vert_min_value < min_metric_simplex->_metric__vert_min_value) {
                min_metric_simplex = sorted_partition[i];
            };
            // Saves unique diameters
            unique_diameter = true;
            for (int j=0; j < diameters.size(); j++) {
                if (diameters[j] == sorted_partition[i]->_diameter) {
                    unique_diameter = false; break;
                };
            };
            if (unique_diameter) {
                diameters.push_back(sorted_partition[i]->_diameter);
            };

            // If this simplex is better then previous with same size swap them.
            found_with_same_size = false;
            for (int j=0; j < best_for_size.size(); j++) {
                if (best_for_size[j]->_diameter == sorted_partition[i]->_diameter){
                    found_with_same_size = true;
                    if (best_for_size[j]->_min_vert_value > sorted_partition[i]->_min_vert_value) {
                        best_for_size.erase(best_for_size.begin()+j);
                        best_for_size.push_back(sorted_partition[i]);
                    };
                };
            };
            if (!found_with_same_size) {
                best_for_size.push_back(sorted_partition[i]);
            };
        };

        vector<Simplex*> selected;
        // Is this OK?  Well compared with examples - its ok.
        if ((best_for_size.size() > 2) ) { // && (min_metric_simplex != best_for_size[best_for_size.size()-1])
            vector<Simplex*> simplexes_below_line;
            double a1 = best_for_size[0]->_diameter;
            double b1 = best_for_size[0]->_min_vert_value;
            // double a1 = min_metric_simplex->_diameter;  // Should be like this based on Direct Matlab implementation
            // double b1 = min_metric_simplex->_min_vert_value;
            double a2 = best_for_size[best_for_size.size()-1]->_diameter;
            double b2 = best_for_size[best_for_size.size()-1]->_min_vert_value;
            
            double slope = (b2 - b1)/(a2 - a1);
            double bias = b1 - slope * a1;

            for (int i=0; i < best_for_size.size(); i++) {
                if (best_for_size[i]->_min_vert_value < slope*best_for_size[i]->_diameter + bias +1e-12) {
                    simplexes_below_line.push_back(best_for_size[i]);
                };
            };
            selected = convex_hull(simplexes_below_line);  // Messes up simplexes_below_line
        } else {
            selected = best_for_size;    // TODO: Why we divide all of them? Could divide only min_metrc_simplex.
                                         // Because practiacally this case does not occur ever.
        };


        for (int i=0; i < selected.size(); i++) {
            selected[i]->_should_be_divided = true;
        };

        // Remove simplexes which do not satisfy condition:   f - slope*d > f_min - epsilon*abs(f_min)
        for (int i=0; i < selected.size() -1; i++) {
            double a1 = selected[selected.size() - i -1]->_diameter;
            double b1 = selected[selected.size() - i -1]->_min_vert_value;
            double a2 = selected[selected.size() - i -2]->_diameter;
            double b2 = selected[selected.size() - i -2]->_min_vert_value;
            double slope = (b2 - double(b1))/(a2 - a1);
            double bias = b1 - slope * a1;

            if (bias > f_min - 0.0001*fabs(f_min)) {   // epsilon
                selected[selected.size() - i -2]->_should_be_divided = false;
            };
        };

        // Remove simplexes which should not be divided
        selected.erase(remove_if(selected.begin(), selected.end(), Simplex::wont_be_divided), selected.end());

        // Select all simplexes which have best _min_vert_value for its size 
        for (int i=0; i < sorted_partition.size(); i++) {
            for (int j=0; j < selected.size(); j++) {
                if ((sorted_partition[i]->_diameter == selected[j]->_diameter) && 
                    (sorted_partition[i]->_min_vert_value == selected[j]->_min_vert_value)) {
                    selected_simplexes.push_back(sorted_partition[i]);
                };
            };
        };

        return selected_simplexes;
    };

    vector<Simplex*> select_simplexes_by_longest_edge_lb(int iteration){
        // Update this method. Should not use convex hull, instead select simplex with lowest lower bound.
        vector<Simplex*> selected_simplexes;
        vector<Simplex*> simplexes = _partition;

        // *****   Strategy:  select simplexes with minimal lower bound value  *****
        // double min_lower_bound_value = numeric_limits<double>::max();
        // for (int sid=0; sid < simplexes.size(); sid++) {
        //     if (simplexes[sid]->_longest_edge_lb_value == min_lower_bound_value) {
        //         selected_simplexes.push_back(simplexes[sid]);
        //     } else {
        //         if (simplexes[sid]->_longest_edge_lb_value < min_lower_bound_value) {
        //             selected_simplexes.clear(); // Test if clear works properly 
        //             min_lower_bound_value = simplexes[sid]->_longest_edge_lb_value;
        //             selected_simplexes.push_back(simplexes[sid]);
        //         };
        //     };
        // };



        // *****   Strategy:  convex hull from smallest lower_bound to largest diameter *****
        // Sort simplexes by their diameter
        vector<Simplex*> sorted_partition = _partition;   // Note: Could sort globally, resorting would take less time
        sort(sorted_partition.begin(), sorted_partition.end(), Simplex::compare_diameter);
        double f_min = _func->_f_min;

        // Find simplex with  minimum metric  and  unique diameters
        // Simplex* min_value_simplex = sorted_partition[0];  // Initial value
        int min_value_simplex_id = 0;
        vector<double> diameters;
        vector<Simplex*> best_for_size;

        bool unique_diameter;
        bool found_with_same_size;
        for (int i=0; i < sorted_partition.size(); i++) {
                // Saves unique diameters
                unique_diameter = true;
                for (int j=0; j < diameters.size(); j++) {
                    if (diameters[j] == sorted_partition[i]->_diameter) {
                        unique_diameter = false; break;
                    };
                };
                if (unique_diameter) {
                    diameters.push_back(sorted_partition[i]->_diameter);
                };

                // If this simplex is better then previous with same size swap them.
                found_with_same_size = false;
                for (int j=0; j < best_for_size.size(); j++) {
                    if (best_for_size[j]->_diameter == sorted_partition[i]->_diameter){
                        found_with_same_size = true;
                        if (best_for_size[j]->_longest_edge_lb_value > sorted_partition[i]->_longest_edge_lb_value) {
                            best_for_size.erase(best_for_size.begin()+j);
                            best_for_size.push_back(sorted_partition[i]);
                        };
                    };
                };
                if (!found_with_same_size) {
                    best_for_size.push_back(sorted_partition[i]);
                };

                if (sorted_partition[i]->_longest_edge_lb_value < best_for_size[min_value_simplex_id]->_longest_edge_lb_value) {
                    // if (iteration == 9) { cout << "Greater" << endl; };
                    min_value_simplex_id = best_for_size.size() - 1;
                };
            // if (iteration == 9) { cout << "Comparing: (" << i<<")"<< sorted_partition[i]->_longest_edge_lb_value << " and ("<< min_value_simplex_id << ")" << best_for_size[min_value_simplex_id]->_longest_edge_lb_value << endl; };
            };



            vector<Simplex*> selected;
            // Select from minimum value simplex to bigest simplex using convex hull
            // Get minimum simplex id in best_for_size:  min_value_simplex_id
            // Find out how minimum value simplex id in best_for_size, not in sorted_partition 
            if ((best_for_size.size() - min_value_simplex_id) > 2) {
                vector<Simplex*> simplexes_below_line;
                // double a1 = best_for_size[0]->_diameter;
                // double b1 = best_for_size[0]->_longest_edge_lb_value;
                double a1 = best_for_size[min_value_simplex_id]->_diameter;  // Should be like this based on Direct Matlab implementation
                double b1 = best_for_size[min_value_simplex_id]->_longest_edge_lb_value;

                double a2 = best_for_size[best_for_size.size()-1]->_diameter;
                double b2 = best_for_size[best_for_size.size()-1]->_longest_edge_lb_value;
                
                double slope = (b2 - b1)/(a2 - a1);
                double bias = b1 - slope * a1;

                for (int i=min_value_simplex_id; i < best_for_size.size(); i++) {
                    if (best_for_size[i]->_longest_edge_lb_value < slope*best_for_size[i]->_diameter + bias +1e-12) {
                        simplexes_below_line.push_back(best_for_size[i]);
                    };
                };
                selected = convex_hull(simplexes_below_line);  // Messes up simplexes_below_line

            } else {
                for (int i=min_value_simplex_id; i < best_for_size.size(); i++) {
                    selected.push_back(best_for_size[i]);
                };
                // if (iteration >= 9) {
                //      cout << iteration << ". Min value simplex: " << min_value_simplex_id << " total simplexes: "<< best_for_size.size() << " selected: " << selected.size() << endl;
                //      cout << "Values: " << best_for_size[0]->_longest_edge_lb_value << ", " << best_for_size[1]->_longest_edge_lb_value << endl;
                //      exit(0);
                // };
            };

            


            // // if (iteration == 24) {cout << (best_for_size.size() > 2) << ", " << (min_metric_simplex != best_for_size[best_for_size.size()-1]);};
            // if ((best_for_size.size() > 2) ) { // && (min_metric_simplex != best_for_size[best_for_size.size()-1])
            //     vector<Simplex*> simplexes_below_line;
            //     double a1 = best_for_size[0]->_diameter;
            //     double b1 = best_for_size[0]->_longest_edge_lb_value;
            //     // double a1 = min_metric_simplex->_diameter;  // Should be like this based on Direct Matlab implementation
            //     // double b1 = min_metric_simplex->_min_value;
            //     double a2 = best_for_size[best_for_size.size()-1]->_diameter;
            //     double b2 = best_for_size[best_for_size.size()-1]->_longest_edge_lb_value;
            //
            //     double slope = (b2 - b1)/(a2 - a1);
            //     double bias = b1 - slope * a1;
            //
            //     for (int i=0; i < best_for_size.size(); i++) {
            //         if (best_for_size[i]->_longest_edge_lb_value < slope*best_for_size[i]->_diameter + bias +1e-12) {
            //             simplexes_below_line.push_back(best_for_size[i]);
            //         };
            //     };
            //     selected = convex_hull(simplexes_below_line);  // Messes up simplexes_below_line
            // } else {
            //     selected = best_for_size;    // TODO: Why we divide all of them? Could divide only min_metrc_simplex.
            //                                  // Because practiacally this case does not occur ever.
            // };

            for (int i=0; i < selected.size(); i++) {
                selected[i]->_should_be_divided = true;
            };

            // Remove simplexes which do not satisfy condition:   f - slope*d > f_min - epsilon*abs(f_min)
            for (int i=0; i < selected.size() -1; i++) {
                double a1 = selected[selected.size() - i -1]->_diameter;
                double b1 = selected[selected.size() - i -1]->_longest_edge_lb_value;
                double a2 = selected[selected.size() - i -2]->_diameter;
                double b2 = selected[selected.size() - i -2]->_longest_edge_lb_value;
                double slope = (b2 - double(b1))/(a2 - a1);
                double bias = b1 - slope * a1;

                if (bias > f_min - _epsilon * fabs(f_min)) {   // epsilon
                    selected[selected.size() - i -2]->_should_be_divided = false;
                };
            };

            // Remove simplexes which should not be divided
            selected.erase(remove_if(selected.begin(), selected.end(), Simplex::wont_be_divided), selected.end());

            // Select all simplexes which have best longest_edge_lb_value for its size 
            for (int i=0; i < sorted_partition.size(); i++) {
                for (int j=0; j < selected.size(); j++) {
                    if ((sorted_partition[i]->_diameter == selected[j]->_diameter) && 
                        (sorted_partition[i]->_longest_edge_lb_value == selected[j]->_longest_edge_lb_value)) {
                        selected_simplexes.push_back(sorted_partition[i]);
                    };
                };
            };




        // *****   Strategy:  convex hull from smallest diameter till largest diameter *****
        // // Sort simplexes by their diameter
        // vector<Simplex*> sorted_partition = _partition;   // Note: Could sort globally, resorting would take less time
        // sort(sorted_partition.begin(), sorted_partition.end(), Simplex::compare_diameter);
        // double f_min = _func->_f_min;
        //
        // // Find simplex with  minimum metric  and  unique diameters
        // Simplex* min_metric_simplex = sorted_partition[0];  // Initial value
        // vector<double> diameters;
        // vector<Simplex*> best_for_size;
        //
        // bool unique_diameter;
        // bool found_with_same_size;
        // for (int i=0; i < sorted_partition.size(); i++) {
        //     if (sorted_partition[i]->_metric__longest_edge_lb < min_metric_simplex->_metric__longest_edge_lb) {
        //         min_metric_simplex = sorted_partition[i];
        //     };
        //         // Saves unique diameters
        //         unique_diameter = true;
        //         for (int j=0; j < diameters.size(); j++) {
        //             if (diameters[j] == sorted_partition[i]->_diameter) {
        //                 unique_diameter = false; break;
        //             };
        //         };
        //         if (unique_diameter) {
        //             diameters.push_back(sorted_partition[i]->_diameter);
        //         };
        //
        //         // If this simplex is better then previous with same size swap them.
        //         found_with_same_size = false;
        //         for (int j=0; j < best_for_size.size(); j++) {
        //             if (best_for_size[j]->_diameter == sorted_partition[i]->_diameter){
        //                 found_with_same_size = true;
        //                 if (best_for_size[j]->_longest_edge_lb_value > sorted_partition[i]->_longest_edge_lb_value) {
        //                     best_for_size.erase(best_for_size.begin()+j);
        //                     best_for_size.push_back(sorted_partition[i]);
        //                 };
        //             };
        //         };
        //         if (!found_with_same_size) {
        //             best_for_size.push_back(sorted_partition[i]);
        //         };
        //     };
        //
        //     vector<Simplex*> selected;
        //     // Is this OK?
        //     // if (iteration == 24) {cout << (best_for_size.size() > 2) << ", " << (min_metric_simplex != best_for_size[best_for_size.size()-1]);};
        //     if ((best_for_size.size() > 2) ) { // && (min_metric_simplex != best_for_size[best_for_size.size()-1])
        //         vector<Simplex*> simplexes_below_line;
        //         double a1 = best_for_size[0]->_diameter;
        //         double b1 = best_for_size[0]->_longest_edge_lb_value;
        //         // double a1 = min_metric_simplex->_diameter;  // Should be like this based on Direct Matlab implementation
        //         // double b1 = min_metric_simplex->_min_value;
        //         double a2 = best_for_size[best_for_size.size()-1]->_diameter;
        //         double b2 = best_for_size[best_for_size.size()-1]->_longest_edge_lb_value;
        //         
        //         double slope = (b2 - b1)/(a2 - a1);
        //         double bias = b1 - slope * a1;
        //
        //         for (int i=0; i < best_for_size.size(); i++) {
        //             if (best_for_size[i]->_longest_edge_lb_value < slope*best_for_size[i]->_diameter + bias +1e-12) {
        //                 simplexes_below_line.push_back(best_for_size[i]);
        //             };
        //         };
        //         selected = convex_hull(simplexes_below_line);  // Messes up simplexes_below_line
        //     } else {
        //         selected = best_for_size;    // TODO: Why we divide all of them? Could divide only min_metrc_simplex.
        //                                      // Because practiacally this case does not occur ever.
        //     };
        //
        //     for (int i=0; i < selected.size(); i++) {
        //         selected[i]->_should_be_divided = true;
        //     };
        //
        //     // Remove simplexes which do not satisfy condition:   f - slope*d > f_min - epsilon*abs(f_min)
        //     for (int i=0; i < selected.size() -1; i++) {
        //         double a1 = selected[selected.size() - i -1]->_diameter;
        //         double b1 = selected[selected.size() - i -1]->_longest_edge_lb_value;
        //         double a2 = selected[selected.size() - i -2]->_diameter;
        //         double b2 = selected[selected.size() - i -2]->_longest_edge_lb_value;
        //         double slope = (b2 - double(b1))/(a2 - a1);
        //         double bias = b1 - slope * a1;
        //
        //         if (bias > f_min - _epsilon * fabs(f_min)) {   // epsilon
        //             selected[selected.size() - i -2]->_should_be_divided = false;
        //         };
        //     };
        //
        //     // Remove simplexes which should not be divided
        //     selected.erase(remove_if(selected.begin(), selected.end(), Simplex::wont_be_divided), selected.end());
        //
        //     // Select all simplexes which have best _min_value for its size 
        //     for (int i=0; i < sorted_partition.size(); i++) {
        //         for (int j=0; j < selected.size(); j++) {
        //             if ((sorted_partition[i]->_diameter == selected[j]->_diameter) && 
        //                 (sorted_partition[i]->_longest_edge_lb_value == selected[j]->_longest_edge_lb_value)) {
        //                 selected_simplexes.push_back(sorted_partition[i]);
        //             };
        //         };
        //     };
        return selected_simplexes;
    };

    vector<Simplex*> select_simplexes_by_lowest_edge_lb(){};

    virtual vector<Simplex*> select_simplexes_to_divide(int iteration=0){
        vector<Simplex*> selected_simplexes;
        if (_lower_bound_strategy == MinVert) {
            selected_simplexes = select_simplexes_by_min_vert();
        };
        if (_lower_bound_strategy == LongestEdgeLB) {
            selected_simplexes = select_simplexes_by_longest_edge_lb(iteration);
        };
        if (_lower_bound_strategy == LowestEdgeLB){
            selected_simplexes = select_simplexes_by_lowest_edge_lb();
        };
        if (_lower_bound_strategy == All) {
            selected_simplexes = _partition;
        };
        for (int i=0; i < selected_simplexes.size(); i++){
            selected_simplexes[i]->_is_in_partition = false;
        };
        return selected_simplexes;
    };

    void minimize(Function* func){
        _func = func;
        timestamp_t start = get_timestamp();
        partition_feasable_region_combinatoricly();
        if (_L_strategy == Neighbours) { Simplex::update_estimates(_partition, _func); };

        int iteration = 0;
        while (_func->_calls <= _max_calls && !_func->is_accurate_enougth()) { // _func->pe() > _min_pe){
            // cout << "Iteration: " << iteration << endl;

            // Selects simplexes to divide
            vector<Simplex*> simplexes_to_divide;
            if (iteration == 0) {
                simplexes_to_divide = _partition;
            } else {
                simplexes_to_divide = select_simplexes_to_divide(iteration);
            };
            Simplex::log_partition(_partition, simplexes_to_divide, "\nIteration ", iteration);
            // test_unique_simplexes();

            // Divides selected simplexes
            vector<Simplex*> new_simplexes;
            for (int i=0; i < simplexes_to_divide.size(); i++) {
                vector<Simplex*> divided_simplexes = divide_simplex(simplexes_to_divide[i]);

                for (int j=0; j < divided_simplexes.size(); j++) {
                    new_simplexes.push_back(divided_simplexes[j]);
                };
            };

            // Remove partitioned simplexes from _partition
            _partition.erase(remove_if(_partition.begin(), _partition.end(), Simplex::not_in_partition), _partition.end());

            // Add new simplexes to _partition and _all_simplexes
            for (int i=0; i < new_simplexes.size(); i++) {
                _partition.push_back(new_simplexes[i]);
                _all_simplexes.push_back(new_simplexes[i]);
            };
            if (_L_strategy == Neighbours) { Simplex::update_estimates(_partition, _func); };

            // Update counters and log the status
            iteration += 1;
            cout << iteration << ". Simplexes: " << _partition.size() << "  calls: " << _func->_calls << endl;

            // if (iteration >= 1) {
            //     break;
            // };
        };
        timestamp_t end = get_timestamp();
        _duration = (end - start) / 1000000.0L;

        // Draw partitioning: output simplex coordinates to file and draw it with Python
    };

    virtual ~Disimplv(){};
};

#endif
