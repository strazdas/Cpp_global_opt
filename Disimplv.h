#ifndef DISIMPLV_H
#define DISIMPLV_H_ 
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
            // Self, ParentSelf, Neighbours
        if (_L_strategy == ParentSelf) {
            if (_parent != 0) {
                _L = _parent_L_part * _parent->_grad_norm + (1 - _parent_L_part) * _grad_norm;
            } else {
                _L = _grad_norm;
            };
        };
        if (_L_strategy == Neighbours) {
            _L = _grad_norm;
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
    };

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

    static void log_partition(vector<Simplex*> simplexes, vector<Simplex*> selected, string label="Partition:", int iteration=0){
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
           log_file << " ("<< simplexes[i]->_diameter << "," << simplexes[i]->_min_vert_value << ")" << endl;
       };
       log_file << "Selected:" << endl;
       for (int i=0; i < selected.size(); i++) {
           for (int j=0; j < selected[i]->_verts.size(); j++) {
               for (int k=0; k < selected[i]->_verts[j]->size(); k++){
                    log_file << selected[i]->_verts[j]->_X[k] << " ";
               };
               log_file << " (" << selected[i]->_verts[j]->_values[0]<<"); ";
           };
           log_file << " ("<< selected[i]->_diameter << "," << selected[i]->_min_vert_value << ")" << endl;
       };

       log_file.close();
    };

    virtual ~Simplex(){
        _verts.clear();
    };  
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

    vector<Simplex*> select_simplexes_by_longest_edge_lb(){
        vector<Simplex*> selected_simplexes;

        // Sort simplexes by their diameter
        vector<Simplex*> sorted_partition = _partition;   // Note: Could sort globally, resorting would take less time
        sort(sorted_partition.begin(), sorted_partition.end(), Simplex::compare_diameter);
        double f_min = _func->_f_min;

        // Find simplex with  minimum metric  and  unique diameters
        Simplex* min_metric_simplex = sorted_partition[0];  // Initial value
        vector<double> diameters;
        vector<Simplex*> best_for_size;

        bool unique_diameter;
        bool found_with_same_size;
        for (int i=0; i < sorted_partition.size(); i++) {
            if (sorted_partition[i]->_metric__longest_edge_lb < min_metric_simplex->_metric__longest_edge_lb) {
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
                        if (best_for_size[j]->_longest_edge_lb_value > sorted_partition[i]->_longest_edge_lb_value) {
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
            // Is this OK?
            // if (iteration == 24) {cout << (best_for_size.size() > 2) << ", " << (min_metric_simplex != best_for_size[best_for_size.size()-1]);};
            if ((best_for_size.size() > 2) ) { // && (min_metric_simplex != best_for_size[best_for_size.size()-1])
                vector<Simplex*> simplexes_below_line;
                double a1 = best_for_size[0]->_diameter;
                double b1 = best_for_size[0]->_longest_edge_lb_value;
                // double a1 = min_metric_simplex->_diameter;  // Should be like this based on Direct Matlab implementation
                // double b1 = min_metric_simplex->_min_value;
                double a2 = best_for_size[best_for_size.size()-1]->_diameter;
                double b2 = best_for_size[best_for_size.size()-1]->_longest_edge_lb_value;
                
                double slope = (b2 - b1)/(a2 - a1);
                double bias = b1 - slope * a1;

                for (int i=0; i < best_for_size.size(); i++) {
                    if (best_for_size[i]->_longest_edge_lb_value < slope*best_for_size[i]->_diameter + bias +1e-12) {
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

            // Select all simplexes which have best _min_value for its size 
            for (int i=0; i < sorted_partition.size(); i++) {
                for (int j=0; j < selected.size(); j++) {
                    if ((sorted_partition[i]->_diameter == selected[j]->_diameter) && 
                        (sorted_partition[i]->_longest_edge_lb_value == selected[j]->_longest_edge_lb_value)) {
                        selected_simplexes.push_back(sorted_partition[i]);
                    };
                };
            };
        return selected_simplexes;
    };

    vector<Simplex*> select_simplexes_by_lowest_edge_lb(){};

    virtual vector<Simplex*> select_simplexes_to_divide(int iteration=0){
        vector<Simplex*> selected_simplexes;
        if (_lower_bound_strategy == MinVert) {
            selected_simplexes = select_simplexes_by_min_vert();
        };
        if (_lower_bound_strategy == LongestEdgeLB) {
            selected_simplexes = select_simplexes_by_longest_edge_lb();
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

        int iteration = 0;
        while (_func->_calls <= _max_calls && !_func->is_accurate_enougth()) { // _func->pe() > _min_pe){
            // Selects simplexes to divide
            vector<Simplex*> simplexes_to_divide;
            if (iteration == 0) {
                simplexes_to_divide = _partition;
            } else {
                simplexes_to_divide = select_simplexes_to_divide(iteration);
            };
            // Simplex::log_partition(_partition, simplexes_to_divide, "\nIteration ", iteration);
            // test_unique_simplexes();

            // Divides selected simplexes
            vector<Simplex*> new_simplexes;
            for (int i=0; i < simplexes_to_divide.size(); i++){
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

            // Update counters and log the status
            iteration += 1;
            // cout << iteration << ". Simplexes: " << _partition.size() << "  calls: " << _func->_calls << endl;

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
