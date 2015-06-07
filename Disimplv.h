#ifndef DISIMPLV_H
#define DISIMPLV_H_ 
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h> 
#include <limits>
#include "functions.h"

using namespace std;


class Simplex {
    Simplex(const Simplex& other){}
    Simplex& operator=(const Simplex& other){}
public:
    Simplex(){
        _is_in_partition = true;
        _parent = 0;
        _diameter = 0;
        _le_v1 = 0;
        _le_v2 = 0;
        _min_vert = 0;
        _min_value = numeric_limits<double>::max();
        _should_be_divided = false;
    };

    vector<Point*> _verts;
    bool _is_in_partition;
    bool _should_be_divided;  // Should be divided in next iteration
    Simplex* _parent;

    Point* _le_v1;      // Longest edge vertex1
    Point* _le_v2; 
    double _diameter;   // Longest edge length
    Point* _min_vert;   // Pointer to vertex with lowest function value 
    double _min_value;  // _min_vert function value 
    double _metric;     // _f_min - glob_f / _diameter
    
    void init_parameters(Function* func){   // Called when all verts have been added
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
            // Finds _min_vert 
            if (_min_value > _verts[a]->_values[0]) {
                _min_value = _verts[a]->_values[0];
                _min_vert = _verts[a];
            };
        }; 
        // epsilon = 1e-4
        double E;
        if (1e-4 * fabs(func->_glob_f) > 1e-8) {
            E = 1e-4 * fabs(func->_glob_f);
        } else {
            E = 1e-8;
        };
        _metric = (_min_value - (func->_glob_f + E)) / _diameter;

        //// Implement after Disimpl-v experiments are repeated
        // double _accurate_lower_bound;  // Find it by solving equation system. 
        // double _approx_lower_bound;  // 
        // double _approx_lower_bound;  // 
    };

    double l2norm(Point* p1, Point* p2) {
        double squared_sum = 0;
        for (int i=0; i < p1->size(); i++){
            squared_sum += pow(p1->_X[i] - p2->_X[i], 2);
        };
        return sqrt(squared_sum);
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
    //     return s1->_metric < s2->_metric;
    // };

    void add_vertex(Point* vertex){
        _verts.push_back(vertex);
    };

    int size() {
        return _verts.size();
    };

    void print(){
        cout << " Simplex  (" << _diameter << ", " << _min_value << "):" << endl;
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
           log_file << " ("<< simplexes[i]->_diameter << "," << simplexes[i]->_min_value << ")" << endl;
       };
       log_file << "Selected:" << endl;
       for (int i=0; i < selected.size(); i++) {
           for (int j=0; j < selected[i]->_verts.size(); j++) {
               for (int k=0; k < selected[i]->_verts[j]->size(); k++){
                    log_file << selected[i]->_verts[j]->_X[k] << " ";
               };
               log_file << " (" << selected[i]->_verts[j]->_values[0]<<"); ";
           };
           log_file << " ("<< selected[i]->_diameter << "," << selected[i]->_min_value << ")" << endl;
       };

       log_file.close();
    };

    virtual ~Simplex(){
        _verts.clear();
    };  
};
// typedef vector<Simplex*> Simplexes;


class Disimplv {
    Disimplv(const Disimplv& other){};
    Disimplv& operator=(const Disimplv& other){};
public:
    // Disimplv(){};
    Disimplv(double min_pe, int max_calls){
        _min_pe = min_pe;
        _max_calls = max_calls;
        ofstream log_file; 
        log_file.open("log/partition.txt");
        log_file.close();
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
            simpl->init_parameters(_func);
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
            c = nextv(nextv(v, m), m);   // d = x = _diameter;  f = y = _min_value;

            double* matrix[3];
            double line1[3] = {simplexes[a]->_diameter, simplexes[a]->_min_value, 1.};
            double line2[3] = {simplexes[b]->_diameter, simplexes[b]->_min_value, 1.};
            double line3[3] = {simplexes[c]->_diameter, simplexes[c]->_min_value, 1.};
            matrix[0] = line1;
            matrix[1] = line2;
            matrix[2] = line3;
            // (double**) 
            // double matrix[3][3] = {{simplexes[a]->_diameter, simplexes[a]->_min_value, 1.},
            //                      {simplexes[b]->_diameter, simplexes[b]->_min_value, 1.},
            //                      {simplexes[c]->_diameter, simplexes[c]->_min_value, 1.}};
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

    vector<Simplex*> select_simplexes_to_divide(int iteration=0, string strategy="min_vert"){

        vector<Simplex*> selected_simplexes;
        if (strategy == "min_vert") {
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
                if (sorted_partition[i]->_metric < min_metric_simplex->_metric) {
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
                        if (best_for_size[j]->_min_value > sorted_partition[i]->_min_value) {
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
                double b1 = best_for_size[0]->_min_value;
                // double a1 = min_metric_simplex->_diameter;  // Should be like this based on Direct Matlab implementation
                // double b1 = min_metric_simplex->_min_value;
                double a2 = best_for_size[best_for_size.size()-1]->_diameter;
                double b2 = best_for_size[best_for_size.size()-1]->_min_value;
                
                double slope = (b2 - b1)/(a2 - a1);
                double bias = b1 - slope * a1;

                for (int i=0; i < best_for_size.size(); i++) {
                    if (best_for_size[i]->_min_value < slope*best_for_size[i]->_diameter + bias +1e-12) {
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

            // Error: wrong stopping condition
            // Remove simplexes which do not satisfy condition:   f - slope*d > f_min - epsilon*abs(f_min)
            for (int i=0; i < selected.size() -1; i++) {
                double a1 = selected[selected.size() - i -1]->_diameter;
                double b1 = selected[selected.size() - i -1]->_min_value;
                double a2 = selected[selected.size() - i -2]->_diameter;
                double b2 = selected[selected.size() - i -2]->_min_value;
                double slope = (b2 - double(b1))/(a2 - a1);
                double bias = b1 - slope * a1;

                if (bias > f_min - 0.01*fabs(f_min)) {
                    selected[selected.size() - i -2]->_should_be_divided = false;
                };
            };

            // Remove simplexes which should not be divided
            selected.erase(remove_if(selected.begin(), selected.end(), Simplex::wont_be_divided), selected.end());

            // Select all simplexes which have best _min_value for its size 
            // vector<Simplex*> selected_simplexes;
            for (int i=0; i < sorted_partition.size(); i++) {
                for (int j=0; j < selected.size(); j++) {
                    if ((sorted_partition[i]->_diameter == selected[j]->_diameter) && 
                        (sorted_partition[i]->_min_value == selected[j]->_min_value)) {
                        selected_simplexes.push_back(sorted_partition[i]);
                    };
                };
            };


        };
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
            // Find middle point
            int n = _func->_D;
            double c[n];
            for (int i=0; i < n; i++) {
                c[i] = (simplex->_le_v1->_X[i] + simplex->_le_v2->_X[i]) / 2.;
            };
            Point* middle_point = _func->get(c, n);

            // Construct two new simplexes using this middle point.
            Simplex* left_simplex = new Simplex();
            Simplex* right_simplex = new Simplex();
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
            left_simplex->init_parameters(_func);
            right_simplex->init_parameters(_func);
            left_simplex->_parent = simplex;
            right_simplex->_parent = simplex;
            simplex->_is_in_partition = false;

            divided_simplexes.push_back(left_simplex);
            divided_simplexes.push_back(right_simplex);
            return divided_simplexes;
        };
    };

    void minimize(Function* func){
        _func = func;
        partition_feasable_region();

        int iteration = 0;
        while (_func->_calls <= _max_calls && !_func->is_accurate_enougth()) { // _func->pe() > _min_pe){
            // Selects simplexes to divide
            vector<Simplex*> simplexes_to_divide;
            if (iteration == 0) {
                simplexes_to_divide = _partition;
            } else {
                simplexes_to_divide = select_simplexes_to_divide(iteration);
            };
            Simplex::log_partition(_partition, simplexes_to_divide, "\nIteration ", iteration);

            // test_unique_simplexes();
            // Simplex::print(simplexes_to_divide, "\n  Simplexes to divide: >> ");
            // cout << "<<<<<\n " << endl;

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
            cout << iteration << ". Simplexes: " << _partition.size() << "  calls: " << _func->_calls << endl;

            // if (iteration >= 9) {
            //     break;
            // };
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
