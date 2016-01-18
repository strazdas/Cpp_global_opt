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
#include "utils.h"
#include "simplex.h"
#include "functions.h"
#include "algorithm.h"

using namespace std;


class Algorithm {
    Algorithm(const Algorithm& other){};
    Algorithm& operator=(const Algorithm& other){};
public:
    Algorithm(){};
    string _name;
    string _status;
    int _max_calls;
    string _stop_criteria;
    double _min_pe;
    double _duration = 0;   // Duration in seconds
    double _max_duration;   // Maximum allowed duration of minimization in seconds
    double _epsilon;    // How small simplexes should still be partitioned 
    vector<Simplex*> _partition;
    vector<Simplex*> _all_simplexes;
    vector<Function*> _funcs;
    vector<Point*> _pareto_front;
 
    LowerBoundStrategy _lower_bound_strategy;   // Simplex lower bound estimation strategy 
    LStrategy _L_strategy;                      // Lipschitz constant estimation strategy
    DivisionStrategy _division_strategy;        // Simplex edge division strategy
    SimplexGradientStrategy _simplex_gradient_strategy;
    double _parent_L_part;


    /* Global optimization strategies */
    // void partition_feasable_region_combinatoricly(){
    //     int n = _funcs[0]->_D;
    //     int number_of_simpleces = 1;
    //     for (int i = 1; i <= n; i++) {
    //         number_of_simpleces *= i;
    //     };
    //
    //     int teta[n];
    //     for (int i=0; i < n; i++){
    //         teta[i] = i;
    //     };
    //
    //     do {
    //         int triangle[n+1][n];
    //
    //         for (int k = 0; k < n; k++) {
    //             triangle[0][k] = 0;
    //         };
    //
    //         for (int vertex=0; vertex < n; vertex++) {
    //             for (int j = 0; j < n + 1; j++) {
    //                 triangle[vertex + 1][j] = triangle[vertex][j];
    //             };
    //             triangle[vertex + 1][teta[vertex]] = 1;
    //         }
    //
    //         Simplex* simpl = new Simplex(_lower_bound_strategy, _L_strategy, _parent_L_part, _simplex_gradient_strategy);
    //         for (int i=0; i < n + 1; i++){
    //             Point* tmp_point = new Point(triangle[i], n);
    //             Point* point = _funcs[0]->get(tmp_point); 
    //             if (tmp_point != point) {
    //                 delete tmp_point;
    //             };
    //             simpl->add_vertex(point);
    //         };
    //         simpl->init_parameters(_funcs);
    //         _partition.push_back(simpl);
    //         _all_simplexes.push_back(simpl);
    //
    //     } while (next_permutation(teta, teta+n));
    // };

    vector<Simplex*> divide_simplex(Simplex* simplex, string strategy="longest_half") {
        vector<Simplex*> divided_simplexes;
        if (strategy== "longest_half") {
            // Find middle point
            int n = _funcs[0]->_D;
            double c[n];
            for (int i=0; i < n; i++) {
                c[i] = (simplex->_le_v1->_X[i] + simplex->_le_v2->_X[i]) / 2.;
            };
            Point* middle_point = _funcs[0]->get(c, n);
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
                simplex->_verts[i]->_neighbours_estimates_should_be_updated();
            };
            middle_point->_neighbours_estimates_should_be_updated();

            left_simplex->_parent = simplex;
            right_simplex->_parent = simplex;
            left_simplex->init_parameters(_funcs);
            right_simplex->init_parameters(_funcs);
            simplex->_is_in_partition = false;

            divided_simplexes.push_back(left_simplex);
            divided_simplexes.push_back(right_simplex);
            return divided_simplexes;
        };
    };

    friend ostream& operator<<(ostream& o, const Algorithm& a){
        o << "alg: '" << a._name << "', func: '" << a._funcs[0]->_name << 
           "', calls: " << a._funcs[0]->_calls << ", subregions: " << a._partition.size() <<
           ", duration: " << a._duration << ", stop_criteria: '" << a._stop_criteria << 
           "', f_min: " << a._funcs[0]->_f_min << ", x_min: " << (*a._funcs[0]->_x_min);
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

#endif
