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

using namespace std;


/* Class definitions */
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
                simplex->_verts[i]->_neighbours_estimates_should_be_updated();
            };
            middle_point->_neighbours_estimates_should_be_updated();

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
        //     if (simplexes[sid]->_min_lb_value == min_lower_bound_value) {
        //         selected_simplexes.push_back(simplexes[sid]);
        //     } else {
        //         if (simplexes[sid]->_min_lb_value < min_lower_bound_value) {
        //             selected_simplexes.clear(); // Test if clear works properly 
        //             min_lower_bound_value = simplexes[sid]->_min_lb_value;
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
                        if (best_for_size[j]->_min_lb_value > sorted_partition[i]->_min_lb_value) {
                            best_for_size.erase(best_for_size.begin()+j);
                            best_for_size.push_back(sorted_partition[i]);
                        };
                    };
                };
                if (!found_with_same_size) {
                    best_for_size.push_back(sorted_partition[i]);
                };

                if (sorted_partition[i]->_min_lb_value < best_for_size[min_value_simplex_id]->_min_lb_value) {
                    // if (iteration == 9) { cout << "Greater" << endl; };
                    min_value_simplex_id = best_for_size.size() - 1;
                };
            // if (iteration == 9) { cout << "Comparing: (" << i<<")"<< sorted_partition[i]->_min_lb_value << " and ("<< min_value_simplex_id << ")" << best_for_size[min_value_simplex_id]->_min_lb_value << endl; };
            };



            vector<Simplex*> selected;
            // Select from minimum value simplex to bigest simplex using convex hull
            // Get minimum simplex id in best_for_size:  min_value_simplex_id
            // Find out how minimum value simplex id in best_for_size, not in sorted_partition 
            if ((best_for_size.size() - min_value_simplex_id) > 2) {
                vector<Simplex*> simplexes_below_line;
                // double a1 = best_for_size[0]->_diameter;
                // double b1 = best_for_size[0]->_min_lb_value;
                double a1 = best_for_size[min_value_simplex_id]->_diameter;  // Should be like this based on Direct Matlab implementation
                double b1 = best_for_size[min_value_simplex_id]->_min_lb_value;

                double a2 = best_for_size[best_for_size.size()-1]->_diameter;
                double b2 = best_for_size[best_for_size.size()-1]->_min_lb_value;
                
                double slope = (b2 - b1)/(a2 - a1);
                double bias = b1 - slope * a1;

                for (int i=min_value_simplex_id; i < best_for_size.size(); i++) {
                    if (best_for_size[i]->_min_lb_value < slope*best_for_size[i]->_diameter + bias +1e-12) {
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
                //      cout << "Values: " << best_for_size[0]->_min_lb_value << ", " << best_for_size[1]->_min_lb_value << endl;
                //      exit(0);
                // };
            };

            


            // // if (iteration == 24) {cout << (best_for_size.size() > 2) << ", " << (min_metric_simplex != best_for_size[best_for_size.size()-1]);};
            // if ((best_for_size.size() > 2) ) { // && (min_metric_simplex != best_for_size[best_for_size.size()-1])
            //     vector<Simplex*> simplexes_below_line;
            //     double a1 = best_for_size[0]->_diameter;
            //     double b1 = best_for_size[0]->_min_lb_value;
            //     // double a1 = min_metric_simplex->_diameter;  // Should be like this based on Direct Matlab implementation
            //     // double b1 = min_metric_simplex->_min_value;
            //     double a2 = best_for_size[best_for_size.size()-1]->_diameter;
            //     double b2 = best_for_size[best_for_size.size()-1]->_min_lb_value;
            //
            //     double slope = (b2 - b1)/(a2 - a1);
            //     double bias = b1 - slope * a1;
            //
            //     for (int i=0; i < best_for_size.size(); i++) {
            //         if (best_for_size[i]->_min_lb_value < slope*best_for_size[i]->_diameter + bias +1e-12) {
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
                double b1 = selected[selected.size() - i -1]->_min_lb_value;
                double a2 = selected[selected.size() - i -2]->_diameter;
                double b2 = selected[selected.size() - i -2]->_min_lb_value;
                double slope = (b2 - double(b1))/(a2 - a1);
                double bias = b1 - slope * a1;

                if (bias > f_min - _epsilon * fabs(f_min)) {   // epsilon
                    selected[selected.size() - i -2]->_should_be_divided = false;
                };
            };

            // Remove simplexes which should not be divided
            selected.erase(remove_if(selected.begin(), selected.end(), Simplex::wont_be_divided), selected.end());

            // Select all simplexes which have best _min_lb_value for its size 
            for (int i=0; i < sorted_partition.size(); i++) {
                for (int j=0; j < selected.size(); j++) {
                    if ((sorted_partition[i]->_diameter == selected[j]->_diameter) && 
                        (sorted_partition[i]->_min_lb_value == selected[j]->_min_lb_value)) {
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
        //     if (sorted_partition[i]->_metric__min_lb < min_metric_simplex->_metric__min_lb) {
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
        //                 if (best_for_size[j]->_min_lb_value > sorted_partition[i]->_min_lb_value) {
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
        //         double b1 = best_for_size[0]->_min_lb_value;
        //         // double a1 = min_metric_simplex->_diameter;  // Should be like this based on Direct Matlab implementation
        //         // double b1 = min_metric_simplex->_min_value;
        //         double a2 = best_for_size[best_for_size.size()-1]->_diameter;
        //         double b2 = best_for_size[best_for_size.size()-1]->_min_lb_value;
        //         
        //         double slope = (b2 - b1)/(a2 - a1);
        //         double bias = b1 - slope * a1;
        //
        //         for (int i=0; i < best_for_size.size(); i++) {
        //             if (best_for_size[i]->_min_lb_value < slope*best_for_size[i]->_diameter + bias +1e-12) {
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
        //         double b1 = selected[selected.size() - i -1]->_min_lb_value;
        //         double a2 = selected[selected.size() - i -2]->_diameter;
        //         double b2 = selected[selected.size() - i -2]->_min_lb_value;
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
        //                 (sorted_partition[i]->_min_lb_value == selected[j]->_min_lb_value)) {
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
        // cout << "Searching estimates for simplexes: " <<  _partition.size() << endl;
        if (_L_strategy == Neighbours) { Simplex::update_estimates(_partition, _func); };

        int iteration = 0;
        while (_func->_calls <= _max_calls && !_func->is_accurate_enougth()) { // _func->pe() > _min_pe){
            // cout << "Iteration: " << iteration << endl;

            // Selects simplexes to divide
            vector<Simplex*> simplexes_to_divide;
            if (iteration == 0) {
                simplexes_to_divide = _partition;
            } else {
                // Note: sort partition here
                simplexes_to_divide = select_simplexes_to_divide(iteration);
            };
            // Simplex::log_partition(_partition, simplexes_to_divide, "\nIteration ", iteration);
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
            // cout << "Searching estimates for simplexes: " <<  _partition.size() << endl;
            if (_L_strategy == Neighbours) { Simplex::update_estimates(_partition, _func); };

            // Update counters and log the status
            iteration += 1;
            cout << iteration << ". Simplexes: " << _partition.size() << "  calls: " << _func->_calls << endl;

        };
        timestamp_t end = get_timestamp();
        _duration = (end - start) / 1000000.0L;

        // Draw partitioning: output simplex coordinates to file and draw it with Python
    };

    virtual ~Disimplv(){};
};

#endif
