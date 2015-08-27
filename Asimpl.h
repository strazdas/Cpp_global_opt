#ifndef ASIMPL_H
#define ASIMPL_H 
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <math.h> 
#include <string>
#include <sys/time.h>
#include <vector>
// #include <forward_list>
#include "utils.h"
#include "simplex.h"
#include "functions.h"
#include "Disimplv.h"

using namespace std;


class Asimpl : public Algorithm {
    Asimpl(const Asimpl& other) {};
    Asimpl& operator=(const Asimpl& other) {};
public:
    Asimpl(double epsilon=0.0001, int max_calls=15000, double max_duration=3600) {
        _lower_bound_strategy = LowestEdgeLB;    // Lowest edge is determined by optimising
        _L_strategy = Neighbours;                // Simplex region to get max L from 
        _division_strategy = LongestHalf;        // Simplex division strategy - longest into two parts
        _simplex_gradient_strategy = FFMinVert;  // Single simplex L determination strategy (grad norm) 
        _stop_criteria = "x_dist_Serg";          // Stopping criteria

        _epsilon = epsilon;                      // Solution accuracy
        _max_calls = max_calls;
        _max_duration = max_duration;

        // Construct algorithm name
        _name = "Asimpl";
        stringstream alg_name; 
        alg_name << _name << "_e" << epsilon;  
        _name =  alg_name.str();

        // Clean partition log file
        ofstream log_file; 
        log_file.open("log/partition.txt");
        log_file.close();
    };
    // forward_list<Simplex*> _simpls;

    // LongestEdgeLB, Neighbours
    vector<Simplex*> select_simplexes_by_lowest_edge_lb() {
        vector<Simplex*> selected_simplexes;
        selected_simplexes.push_back(_partition[0]);  // _partition should be sorted ascending by lb min
        return selected_simplexes;                    
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
            c = nextv(nextv(v, m), m);   // d = x = _diameter;  f = y = _min_lb_value;

            double* matrix[3];
            double line1[3] = {simplexes[a]->_diameter, simplexes[a]->_min_lb_value, 1.};
            double line2[3] = {simplexes[b]->_diameter, simplexes[b]->_min_lb_value, 1.};
            double line3[3] = {simplexes[c]->_diameter, simplexes[c]->_min_lb_value, 1.};
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

    vector<Simplex*> select_simplexes_by_lowest_edge_lb_and_diameter_convex_hull() {
        vector<Simplex*> selected_simplexes;

        // Sort simplexes by their diameter
        vector<Simplex*> sorted_partition = _partition;   // Note: Could sort globally, resorting would take less time

        // Simplex::print(sorted_partition, "Selecting for division from: ");
        // sort(sorted_partition.begin(), sorted_partition.end(), Simplex::compare_diameter);
        double f_min = _func->_f_min;

        // Find simplex with  minimum metric  and  unique diameters
        Simplex* min_metric_simplex = sorted_partition[0]; // Initial value
        vector<double> diameters;
        vector<Simplex*> best_for_size;

        bool unique_diameter;
        bool found_with_same_size;
        for (int i=0; i < sorted_partition.size(); i++) {
            if (sorted_partition[i]->_metric__min_lb < min_metric_simplex->_metric__min_lb) {
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
                    if (best_for_size[j]->_min_lb_value > sorted_partition[i]->_min_lb_value) {
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
            double b1 = best_for_size[0]->_min_lb_value;
            // double a1 = min_metric_simplex->_diameter;  // Should be like this based on Direct Matlab implementation
            // double b1 = min_metric_simplex->_min_lb_value;
            double a2 = best_for_size[best_for_size.size()-1]->_diameter;
            double b2 = best_for_size[best_for_size.size()-1]->_min_lb_value;

            double slope = (b2 - b1)/(a2 - a1);
            double bias = b1 - slope * a1;

            for (int i=0; i < best_for_size.size(); i++) {
                if (best_for_size[i]->_min_lb_value < slope*best_for_size[i]->_diameter + bias +1e-12) {
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
            double b1 = selected[selected.size() - i -1]->_min_lb_value;
            double a2 = selected[selected.size() - i -2]->_diameter;
            double b2 = selected[selected.size() - i -2]->_min_lb_value;
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
                    (sorted_partition[i]->_min_lb_value == selected[j]->_min_lb_value)) {
                    selected_simplexes.push_back(sorted_partition[i]);
                };
            };
        };

        return selected_simplexes;
    };


    virtual vector<Simplex*> select_simplexes_to_divide() {
        // vector<Simplex*> selected_simplexes = select_simplexes_by_lowest_edge_lb();
        vector<Simplex*> selected_simplexes = select_simplexes_by_lowest_edge_lb_and_diameter_convex_hull();

        for (int i=0; i < selected_simplexes.size(); i++){
            selected_simplexes[i]->_is_in_partition = false;
        };
        return selected_simplexes;
    };

    void minimize(Function* func){
        _func = func;
        timestamp_t start = get_timestamp();
        partition_feasable_region_combinatoricly();     // Note: Should not use global variables
        Simplex::update_estimates(_partition, _func);
        // sort(_partition.begin(), _partition.end(), Simplex::ascending_min_lb_value);
        sort(_partition.begin(), _partition.end(), Simplex::ascending_diameter);

        int iteration = 0;
        while (_func->_calls <= _max_calls && _duration <= _max_duration && !_func->is_accurate_enougth()) { // _func->pe() > _min_pe){
            // Selects simplexes to divide
            vector<Simplex*> simplexes_to_divide;
            if (iteration == 0) {
                simplexes_to_divide = _partition;
            } else {
                simplexes_to_divide = select_simplexes_to_divide();
            };

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
            // Delete simplexes?

            // Add new simplexes to _partition and _all_simplexes
            for (int i=0; i < new_simplexes.size(); i++) {
                _partition.push_back(new_simplexes[i]);
                _all_simplexes.push_back(new_simplexes[i]);
            };
            // cout << "Searching estimates for simplexes: " <<  _partition.size() << endl;
            Simplex::update_estimates(_partition, _func);
            // sort(_partition.begin(), _partition.end(), Simplex::ascending_min_lb_value);
            sort(_partition.begin(), _partition.end(), Simplex::ascending_diameter);

            // Update counters and log the status
            iteration += 1;
            // cout << iteration << ". Simplexes: " << _partition.size() << "  calls: " << _func->_calls << endl;
            timestamp_t end = get_timestamp();
            _duration = (end - start) / 1000000.0L;

        };

        if ((_func->_calls <= _max_calls) && (_duration <= _max_duration)) {
            _status = "D"; 
        } else {
            _status = "S";
        };

        // Draw partitioning: output simplex coordinates to file and draw it with Python
    };

    virtual ~Asimpl(){};
};

#endif
