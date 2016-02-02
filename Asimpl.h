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
#include "utils.h"
#include "simplex.h"
#include "functions.h"
#include "algorithm.h"

using namespace std;


class Asimpl : public Algorithm {
    Asimpl(const Asimpl& other) {};
    Asimpl& operator=(const Asimpl& other) {};
public:
    Asimpl(int max_calls=15000, double max_duration=3600, double epsilon=0.0001) {
        _lower_bound_strategy = LowestEdgeLB;    // Lowest edge is determined by optimising
        _L_strategy = Neighbours;                // Simplex region to get max L from 
        // _max_diff_verts_to_be_neighbour = 2;     // Max number of different verts to still be a neighbour
        _min_same_verts_to_be_neighbour = 1;     // Max number of different verts to still be a neighbour
        _division_strategy = LongestHalf;        // Simplex division strategy - longest into two parts
        _simplex_gradient_strategy = FFMinVert;  // Single simplex L determination strategy (grad norm) 
        _stop_criteria = "x_dist_Serg";          // Stopping criteria

        _epsilon = epsilon;                      // Solution accuracy
        _max_calls = max_calls;
        _max_duration = max_duration;

        _iteration = 0;         // Algorithm iteration for debuging purposes

        // Construct algorithm name
        _name = "Asimpl";
        stringstream alg_name; 
        alg_name << _name << "_e" << epsilon;  
        _name =  alg_name.str();

        // Clean partition log file
        ofstream log_file; 
        log_file.open("log/partition.txt");
        log_file.close();
        log_file.open("log/front.txt");
        log_file.close();
    };

    int _max_diff_verts_to_be_neighbour;
    int _min_same_verts_to_be_neighbour;
    int _iteration;

    void partition_feasable_region_combinatoricly() {
        int n = _funcs[0]->_D;
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

            Simplex* simpl = new Simplex(_lower_bound_strategy, _L_strategy, _simplex_gradient_strategy);
            for (int i=0; i < n + 1; i++){
                Point* tmp_point = new Point(triangle[i], n);
                
                Point* point = _funcs[0]->get(tmp_point); 
                if (tmp_point != point) {
                    delete tmp_point;
                } else {
                    for (int j=1; j < _funcs.size(); j++) {
                        _funcs[j]->update_meta(point);     // Note: update_meta - is not very verbose name
                    };                                     // update would be better name
                    update_pareto_front(point);
                };
                simpl->add_vertex(point);
            };
            simpl->init_parameters(_funcs);
            _partition.push_back(simpl);

        } while (next_permutation(teta, teta+n));

        for (int i=0; i < _partition.size()-1; i++) {
            for (int j=i+1; j < _partition.size(); j++) {
                if (are_neighbours(_partition[i], _partition[j])) {
                    _partition[i]->_neighbours.push_back(_partition[j]);
                    _partition[j]->_neighbours.push_back(_partition[i]);
                };
                // Find out if they are neighbours if yes add each other to each Simplex::neighbours list
            };
        };
    };

    void print_neighbours() {  // helper function
        cout << "Total: " << _partition.size() << endl;    
        for (int sid=0; sid < _partition.size(); sid++) {
            _partition[sid]->print(); 
            cout << " Simplex: ";
            cout << _partition[sid];
            cout << "("; 
            for (int i=0; i < _partition[sid]->_Ls.size(); i++) {
                cout << _partition[sid]->_Ls[i];
                if (i != _partition[sid]->_Ls.size() - 1) {
                    cout << ",";
                };
            };
            cout << ")";
            cout << "      Neighbours are:";
            cout << "(total: " << _partition[sid]->_neighbours.size() << ")";

            // for (int nid=0; nid < _partition[sid]->_neighbours.size(); nid++) {
            Simplex* neighbour;
            for (list<Simplex*>::iterator it=_partition[sid]->_neighbours.begin(); it != _partition[sid]->_neighbours.end(); ++it) {
                neighbour = *it;
                cout << ' ' << (*it);
            };
            cout << endl << endl;
        };
    };

    bool are_neighbours(Simplex* s1, Simplex* s2) {
        // vector<Point*> same_verts(s1->_verts.size() + s2->_verts.size());
        // vector<Point*>::iterator it;
        //
        // vector<Point*> s1_verts = s1->_verts;
        // vector<Point*> s2_verts = s2->_verts;
        //
        //
        // sort(s1_verts.begin(), s1_verts.end());
        // sort(s2_verts.begin(), s2_verts.end());
        // it = set_intersection(s1_verts.begin(), s1_verts.end(), s2_verts.begin(), s2_verts.end(), same_verts.begin());
        //
        // int same_verts_count = it - same_verts.begin();
        // int diff_verts_count = s1_verts.size() - same_verts_count;
        //
        // same_verts.clear();
        // s1_verts.clear();
        // s2_verts.clear();
        //
        // if (same_verts_count >= _min_same_verts_to_be_neighbour) {
        //     return true;
        // };
 
        //// if (diff_verts_count <= _max_diff_verts_to_be_neighbour) {
        ////     return true;
        //// };
        return false;
    };

    // LongestEdgeLB, Neighbours    // Worning: this method could be outdated, _partition could not be sorted
    // vector<Simplex*> select_simplexes_by_lowest_edge_lb() {
    //     vector<Simplex*> selected_simplexes;
    //     selected_simplexes.push_back(_partition[0]);  // _partition should be sorted ascending by lb min
    //     return selected_simplexes;                    
    // };

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
            c = nextv(nextv(v, m), m);   // d = x = _diameter;  f = y = _tolerance;

            double* matrix[3];
            double line1[3] = {simplexes[a]->_diameter, simplexes[a]->_tolerance, 1.};
            double line2[3] = {simplexes[b]->_diameter, simplexes[b]->_tolerance, 1.};
            double line3[3] = {simplexes[c]->_diameter, simplexes[c]->_tolerance, 1.};
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

    vector<Simplex*> select_simplexes_by_lb_estimate_and_diameter_convex_hull() {   // Selects from best
        vector<Simplex*> selected_simplexes;

        // Sort simplexes by their diameter
        vector<Simplex*> sorted_partition = _partition;   // Note: Could sort globally, resorting would take less time

        // Simplex::print(sorted_partition, "Selecting for division from: ");
        sort(sorted_partition.begin(), sorted_partition.end(), Simplex::ascending_diameter);
        double f_min = _funcs[0]->_f_min;

        // Find simplex with  minimum metric  and  unique diameters
        Simplex* min_metric_simplex = sorted_partition[0];  // Initial value
        vector<double> diameters;
        vector<Simplex*> best_for_size;

        bool unique_diameter;
        bool found_with_same_size;
        for (int i=0; i < sorted_partition.size(); i++) {
            if (sorted_partition[i]->_tolerance < min_metric_simplex->_tolerance) {
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
                    if (best_for_size[j]->_tolerance > sorted_partition[i]->_tolerance) {
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
        if (min_metric_simplex == best_for_size[best_for_size.size()-1]) {
            selected.push_back(min_metric_simplex);
        } else {
            // Is this OK?  Well compared with examples - its ok.
            if ((best_for_size.size() > 2) && (min_metric_simplex != best_for_size[best_for_size.size()-1])) {
                vector<Simplex*> simplexes_below_line;
                // double a1 = best_for_size[0]->_diameter;
                // double b1 = best_for_size[0]->_tolerance;
                double a1 = min_metric_simplex->_diameter;  // Should be like this based on Direct Matlab implementation
                double b1 = min_metric_simplex->_tolerance;
                double a2 = best_for_size[best_for_size.size()-1]->_diameter;
                double b2 = best_for_size[best_for_size.size()-1]->_tolerance;

                double slope = (b2 - b1)/(a2 - a1);
                double bias = b1 - slope * a1;

                for (int i=0; i < best_for_size.size(); i++) {
                    if (best_for_size[i]->_tolerance < slope*best_for_size[i]->_diameter + bias +1e-12) {
                        simplexes_below_line.push_back(best_for_size[i]);
                    };
                };
                selected = convex_hull(simplexes_below_line);  // Messes up simplexes_below_line
            } else {
                selected = best_for_size;    // TODO: Why we divide all of them? Could divide only min_metrc_simplex.
                                             // Because practiacally this case does not occur ever.
            };
        };

        for (int i=0; i < selected.size(); i++) {
            selected[i]->_should_be_divided = true;
        };

        //// Should check all criterias, not only first
        // This part is very irational, because tolerance line and f_min point
        // are compared.
        //
        // Should f1, f2 lines be constructed, but in this case lbm should
        // known.
        //
        // Should use min_lbs values instead of tolerance.
        //
        // 




        // Remove simplexes which do not satisfy condition:   f - slope*d > f_min - epsilon*abs(f_min)
        for (int i=0; i < (signed) selected.size() -1; i++) {  // I gess error here - bias is incorrect
            //// Version2: All functions should be improved
            // int improvable = 0;
            // for (int j; j < _funcs.size(); j++) {
            //     if (selected[i]->_min_lbs[j]->_values[0] < _funcs[j]->_f_min - 0.0001 * fabs(_funcs[j]->_f_min)) {
            //         improvable += 1;
            //     };
            // };
            // if (improvable == 0) {
            //     selected[i]->_should_be_divided = false;
            // };

            //// Version3: Too small simplexes should not be divided (but tolerance should ensure this)
            // if (selected[i]->_diameter < 1e-5) {
            //     selected[i]->_should_be_divided = false;
            // };

            //// Version4: tolerance should ensure this.

            //// VersionOriginal: for single criteria
            // double a1 = selected[selected.size() - i -1]->_diameter;
            // double b1 = selected[selected.size() - i -1]->_tolerance;
            // double a2 = selected[selected.size() - i -2]->_diameter;
            // double b2 = selected[selected.size() - i -2]->_tolerance;
            // double slope = (b2 - double(b1))/(a2 - a1);
            // double bias = b1 - slope * a1;
            // // cout << "slope = (b2 - double(b1))/(a2 - a1);   bias = b1 - slope * a1;" << endl;
            // // cout << "b-remove if    bias > f_min - 0.0001*fabs(f_min)" << endl;
            // // cout << "a1 " << a1 << " b1 " << b1 << " a2 " << a2 << " b2 " << b2 << endl;
            // // cout << "slope " << slope << " bias " << bias << endl;
            // // cout << "bias: " << bias << " fmin " << f_min << endl;
            // // cout << endl;
            //
            // if (bias > f_min - 0.0001*fabs(f_min)) {   // epsilon
            //     selected[selected.size() - i -2]->_should_be_divided = false;
            // };
        };



        // Remove simplexes which should not be divided
        selected.erase(remove_if(selected.begin(), selected.end(), Simplex::wont_be_divided), selected.end());

        // Select all simplexes which have best _tolerance for its size 
        for (int i=0; i < sorted_partition.size(); i++) {
            for (int j=0; j < selected.size(); j++) {
                if ((sorted_partition[i]->_diameter == selected[j]->_diameter) && 
                    (sorted_partition[i]->_tolerance == selected[j]->_tolerance)) {
                    selected_simplexes.push_back(sorted_partition[i]);
                };
            };
        };

        return selected_simplexes;
    };


    virtual vector<Simplex*> select_simplexes_to_divide() {
        vector<Simplex*> selected_simplexes = select_simplexes_by_lb_estimate_and_diameter_convex_hull();
        return selected_simplexes;
    };

    vector<Simplex*> divide_simplex(Simplex* simplex, string strategy="longest_half") {
        vector<Simplex*> divided_simplexes;
        if (strategy== "longest_half") {
            // Find middle point

            // if (_iteration == 199) {
            //     cout << "It" << _iteration << endl;
            //     simplex->_le_v1->print();
            //     simplex->_le_v2->print();
            //     simplex->print();
            //     cout << " ---- " << endl;
            // };
                
            int n = _funcs[0]->_D;
            double c[n];
            for (int i=0; i < n; i++) {
                c[i] = (simplex->_le_v1->_X[i] + simplex->_le_v2->_X[i]) / 2.;
            };

            // Construct middle point
            Point* tmp_point = new Point(c, n);
            Point* middle_point = _funcs[0]->get(tmp_point); 
            if (tmp_point != middle_point) {
                delete tmp_point;
            } else {
                for (int j=1; j < _funcs.size(); j++) {
                    _funcs[j]->update_meta(middle_point); // Note: update_meta - is not very verbose name
                };                                        // update would be better name
                update_pareto_front(middle_point);
            };

            // Construct two new simplexes using this middle point.
            Simplex* left_simplex = new Simplex(_lower_bound_strategy, _L_strategy, _simplex_gradient_strategy);
            Simplex* right_simplex = new Simplex(_lower_bound_strategy, _L_strategy, _simplex_gradient_strategy);

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
            middle_point->_neighbours_estimates_should_be_updated();   // Note: this method should also be updated to take into account Algorithm._max_diff_verts_to_be_neighbour

            left_simplex->init_parameters(_funcs);
            right_simplex->init_parameters(_funcs);

            // Update new simplexes neighbours based on their parent neighbours
            // Divided simplexes are neighbours.
            left_simplex->_neighbours.push_back(right_simplex);
            right_simplex->_neighbours.push_back(left_simplex);

            Simplex* neighbour;
            // Iterate through parents neighbours and check if any of them are  child neighbours
            for (list<Simplex*>::iterator it=simplex->_neighbours.begin(); it != simplex->_neighbours.end(); ++it) {
                neighbour = *it;
                neighbour->_neighbours.remove(simplex);
                if (are_neighbours(left_simplex, neighbour)) {
                    left_simplex->_neighbours.push_back(neighbour);
                    neighbour->_neighbours.push_back(left_simplex);
                };
                if (are_neighbours(right_simplex, neighbour)) {
                    right_simplex->_neighbours.push_back(neighbour);
                    neighbour->_neighbours.push_back(right_simplex);
                };
            };

            simplex->_is_in_partition = false;

            //     cout << "It" << _iteration << endl;
            //     simplex->_le_v1->print();
            //     simplex->_le_v2->print();
            //     simplex->print();
            //     cout << " -- " << endl;
            //     left_simplex->print();
            //     right_simplex->print();
            //     cout << " ---- " << endl;
            // };

            divided_simplexes.push_back(left_simplex);
            divided_simplexes.push_back(right_simplex);
            return divided_simplexes;
        };
    };

    vector<Simplex*> divide_simplexes(vector<Simplex*> simplexes) {
        // Divides selected simplexes
        vector<Simplex*> new_simplexes;
        for (int i=0; i < simplexes.size(); i++) {
            vector<Simplex*> divided_simplexes = divide_simplex(simplexes[i]);

            for (int j=0; j < divided_simplexes.size(); j++) {
                new_simplexes.push_back(divided_simplexes[j]);
            };
        };

        // Update neighbours for new_simplexes: remove neighbours, who are not in partition
        for (int i=0; i < new_simplexes.size(); i++) {
            list<Simplex*> neighbours = new_simplexes[i]->_neighbours;
            neighbours.erase(remove_if(neighbours.begin(), neighbours.end(), Simplex::not_in_partition), neighbours.end());
        };

        // Update neighbours for each new_simplex:  add newly created neighbours
        list<Simplex*> neighbours;
        for (int i=0; i < new_simplexes.size(); i++) {
            neighbours = new_simplexes[i]->_neighbours;
            for (int j=i+1; j < new_simplexes.size(); j++) {
                if (are_neighbours(new_simplexes[i], new_simplexes[j])) {
                    bool found = (find(neighbours.begin(), neighbours.end(), new_simplexes[j]) != neighbours.end());
                    if (!found) {
                        new_simplexes[i]->_neighbours.push_back(new_simplexes[j]);
                        new_simplexes[j]->_neighbours.push_back(new_simplexes[i]);
                    };
                };
            };
        };
        return new_simplexes;
    };

    bool is_accurate_enough() {
        for (int i=0; i < _funcs.size(); i++) {
            if (!_funcs[i]->is_accurate_enougth()) {
                return false;
            }
        };
        return true;
    };

    // domination returns 0 > , 1 = , 2 <
    int domination(Point* q, Point* p) {
        // 0 - q domiates p
        // 1 - none of them dominates
        // 2 - p dominates q
        int q_better = 0;
        int p_better = 0;
        for (int i=0; i < q->_values.size(); i++) {
            if (q->_values[i] < p->_values[i]) { q_better++; };
            if (q->_values[i] > p->_values[i]) { p_better++; };
        };
        if (q_better == q->_values.size()) { return 0; }
        if (p_better == q->_values.size()) { return 2; }
        return 1;
    };

    bool update_pareto_front(Point* p) { // Returns if value if point was added to pareto_front
        // if any(q > p for q in S):
        //     return
        // for q in [q for q in S if p > q]:
        //     S.remove(q)
        //     S.add(p)

        vector<int> dominated; 
        int relation;
        for (int i=0; i < _pareto_front.size(); i++) {
            relation = domination(_pareto_front[i], p);  // domination returns 0 > , 1 = , 2 <
            if (relation == 0) { return false; };
            if (relation == 2) { dominated.push_back(i); };
        };
        for (int i = dominated.size() - 1; i >= 0; i--) {
            _pareto_front.erase(_pareto_front.begin() + dominated[i], _pareto_front.begin() + dominated[i] +1);
        };
        _pareto_front.push_back(p);
        return true;
    };

    void minimize(vector<Function*> funcs){
        _funcs = funcs;
        timestamp_t start = get_timestamp();
        partition_feasable_region_combinatoricly();     // Note: Should not use global variables
        sort(_partition.begin(), _partition.end(), Simplex::ascending_diameter);
        Simplex::max_diameter = _partition[_partition.size()-1]->_diameter;
        Simplex::update_estimates(_partition, _funcs, _pareto_front, 0);

        while (_funcs[0]->_calls <= _max_calls && _duration <= _max_duration && !is_accurate_enough()) { // _func->pe() > _min_pe){
            // Selects simplexes to divide
            vector<Simplex*> simplexes_to_divide;
            if (_iteration == 0) {
                simplexes_to_divide = _partition;
            } else {
                simplexes_to_divide = select_simplexes_to_divide();
            };

            // if (_iteration == 200) {
            //     cout << "Selected simpls: " << _iteration << endl;
            //     for (int i=0; i < simplexes_to_divide.size(); i++) {
            //         // simplexes_to_divide[i]->print();
            //         cout << "(" << simplexes_to_divide[i]->_diameter << ", " << simplexes_to_divide[i]->_tolerance<< "), ";
            //     };
            //     cout << endl;
            // }

            // Divide seletected simplexes
            vector<Simplex*> new_simplexes = divide_simplexes(simplexes_to_divide);

            // if (_iteration == 10) {
            //     Simplex::log_partition(_partition, simplexes_to_divide);
            //     cout << "Loged successfuly" << endl;
            //     exit(0);
            // };

            //// Draw partition in each iteration:
            // Simplex::log_partition(_partition, simplexes_to_divide);
            // FILE* testp = popen("python log/show_partition.py log/partition.txt", "r");
            // pclose(testp);

            //// Draw pareto front in each iteration:
            // Simplex::log_front(_pareto_front, simplexes_to_divide);
            // FILE* test = popen("python log/show_front.py log/front.txt", "r");
            // pclose(test);

            // Remove partitioned simplexes from _partition
            _partition.erase(remove_if(_partition.begin(), _partition.end(), Simplex::not_in_partition), _partition.end());
            // Delete simplexes
            for (int i=0; i < simplexes_to_divide.size(); i++) {
                delete simplexes_to_divide[i];
            };
            simplexes_to_divide.clear();

            // Add new simplexes to _partition
            for (int i=0; i < new_simplexes.size(); i++) {
                _partition.push_back(new_simplexes[i]);
            };

            // Update estimates
            sort(_partition.begin(), _partition.end(), Simplex::ascending_diameter);
            Simplex::max_diameter = _partition[_partition.size()-1]->_diameter;
            Simplex::update_estimates(_partition, _funcs, _pareto_front, _iteration);

            // Update counters and log the status
            _iteration += 1;
            cout << _iteration << ". Simplexes: " << _partition.size() << "  calls: " << _funcs[0]->_calls << "  f_min:" << _funcs[0]->_f_min << endl;

            timestamp_t end = get_timestamp();
            _duration = (end - start) / 1000000.0L;
        };

        if ((_funcs[0]->_calls <= _max_calls) && (_duration <= _max_duration)) {
            _status = "D"; 
        } else {
            _status = "S";
        };

        // Draw partitioning: output simplex coordinates to file and draw it with Python
    };

    virtual ~Asimpl(){};
};

#endif
