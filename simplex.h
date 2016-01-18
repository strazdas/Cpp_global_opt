#ifndef SIMPLEX_H
#define SIMPLEX_H 
/* Simplex is a simplex used for solving main problem, e.g. minimizing GKLS function */
#include <list>
#include "Eigen/Dense"
#include "Elbme.h"
#include "Conte.h"

using namespace std;


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


class Simplex {  // Designed for outer problems
    Simplex(const Simplex& other){}
    Simplex& operator=(const Simplex& other){}
public:
    Simplex(LowerBoundStrategy lower_bound_strategy,
            LStrategy L_strategy,
            double parent_L_part,
            SimplexGradientStrategy simplex_gradient_strategy
        ) {
        _lower_bound_strategy = lower_bound_strategy;
        _L_strategy = L_strategy;
        _parent_L_part = parent_L_part;
        _simplex_gradient_strategy = simplex_gradient_strategy;
        _is_in_partition = true;
        _parent = 0;
        _diameter = 0;
        _tolerance = 0;
        _le_v1 = 0;
        _le_v2 = 0;
        // _min_vert = 0;
        // _max_vert = 0;
        // _max_vert_value = -numeric_limits<double>::max();
        // _min_vert_value = numeric_limits<double>::max();
        _should_be_divided = false;
        _should_estimates_be_updated = true;
        // _min_lbs = 0;
        // _min_lb_value = 0;
        _D = 0;
        //
    };

    LowerBoundStrategy _lower_bound_strategy;
    LStrategy _L_strategy;
    SimplexGradientStrategy _simplex_gradient_strategy;

    int _D;                       // Variable space dimension
    int _C;                       // Criteria (objective function) space dimension
    double _tolerance;            // Distance between _min_lb and _min_vert  or _min_lb_value and _min_vert_value

    vector<Point*> _verts;        // Simplex vertexes (points with coordinates and values)
    bool _is_in_partition;
    bool _should_be_divided;  // Should be divided in next iteration
    bool _should_estimates_be_updated;   // Should Lipschitz constant estimate and its lower bound be updated
    Simplex* _parent;
    list<Simplex*> _neighbours;

    Point* _le_v1;      // Longest edge vertex1
    Point* _le_v2; 
    double _diameter;   // Longest edge length

    vector<double> _Ls;          // Cumulative estimates of Lipschitz constants for each criteria
    vector<double> _grad_norms;  // Lipschitz constant estimate calculated by Simplex Gradient Euclidean norm.
    double _parent_L_part;

    // Point* _min_vert;   // Pointer to vertex with lowest function value 
    // double _min_vert_value;  // _min_vert function value 
    // Point* _max_vert;
    // double _max_vert_value;
    // double _metric__vert_min_value;     // _f_min - glob_f / _diameter

    vector<Point*> _min_lbs;
    // double _min_lb_value;
    // double _metric__min_lb;       // _f_min - glob_f / _diameter
    
    void init_parameters(vector<Function*> funcs) {   // Called when all verts have been added
        _D = _verts.size() - 1;
        _C = funcs.size();
        for (int i=0; i < _C; i++) {
            _Ls.push_back(0);
            _grad_norms.push_back(0);
        };
        
        // Note: claculating metrics needed by algorithm here would reduce calculations

        // Sorts vertexes using first criteria values
        sort(_verts.begin(), _verts.end(), Point::compare_by_value);  // Is there calculation of intersection anywhere in the algorithm, in this case sorting by adress is needed?

        // Find longest edge length and verts
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

        // Find adaptive Lipschitz constant
        // double E;
        // if (1e-4 * fabs(funcs[0]->_glob_f) > 1e-8) {
        //     E = 1e-4 * fabs(funcs[0]->_glob_f);
        // } else {
        //     E = 1e-8;
        // };

        for (int i=0; i < funcs.size(); i++) {
            _grad_norms[i] = find_simplex_gradient_norm(i, _simplex_gradient_strategy);      // Check in the article if global Lipschitz constant is defined
        };
        // ToDo: Should rename this accordingly to the article?
    };

    // Need a scenario where a single simplex is created and I can test with it  
    vector<Point*> find_accurate_lb_min_estimates(vector<Point*> verts, vector<double> Ls) {
        vector<Point*> estimates_of_accurate_lb_min;
        for (int i=0; i < Ls.size(); i++) {
            // Elbme* alg = new Elbme(verts, Ls, i);
            Conte* alg = new Conte(verts, Ls, i);
            Point* estimate_of_accurate_lb_min = alg->minimize(); // ->copy();
            delete alg;
            estimates_of_accurate_lb_min.push_back(estimate_of_accurate_lb_min);
        };
        return estimates_of_accurate_lb_min;
    };

    // static void extend_region_with_vertex_neighbours(Point* vertex, SimplexTree* region, int depth);

    static void update_estimates(vector<Simplex*> simpls, vector<Function*> funcs, vector<Point*> pareto_front, int iteration);

    double find_simplex_gradient_norm(int crit_id, SimplexGradientStrategy simplex_gradient_strategy){ 
        double L_estimate = 0;
        // Eigen::VectorXd f_diff(_D);
        // Eigen::MatrixXd x_diff(_D, _D);
        // Eigen::MatrixXd x_diff_inv_T(_D, _D);
        // Eigen::VectorXd grad(_D);
        //
        // if (simplex_gradient_strategy == FFMinVert) {  // Gradient at min vertex
        //     // throw "FFMaxVert gradient strategy not implemented yet";
        // //     for (int i=1; i < D+1; i++) { 
        // //         f_diff(i - 1) = _verts[i]->_values[0] - _min_vert->_values[0];
        // //     }; 
        // //
        // //     for (int i=1; i < D+1; i++) {
        // //         for (int j=0; j < D; j++) {
        // //             x_diff(i-1, j) = _verts[i]->_X[j] - _min_vert->_X[j];
        // //         };
        // //     };
        // };
        // if (simplex_gradient_strategy == FFMaxVert) {  // Gradient at max vertex
        //     // throw "FFMaxVert gradient strategy not implemented yet";
        // //     for (int i=0; i < D; i++) { 
        // //         f_diff(i) = _verts[i]->_values[0] - _max_vert->_values[0];
        // //     }; 
        // //
        // //     for (int i=0; i < D; i++) {
        // //         for (int j=0; j < D; j++) {
        // //             x_diff(i, j) = _verts[i]->_X[j] - _max_vert->_X[j];
        // //         };
        // //     };
        // };
        // if (simplex_gradient_strategy == FFMaxVert) {
        //     // throw "FFMaxVert gradient strategy not implemented yet";
        //     // FFAllVertMean
        // };
        //
        // // Find gradient at lowest point
        // x_diff_inv_T = x_diff.inverse().transpose();
        // for (int i=0; i < _D; i++) {
        //     grad[i] = x_diff_inv_T.row(i).dot(f_diff);
        // };
        //
        // // cout << grad << endl;
        //
        // // Find norm of gradient at _min_vert 
        // for (int i=0; i < grad.size(); i++){
        //     L_estimate += pow(grad[i], 2);
        // };
        // L_estimate = sqrt(L_estimate);

        // Find minimum simplex L and use it if its greater then the estimate
        double simplex_min_L = find_simplex_min_L_l1norm(crit_id);
        if (simplex_min_L > L_estimate) {
            L_estimate = simplex_min_L;
        };
        return L_estimate;
    };

    double find_simplex_min_L(int crit_id) {  // Finds L using Euclidean (l2norm)
        double dist;
        double f_diff;
        double edge_L;
        double max_edge_L = -numeric_limits<double>::max();
        for (int i=0; i < _verts.size(); i++) {
            for (int j=i+1; j < _verts.size(); j++) {
                f_diff = fabs(_verts[i]->_values[crit_id] - _verts[j]->_values[crit_id]);
                dist = l2norm(_verts[i], _verts[j]);
                edge_L = f_diff / dist;   // Note: maybe dist (division by zero) protection is needed?
                if (edge_L > max_edge_L) {
                    max_edge_L = edge_L;
                };
            };
        };
        return max_edge_L;
    };

    double find_simplex_min_L_l1norm(int crit_id) {  // Finds L using City block (l1norm)
        double dist;
        double f_diff;
        double edge_L;
        double max_edge_L = -numeric_limits<double>::max();
        for (int i=0; i < _verts.size(); i++) {
            for (int j=i+1; j < _verts.size(); j++) {
                f_diff = fabs(_verts[i]->_values[crit_id] - _verts[j]->_values[crit_id]);
                dist = l1norm(_verts[i], _verts[j]);
                edge_L = f_diff / dist;   // Note: maybe dist (division by zero) protection is needed?
                if (edge_L > max_edge_L) {
                    max_edge_L = edge_L;
                };
            };
        };
        return max_edge_L;
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

    double find_tolerance(vector<Point*> pareto_front) {
        // _min_lbs - is initialized and contains Ms for each criteria

        // Sekti pareto frontą:  naudojant tentą paprastą algoritmą.
        // rasti atstumą iki pareto fronto.
        // For now find lowest distance to vertex.
        vector<double> M;
        for (int i=0; i < _min_lbs.size(); i++) {
            M.push_back(_min_lbs[i]->_values[0]);
            // cout << M[i] << ", ";
        };
        // cout << pareto_front[0]->_values[0] << " , " << pareto_front[0]->_values[0] << endl;
        // cout << pareto_front[0]->_values[0] << " , " << pareto_front[0]->_values[1] << endl;
        // cout << pareto_front[1]->_values[0] << " , " << pareto_front[1]->_values[0] << endl;
        // cout << pareto_front[1]->_values[0] << " , " << pareto_front[1]->_values[1] << endl;
        // cout << endl;
        // exit(0);

        double min_dist = numeric_limits<double>::max();

        for (int i=0; i < pareto_front.size(); i++) {
            if ((M[0] < pareto_front[i]->_values[0]) && (M[1] < pareto_front[i]->_values[1])) {
                double dist = gtl1norm(pareto_front[i]->_values, M);
                if (dist < min_dist) {
                    min_dist = dist;
                };
            }
        };
        // If M is dominated, than do not divide this simplex
        if (min_dist == numeric_limits<double>::max()) {
            return 0;
        };
        return -min_dist;  // minus because in convex_hull we need to find max tolerance
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

    // static double ascending_min_lb_value(Simplex* s1, Simplex* s2) {
    //     return s1->_min_lb_value < s2->_min_lb_value;
    // };

    static double ascending_tolerance(Simplex* s1, Simplex* s2) {
        return s1->_tolerance < s2->_tolerance;
    };

    static double ascending_diameter(Simplex* s1, Simplex* s2) {
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
        cout << " Simplex   diam:  " << _diameter << "   tol:  " << _tolerance;
        cout << "   L:  ";
        for (int i=0; i < _D; i++) {
            cout << _Ls[i] << " ";
        };
        cout << "   grad_norm:  ";
        for (int i=0; i < _D; i++) {
            cout << _grad_norms[i] << " ";
        };
        cout << endl;
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

    static void log_front(vector<Point*> pareto_front, vector<Simplex*> simplexes_to_divide) {
       ofstream log_file; 
       log_file.open("log/front.txt");
       log_file.close();
       log_file.open("log/front.txt", ios::app);
       for (int j=0; j < pareto_front.size(); j++) {
           for (int i=0; i < pareto_front[j]->size(); i++) {
                log_file << pareto_front[j]->_X[i] << " ";
           };
           log_file << " -> ";
           for (int i=0; i < pareto_front[j]->_values.size(); i++) {
               log_file << pareto_front[j]->_values[i] << " ";
           };
           log_file << endl;
       };

       log_file << "Tolerances" << endl;
       for (int j=0; j < simplexes_to_divide.size(); j++) {
           for (int i=0; i < simplexes_to_divide[j]->_min_lbs.size(); i++) {
              log_file <<  simplexes_to_divide[j]->_min_lbs[i]->_values[0] << " ";
           };
           log_file << endl;
       };
       log_file.close();
    };

    static void log_partition(vector<Simplex*> simplexes,
                              vector<Simplex*> selected,
                              string label="Partition:",
                              int iteration=0) {
       ofstream log_file; 
       log_file.open("log/partition.txt");
       log_file.close();
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
                log_file << " ("<< simplexes[i]->_diameter << "," << simplexes[i]->_tolerance << ")" << endl;
           };
           if (simplexes[i]->_L_strategy == Neighbours) {
                log_file << " ("<< simplexes[i]->_diameter << "," << simplexes[i]->_tolerance << ")" << endl;
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
               log_file << " ("<< selected[i]->_diameter << "," << selected[i]->_tolerance << ")" << endl;
           };
           if (selected[i]->_L_strategy == Neighbours) {
               log_file << " ("<< selected[i]->_diameter << "," << selected[i]->_tolerance << ")" << endl;
           };
       };

       log_file.close();
    };

    virtual ~Simplex(){
        for (int i=0; i < _min_lbs.size(); i++) {
            delete _min_lbs[i];
        };
        _min_lbs.clear();
        // delete _min_lb;
        _verts.clear();
        _neighbours.clear();
    };  
};

void Point::_neighbours_estimates_should_be_updated() {
    for (int sid=0; sid < _simplexes.size(); sid++) {
        _simplexes[sid]->_should_estimates_be_updated = true;
    };
};


// class SimplexTreeNode {
//     SimplexTreeNode(const SimplexTreeNode& other){}
//     SimplexTreeNode& operator=(const SimplexTreeNode& other){}
// public:                
//     SimplexTreeNode(Simplex* value){
//         _height = 1;
//         _value = value;
//         _parent = 0;
//         _left = 0;
//         _right = 0;
//     };
//     int _height;
//     Simplex* _value;
//     SimplexTreeNode* _parent;
//     SimplexTreeNode* _left;
//     SimplexTreeNode* _right;
//
//     void print(){
//         if (_left != 0) { cout << "l"; _left->print(); };
//         cout << _value << "("<< _height << ")";
//         if (_right != 0) { cout << "r"; _right->print(); };
//     };
//
//     virtual ~SimplexTreeNode();
// };
//
// class SimplexTree {  // Binary balancing Simplex tree for storing simplex neighbours 
//                      // Simplex adresses are stored in this tree
//     SimplexTree(const SimplexTree& other){}
//     SimplexTree& operator=(const SimplexTree& other){}
// public:
//     SimplexTree(){
//          _max_grad_norm = -numeric_limits<double>::max();
//          _tree_root = 0;
//     };
//     double _max_grad_norm;
//     SimplexTreeNode* _tree_root;
//
//     void update_height(SimplexTreeNode* node) {
//         int lh = 0;
//         int rh = 0;
//         if (node->_left != 0) { lh = node->_left->_height; }; 
//         if (node->_right != 0) { rh = node->_right->_height; };
//         if (lh > rh) {
//             node->_height = lh + 1;
//         } else {
//             node->_height = rh + 1;
//         };
//         // Also update all ancestors heights
//         if (node->_parent != 0) {
//             update_height(node->_parent);
//         };
//     };
//     void left_right_rebalance(SimplexTreeNode* node) {
//         SimplexTreeNode* diatteched_node;
//         // node left right  <-  node left right left 
//         diatteched_node = node->_left->_right;
//         node->_left->_right = node->_left->_right->_left;
//         if (node->_left->_right != 0) { node->_left->_right->_parent = node->_left; };
//         // Diatteched left = node->_left
//         diatteched_node->_left = node->_left;
//         node->_left->_parent = diatteched_node;
//         // node left  <-  node left right
//         node->_left = diatteched_node;
//         diatteched_node->_parent = node;
//         // Update heights
//         update_height(node);
//         update_height(diatteched_node);
//         update_height(diatteched_node->_left);
//     };
//     void left_left_rebalance(SimplexTreeNode* node) {
//         SimplexTreeNode* diatteched;
//         diatteched = node->_left;
//         node->_left = node->_left->_right;
//         if (node->_left != 0) { node->_left->_parent = node; };  
//         diatteched->_parent = node->_parent;
//         if (node->_parent != 0) {
//             if (node->_parent->_left == node) {
//                 node->_parent->_left = diatteched;
//             } else {
//                 node->_parent->_right = diatteched;
//             };
//         } else {
//             _tree_root = diatteched;
//         };
//         diatteched->_right = node;
//         node->_parent = diatteched;
//         // Update heights
//         update_height(node);
//         update_height(diatteched);
//     };
//     void right_left_rebalance(SimplexTreeNode* node) {
//         SimplexTreeNode* diatteched_node;
//         // node left right  <-  node left right left 
//         diatteched_node = node->_right->_left;
//         node->_right->_left = node->_right->_left->_right;
//         if (node->_right->_left != 0) { node->_right->_left->_parent = node->_right; };
//         // Diatteched left = node->_left
//         diatteched_node->_right = node->_right;
//         node->_right->_parent = diatteched_node;
//         // node left  <-  node left right
//         node->_right = diatteched_node;
//         diatteched_node->_parent = node;
//         // Update heights
//         update_height(node);
//         update_height(diatteched_node);
//         update_height(diatteched_node->_right);
//     };
//     void right_right_rebalance(SimplexTreeNode* node) {
//         SimplexTreeNode* diatteched;
//         diatteched = node->_right;
//         node->_right = node->_right->_left;
//         if (node->_right != 0) { node->_right->_parent = node; };
//         diatteched->_parent = node->_parent;                       
//         if (node->_parent != 0) {
//             if (node->_parent->_left == node) {
//                 node->_parent->_left = diatteched;
//             } else {
//                 node->_parent->_right = diatteched;
//             };
//         } else {
//             _tree_root = diatteched;
//         };
//         diatteched->_left = node;
//         node->_parent = diatteched;
//         // Update heights
//         update_height(node);
//         update_height(diatteched);
//     };
//
//     void check_if_balanced(SimplexTreeNode* node) {  // Rebalances tree if its not balanced
//         int lh = 0;
//         int rh = 0;
//         int llh = 0;
//         int lrh = 0;
//         int rlh = 0;
//         int rrh = 0;
//         if (node->_left != 0) {
//             lh = node->_left->_height;
//             if (node->_left->_left != 0) { llh = node->_left->_left->_height; };
//             if (node->_left->_right != 0) { lrh = node->_left->_right->_height; };
//         };
//         if (node->_right != 0) {
//             rh = node->_right->_height;
//             if (node->_right->_left != 0) { rlh = node->_right->_left->_height; };
//             if (node->_right->_right != 0) { rrh = node->_right->_right->_height; };
//         };
//         if (abs(rh - lh) > 1) {
//             // Not balanced, so rebalance
//             if (rh > lh) {
//                 if (rrh > rlh) {
//                     right_right_rebalance(node);
//                 } else {
//                     right_left_rebalance(node);
//                     right_right_rebalance(node);
//                 };
//             };
//             if (rh < lh) {
//                 if (llh > lrh) {
//                     left_left_rebalance(node);
//                 } else {
//                     left_right_rebalance(node);
//                     left_left_rebalance(node);
//                 };
//             };
//         };
//         if (node->_parent != 0) {
//             check_if_balanced(node->_parent);
//         };
//     };
//
//     Simplex* add(Simplex* value) {
//         // Warning: how is this working, if simplex but not its value is used in comparisons?
//
//         // Get same point or insert given (if inserted returns 0)
//         SimplexTreeNode* node = _tree_root;
//         double grad_norm = value->_grad_norm;
//         if (grad_norm > _max_grad_norm) {
//             _max_grad_norm = grad_norm;
//         };
//         if (_tree_root == 0) {  // Create first tree node
//             _tree_root = new SimplexTreeNode(value);
//         } else {
//             while (true) {  // Walk through tree
//                 if (value > node->_value) {
//                     if (node->_right == 0) {
//                         node->_right = new SimplexTreeNode(value);
//                         node->_right->_parent = node;
//                         update_height(node->_right);
//                         check_if_balanced(node->_right);
//                         return 0;
//                     };
//                     node = node->_right;
//                 } else if (value < node->_value) {
//                     if (node->_left == 0) {
//                         node->_left = new SimplexTreeNode(value);
//                         node->_left->_parent = node;
//                         update_height(node->_left);
//                         check_if_balanced(node->_left);
//                         return 0;
//                     };
//                     node = node->_left;
//                 } else {
//                     // Node value matches given simplex value
//                     return value;
//                 };
//             };
//         };
//     };
//
//     void print() {
//         _tree_root->print();
//         cout << endl;
//     };
//
//     virtual ~SimplexTree();
//
// };
//
// SimplexTreeNode::~SimplexTreeNode() {
//     if (_left != 0) { delete _left; };
//     if (_right != 0) { delete _right; };
// };
//
// SimplexTree::~SimplexTree() {
//     delete _tree_root;
// };


// void Simplex::extend_region_with_vertex_neighbours(Point* vertex, SimplexTree* region, int depth) {
//     // Recursively adds vertex neighbours to region
//     for (int sid=0; sid < vertex->_simplexes.size(); sid++) {
//         Simplex* simpl = vertex->_simplexes[sid];
//         Simplex* result = region->add(simpl);
//         if (depth != 0) {
//             if (result == 0 && simpl->_is_in_partition) {
//                 for (int vid=0; vid < simpl->_verts.size(); vid++) {
//                     if (simpl->_verts[vid] != vertex) {
//                         extend_region_with_vertex_neighbours(simpl->_verts[vid], region, depth-1);
//                     };
//                 };
//             };
//         };    
//     };
// };

void Simplex::update_estimates(vector<Simplex*> simpls, vector<Function*> funcs, vector<Point*> pareto_front, int iteration) {   // Neighbours strategy - updates estimates
    for (int sid=0; sid < simpls.size(); sid++) {
        if (simpls[sid]->_should_estimates_be_updated) {
            // Use simplex's \hat{L} as initial max_grad_norms value
            vector<double> max_grad_norms;
            for (int i=0; i < simpls[sid]->_grad_norms.size(); i++) {
                max_grad_norms.push_back(simpls[sid]->_grad_norms[i]);   
            };

            // Find max \hat{L} among neighbours
            for (list<Simplex*>::iterator it=simpls[sid]->_neighbours.begin(); it != simpls[sid]->_neighbours.end(); ++it) {
                for (int i=0; i < funcs.size(); i++) {
                    if ((*it)->_grad_norms[i] > max_grad_norms[i]) {
                        max_grad_norms[i] = (*it)->_grad_norms[i];
                    };
                };
            };

            // Update simplex's L
            for (int i=0; i < funcs.size(); i++) {
                simpls[sid]->_Ls[i] = max_grad_norms[i];
            };

            //// Find accurate lower bound point and value estimates with given precision
            for (int i=0; i < simpls[sid]->_min_lbs.size(); i++) {
                delete simpls[sid]->_min_lbs[i];
            };

            simpls[sid]->_min_lbs = simpls[sid]->find_accurate_lb_min_estimates(simpls[sid]->_verts, simpls[sid]->_Ls);

            simpls[sid]->_tolerance = simpls[sid]->find_tolerance(pareto_front);

            simpls[sid]->_should_estimates_be_updated = false;
        };
    };

    // Note: gali būti, kad slope apibrėžimas pas mane netinkamas atmetant
    // simpleksus su epsilon (potencialiai optimalių simpleksų parinkimo metu).
    //     simpls[sid]->_lb = simpls[sid].find_lb();
    // };
};

#endif
