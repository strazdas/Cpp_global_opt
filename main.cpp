#include <iostream>
#include "Disimplv.h"

using namespace std;


void print_stats(int calls[], int subregions[], int n) {
    sort(calls, calls+n);
    sort(subregions, subregions+n);
    int calls_sum = 0;
    for (int i=0; i < n; i++) {
        calls_sum += calls[i];
    };
    cout << "Calls50: " << calls[49] << " Calls100: " << calls[99] << " Average: " << calls_sum/100. <<
        " Subregions50: " << subregions[49] << " Subregions100: " << subregions[99] << endl; 
};


int main() {
    GKLSFunction* func;
    Disimplv* alg;
    int n = 100;
    for (int cls=4; cls <= 4; cls++) {
        int calls[100];
        int subregions[100];
        for (int fid=1; fid <= n; fid++) {
            //// Laisvės laipsniai:
            // Kurioje viršūnėje skaičiuoti simplekso gradientą (centre?).
            // LowerBoundStrategy { MinVert, LongestEdgeLB, LowestEdgeLB };
            // LStrategy { Self, ParentSelf, NeighboursSelf, ParentNeighboursSelf };
            // parent_L_part
            alg = new Disimplv(LongestEdgeLB, Neighbours);
            // alg = new Disimplv();
            func = new GKLSFunction(cls, fid);

            alg->minimize(func);
            if (fid == 1 && cls == 1) { cout << alg->_name << " " << alg->_parent_L_part << endl; };
            cout << (*alg) << endl;

            calls[fid-1] = func->_calls;
            subregions[fid-1] = alg->_partition.size();

            delete alg;
            delete func;
            // break;
        };
        cout << "Class " << cls << " "; 
        print_stats(calls, subregions, n);
    };
    return 0;

    // GKLSFunction* func;
    // SimplexTree* tree;
    // tree = new SimplexTree();
    //
    // func = new GKLSFunction(1, 1);
    //
    // Simplex* simpl = new Simplex(LongestEdgeLB, Neighbours, 0, FFMinVert);
    // double coords[2] = {0, 0}; // 1
    // Point* point = new Point(coords, 2);  point->add_value(0.5);
    // simpl->add_vertex(point);
    // coords[0] = 1;  coords[1] = 0; // 2
    // point = new Point(coords, 2);  point->add_value(0.3);
    // simpl->add_vertex(point);
    // coords[0] = 1;  coords[1] = 1; // 3
    // point = new Point(coords, 2);  point->add_value(0.6);
    // simpl->add_vertex(point);
    // simpl->init_parameters(func);
    // tree->add(simpl);
    // cout << "1. " << simpl->_grad_norm << endl;
    // // simpl->print();
    //
    // simpl = new Simplex(LongestEdgeLB, Neighbours, 0, FFMinVert);
    // coords[0] = 0;  coords[1] = 0;
    // point = new Point(coords, 2);  point->add_value(0.7);
    // simpl->add_vertex(point);
    // coords[0] = 1;  coords[1] = 0;
    // point = new Point(coords, 2);  point->add_value(0.0);
    // simpl->add_vertex(point);
    // coords[0] = 1;  coords[1] = 1;
    // point = new Point(coords, 2);  point->add_value(0.8);
    // simpl->add_vertex(point);
    // simpl->init_parameters(func);
    // tree->add(simpl);
    // cout << "2. " << simpl->_grad_norm << endl;
    // // simpl->print();
    //
    // simpl = new Simplex(LongestEdgeLB, Neighbours, 0, FFMinVert);
    // coords[0] = 0;  coords[1] = 0;
    // point = new Point(coords, 2);  point->add_value(0.4);
    // simpl->add_vertex(point);
    // coords[0] = 1;  coords[1] = 0;
    // point = new Point(coords, 2);  point->add_value(0.3);
    // simpl->add_vertex(point);
    // coords[0] = 1;  coords[1] = 1;
    // point = new Point(coords, 2);  point->add_value(100.9);
    // simpl->add_vertex(point);
    // simpl->init_parameters(func);
    // tree->add(simpl);
    // cout << "3. " << simpl->_grad_norm << endl;
    // // simpl->print();
    //
    // simpl = new Simplex(LongestEdgeLB, Neighbours, 0, FFMinVert);
    // coords[0] = 0;  coords[1] = 0;
    // point = new Point(coords, 2);  point->add_value(4.7);
    // simpl->add_vertex(point);
    // coords[0] = 1;  coords[1] = 0;
    // point = new Point(coords, 2);  point->add_value(-1.0);
    // simpl->add_vertex(point);
    // coords[0] = 1;  coords[1] = 1;
    // point = new Point(coords, 2);  point->add_value(4.3);
    // simpl->add_vertex(point);
    // simpl->init_parameters(func);
    // cout << "Pirmas pridėjimas " << tree->add(simpl) << endl;
    // cout << "Pakartotinis pridėjimas " << tree->add(simpl) << endl;
    // cout << "4. " << simpl->_grad_norm << endl;
    // // simpl->print();
    //
    // tree->print();
    // cout << "Max grad norm" << tree->_max_grad_norm << endl;
    //
    // // Add four simplexes 
    // // Print tree after adding each simplex
    //
    // exit(0);
};
