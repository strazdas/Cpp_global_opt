#include <iostream>
#include "Elleme.h"

using namespace std;

enum LowerBoundStrategy { MinVert, LongestEdgeLB, LowestEdgeLB, All };
enum LStrategy { Self, ParentSelf, Neighbours };
enum DivisionStrategy { LongestHalf };
enum SimplexGradientStrategy { FFMinVert, FFMaxVert, FFAllVertMean };

void partition_feasable_region_combinatoricly(GKLSFunction func, vector<Simplex*> partition){
    int n = func->_D;
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
            Point* tmp_point = new Point(triangle[i], n);
            Point* point = func->get(tmp_point); 
            if (tmp_point != point) {
                delete tmp_point;
            };
            simpl->add_vertex(point);
        };
        simpl->init_parameters(func);
        partition.push_back(simpl);

    } while (next_permutation(teta, teta+n));
    return partition;
};



int main() {
    GKLSFunction* func = new GKLSFunction(1, 1);
    // partition_feasible_region
    vector<Simplex*> _partition;
    Elleme* alg = new Elleme();
    // solve problem with one simplex from partitioned_feasible_region
    Point* result = alg->minimize(simplexes[0]);
    // print results
    cout << result << endl;
    delete alg;
    delete func;
    return 0;
};
