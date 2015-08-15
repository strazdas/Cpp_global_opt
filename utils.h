#ifndef UTILS_H
#define UTILS_H 
#include "functions.h"
#include <sys/time.h>

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

#endif
