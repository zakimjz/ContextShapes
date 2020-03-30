#include <cmath>
#include <vector>
#include "util.h"

#define SMALL_NUM 0.00000001  // anything that avoids division overflow

using namespace std;

/** 
 * This subroutine compute the area of a triangular region, each vertex of the 
 * triangle are given as vector (of size 3) of double
 */
double compute_triangle_area(vector<double>& v1, vector<double>& v2, vector<double>& v3) {
    double a12, a23, a31, p, rtn_area;
    a12  = (v1[0] - v2[0]) * (v1[0] - v2[0]); 
    a12 += (v1[1] - v2[1]) * (v1[1] - v2[1]); 
    a12 += (v1[2] - v2[2]) * (v1[2] - v2[2]); 
    if (a12 <= SMALL_NUM) return 0;
    a12 = sqrt(a12);

    a23  = (v3[0] - v2[0]) * (v3[0] - v2[0]); 
    a23 += (v3[1] - v2[1]) * (v3[1] - v2[1]); 
    a23 += (v3[2] - v2[2]) * (v3[2] - v2[2]); 
    if (a23 <= SMALL_NUM) return 0;
    a23 = sqrt(a23);

    a31  = (v1[0] - v3[0]) * (v1[0] - v3[0]); 
    a31 += (v1[1] - v3[1]) * (v1[1] - v3[1]); 
    a31 += (v1[2] - v3[2]) * (v1[2] - v3[2]); 
    if (a31 <= SMALL_NUM) return 0;
    a31 = sqrt(a31);

    p = (a12 + a23 + a31) / 2;

    rtn_area = sqrt(p * (p - a12) * (p - a23) * (p - a31));

    return rtn_area;
}

double rdl_vector_ssd(double pm_v1[], double pm_v2[]) {
    double     ssd;
    ssd  = (pm_v1[0] - pm_v2[0]) * (pm_v1[0] - pm_v2[0]);
    ssd += (pm_v1[1] - pm_v2[1]) * (pm_v1[1] - pm_v2[1]);
    ssd += (pm_v1[2] - pm_v2[2]) * (pm_v1[2] - pm_v2[2]);
    ssd = sqrt(ssd);
    return ssd;
}

/**
 * Finding volume of a conic object
 */
double rdl_cone_volume(vector<pair<double, double> >& pm_segments, int pm_num_points) {
    size_t i;
    double  v = 0;
    double  c = (4.0 * M_PI) / (3.0 * pm_num_points);
    double  r2_cubic, r1_cubic;

    for (i = 0; i < pm_segments.size(); i ++){
        r2_cubic = pm_segments[i].second * pm_segments[i].second * pm_segments[i].second;
        r1_cubic = pm_segments[i].first * pm_segments[i].first * pm_segments[i].first;
        v += c * (r2_cubic - r1_cubic);
    }

    return v;
}

/*
 * A matrix-vector multiplication
 */
vector<double> rdl_matrix_application_3d(const vector<vector<double> >& pm_matrix, const vector<double>& pm_vector) {
    vector<double>  rtn_vector(3, 0);
    int i = 0;

    for (i = 0; i < 3; ++ i){
        rtn_vector[i] = pm_matrix[i][0] * pm_vector[0] + pm_matrix[i][1] * pm_vector[1] + pm_matrix[i][2] * pm_vector[2];
    }
    return rtn_vector;
}

/*
 * Transpose a matrix
 */
vector<vector<double> > rdl_matrix_transpose_3d(const vector<vector<double> >& pm_matrix) {
    vector<vector<double> > rtn_matrix(3, vector<double>(3, 0));

    rtn_matrix[0][1] = pm_matrix[1][0];
    rtn_matrix[1][0] = pm_matrix[0][1];
    rtn_matrix[0][2] = pm_matrix[2][0];
    rtn_matrix[2][0] = pm_matrix[0][2];
    rtn_matrix[1][2] = pm_matrix[2][1];
    rtn_matrix[2][1] = pm_matrix[1][2];
    rtn_matrix[0][0] = pm_matrix[0][0];
    rtn_matrix[1][1] = pm_matrix[1][1];
    rtn_matrix[2][2] = pm_matrix[2][2];

    return rtn_matrix;
}

/*
 * Two 3X3 matrix multiplication
 */
vector<vector<double> > rdl_matrix_multiplication_3d(const vector<vector<double> >& pm_matrix_A, const vector<vector<double> >& pm_matrix_B) {
    vector<vector<double> > rtn_matrix(3, vector<double>(3, 0));
    int i, j;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            rtn_matrix[i][j] = pm_matrix_A[i][0] * pm_matrix_B[0][j] + pm_matrix_A[i][1] * pm_matrix_B[1][j] + pm_matrix_A[i][2] * pm_matrix_B[2][j];     
        }
    }

    return rtn_matrix;
}
