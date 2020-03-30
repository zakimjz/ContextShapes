#ifndef _UTIL_H_
#define _UTIL_H_

#include <cmath>
#include <vector>
#include <iostream>

#define SMALL_NUM 0.00000001  // anything that avoids division overflow

using namespace std;

/**
 * 20 uniform points on a unit sphere
 */
static size_t g_uniform_CP_template_size = 20;
static double g_uniform_CP_template[20][3] = {
    {0.000000000000 , 0.000000000000 ,  1.000000000000 },
    {-0.294617976538,   -0.618225921447,  0.728695380767 },
    {0.608925612730 , -0.339257732343,  0.717017286547 },
    {-0.104735221856,   0.691394698346 ,  0.714845370974 },
    {0.684908582035 , 0.439622057675 ,  0.581061684041 },
    {-0.757605356218,   0.312817091680 ,  0.572869611154 },
    {-0.874764920042,   -0.378746461147,  0.302220867633 },
    {0.018257132938 , -0.986007156931,  0.165700222020 },
    {0.740478385735 , -0.661268064208,  0.120067928762 },
    {0.176416014581 , 0.980045628564 ,  0.091585783464 },
    {0.935190188535 , 0.345746250636 ,  -0.076673603271},
    {-0.584115620751,   0.807883650985 ,  -0.078313141080},
    {-0.592148763977,   -0.746071641816,  -0.304527415183},
    {-0.943264692334,   0.120121163696 ,  -0.309552299664},
    {0.837544872724 , -0.246886554412,  -0.487407032596},
    {0.399836533541 , 0.735147635884 ,  -0.547438306935},
    {0.217601482779 , -0.767863216534,  -0.602524418913},
    {-0.363154390171,   0.520840309655 ,  -0.772556962779},
    {-0.410955115499,   -0.275888234942,  -0.868908266082},
    {0.317474868772 , 0.063550278964 ,  -0.946134805269}
};
/** 
 * This subroutine compute the area of a triangular region, each vertex of the 
 * triangle are given as vector (of size 3) of double
 */
double compute_triangle_area(vector<double>& v1, vector<double>& v2, vector<double>& v3); 

/**
 * Rotate a point (template function, so defining in header file)
 */
template <class T>
void rdl_rotate(vector<T> const &pm_axis, vector<T> const &pm_old_point, T pm_angle, vector<T> &pm_new_point) {
    pm_new_point[0] = 0.0;
    pm_new_point[1] = 0.0;
    pm_new_point[2] = 0.0;
    double costheta, sintheta;

    costheta = cos(pm_angle);
    sintheta = sin(pm_angle);

    pm_new_point[0] += (costheta + (1 - costheta) * pm_axis[0] * pm_axis[0]) * pm_old_point[0];
    pm_new_point[0] += ((1 - costheta) * pm_axis[0] * pm_axis[1] - pm_axis[2] * sintheta) * pm_old_point[1];
    pm_new_point[0] += ((1 - costheta) * pm_axis[0] * pm_axis[2] + pm_axis[1] * sintheta) * pm_old_point[2];

    pm_new_point[1] += ((1 - costheta) * pm_axis[0] * pm_axis[1] + pm_axis[2] * sintheta) * pm_old_point[0];
    pm_new_point[1] += (costheta + (1 - costheta) * pm_axis[1] * pm_axis[1]) * pm_old_point[1];
    pm_new_point[1] += ((1 - costheta) * pm_axis[1] * pm_axis[2] - pm_axis[0] * sintheta) * pm_old_point[2];

    pm_new_point[2] += ((1 - costheta) * pm_axis[0] * pm_axis[2] - pm_axis[1] * sintheta) * pm_old_point[0];
    pm_new_point[2] += ((1 - costheta) * pm_axis[1] * pm_axis[2] + pm_axis[0] * sintheta) * pm_old_point[1];
    pm_new_point[2] += (costheta + (1 - costheta) * pm_axis[2] * pm_axis[2]) * pm_old_point[2];

    return;
}

/**
 * Compute dot procut of two vectors
 */
template <class T>
T rdl_dot_3d(const vector<T>& pm_v1, const vector<T>& pm_v2) {

    T     r;
    r  = pm_v1[0] * pm_v2[0];
    r += pm_v1[1] * pm_v2[1];
    r += pm_v1[2] * pm_v2[2];

    return r;
}
/**
 * Compute the cross-product of two vectors
 */
template <class T>
vector<T> rdl_cross_3d(const vector<T>& pm_v1, const vector<T>& pm_v2){ 
    vector<T>     vr(3, 0);

    vr[0] = pm_v1[1] * pm_v2[2] - pm_v1[2] * pm_v2[1];
    vr[1] = pm_v1[2] * pm_v2[0] - pm_v1[0] * pm_v2[2];
    vr[2] = pm_v1[0] * pm_v2[1] - pm_v1[1] * pm_v2[0];

    return vr;
}
/**
 * Compute the euclid disance between two points
 */
template <class T>
T rdl_vector_ssd(const vector<T>& pm_v1, const vector<T>& pm_v2) {
    T     ssd;
    ssd  = (pm_v1[0] - pm_v2[0]) * (pm_v1[0] - pm_v2[0]);
    ssd += (pm_v1[1] - pm_v2[1]) * (pm_v1[1] - pm_v2[1]);
    ssd += (pm_v1[2] - pm_v2[2]) * (pm_v1[2] - pm_v2[2]);
    ssd = sqrt(ssd);
    return ssd;
}
/**
 * Return the euclid dist-sqr between two points
 */
template <class T>
T rdl_vector_ssd_sqr(const vector<T>& pm_v1, const vector<T>& pm_v2) {
    T     ssd;
    ssd  = (pm_v1[0] - pm_v2[0]) * (pm_v1[0] - pm_v2[0]);
    ssd += (pm_v1[1] - pm_v2[1]) * (pm_v1[1] - pm_v2[1]);
    ssd += (pm_v1[2] - pm_v2[2]) * (pm_v1[2] - pm_v2[2]);
    return ssd;
}

double rdl_vector_ssd(double pm_v1[], double pm_v2[]) ;

/**
 * Vector sum operation
 */
template <class T>
vector<T> operator + (const vector<T>& pm_vert_end, const vector<T>& pm_vert_start){
    vector<T> vr(3, 0);
    vr[0] = pm_vert_end[0] + pm_vert_start[0];
    vr[1] = pm_vert_end[1] + pm_vert_start[1];
    vr[2] = pm_vert_end[2] + pm_vert_start[2];

    return vr;
}
/**
 * Vector subtract operation
 */
template <class T>
vector<T> operator - (const vector<T>& pm_vert_end, const vector<T>& pm_vert_start) {
    vector<T> vr(3, 0);
    vr[0] = pm_vert_end[0] - pm_vert_start[0];
    vr[1] = pm_vert_end[1] - pm_vert_start[1];
    vr[2] = pm_vert_end[2] - pm_vert_start[2];

    return vr;
}
/**
 * Scaler scaling operation
 */
template <class T>
vector<T> operator * (T c, const vector<T>& pm_vert_start) {
    vector<T> vr(3, 0);
    vr[0] = c * pm_vert_start[0];
    vr[1] = c * pm_vert_start[1];
    vr[2] = c * pm_vert_start[2];

    return vr;
}
/**
 * Finding volume of a conic object
 */
double rdl_cone_volume(vector<pair<double, double> >& pm_segments, int pm_num_points) ; 
/*
 * Normalize a vector to have unit length
 */
template <class T>
void normalize_unit(T &pm_x, T &pm_y, T &pm_z) { 
    if (fabs(pm_x) < SMALL_NUM && fabs(pm_y) < SMALL_NUM && fabs(pm_z) < SMALL_NUM) return;

    T temp_r;
    temp_r  = pm_x * pm_x;
    temp_r += pm_y * pm_y;
    temp_r += pm_z * pm_z;
    temp_r = sqrt(temp_r);

    pm_x /= temp_r;
    pm_y /= temp_r;
    pm_z /= temp_r;

    return;
}
template <class T>
void normalize_unit_3d(vector<T> &pm_v){ 
    normalize_unit(pm_v[0], pm_v[1], pm_v[2]);
}
/*
 * Return the norm of a vector
 */
template <class T>
T norm_3d(vector<T> &pm_v){ 
    T n;
    n  = pm_v[0] * pm_v[0];
    n += pm_v[1] * pm_v[1];
    n += pm_v[2] * pm_v[2];

    return sqrt(n);
}
/*
 * A matrix-vector multiplication
 */
vector<double> rdl_matrix_application_3d(const vector<vector<double> >& pm_matrix, const vector<double>& pm_vector); 

/*
 * Transpose a matrix
 */
vector<vector<double> > rdl_matrix_transpose_3d(const vector<vector<double> >& pm_matrix); 

/*
 * Two 3X3 matrix multiplication
 */
vector<vector<double> > rdl_matrix_multiplication_3d(const vector<vector<double> >& pm_matrix_A, const vector<vector<double> >& pm_matrix_B); 

/**
 * Print a point
 */
template <class T>
void print_vector_3d(vector<T>& pm_v){
    cout << pm_v[0] << " "
        << pm_v[1] << " "
        << pm_v[2] << " ";

}

#endif
