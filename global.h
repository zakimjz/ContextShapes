#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <vector>
#include <sstream>

#define GET_WEIGHT precomputed32_weight

using namespace std;
#define SMALL_NUM 0.00000001

double g_sample_region_radius;
double g_grid_cell_len;

extern const unsigned int g_solid_volume_num_dynamicR;
int g_bov_filter_selector = 1;

double g_bov_inInner2A_CORE_upper_bound = 25;

int g_solid_volume_radius_selector;
double g_bov_inInner1A_upper_bound   = 100;
double g_bov_inInner2A_upper_bound   = 50;
double g_bov_inInner3A_upper_bound   = 15;
double g_bov_inInner4A_upper_bound   = 0;
double g_bov_inCore_upper_bound      = 0;

string g_dir_pdb_DB;
string g_dir_pdb_Query;
string g_dir_database_DB;
string g_dir_database_Query;
string g_dir_table;
string g_dir_output;
string g_cb_suffix_DB;
string g_cb_suffix_Query;
string g_match_pair_str;

int     g_task = 100;

double g_solid_volume_min_radius;
double g_solid_angle_lower_bound;
double g_solid_angle_upper_bound;
int _solid_volume_radius_selector = 5;
double g_distance_true_pair = 2.8;
double g_rmsd_true_pair = 2.0;

int g_matching_table_selector = 40;
unsigned int g_num_predicates = 4000;
unsigned int g_num_threads = 2;

double g_BSA_ratio_inner = 0;
double  g_BSA_ratio_core_4A = 0;
double  g_BSA_ratio_core_3A = 0;
double  g_BSA_ratio_core_2A = 1;
double  g_BSA_ratio_core_1A = 1;
double  g_BSA_ratio_shell_1A = 1;
double  g_BSA_ratio_shell_2A = 1;
double  g_BSA_ratio_shell_3A = 1;
double  g_BSA_ratio_shell_4A = 0;

unsigned long  g_vmpeak;
unsigned long  g_sizeof_mt;
unsigned long  g_sizeof_grid;

unsigned long  g_sizeof_single_cs_coarse;
unsigned long  g_sizeof_cs_coarse;
unsigned long  g_numberof_cs_coarse;

unsigned long  g_sizeof_single_cs_dense;
unsigned long  g_sizeof_cs_dense;
unsigned long  g_numberof_cs_dense;

vector<pair<double, pair<vector<vector<int> >, vector<vector<int> > > > > g_dynamic_radius; 

         

template <class FirstType, class SecondType>
struct great_than_univ{
    bool operator()(pair<FirstType, SecondType> pm_A, pair<FirstType, SecondType> pm_B)
    {   
        if (pm_A.first > pm_B.first) return true;
        else return false;
    }
};

template <class FirstType, class SecondType>
struct equal_to_univ{
    bool operator()(pair<FirstType, SecondType> pm_A, pair<FirstType, SecondType> pm_B)
    {
        if (fabs(pm_A.first - pm_B.first) < SMALL_NUM) return true;
        else return false;
    }
};

template <class FirstType, class SecondType>
struct great_than_predicate{
    great_than_predicate(FirstType pm_value)
    {
        m_value = pm_value;
    }

    bool operator()(pair<FirstType, SecondType> pm_A)
    {
        if (pm_A.first > m_value) return true;
        else return false;
    }

private:
    FirstType  m_value;
};


stringstream    g_parameters_str;


void compute_dynamic_volume_template_R(double min_r, double max_r, vector<vector<int> >& out_inner, vector<vector<int> >& out_surface);

void compute_dynamic_volume_template()
{
    vector<double>  tmp_radius;
    double r;
    for (r = g_solid_volume_min_radius; r < g_solid_volume_min_radius + g_solid_volume_num_dynamicR; r += 1){
        tmp_radius.push_back(r);
    }

    size_t i;

    g_dynamic_radius.clear();
    //handle the first r
    vector<vector<int> > tmp_inner;
    vector<vector<int> > tmp_surface;
    compute_dynamic_volume_template_R(0, tmp_radius[0], tmp_inner, tmp_surface);
    g_dynamic_radius.push_back(pair<double, pair<vector<vector<int> >, vector<vector<int> > > >(tmp_radius[0], pair<vector<vector<int> >, vector<vector<int> > >(tmp_inner, tmp_surface)));

    //handle the rest r
    for (i = 1; i < tmp_radius.size(); i ++){
        tmp_inner.clear();
        tmp_surface.clear();
        compute_dynamic_volume_template_R(tmp_radius[i-1], tmp_radius[i], tmp_inner, tmp_surface);
        g_dynamic_radius.push_back(pair<double, pair<vector<vector<int> >, vector<vector<int> > > >(tmp_radius[i], pair<vector<vector<int> >, vector<vector<int> > >(tmp_inner, tmp_surface)));
    }
}

void compute_dynamic_volume_template_R(double min_r, double max_r, vector<vector<int> >& out_inner, vector<vector<int> >& out_surface)
{
    out_inner.clear();
    out_surface.clear();

    vector<int>             tmp_cell(3, 0);
    double                  tmp_min_dist_sq = min_r * min_r;
    double                  tmp_max_dist_sq = max_r * max_r;

    int ix, iy, iz;
    int tmp_grid_size = static_cast<int>(max_r / g_grid_cell_len + 1);

    for (ix = 0 - tmp_grid_size; ix <= tmp_grid_size; ix ++){
        for (iy = 0 - tmp_grid_size; iy <= tmp_grid_size; iy ++){
            for (iz = 0 - tmp_grid_size; iz <= tmp_grid_size; iz ++){
                bool has_inner = false;
                bool has_outer_min = false;
                bool has_outer_max = false;
                //bool has_surf_min = false;
                //bool has_surf_max = false;

                double  tmp_dist;
                tmp_cell[0] = ix;
                tmp_cell[1] = iy;
                tmp_cell[2] = iz;

                //000
                tmp_dist  = (ix * g_grid_cell_len) * (ix * g_grid_cell_len);
                tmp_dist += (iy * g_grid_cell_len) * (iy * g_grid_cell_len);
                tmp_dist += (iz * g_grid_cell_len) * (iz * g_grid_cell_len);
                if (tmp_dist < tmp_max_dist_sq && tmp_dist > tmp_min_dist_sq) has_inner = true;
                if (tmp_dist > tmp_max_dist_sq) has_outer_max = true;
                if (tmp_dist < tmp_min_dist_sq) has_outer_min = true;
                //if (tmp_dist == tmp_min_dist_sq) has_surf_min = true;
                //if (tmp_dist == tmp_max_dist_sq) has_surf_max = true;

                //100
                tmp_dist  = ((ix + 1) * g_grid_cell_len) * ((ix + 1) * g_grid_cell_len);
                tmp_dist += (iy * g_grid_cell_len) * (iy * g_grid_cell_len);
                tmp_dist += (iz * g_grid_cell_len) * (iz * g_grid_cell_len);
                if (tmp_dist < tmp_max_dist_sq && tmp_dist > tmp_min_dist_sq) has_inner = true;
                if (tmp_dist > tmp_max_dist_sq) has_outer_max = true;
                if (tmp_dist < tmp_min_dist_sq) has_outer_min = true;
                if (has_inner && has_outer_max) {
                    out_surface.push_back(tmp_cell);
                    continue;
                }

                //010
                tmp_dist  = (ix * g_grid_cell_len) * (ix * g_grid_cell_len);
                tmp_dist += ((iy + 1) * g_grid_cell_len) * ((iy + 1) * g_grid_cell_len);
                tmp_dist += (iz * g_grid_cell_len) * (iz * g_grid_cell_len);
                if (tmp_dist < tmp_max_dist_sq && tmp_dist > tmp_min_dist_sq) has_inner = true;
                if (tmp_dist > tmp_max_dist_sq) has_outer_max = true;
                if (tmp_dist < tmp_min_dist_sq) has_outer_min = true;
                if (has_inner && has_outer_max) {
                    out_surface.push_back(tmp_cell);
                    continue;
                }

                //001
                tmp_dist  = (ix * g_grid_cell_len) * (ix * g_grid_cell_len);
                tmp_dist += (iy * g_grid_cell_len) * (iy * g_grid_cell_len);
                tmp_dist += ((iz + 1) * g_grid_cell_len) * ((iz + 1) * g_grid_cell_len);
                if (tmp_dist < tmp_max_dist_sq && tmp_dist > tmp_min_dist_sq) has_inner = true;
                if (tmp_dist > tmp_max_dist_sq) has_outer_max = true;
                if (tmp_dist < tmp_min_dist_sq) has_outer_min = true;
                if (has_inner && has_outer_max) {
                    out_surface.push_back(tmp_cell);
                    continue;
                }

                //110
                tmp_dist  = ((ix + 1) * g_grid_cell_len) * ((ix + 1) * g_grid_cell_len);
                tmp_dist += ((iy + 1) * g_grid_cell_len) * ((iy + 1) * g_grid_cell_len);
                tmp_dist += (iz * g_grid_cell_len) * (iz * g_grid_cell_len);
                if (tmp_dist < tmp_max_dist_sq && tmp_dist > tmp_min_dist_sq) has_inner = true;
                if (tmp_dist > tmp_max_dist_sq) has_outer_max = true;
                if (tmp_dist < tmp_min_dist_sq) has_outer_min = true;
                if (has_inner && has_outer_max) {
                    out_surface.push_back(tmp_cell);
                    continue;
                }

                //101
                tmp_dist  = ((ix + 1) * g_grid_cell_len) * ((ix + 1) * g_grid_cell_len);
                tmp_dist += (iy * g_grid_cell_len) * (iy * g_grid_cell_len);
                tmp_dist += ((iz + 1) * g_grid_cell_len) * ((iz + 1) * g_grid_cell_len);
                if (tmp_dist < tmp_max_dist_sq && tmp_dist > tmp_min_dist_sq) has_inner = true;
                if (tmp_dist > tmp_max_dist_sq) has_outer_max = true;
                if (tmp_dist < tmp_min_dist_sq) has_outer_min = true;
                if (has_inner && has_outer_max) {
                    out_surface.push_back(tmp_cell);
                    continue;
                }

                //011
                tmp_dist  = (ix * g_grid_cell_len) * (ix * g_grid_cell_len);
                tmp_dist += ((iy + 1) * g_grid_cell_len) * ((iy + 1) * g_grid_cell_len);
                tmp_dist += ((iz + 1) * g_grid_cell_len) * ((iz + 1) * g_grid_cell_len);
                if (tmp_dist < tmp_max_dist_sq && tmp_dist > tmp_min_dist_sq) has_inner = true;
                if (tmp_dist > tmp_max_dist_sq) has_outer_max = true;
                if (tmp_dist < tmp_min_dist_sq) has_outer_min = true;
                if (has_inner && has_outer_max) {
                    out_surface.push_back(tmp_cell);
                    continue;
                }

                //111
                tmp_dist  = ((ix + 1) * g_grid_cell_len) * ((ix + 1) * g_grid_cell_len);
                tmp_dist += ((iy + 1) * g_grid_cell_len) * ((iy + 1) * g_grid_cell_len);
                tmp_dist += ((iz + 1) * g_grid_cell_len) * ((iz + 1) * g_grid_cell_len);
                if (tmp_dist < tmp_max_dist_sq && tmp_dist > tmp_min_dist_sq) has_inner = true;
                if (tmp_dist > tmp_max_dist_sq) has_outer_max = true;
                if (tmp_dist < tmp_min_dist_sq) has_outer_min = true;

                if (has_inner && has_outer_max) {
                    out_surface.push_back(tmp_cell);
                }
                else if (has_inner && (!has_outer_max) && (!has_outer_min)){
                    out_inner.push_back(tmp_cell);
                }
            }//for iz
        }//for iy
    }//for ix
}

template <class X>
string rdl_x2str(X num)
{
  ostringstream myStream; //creates an ostringstream object
  myStream << num << flush;

  /*
 *    * outputs the number into the string stream and then flushes
 *       * the buffer (makes sure the output is put into the stream)
 *          */

  return(myStream.str()); //returns the string form of the stringstream object
}

#endif


