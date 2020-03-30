#ifndef _CONTEXTBALL_H_
#define _CONTEXTBALL_H_

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <ext/algorithm>
#include <numeric>
#include <sys/timeb.h>
#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <unistd.h>

#include "rbt_coarse.h"
#include "rbt_dense.h"
#include "rbt_pose.h"
#include "typedefs.h"
#include "SES.h"
#include "util.h"

using namespace std;

extern int g_matching_table_selector;
extern string g_dir_table;
extern double g_sample_region_radius;

const unsigned int g_solid_volume_num_dynamicR = 7;

struct MyBall
{
  UINT_Range_tree_3_type*      m_index;  // an range tree built with the points on the ball, used
                                         // to query points in the ball that are within some distance
                                         // from some point
  double                      m_nn_r;

  unsigned long               m_size;    // memory footprint of this struct, used for resource consumption statistics

  // Constructor
  MyBall();

  // Destructor
  ~MyBall();

  // find all the points of the ball that are within pm_r distance from pm_point
  // used CGAL routine window_query to find it out.
  void query_local(const vector<double>&  pm_point, double pm_r, vector<unsigned int>& pm_out);

  // find the point in the ball that is the nearest point of the point pm_point.
  unsigned short query_nn(const vector<double>&   pm_point);

};

// A coarse ball, used for receptor
// number of points on the ball is small, default 1256
// The points are taken from the file rbt_coarse.h
struct MyBall_Coarse : public MyBall{

    void init();
};

// A dense ball, used for ligand
// number of points on the ball is large, default 3768
// The points are taken from the file rbt_dense.h
struct MyBall_Dense : public MyBall{

    void init();
};

// A more dense ball for pose generation, used for ligand
// number of points on the ball is very large, default 12560
// The points are taken from the file rbt_pose.h
struct MyBall_Pose : public MyBall{
    vector<vector<vector<double> > >        m_arcs;

    void init();

    void init_pose();
};

// A ball with unit radius, mostly used as template
struct MyUnitBall {
    vector<vector<double> >     m_unitball_coarse;
    vector<vector<double> >     m_unitball_dense;
    vector<vector<double> >     m_unitball_pose;

    MyBall_Coarse           m_coarse;
    MyBall_Dense            m_dense;
    MyBall_Pose             m_pose;

    vector<vector<unsigned short> >             m_table;                
    vector<vector<unsigned short> >             m_table_coarse2dense;   
    vector<vector<unsigned short> >             m_table_dense2coarse;   
    vector<vector<vector<double> > >            m_rmatrix;


    double                  m_farthest_neighbor_dist_coarse;
    double                  m_farthest_neighbor_dist_dense;

    unsigned long size();

    void init();
    
    void sload_matching_table(string pm_filename);

    void sload_matching_table_coarse2dense(string pm_filename);

    void sload_matching_table_dense2coarse(string pm_filename);

    void sload_matching_table_rmatrix(string pm_filename);

    void ssave_matching_table(string pm_filename);

    void ssave_matching_table_coarse2dense(string pm_filename);

    void ssave_matching_table_dense2coarse(string pm_filename);

    void ssave_matching_table_rmatrix(string pm_filename);

    void load_unitball();

    void compute_farthest_neighbor_dist();
};

// Actual context ball data structure
struct MyContextBall{
    double                m_center_area;

    vector<vector<double> >*            m_unitball_template;
    int                                 m_ball_id;  //id of SES face
    size_t                m_ball_id_serialno; //sn of the actual points; 
    //when a convex SES face is too big, it is uniformly partitioned into several pieces
    double                              m_ball_center[3];

    //by MyGridSES::compute_CB_solid_volume(MyContextBall*)
    double                              m_solid_volume_dynamicR[g_solid_volume_num_dynamicR];
    double                              m_solid_vector[g_solid_volume_num_dynamicR][3]; //3D vector

    //by MyGridSES::compute_CB_solid_angle(MyContextBall*)
    double                              m_solid_angle[g_solid_volume_num_dynamicR];

    //by MyGridSES::compute_CB_rotation_matrix(MyContextBall*)
    double                              m_rotation_matrix[3][3];
    //this matrix is used to rotation the normalized pose to the current local reference frame with m_solid_dir as Z

    //by MyGridSES::compute_CB_context_rays(MyContextBall*)  --  "shell" is the solid part between two surfaces
    vector<CONTEXT_RAY_INT_TYPE>                m_context_rays;
    vector<CONTEXT_RAY_INT_TYPE>                m_context_shell_L2; //always dense
    vector<CONTEXT_RAY_INT_TYPE>                m_context_core_L2;  //always dense
    vector<CONTEXT_RAY_INT_TYPE>                m_context_inner;    //always dense

    vector<CONTEXT_RAY_INT_TYPE>                m_context_shell_1A; //always dense
    vector<CONTEXT_RAY_INT_TYPE>                m_context_shell_2A; //always dense
    vector<CONTEXT_RAY_INT_TYPE>                m_context_shell_3A; //always dense
    vector<CONTEXT_RAY_INT_TYPE>                m_context_shell_4A; //always dense

    vector<CONTEXT_RAY_INT_TYPE>                m_context_core_1A;  //always dense
    vector<CONTEXT_RAY_INT_TYPE>                m_context_core_2A;  //always dense
    vector<CONTEXT_RAY_INT_TYPE>                m_context_core_3A;  //always dense
    vector<CONTEXT_RAY_INT_TYPE>                m_context_core_4A;  //always dense

    vector<CONTEXT_RAY_INT_TYPE>                m_context_inner_2A;  //always dense

    vector<CONTEXT_RAY_INT_TYPE>                m_context_inner_core3A;  //always dense
    vector<CONTEXT_RAY_INT_TYPE>                m_context_core1A_shell3A;  //always dense
    //core_3A_inDB;
    //core_1A_inDB;
    //shell_1A_inDB;
    //shell_2A_inDB;
    //shell_3A_inDB;

    vector<CONTEXT_RAY_INT_TYPE>                m_context_core1Ashell3A;  //always dense


    //by MyGridSES::compute_CB_buried_vert(MyContextBall*)
    vector<MyBuriedVert*>               m_buried_verts; 

    unsigned long size();


    vector<double>  get_center_point();

    void sload(boost::archive::text_iarchive &la);

    void ssave(boost::archive::text_oarchive &la);

    MyContextBall(int pm_ball_id, vector<double>& pm_ball_center, vector<vector<double> >* pm_unitball_template);

    MyContextBall();

    void write_shell_vtk(vector<unsigned int>& pm_rays, vector<vector<double> >& pm_unitball, string out_fn);

    void write_vtk(string out_fn);


    void write_vtk_dockinfo(vector<vector<double> >& pm_unitball, string out_fn);

};

#endif
