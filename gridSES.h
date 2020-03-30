#ifndef _GRID_SES_H_
#define _GRID_SES_H_

#include <iostream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <list>
#include <queue>
#include <ext/algorithm>
#include <numeric>
#include <sys/timeb.h>
#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <unistd.h>

#include "contextBall.h"
#include "util.h"
#include "memusage.h"

using namespace std;

extern double g_grid_cell_len;
extern string g_dir_output;
extern int g_solid_volume_radius_selector;
extern vector<pair<double, pair<vector<vector<int> >, vector<vector<int> > > > > g_dynamic_radius;
extern string g_dir_database_DB;
extern string g_dir_database_Query;
extern string g_cb_suffix_DB;
extern string g_cb_suffix_Query;
extern unsigned long g_sizeof_grid;

const unsigned int RDL_SES              = 0x10000000u; //0001
const unsigned int RDL_INNER            = 0x20000000u; //0010
const unsigned int RDL_OUTER            = 0x30000000u; //0011 
const unsigned int RDL_SHELL_L1         = 0x40000000u; //0100
const unsigned int RDL_SHELL_L1_SURF    = 0x50000000u; //0101
const unsigned int RDL_SHELL_L2         = 0x60000000u; //0110
const unsigned int RDL_SHELL_L2_SURF    = 0x70000000u; //0111
const unsigned int RDL_CORE_L1          = 0x80000000u; //1000
const unsigned int RDL_CORE_L1_SURF     = 0x90000000u; //1001
const unsigned int RDL_CORE_L2          = 0xA0000000u; //1010
const unsigned int RDL_CORE_L2_SURF     = 0xB0000000u; //1011
const unsigned int RDL_MASK             = 0xF0000000u; //1111




typedef unsigned int    TID_VERT;
typedef unsigned int    TID_VERT_TYPE;
typedef unsigned int    TID_ATOM;


/*********** struct MyGridSES_SPD ***********
 * This class allocates a big grid, put the protein with SES in the grid, 
 * compute the layers, compute the critical points, computer context ball for each critical point.
 */
struct MyGridSES_SPD // SPD stands for sharp penetration detect
{ 

  // The follwoing 4 variables are initialized by constructor
    MySES*      m_surface;                          // The SES surface for which the grid is generated
    MyUnitBall* m_unitball;                         // The templates with uniform points on Unit sphere
    double      m_step_len;                         // grid side length, same length in x, y and z axis
    double      m_step_multiplier;                  // inverse of grid side length

    // the following 10 vaiables are initialized in the allocate_grid routine
    double      m_min_x;                             // Minimum x value of the protein
    double      m_min_y;                             // Minimum y value of the protein
    double      m_min_z;                             // Minimum z value of the protein
    double      m_max_x;                             // Maximum x value of the protein
    double      m_max_y;                             // Maximum y value of the protein
    double      m_max_z;                             // Maximum z value of the protein

    int         m_size_x;                            // number of grid cells along the x-axis
    int         m_size_y;                            // number of grid cells along the y-axis
    int         m_size_z;                            // number of grid cells along the z-axis

    vector<vector<vector<unsigned int> > >  m_grid;  // storing all the grid cell's value in marching cube algorith. (three vector for
                                                     // x, y and z direction. Now,
                                                     // An integer has 4 bytes, we use the 4 bytes to represent 4 numbers
                                                     // one bytes for the type of surface, and 3 bytes to represent the 
						     // distance to this cell's nearest SES grid cell
    
    // A context ball shall be placed at the critical point on the SES surface. The following
    // variables holds those context balls.
    vector<MyContextBall*>                  m_context_balls;
    // refined context balls for those context balls that are located on the center of big contact faces (convex faces)
    // some contact faces are very big. Using the center of the faces can not sample the local region 
    // very well. The approach adopted here to overcome this problem is to use more surface points on the 
    // same contact face as critical points. Then for each of such "refined" critical point, we generate
    // a context ball.
    vector<vector<MyContextBall*> >         m_context_balls_refine; // outer vector represent each face, inner vector represent
                                                                    // all the contexball generated for this face. Note that for
								    // concave only one ball is generated, but for convex, if the 
								    // face is big, multiple is generated.

    /*	Critical points are defined as the center of the SES faces. There are 3 kinds of SES faces:
     *  1) Contact faces (convex); 2) re-entrant faces (concave); and 3) toroidal faces (belt).
     *  For each SES face, its center is computed and a critical point is derived from it. 
     *  If the area of a contact face is very big, more points on that contact face will be used as
     *  "refined" critical points.
     */
    vector<vector<double> >                 m_critical_belt;       // belt critical points (inner vector fox x,y,z)
    vector<vector<double> >                 m_critical_concave;    // concave critical points
    vector<vector<double> >                 m_critical_convex;     // convex critical points
    vector<vector<double> >                 m_critical_points;     // all critical points

   // corresponding to m_context_balls_refine, each refined critical point has an entry here.
   // first number is the serialno of m_context_balls, second number is the serialno of refined.
   // For example, if we have 1000 contact faces, with ID 1,2,...1000. Say, face 1 has only 1
   // context ball, but face 2 has 4, then there will be the following entry:(1, 1), (2,1),(2,2)
   // (2,3) and (2,4) and so on.
    vector<pair<size_t, size_t> >           m_critical_convex_ids; // Read the comment above

    vector<double>                 m_critical_belt_area;           // area's of a critical belt
    vector<double>                 m_critical_concave_area;        // area's of a critical concave
    vector<double>                 m_critical_convex_area;         // area's of a critical convex
    vector<double>                 m_critical_points_area;         // area of all critical point

    vector<vector<vector<double> > >        m_mc_table;         // holds distance of each grid point from (0,0,0). 
                                                                // used for speeding up marching cube computation
                                                               
    size_t                                  m_mc_size;          // How many grid point along an axis. this value is
                                                                // maximum index of the above vector m_mc_table

    double                                  m_mc_localR; //=2.8


    // constructor routine
    MyGridSES_SPD(MySES* pm_ses, MyUnitBall* pm_ball);

    // destructor routine
    ~MyGridSES_SPD();

    // for a given point, pm_point, we like to know it's cell type. First, we check whether the point
    // is less-than min_axis_val or higher-than max_axis_val for all 3 axis. then it is an OUTER point,
    // otherwise, just return its value. However, the point's actual co-ordinates needs to be converted
    // to its grid co-ordinates.
    // Precondition: All the grid value is already filled up in m_grid[][][]
    unsigned int check_cell_type(vector<double>& pm_point);
    
    // the following mc_* functions are used to calculate the layers in a Marching-Cubes-like algorithm.
    // each cell in the grid is assigned an 4-byte unsigned integer. The first byte from MSB is used as Type marker,
    // the rest 3 bytes are used to record the x, y, z distance to its nearest SES cell, in terms of number
    // of cells along the x, y and z axis. 
    //inline size_t mc_get_x(unsigned int pm_label);
    size_t mc_get_x(unsigned int pm_label);

    //inline size_t mc_get_y(unsigned int pm_label);
    size_t mc_get_y(unsigned int pm_label);

    //inline size_t mc_get_z(unsigned int pm_label);
    size_t mc_get_z(unsigned int pm_label);

    // init the lookup table that is used by the MC algorithm. Basically this table is a cube. 
    // the center of the cube is (0, 0, 0). Each cell of this table is assign a double number which
    // is the real distance to the center.
    // Now, note that out cube size is a 4A X 4A X 4A. This is because we used MC algorithm to find
    // the CORE_L2 and SHELL_L2 layer. They are 4A distance from the SES. So, we shall start our marching
    // from the SES and shall march only atmost 4A distance along any direction. But, actual grid point
    // along each axis need to be calculated by the eq: 4/(grid_step_size) + 1
    void mc_init();

    // load the context balls from disk file archive. Actually in the first stage we make all different
    // context ball and save them in disk. Matching algorithm retrieve those from disk. Boost archiving
    // algorithm is used here
    void sload_context_balls(boost::archive::text_iarchive &la);

    // Save the context balls from disk file archive. Actually in the first stage we make all different
    // context ball and save them in disk. Matching algorithm retrieve those from disk. Boost archiving
    // algorithm is used here
    void ssave_context_balls(boost::archive::text_oarchive &la);

	
    // for debug purpose only	
    void write_context_ball_vtk_criticalPointBased(string pm_protein_name);

    // for debug purpose only	
    void write_context_ball_vtk_localverts(string pm_protein_name);

    // for debug purpose only	
    void write_grid_vtk();

    // for debug purpose only	
    void write_grid_vtk(unsigned int pm_type, string pm_fn);

    // for debug purpose only	
    void write_critical_points(string pm_fn);

    // for debug purpose only	
    void write_critical_points(vector<vector<double> >& pm_points, string pm_fn);

    // For each of the context ball, from the variable m_context_ball and m_context_ball_refine
    // it computes, solid_volume, solid_angle, rotation_matrix, context rays, context shell and
    // burried vertices. Just simply called the corresponding routines.
    void compute_context_balls_criticalPointBased();

    // this is the dockinfo version. this function is deprecated.
    // dockinfo is now calculated by a seperate class MyDockinfo
    void compute_context_balls_criticalPointBased_dockinfo();

    // Computing the critical points; note that, they only constitute of concave and convex regions
    // before this routine is called, the vector m_critical_convex and m_critical_concave neeeds to
    // be filled, it just combined those two type of critical points in one vector
    void compute_critical_points();

    // This routine fills the vector m_critical_convex, thus putting all the convex critical points
    // in the desired vector (convex critical points are computed from the contact faces
    void compute_critical_points_convex();

    // deprecated
    void compute_critical_points_convex_as_concave();

    // very similar to as the case of convex, but no refined points here, just the center of the
    // face is taken as critical point.
    void compute_critical_points_concave();

    // deprecated
    void print_critical_concave_vs_ses_vertex();

    // deprecated
    void print_critical_convex_vs_ses_vertex();

    // allocating memory for context balls
    void create_context_balls(vector<vector<double> >* pm_unitball_template);

    //compute solid volume and solid vector, solid volume is calculated by using a template spheric
    //grid and summing up the volume of grid cells that belongs to the surface.
    void compute_CB_solid_volume(MyContextBall* pm_cb);

    //compute solid angle
    void compute_CB_solid_angle(MyContextBall* pm_cb);

    //compute Context Ball's m_rotation_matrix
    void compute_CB_rotation_matrix(MyContextBall* pm_cb, int pm_sv);

    // compute Context Ball's m_context_rays 
    // and m_context_shell_L2 based on the same unitball_template
    void compute_CB_context_rays(MyContextBall* pm_cb);

    // deprecated
    void compute_CB_context_rays_dockinfo(MyContextBall* pm_cb);

    //like buried verts, shell_L2 should be always densely sampled
    void compute_CB_context_shell_L2(MyContextBall* pm_cb);

    void compute_CB_context_shell_L2_dockinfo(MyContextBall* pm_cb);

    // It computes the buried vertices
    void compute_CB_buried_verts(MyContextBall* pm_cb);

    void compute_CB_buried_verts_dockinfo(MyContextBall* pm_cb);

    void init(double pm_min_x, double pm_min_y, double pm_min_z, double pm_max_x, double pm_max_y, double pm_max_z);

private:
    void allocate_grid(double pm_min_x, double pm_min_y, double pm_min_z, double pm_max_x, double pm_max_y, double pm_max_z);

    void set_surface_cells(MySES* pm_surface, unsigned int pm_mark);

    void mc_cells(unsigned int pm_mark_init, unsigned int pm_mark, size_t pm_i, size_t pm_j, size_t pm_k);

    void mc_looper();

    //go from any atom center as seed, mark pm_overriden with pm_mark
    void set_inner_cells(unsigned int pm_mark, unsigned int pm_overriden);

    //check pm_check's neighbors, if pm_overriden then mark with pm_mark
    void set_shell_cells(unsigned int pm_check, unsigned int pm_mark, unsigned int pm_overriden);

    void mark_cells(int pm_ix, int pm_iy, int pm_iz, unsigned int pm_mark, unsigned int pm_overriden);

    void compute_trianglet(vector<vector<double> >& pm_triangle, double pm_threshold, set<vector<double> >& rtn_points);
};

struct MyGridSES_Dense_SPD: public MyGridSES_SPD {
    MyGridSES_Dense_SPD(MySES* pm_ses, MyUnitBall* pm_ball) : MyGridSES_SPD(pm_ses, pm_ball) {}

    void sload_CB(string pm_PDB_code);

    void ssave_CB(string pm_PDB_code);

    void compute_context_balls_criticalPointBased();

    void compute_context_balls_criticalPointBased_dockinfo();

    void init(double pm_min_x, double pm_min_y, double pm_min_z, double pm_max_x, double pm_max_y, double pm_max_z);
};



struct MyGridSES_Coarse_SPD: public MyGridSES_SPD {
    MyGridSES_Coarse_SPD(MySES* pm_ses, MyUnitBall* pm_ball) : MyGridSES_SPD(pm_ses, pm_ball) {}

    void sload_CB(string pm_PDB_code);

    void ssave_CB(string pm_PDB_code);

    void compute_context_balls_criticalPointBased();

    void compute_context_balls_criticalPointBased_dockinfo();

    void init(double pm_min_x, double pm_min_y, double pm_min_z, double pm_max_x, double pm_max_y, double pm_max_z);
};

#endif
