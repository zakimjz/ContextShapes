#ifndef _SES_H_
#define _SES_H_

#include <set>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

/*
//index bbox of surface triangles
#include "CGAL/Cartesian.h"
#include "CGAL/Segment_tree_k.h"
#include "CGAL/Range_segment_tree_traits.h"

//used to index triangle, to find closest pairs
#include "CGAL/basic.h"
#include "CGAL/Point_3.h"
#include "CGAL/Range_tree_k.h"
*/

#include "typedefs.h"
using namespace std;

/*
typedef unsigned int    TID_VERT;
typedef unsigned int    TID_VERT_TYPE;
typedef unsigned int    TID_ATOM;
typedef unsigned int  CONTEXT_RAY_INT_TYPE;
*/

const unsigned int CONTEXT_RAY_INT_SIZE = sizeof(CONTEXT_RAY_INT_TYPE) * 8;
/*

typedef CGAL::Cartesian<double> UINT_Representation;
typedef CGAL::Range_tree_map_traits_3<UINT_Representation, unsigned int> UINT_Traits;
typedef CGAL::Range_tree_3<UINT_Traits> UINT_Range_tree_3_type;
typedef UINT_Traits::Key UINT_Key;
typedef UINT_Traits::Pure_key UINT_Pure_key;
typedef UINT_Traits::Interval UINT_Interval;
*/

typedef CGAL::Cartesian<double> PTR_Representation;
typedef CGAL::Range_tree_map_traits_3<PTR_Representation, void*> PTR_Traits;
typedef CGAL::Range_tree_3<PTR_Traits> PTR_Range_tree_3_type;
typedef PTR_Traits::Key PTR_Key;
typedef PTR_Traits::Pure_key PTR_Pure_key;
typedef PTR_Traits::Interval PTR_Interval;


extern std::string g_dir_interface_ca;

extern double compute_triangle_area(vector<double>&, vector<double>&, vector<double>&);

template <class FirstType, class SecondType>
struct less_than_univ{
    bool operator()(pair<FirstType, SecondType> pm_A, pair<FirstType, SecondType> pm_B)
    {   
        if (pm_A.first < pm_B.first) return true;
        else return false;
    }   
}; 
/**
 * This class represents an atom of a molecule for the molecular docking
 * 
 */
struct MyAtom {
    TID_ATOM        m_id;    // ID for the atom
    double          m_x;     // X co-ordinate of the center
    double          m_y;     // Y co-ordinate of the center
    double          m_z;     // Z co-ordinate of the center
    double          m_r;     // radius of the molecule

    bool      m_is_ca;       // Is it a C-Alpha atom? applicable only for protein molecule

    //TID_ATOM        m_id_in_pdb;   // molecule ID in pdb file

    //vector<unsigned int>   m_local_verts;

    // A constructor, that initializes the above variables, by default an
    // atom is NOT a C-alpha atom.
    MyAtom(TID_ATOM pm_id, double pm_x, double pm_y, double pm_z, double pm_r) 
    {
        m_id = pm_id;
        m_x = pm_x;
        m_y = pm_y;
        m_z = pm_z;
        m_r = pm_r;
        m_is_ca = false;
        //m_local_verts.clear();
    }

    // This returns the center as a vector<double>
    vector<double> point()
    {
        vector<double> rtn_center(3, 0);
        rtn_center[0] = m_x;
        rtn_center[1] = m_y;
        rtn_center[2] = m_z;
        return rtn_center;
    }

    inline double  x() {return m_x;}
    inline double  y() {return m_y;}
    inline double  z() {return m_z;}
    inline double  r() {return m_r;}
};

/** 
 * Represent a point on the molecular surface
 */
struct MyVert {
    TID_VERT        m_id;         // An id for the vertex point.
    TID_VERT_TYPE   m_vert_type;  // Type of the vertex point, 
    TID_ATOM        m_atom_id;    // This vertex point belongs to some
                                  // atom in the molecule, this variable 
				  // is refereing to that atom.
				  //
    // X, Y and Z co-ordinate of the point
    double      m_x;
    double      m_y;
    double      m_z;

    // scattered points are used to represent a continious surface. In
    // this representation every point has an area, which is like a voronoi cell
    // Although, geometrically a point has a zero area, our vertex point will have
    // a finite non-zero area.
    double      m_area;

    vector<double> point()   // returning the vertex point as a vector<double>
    {
        vector<double> p(3, 0);
        p[0] = m_x;
        p[1] = m_y;
        p[2] = m_z;

        return p;
    }

    inline double  x() 
    {
        return m_x;
    }

    inline double  y() 
    {
        return m_y;
    }

    inline double  z() 
    {
        return m_z;
    }

    // increase the area by pm_delta amount
    inline void increase_area(double pm_delta)
    {
        m_area += pm_delta;
    }

    // A typical constructor, area value is initialized to zero.
    MyVert(TID_VERT id, TID_VERT_TYPE type, TID_ATOM atom_id, double x, double y, double z) {
        m_id = id;
        m_vert_type = type;
        m_atom_id = atom_id;
        m_x = x;
        m_y = y;
        m_z = z;

        m_area = 0;
    }
};

/** This struct represent a surface patch as obtained from the MSMS program.
 *  Surface patches are triangulated by the MSMS program, so a patch is basically
 *  a set of triangles that tile the patch.
 *  Used in conjunction with MySES_Raw struct to represent a molecular surface
 *
 */

struct MySES_Patch{
    int                         m_patch_type;    // Type of patch, belt, convex, concave
    vector<vector<TID_VERT> >   m_triangles;     // Inner vector is the ID of the vertex
                                                 // that makes up a triangle, outer vector
						 // is for the all the triangles that made 
						 // this patch.
    set<TID_VERT>               m_verts;         // Set of all the vertices, that belongs to
    			                         // this patch.
    double                      m_area;          // Area of this patch

    vector<MySES_Patch*>        m_nb_belt;       // vector of all neighboring belt patches of this patch
    vector<MySES_Patch*>        m_nb_concave;    // vector of all neighboring concave patches of this patch.
    vector<MySES_Patch*>        m_nb_convex;     // vector of all neighboring convex patches of this patch.

    // A constructor that basically initializing all the values to 0
    MySES_Patch()            
    {
        m_patch_type = 0;          // we don't know what kind of patch is this, at this moment, so initialized with 0
	                           // but note that is will use type = 1 (belt), 2(concave) and 3(convex)
        m_area = 0;                // Area of this patch

        m_triangles.clear();       // All the triangles that made up the patch.
        m_verts.clear();             
        m_nb_belt.clear();
        m_nb_concave.clear();
        m_nb_convex.clear();
    }
};

/** This struct represents the data stores of a molecular surface that will be 
 *  extended by the class MySES. It loads the data from the files obtained from 
 *  the MSMS program  It Uses the MySES_Patch class to represent a surface patch. 
 *  An entire surface is basically a collection of surface patches.
 *
 *  Input: 
 *  Besides the pdb file of the protein, all the files that are generated 
 *  from the MSMS program, like vertex file, face file, xyzrn file will be used as 
 *  input.
 *
 *  Significant Output:
 *  m_patch_all, which is a vector of MySES_Patch*. It stores all the patches, that
 *  constitute this molecular surface.
 *
 *  Other Important Methods:
 *  An surface can be viewed using vtk program; this class also has routines to write 
 *  a vtk file for the surface.
 *  Author: Zujun Shentu
 *  Commented by: Mohammad Hasan
*/

struct MySES_Raw {
public:
    string                  m_dir_pdb;         // Directory that contain the pdb file
    string                  m_filename_pdb;    // Name of the pdb file
    string                  m_filename_xyzrn;  // Name of the xyzrn file (MSMS generate this file from pdb file)
    string                  m_filename_vert;   // Name of the vert file (MSMS generate this file from pdb file)
    string                  m_filename_face;   // Name of the face file (MSMS generate this file from pdb file)

    double                  m_density;         // Density of the triangulation, MSMS argument (my guess)
    double                  m_probe_radius;    // Radius of the probe used to generate the molecular surface, MSMS
                                               // use this as arguments
    unsigned int            m_num_vertex;      // Total number of vertices in the surfaces (a vertex is a corner
                                               // of the triangle)
    unsigned int            m_num_triangle;    // Total number of triangles.
    unsigned int            m_num_atom;        // Total number of atom in the molecular surface (my guess)

    //the following fields are from .pdb file
    vector<vector<double> > m_atom_xyz_ca;       // xyz co-ordinates of the C-alpha atom, inner vector
    				                 // is for a point's x,y,z co-ordinate, outer vector
						 // for many C-alpha atoms.
    vector<unsigned int >   m_atom_index_ca;     // index number of all the C-alpha atoms

    //the following fields are from .xyzrn file
    vector<vector<double> > m_atom_xyz;         //c1 ~ c3, coordinates xyz
    vector<double>          m_atom_radius;      //c4
    vector<int>             m_atom_c5;          //c5    
    vector<string>          m_atom_name;        //c6

    //the following fields are from .vert file
    vector<vector<double> > m_vert_xyz;         //c1 ~ c3, coordinates xyz
    vector<vector<double> > m_vert_normal;      //c4 ~ c6
    vector<int>             m_vert_c7;          //c7        
    vector<int>             m_vert_atom_id;     //c8, closest atom id
    vector<int>             m_vert_type;        //c9, 1 for belt, 2 for concave, 3 for convex
    vector<string>          m_vert_atom_name;   //c10, closest atom name

    //the following fields are from .face file
    vector<vector<TID_VERT> >   m_face_triangle; //c1 ~ c3
    vector<int>                 m_face_type;     //c4, 1 for belt, 2 for concave, 3 for convex
    vector<int>                 m_face_patch_id; //c5 (id of the triangular face)


    //derived from the coordinates for vertices (should be used to get a rectangular bounding box around the surface)
    double               m_min_x;                // minimum x co-ordinate
    double               m_max_x;                // maximum x co-ordinates
    double               m_min_y;                // minimum y co-ordinate
    double               m_max_y;                // maximum y co-ordinates
    double               m_min_z;                // minimum z co-ordinate
    double               m_max_z;                // maximum z co-ordinates

    //raw data from MSMS .face file
    map<unsigned int, set<TID_VERT> >   m_patch_belt;    // triangle id --> set of vertex-id of all belt type patches
    map<unsigned int, set<TID_VERT> >   m_patch_concave; // triangle id --> set of vertex-id of all concave type patches
    map<unsigned int, set<TID_VERT> >   m_patch_convex;  // triangle id --> set of vertex-id of all convex type patches

    vector<MySES_Patch*>                m_patch_all;     // All the patches that made up this surface

    // This routine finds a set of belt patch with two (or one) concave patches, all has areas below the min_area as defined
    // below. then it merge those patches and mark it as a concave patch. all the neighbor lists are updated as 
    // required.
    void melt_patch_belt();
    void melt_patch_init();

    // This routine writes the patch_info of a surface, in 3 files one for each type of surface
    void write_patch_vtk(string pm_filename);

    void write_patch_vtk(string pm_filename, int pm_patch_type);

    void write_patch_vtk_single(MySES_Patch* pm_patch, string pm_filename);


    void write_vertex_vtk(string pm_filename);


    void write_vertex_vtk(string pm_filename, int pm_vert_type);

    void init(string pm_dir_pdb, string pm_filename_pdb, string pm_filename_atom, string pm_filename_vert, string pm_filename_face);

private:

    // based on triangle info
    void read_patch();



    void read_face();


    void read_vert();

    void read_pdb_ca();

    void read_xyzrn();

    void get_bbox();

};

/**
 * This class represent a molecular surface. It is extended from
 * MySES_Raw. Note that MySES_Raw will remember all low-level-details of the
 * surface, like all its faces, triangles, vertices etc.
 */

struct MySES: public MySES_Raw
{
    vector<MyAtom*>             m_atom;             // All the atoms that constitute the molecule
    vector<MyVert*>             m_vert;             // All the verticess
    vector<vector<TID_VERT> >   m_triangle;         // All the triangles, inner vector represent
                                                    // 3 vertex-id that made up that triangle

    vector<vector<double> >   m_atom_interface_ca;  // C-Alpha atoms that make the interface region while docking.
                                                    // inner vector represent the center co-ordinate of those atoms
						    // (my guess)
    vector<double>        m_atom_radius_ca;         // Radius of the C-Alpha atoms

    //used to find NN for points between SIS and SES
    UINT_Range_tree_3_type*     m_index_vert;       // This type is taken from CGAL for solving the range query
    UINT_Range_tree_3_type*     m_index_atom;


    // Given another MySES*, this routine defines the interface C-Alpha
    void populate_atom_interface_ca(MySES* pm_other_ses);

    // If you don't want to calculate interface C-Alpha all the time, you can load them from a file
    // This routine load interface C-Alpha from a file
    void populate_atom_interface_ca();

    // This routine populate the vertices from the variable m_vert_xyz, m_vert_type defined in MySES_Raw
    void populate_vert();

    // This routine populate the atoms from the variable m_atom_xyz[][], m_atom_radius, defined in MySES_Raw
    void populate_atom();

    // This routine populate the triangles from the variable m_face_triangle defined in MySES_Raw
    void populate_triangle();

    // A destructor
    ~MySES();

    // I guess, this is indexing all the vertex point, using some fast data structure
    void build_index_vert();

    // I guess, this is indexing all the center of the atoms, using some fast data structure
    void build_index_atom();

    // For a given point(1st arg, pm_point) return all the vertex_id as a vector (3rd arg, pm_out),
    // which are within the distance of pm_local_r, from pm_point
    // CGAL, window_query is used. but window query works for rectangular region, so need to do a
    // pruning at the end. (used for atom center)
    void query_local_atom(const vector<double>&  pm_point, double pm_local_r, vector<unsigned int>& pm_out);

    // For a given point(1st arg, pm_point) return all the vertex_id as a vector (3rd arg, pm_out),
    // which are within the distance of pm_local_r, from pm_point
    // CGAL, window_query is used. but window query works for rectangular region, so need to do a
    // pruning at the end. (used for vertices)
    void query_local(const vector<double>&  pm_point, double pm_local_r, vector<unsigned int>& pm_out);

    // Retrun the nearest neighbor of a point pm_point (used for vertices)
    unsigned int query_nn(const vector<double>& pm_point);

    // Initializing routine for MySES. It first calls MySES_Raw::init, then it calls all the populate
    // routines. Take all the input file names as argument
    void init(string pm_dir_pdb, string pm_filename_pdb, string pm_filename_atom, string pm_filename_vert, string pm_filename_face);

};

/**
 * This class represets a burried vertex
 */
struct MyBuriedVert{
    vector<double>  m_xyz;
    double          m_dist;     //distance to the center
    unsigned int  m_ray_id;   //dense or coarse ray_id
    CONTEXT_RAY_INT_TYPE    m_vert_mask;
    TID_VERT        m_vert_id;
    double          m_area;

    vector<double>  get_point()
    {
        return m_xyz;
    }

    void sload(boost::archive::text_iarchive &la)
    {
        la & m_xyz;
        la & m_dist;
        la & m_ray_id;
        la & m_vert_mask; 
        la & m_vert_id;
        la & m_area;
    }

    void ssave(boost::archive::text_oarchive &la)
    {
        la & m_xyz;
        la & m_dist;
        la & m_ray_id;
        la & m_vert_mask; 
        la & m_vert_id;
        la & m_area;
    }

};

#endif
