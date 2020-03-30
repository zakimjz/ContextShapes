/* 
*
	v1: 1) only compute the buried area in the layers that are used in the scoring function. 
		2) optimized by binary search.
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>
#include <ext/algorithm>
#include <numeric>
#include <sys/timeb.h>
#include <time.h>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <pthread.h>
#include <unistd.h>

#include "rbt_coarse.h"
#include "rbt_dense.h"
#include "rbt_pose.h"
#include "global.h"
#include "SES.h"
#include "contextBall.h"
#include "gridSES.h"
#include "memusage.h"
#include "MyDockinfo.h"
#include "MyMatchPair.h"
#include "MyCompare_SPD.h"
#include "typedefs.h"

using namespace std;

double weight_in_1_16bits [0x1u << 16] ;
double weight_in_2_16bits [0x1u << 16] ;


std::string g_dir_interface_ca; 


std::string rdl_trim(std::string& s, const std::string& drop = " ")
{
	std::string r=s.erase(s.find_last_not_of(drop)+1);
	return r.erase(0,r.find_first_not_of(drop));
}


void* mt_subtask_match(void* pm_subtask)
{
	MySubtask* tmp_subtask = (MySubtask*) pm_subtask;
	//tmp_subtask->run();
	tmp_subtask->m_compare->subtask_match(tmp_subtask);
	return NULL;
}


void task_build_database_SPD(int argc, char *argv[])
{
    string pdb_code_DB = argv[3];
    int    type_selector = atoi(argv[4]);

    string filename_pdb_DB      = pdb_code_DB + ".pdb";
    string filename_xyzrn_DB    = pdb_code_DB + ".xyzrn";
    string filename_vert_DB     = pdb_code_DB + ".vert";
    string filename_face_DB     = pdb_code_DB + ".face";

    //string filename_L1_xyzrn_DB    = pdb_code_DB + "_L1.xyzrn";
    //string filename_L1_vert_DB     = pdb_code_DB + "_L1.vert";
    //string filename_L1_face_DB     = pdb_code_DB + "_L1.face";

    double time_begin = 0, time_end = 0;
    struct timeb tp;

    ftime(&tp);
    time_begin = tp.time + (tp.millitm * 1.0) / 1000;

    compute_weight_in_16bits();
    //compute_bits_in_16bits ();

    cout << "*****init: unitball..." << endl;
    MyUnitBall myBall;
    myBall.init();

    g_sizeof_mt = myBall.size();

    ftime(&tp);
    time_end = tp.time + (tp.millitm * 1.0) / 1000;
    cout << "INFO: time to init unitball & weight_table " << time_end - time_begin << " seconds." << endl; 

    double minx, miny, minz, maxx, maxy, maxz;

    //***************************************
    //********compute protein_DB
    //***************************************
    ftime(&tp);
    time_begin = tp.time + (tp.millitm * 1.0) / 1000;

    cout << "*****init: SES_DB..." << endl;
    MySES      mySES_DB;
    mySES_DB.init(g_dir_pdb_DB, filename_pdb_DB, filename_xyzrn_DB, filename_vert_DB, filename_face_DB);
    mySES_DB.melt_patch_init();


    minx = mySES_DB.m_min_x - (4 + g_grid_cell_len);
    miny = mySES_DB.m_min_y - (4 + g_grid_cell_len);
    minz = mySES_DB.m_min_z - (4 + g_grid_cell_len);
    maxx = mySES_DB.m_max_x + (4 + g_grid_cell_len);
    maxy = mySES_DB.m_max_y + (4 + g_grid_cell_len);
    maxz = mySES_DB.m_max_z + (4 + g_grid_cell_len);


    cout << "*****GridSES_DB.init(...)..." << endl;

    if (type_selector == 1){
        MyGridSES_Dense_SPD myGridSES_DB(&mySES_DB, &myBall);
        myGridSES_DB.init(minx, miny, minz, maxx, maxy, maxz);
        myGridSES_DB.compute_critical_points_concave();
        myGridSES_DB.compute_context_balls_criticalPointBased();
        myGridSES_DB.ssave_CB(pdb_code_DB);

        g_numberof_cs_dense = myGridSES_DB.m_context_balls.size();
        g_sizeof_single_cs_dense  = myGridSES_DB.m_context_balls[0]->size();
        g_sizeof_cs_dense = g_numberof_cs_dense * g_sizeof_single_cs_dense;
    }
    else if (type_selector == 2){
        MyGridSES_Coarse_SPD myGridSES_DB(&mySES_DB, &myBall);
        myGridSES_DB.init(minx, miny, minz, maxx, maxy, maxz);
        myGridSES_DB.compute_critical_points_convex();
	myGridSES_DB.compute_context_balls_criticalPointBased();

        myGridSES_DB.ssave_CB(pdb_code_DB);

        g_numberof_cs_coarse = myGridSES_DB.m_context_balls.size();
        g_sizeof_single_cs_coarse  = myGridSES_DB.m_context_balls[0]->size();
        g_sizeof_cs_coarse = g_numberof_cs_coarse * g_sizeof_single_cs_coarse;
    }

    ftime(&tp);
    time_end = tp.time + (tp.millitm * 1.0) / 1000;
    cout << "INFO: time to build DB " << time_end - time_begin << " seconds." << endl; 
}

//keep
void task_dockinfo(int argc, char *argv[])
{
    string pdb_code_DB = argv[3];
    string pdb_code_Query = argv[4];

    string filename_pdb_DB      = pdb_code_DB + ".pdb";
    string filename_xyzrn_DB    = pdb_code_DB + ".xyzrn";
    string filename_vert_DB     = pdb_code_DB + ".vert";
    string filename_face_DB     = pdb_code_DB + ".face";

    string filename_pdb_Query   = pdb_code_Query + ".pdb";
    string filename_xyzrn_Query = pdb_code_Query + ".xyzrn";
    string filename_vert_Query  = pdb_code_Query + ".vert";
    string filename_face_Query  = pdb_code_Query + ".face";

    double time_begin = 0, time_end = 0;
    struct timeb tp;
	double minx, miny, minz, maxx, maxy, maxz;

    ftime(&tp);
    time_begin = tp.time + (tp.millitm * 1.0) / 1000;


    compute_weight_in_16bits();
    //test_weight_table_32();
    //return;
    //compute_bits_in_16bits ();

    cout << "*****init: unitball..." << endl;
    MyUnitBall myBall;
    myBall.init();

    ftime(&tp);
    time_end = tp.time + (tp.millitm * 1.0) / 1000;
    cout << "INFO: time to init unitball & weight_table " << time_end - time_begin << " seconds." << endl; 

    //***************************************
    //********compute protein_DB
    //***************************************
    ftime(&tp);
    time_begin = tp.time + (tp.millitm * 1.0) / 1000;

    cout << "*****init: SES_DB..." << endl;
    MySES      mySES_DB;
    mySES_DB.init(g_dir_pdb_DB, filename_pdb_DB, filename_xyzrn_DB, filename_vert_DB, filename_face_DB);
    mySES_DB.melt_patch_init();

    minx = mySES_DB.m_min_x - (4 + g_grid_cell_len);
    miny = mySES_DB.m_min_y - (4 + g_grid_cell_len);
    minz = mySES_DB.m_min_z - (4 + g_grid_cell_len);
    maxx = mySES_DB.m_max_x + (4 + g_grid_cell_len);
    maxy = mySES_DB.m_max_y + (4 + g_grid_cell_len);
    maxz = mySES_DB.m_max_z + (4 + g_grid_cell_len);

    cout << "*****GridSES_DB.init(...)..." << endl;
    MyGridSES_Dense_SPD myGridSES_DB(&mySES_DB, &myBall);
        myGridSES_DB.init(minx, miny, minz, maxx, maxy, maxz);
        myGridSES_DB.compute_critical_points_concave();
        myGridSES_DB.compute_context_balls_criticalPointBased_dockinfo();
    //myGridSES_DB.sload_CB(pdb_code_DB);


    ftime(&tp);
    time_end = tp.time + (tp.millitm * 1.0) / 1000;
    cout << "INFO: time to build DB " << time_end - time_begin << " seconds." << endl; 


    //***************************************
    //********compute protein_QUERY
    //***************************************
    ftime(&tp);
    time_begin = tp.time + (tp.millitm * 1.0) / 1000;

    cout << "*****init: SES_Query..." << endl;
    MySES      mySES_Query;
    mySES_Query.init(g_dir_pdb_Query, filename_pdb_Query, filename_xyzrn_Query, filename_vert_Query, filename_face_Query);
    mySES_Query.melt_patch_init();  

    minx = mySES_Query.m_min_x - (4 + g_grid_cell_len);
    miny = mySES_Query.m_min_y - (4 + g_grid_cell_len);
    minz = mySES_Query.m_min_z - (4 + g_grid_cell_len);
    maxx = mySES_Query.m_max_x + (4 + g_grid_cell_len);
    maxy = mySES_Query.m_max_y + (4 + g_grid_cell_len);
    maxz = mySES_Query.m_max_z + (4 + g_grid_cell_len);


    mySES_Query.populate_atom_interface_ca(&mySES_DB);
    mySES_DB.populate_atom_interface_ca(&mySES_Query);

    //mySES_Query.write_patch_vtk(g_dir_output + pdb_code_Query);

    cout << "*****GridSES_Query.init(...)..." << endl;
    MyGridSES_Coarse_SPD    myGridSES_Query(&mySES_Query, &myBall);
        myGridSES_Query.init(minx, miny, minz, maxx, maxy, maxz);
        myGridSES_Query.compute_critical_points_convex();
        myGridSES_Query.compute_context_balls_criticalPointBased_dockinfo();
	//myGridSES_Query.sload_CB(pdb_code_Query);

    ftime(&tp);
    time_end = tp.time + (tp.millitm * 1.0) / 1000;
    cout << "INFO: time to build Query " << time_end - time_begin << " seconds." << endl; 

    cout << "*****OK: dockinfo matching..." << endl;
    ftime(&tp);
    time_begin = tp.time + (tp.millitm * 1.0) / 1000;

    MyDockinfo   myDockinfo(&myGridSES_DB, &myGridSES_Query, &myBall);
    myDockinfo.match_dockinfo();

    ftime(&tp);
    time_end = tp.time + (tp.millitm * 1.0) / 1000;
    cout << "INFO: TOTAL time to match " << time_end - time_begin << " seconds." << endl; 

}


void task_match_sload_SPD(int argc, char *argv[])
{
    string pdb_code_DB = argv[3];
    string pdb_code_Query = argv[4];

    string filename_pdb_DB      = pdb_code_DB + ".pdb";
    string filename_xyzrn_DB    = pdb_code_DB + ".xyzrn";
    string filename_vert_DB     = pdb_code_DB + ".vert";
    string filename_face_DB     = pdb_code_DB + ".face";

    string filename_pdb_Query   = pdb_code_Query + ".pdb";
    string filename_xyzrn_Query = pdb_code_Query + ".xyzrn";
    string filename_vert_Query  = pdb_code_Query + ".vert";
    string filename_face_Query  = pdb_code_Query + ".face";

    double time_begin = 0, time_end = 0;
    struct timeb tp;

    ftime(&tp);
    time_begin = tp.time + (tp.millitm * 1.0) / 1000;

    compute_weight_in_16bits();
    //test_weight_table_32();
    //return;
    //compute_bits_in_16bits ();

    cout << "*****init: unitball..." << endl;
    MyUnitBall myBall;
    myBall.init();

    g_sizeof_mt = myBall.size();

    ftime(&tp);
    time_end = tp.time + (tp.millitm * 1.0) / 1000;
    cout << "INFO: time to init unitball & weight_table " << time_end - time_begin << " seconds." << endl; 

    //***************************************
    //********compute protein_DB
    //***************************************
    ftime(&tp);
    time_begin = tp.time + (tp.millitm * 1.0) / 1000;

    cout << "*****init: SES_DB..." << endl;
    MySES      mySES_DB;
    mySES_DB.init(g_dir_pdb_DB, filename_pdb_DB, filename_xyzrn_DB, filename_vert_DB, filename_face_DB);
    mySES_DB.melt_patch_init();

    cout << "*****GridSES_DB.init(...)..." << endl;
    MyGridSES_Dense_SPD myGridSES_DB(&mySES_DB, &myBall);
    //myGridSES_DB.init(minx, miny, minz, maxx, maxy, maxz);
    //myGridSES_DB.compute_critical_points_concave();
    //myGridSES_DB.compute_critical_points();
    //myGridSES_DB.compute_context_balls_criticalPointBased();
    myGridSES_DB.sload_CB(pdb_code_DB);

    g_numberof_cs_dense = myGridSES_DB.m_context_balls.size();
    g_sizeof_single_cs_dense  = myGridSES_DB.m_context_balls[0]->size();
    g_sizeof_cs_dense = g_numberof_cs_dense * g_sizeof_single_cs_dense;


    ftime(&tp);
    time_end = tp.time + (tp.millitm * 1.0) / 1000;
    cout << "INFO: time to build DB " << time_end - time_begin << " seconds." << endl; 


    //***************************************
    //********compute protein_QUERY
    //***************************************
    ftime(&tp);
    time_begin = tp.time + (tp.millitm * 1.0) / 1000;

    cout << "*****init: SES_Query..." << endl;
    MySES      mySES_Query;
    mySES_Query.init(g_dir_pdb_Query, filename_pdb_Query, filename_xyzrn_Query, filename_vert_Query, filename_face_Query);
    mySES_Query.melt_patch_init();  


    mySES_Query.populate_atom_interface_ca(&mySES_DB);
    mySES_DB.populate_atom_interface_ca(&mySES_Query);

    //mySES_Query.write_patch_vtk(g_dir_output + pdb_code_Query);

    cout << "*****GridSES_Query.init(...)..." << endl;
    MyGridSES_Coarse_SPD    myGridSES_Query(&mySES_Query, &myBall);
    //myGridSES_Query.init(minx, miny, minz, maxx, maxy, maxz);
    //myGridSES_Query.compute_critical_points_convex();
    //myGridSES_Query.compute_critical_points();
    //myGridSES_Query.compute_context_balls_criticalPointBased();
    myGridSES_Query.sload_CB(pdb_code_Query);

    g_numberof_cs_coarse = myGridSES_Query.m_context_balls.size();
    g_sizeof_single_cs_coarse  = myGridSES_Query.m_context_balls[0]->size();
    g_sizeof_cs_coarse = g_numberof_cs_coarse * g_sizeof_single_cs_coarse;



    ftime(&tp);
    time_end = tp.time + (tp.millitm * 1.0) / 1000;
    cout << "INFO: time to build Query " << time_end - time_begin << " seconds." << endl; 

    cout << "*****OK: matching..." << endl;
    ftime(&tp);
    time_begin = tp.time + (tp.millitm * 1.0) / 1000;

    MyCompare_SPD   myCompare(&myGridSES_DB, &myGridSES_Query, &myBall);
	if (g_task == 510){
		if (g_num_threads > 1){
			myCompare.match_SPD_mt();
		}
		else{
			myCompare.match_SPD();
		}
	}
	else if (g_task == 500){
		myCompare.match_SPD_v1_mt();
	}
	else if (g_task == 501){
		myCompare.match_SPD_1pair(atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));
	}
	else if (g_task == 509){
		myCompare.match_SPD_test();
	}

    //cout << "hit_non = " << myCompare.m_non_hit << endl;
    //cout << "hit_core_3A = " << myCompare.m_hit_core_3A << endl;
    //cout << "hit_core_1A = " << myCompare.m_hit_core_1A << endl;
    //cout << "hit_shell_1A = " << myCompare.m_hit_shell_1A << endl;
    //cout << "hit_shell_2A = " << myCompare.m_hit_shell_2A << endl;
    //cout << "hit_shell_3A = " << myCompare.m_hit_shell_3A << endl;

    ftime(&tp);
    time_end = tp.time + (tp.millitm * 1.0) / 1000;
    cout << "INFO: TOTAL time to match " << time_end - time_begin << " seconds." << endl; 

	ofstream out_time;
	out_time.open ("smatch_timings.log", ofstream::out | ofstream::app);
	out_time << g_match_pair_str << endl;
	out_time << "INFO: TOTAL time to match " << time_end - time_begin << " seconds." << endl;
	out_time.close();


}


void parse_param(string pm_param_filename)
{
    ifstream in_file(pm_param_filename.c_str());
    if (!in_file.is_open()){
        cerr << "ERROR: cannot open file " << pm_param_filename << endl;
        exit(1);
    }
    //g_msms_probe_radius = 1.4;
    //g_msms_density = 2.0;
    g_sample_region_radius = 10;
    g_grid_cell_len = 0.2;

	g_bov_filter_selector = 1;
//bov_filter_selector = 1
	g_bov_inInner2A_CORE_upper_bound = 25;
//bov_filter_selector = 2
	g_bov_inInner1A_upper_bound   = 100;
	g_bov_inInner2A_upper_bound   = 25;
	g_bov_inInner3A_upper_bound   = 10;
	g_bov_inInner4A_upper_bound   = 0;
	g_bov_inCore_upper_bound      = 0;


    g_dir_pdb_DB = "./dir_PDB/";
    g_dir_pdb_Query = "./dir_PDB/";
    g_dir_database_DB = "./dir_database/";
    g_dir_database_Query = "./dir_database/";

    g_dir_table = "./dir_table/";
    g_dir_output = "./dir_output/";

    g_solid_volume_min_radius = 4;
    //g_solid_volume_lower_bound = 0.7;
    //g_solid_volume_upper_bound = 1.0;
    g_solid_angle_lower_bound = 0.75;
    g_solid_angle_upper_bound = 1.05;
    g_solid_volume_radius_selector = 5;

    g_distance_true_pair = 2.8;
    g_rmsd_true_pair = 2.0;

    g_matching_table_selector = 40;
//    g_use_table_sort = 1;
    g_num_predicates = 4000;

    g_cb_suffix_DB = "_DB";
    g_cb_suffix_Query = "_Query";



    string str_line;
    //size_t i;

    //char line[1024];
    //line[1023] = '\0';

    while(!in_file.eof()){
        //in_file.getline(line, 1023);
        //str_line = line;
        getline(in_file, str_line);
        if (str_line.length() == 0) continue;
        if (str_line[0] == '#') continue;

        //if (str_line.substr(0, strlen("msms_probe_radius")) == "msms_probe_radius"){
        //    g_msms_probe_radius = atof(str_line.substr(str_line.find("=", strlen("msms_probe_radius")) + 1, str_line.length()).c_str());
        //    g_parameters_str << "#msms_probe_radius = " << g_msms_probe_radius << endl;
        //}
        //else if (str_line.substr(0, strlen("msms_density")) == "msms_density"){
        //    g_msms_density = atof(str_line.substr(str_line.find("=", strlen("msms_density")) + 1, str_line.length()).c_str());
        //    g_parameters_str << "#msms_density = " << g_msms_density << endl;
        //}
        //else 
		if (str_line.substr(0, strlen("sample_region_radius")) == "sample_region_radius"){
            g_sample_region_radius = atof(str_line.substr(str_line.find("=", strlen("sample_region_radius")) + 1, str_line.length()).c_str());
            g_parameters_str << "#sample_region_radius = " << g_sample_region_radius << endl;
        }
        else if (str_line.substr(0, strlen("grid_cell_len")) == "grid_cell_len"){
            g_grid_cell_len = atof(str_line.substr(str_line.find("=", strlen("grid_cell_len")) + 1, str_line.length()).c_str());
            g_parameters_str << "#grid_cell_len = " << g_grid_cell_len << endl;
        }
        else if (str_line.substr(0, strlen("solid_angle_lower_bound")) == "solid_angle_lower_bound"){
            g_solid_angle_lower_bound = atof(str_line.substr(str_line.find("=", strlen("solid_angle_lower_bound")) + 1, str_line.length()).c_str());
            g_parameters_str << "#solid_angle_lower_bound = " << g_solid_angle_lower_bound << endl;
        }
        else if (str_line.substr(0, strlen("solid_angle_upper_bound")) == "solid_angle_upper_bound"){
            g_solid_angle_upper_bound = atof(str_line.substr(str_line.find("=", strlen("solid_angle_upper_bound")) + 1, str_line.length()).c_str());
            g_parameters_str << "#solid_angle_upper_bound = " << g_solid_angle_upper_bound << endl;
        }
		else if (str_line.substr(0, strlen("bov_filter_selector")) == "bov_filter_selector"){
            g_bov_filter_selector = atoi(str_line.substr(str_line.find("=", strlen("bov_filter_selector")) + 1, str_line.length()).c_str());
            g_parameters_str << "#bov_filter_selector = " << g_bov_filter_selector << endl;
        }
		else if (str_line.substr(0, strlen("bov_inInner2A_CORE_upper_bound")) == "bov_inInner2A_CORE_upper_bound"){
            g_bov_inInner2A_CORE_upper_bound = atof(str_line.substr(str_line.find("=", strlen("bov_inInner2A_CORE_upper_bound")) + 1, str_line.length()).c_str());
            g_parameters_str << "#bov_inInner2A_CORE_upper_bound = " << g_bov_inInner2A_CORE_upper_bound << endl;
        }
		else if (str_line.substr(0, strlen("bov_inInner1A_upper_bound")) == "bov_inInner1A_upper_bound"){
            g_bov_inInner1A_upper_bound = atof(str_line.substr(str_line.find("=", strlen("bov_inInner1A_upper_bound")) + 1, str_line.length()).c_str());
            g_parameters_str << "#bov_inInner1A_upper_bound = " << g_bov_inInner1A_upper_bound << endl;
        }
		else if (str_line.substr(0, strlen("bov_inInner2A_upper_bound")) == "bov_inInner2A_upper_bound"){
            g_bov_inInner2A_upper_bound = atof(str_line.substr(str_line.find("=", strlen("bov_inInner2A_upper_bound")) + 1, str_line.length()).c_str());
            g_parameters_str << "#bov_inInner2A_upper_bound = " << g_bov_inInner2A_upper_bound << endl;
        }
		else if (str_line.substr(0, strlen("bov_inInner3A_upper_bound")) == "bov_inInner3A_upper_bound"){
            g_bov_inInner3A_upper_bound = atof(str_line.substr(str_line.find("=", strlen("bov_inInner3A_upper_bound")) + 1, str_line.length()).c_str());
            g_parameters_str << "#bov_inInner3A_upper_bound = " << g_bov_inInner3A_upper_bound << endl;
        }
		else if (str_line.substr(0, strlen("bov_inInner4A_upper_bound")) == "bov_inInner4A_upper_bound"){
            g_bov_inInner4A_upper_bound = atof(str_line.substr(str_line.find("=", strlen("bov_inInner4A_upper_bound")) + 1, str_line.length()).c_str());
            g_parameters_str << "#bov_inInner4A_upper_bound = " << g_bov_inInner4A_upper_bound << endl;
        }

		else if (str_line.substr(0, strlen("bov_inCore_upper_bound")) == "bov_inCore_upper_bound"){
            g_bov_inCore_upper_bound = atof(str_line.substr(str_line.find("=", strlen("bov_inCore_upper_bound")) + 1, str_line.length()).c_str());
            g_parameters_str << "#bov_inCore_upper_bound = " << g_bov_inCore_upper_bound << endl;
        }

        else if (str_line.substr(0, strlen("solid_volume_min_radius")) == "solid_volume_min_radius"){
            g_solid_volume_min_radius = atof(str_line.substr(str_line.find("=", strlen("solid_volume_min_radius")) + 1, str_line.length()).c_str());
            g_parameters_str << "#solid_volume_min_radius = " << g_solid_volume_min_radius << endl;
        }
        //else if (str_line.substr(0, strlen("solid_volume_lower_bound")) == "solid_volume_lower_bound"){
        //    g_solid_volume_lower_bound = atof(str_line.substr(str_line.find("=", strlen("solid_volume_lower_bound")) + 1, str_line.length()).c_str());
        //    g_parameters_str << "#solid_volume_lower_bound = " << g_solid_volume_lower_bound << endl;
        //}
        //else if (str_line.substr(0, strlen("solid_volume_upper_bound")) == "solid_volume_upper_bound"){
        //    g_solid_volume_upper_bound = atof(str_line.substr(str_line.find("=", strlen("solid_volume_upper_bound")) + 1, str_line.length()).c_str());
        //    g_parameters_str << "#solid_volume_upper_bound = " << g_solid_volume_upper_bound << endl;
        //}
        else if (str_line.substr(0, strlen("distance_true_pair")) == "distance_true_pair"){
            g_distance_true_pair = atof(str_line.substr(str_line.find("=", strlen("distance_true_pair")) + 1, str_line.length()).c_str());
            g_parameters_str << "#distance_true_pair = " << g_distance_true_pair << endl;
        }
        else if (str_line.substr(0, strlen("rmsd_true_pair")) == "rmsd_true_pair"){
            g_rmsd_true_pair = atof(str_line.substr(str_line.find("=", strlen("rmsd_true_pair")) + 1, str_line.length()).c_str());
            g_parameters_str << "#rmsd_true_pair = " << g_rmsd_true_pair << endl;
        }
        //else if (str_line.substr(0, strlen("task")) == "task"){
        //    g_task = atoi(str_line.substr(str_line.find("=", strlen("task")) + 1, str_line.length()).c_str());
        //    g_parameters_str << "#task = " << g_task << endl;
        //}
        //else if (str_line.substr(0, strlen("use_table_sort")) == "use_table_sort"){
        //    g_use_table_sort = atoi(str_line.substr(str_line.find("=", strlen("use_table_sort")) + 1, str_line.length()).c_str());
        //    g_parameters_str << "#use_table_sort = " << g_use_table_sort << endl;
        //}
        else if (str_line.substr(0, strlen("matching_table_selector")) == "matching_table_selector"){
            g_matching_table_selector = atoi(str_line.substr(str_line.find("=", strlen("matching_table_selector")) + 1, str_line.length()).c_str());
            g_parameters_str << "#matching_table_selector = " << g_matching_table_selector << endl;
        }
        else if (str_line.substr(0, strlen("solid_volume_radius_selector")) == "solid_volume_radius_selector"){
            g_solid_volume_radius_selector = atoi(str_line.substr(str_line.find("=", strlen("solid_volume_radius_selector")) + 1, str_line.length()).c_str());
            g_parameters_str << "#solid_volume_radius_selector = " << g_solid_volume_radius_selector << endl;
        }
        else if (str_line.substr(0, strlen("num_predicates")) == "num_predicates"){
            g_num_predicates = atoi(str_line.substr(str_line.find("=", strlen("num_predicates")) + 1, str_line.length()).c_str());
            g_parameters_str << "#num_predicates = " << g_num_predicates << endl;
        }
        else if (str_line.substr(0, strlen("num_threads")) == "num_threads"){
            g_num_threads = atoi(str_line.substr(str_line.find("=", strlen("num_threads")) + 1, str_line.length()).c_str());
            g_parameters_str << "#num_threads = " << g_num_threads << endl;
        }

		else if (str_line.substr(0, strlen("cb_suffix_Query")) == "cb_suffix_Query"){
            g_cb_suffix_Query = str_line.substr(str_line.find("=", strlen("cb_suffix_Query")) + 1);
			g_cb_suffix_Query = rdl_trim(g_cb_suffix_Query);
			g_parameters_str << "#cb_suffix_Query = " << g_cb_suffix_Query << endl;
        } 
        else if (str_line.substr(0, strlen("cb_suffix_DB")) == "cb_suffix_DB"){
            g_cb_suffix_DB = str_line.substr(str_line.find("=", strlen("cb_suffix_DB")) + 1);
			g_cb_suffix_DB = rdl_trim(g_cb_suffix_DB);
            g_parameters_str << "#cb_suffix_DB = " << g_cb_suffix_DB << endl;
        } 

		//else if (str_line.substr(0, strlen("dir_input_other")) == "dir_input_other"){
		//	g_dir_input_other = str_line.substr(str_line.find("=", strlen("dir_input_other")) + 1);
		//	g_dir_input_other = rdl_trim(g_dir_input_other);
		//	if (g_dir_input_other.length() > 0){
		//		if (g_dir_input_other[g_dir_input_other.length() - 1] != '/')
		//			g_dir_input_other += "/";
		//	}
		//	g_parameters_str << "#dir_input_other = " << g_dir_input_other << endl;
		//} 
        else if (str_line.substr(0, strlen("dir_pdb_DB")) == "dir_pdb_DB"){
            g_dir_pdb_DB = str_line.substr(str_line.find("=", strlen("dir_pdb_DB")) + 1);
            g_dir_pdb_DB = rdl_trim(g_dir_pdb_DB);
			if (g_dir_pdb_DB.length() > 0){
				if (g_dir_pdb_DB[g_dir_pdb_DB.length() - 1] != '/')
					g_dir_pdb_DB += "/";
			}

			g_parameters_str << "#dir_pdb_DB = " << g_dir_pdb_DB << endl;
        } 
        else if (str_line.substr(0, strlen("dir_pdb_Query")) == "dir_pdb_Query"){
            g_dir_pdb_Query = str_line.substr(str_line.find("=", strlen("dir_pdb_Query")) + 1);
			g_dir_pdb_Query = rdl_trim(g_dir_pdb_Query);
			if (g_dir_pdb_Query.length() > 0){
				if (g_dir_pdb_Query[g_dir_pdb_Query.length() - 1] != '/')
					g_dir_pdb_Query += "/";
			}

			g_parameters_str << "#dir_pdb_Query = " << g_dir_pdb_Query << endl;
        } 
        else if (str_line.substr(0, strlen("dir_table")) == "dir_table"){
            g_dir_table = str_line.substr(str_line.find("=", strlen("dir_table")) + 1);
			g_dir_table = rdl_trim(g_dir_table);
			if (g_dir_table.length() > 0){
				if (g_dir_table[g_dir_table.length() - 1] != '/')
					g_dir_table += "/";
			}

			g_parameters_str << "#dir_table = " << g_dir_table << endl;
        } 
        else if (str_line.substr(0, strlen("dir_database_DB")) == "dir_database_DB"){
            g_dir_database_DB = str_line.substr(str_line.find("=", strlen("dir_database_DB")) + 1);
			g_dir_database_DB = rdl_trim(g_dir_database_DB);
			if (g_dir_database_DB.length() > 0){
				if (g_dir_database_DB[g_dir_database_DB.length() - 1] != '/')
					g_dir_database_DB += "/";
			}

			g_parameters_str << "#dir_database_DB = " << g_dir_database_DB << endl;
        } 
        else if (str_line.substr(0, strlen("dir_database_Query")) == "dir_database_Query"){
            g_dir_database_Query = str_line.substr(str_line.find("=", strlen("dir_database_Query")) + 1);
			g_dir_database_Query = rdl_trim(g_dir_database_Query);
			if (g_dir_database_Query.length() > 0){
				if (g_dir_database_Query[g_dir_database_Query.length() - 1] != '/')
					g_dir_database_Query += "/";
			}

			g_parameters_str << "#dir_database_Query = " << g_dir_database_Query << endl;
        } 

        else if (str_line.substr(0, strlen("dir_output")) == "dir_output"){
            g_dir_output = str_line.substr(str_line.find("=", strlen("dir_output")) + 1);
			g_dir_output = rdl_trim(g_dir_output);
			if (g_dir_output.length() > 0){
				if (g_dir_output[g_dir_output.length() - 1] != '/')
					g_dir_output += "/";
			}

			g_parameters_str << "#dir_output = " << g_dir_output << endl;
        } 
        else if (str_line.substr(0, strlen("dir_interface_ca")) == "dir_interface_ca"){
            g_dir_interface_ca = str_line.substr(str_line.find("=", strlen("dir_interface_ca")) + 1);
			g_dir_interface_ca = rdl_trim(g_dir_interface_ca);
			if (g_dir_interface_ca.length() > 0){
				if (g_dir_interface_ca[g_dir_interface_ca.length() - 1] != '/')
					g_dir_interface_ca += "/";
			}
            
			g_parameters_str << "#dir_interface_ca = " << g_dir_interface_ca << endl;
        } 
    }//end of while
    //cout << g_parameters_str.str();
}//end of function

int main(int argc, char *argv[])
{
	if (argc < 5){
		cerr << "USAGE: " << argv[0] << " <param_filename> <taskid> <PDB_codename> <1 or 2 to compute CB | PDB_codename for ligand>" << endl;
		return -1;
	}

    string param_fn = argv[1];
    parse_param(param_fn);
    compute_dynamic_volume_template();

	g_task = atoi(argv[2]);
    cout << "===^===" << argv[3] << " vs " << argv[4] << "===^===" << endl;

    g_match_pair_str = argv[3];
    g_match_pair_str += "_vs_";
    g_match_pair_str += argv[4];

    if (g_task == 100){
        task_dockinfo(argc, argv);
    }
    else if (g_task == 410){
        task_build_database_SPD(argc, argv);
    }
    else if (g_task >= 500){
        task_match_sload_SPD(argc, argv);
    }
    g_vmpeak = get_memusage_vmpeak();

    //cout << "g_vmpeak = " << (g_vmpeak * 1.0)/1024.0 << " MB"<< endl;
    //cout << "g_sizeof_mt = " << g_sizeof_mt/(1024.0*1024.0) << " MB" << endl;
    //cout << "g_sizeof_grid = " << g_sizeof_grid/(1024.0*1024.0) << " MB" << endl;
    //cout << "g_sizeof_single_cs_coarse = " << g_sizeof_single_cs_coarse/(1024.0*1024.0) << " MB" << endl;
    //cout << "g_sizeof_cs_coarse = " << g_sizeof_cs_coarse/(1024.0*1024.0) << " MB" << endl;
    //cout << "g_numberof_cs_coarse = " << g_numberof_cs_coarse << endl;
    //cout << "g_sizeof_single_cs_dense = " << g_sizeof_single_cs_dense/(1024.0*1024.0) << " MB"  << endl;
    //cout << "g_sizeof_cs_dense = " << g_sizeof_cs_dense/(1024.0*1024.0) << " MB"  << endl;
    //cout << "g_numberof_cs_dense = " << g_numberof_cs_dense << endl;

    return 0;
}
