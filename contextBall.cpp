#include "contextBall.h"


using namespace std;

MyBall :: ~MyBall() {
  delete m_index;
}

MyBall :: MyBall() {
  m_index = NULL;
  m_size = sizeof(unsigned long) + sizeof(double) + sizeof(UINT_Range_tree_3_type*);
}

/*
 * Finds all the points of the ball that are within pm_r distance from pm_point
 * The 1st parameter is a point (vector is used to store x, y and z co-ordinate of the point
 * The second parameter is the distance threshold
 * The third parameter is output parameter, it stores the ID of all the points on ball that
 * is within pm_r distance from the point pm_point
 * CGAL range query routine is used.
 */
void
MyBall :: query_local(const vector<double>&  pm_point, double pm_r, vector<unsigned int>& pm_out) {
  pm_out.clear();

  vector<UINT_Key> OutputList;

  // defining a bounding box around the point pm_point, where any point that is within
  // pm_r distant from pm_point is within the box
  UINT_Pure_key a = UINT_Pure_key(pm_point[0] - pm_r, pm_point[1] - pm_r, pm_point[2] - pm_r);
  UINT_Pure_key b = UINT_Pure_key(pm_point[0] + pm_r, pm_point[1] + pm_r, pm_point[2] + pm_r);
  UINT_Interval win = UINT_Interval(UINT_Key(a, 0), UINT_Key(b, 0));

  // using the range-tree of the ball (m_index) to find those points that are within the bounding box
  m_index->window_query(win, back_inserter(OutputList));

  if (OutputList.size() == 0){
    cerr << "ERROR: OutputList.size() == 0 in MyBall::query_local(...)! " << endl;
    cerr << "\twhen query (" << pm_point[0] << ", " << pm_point[1] << ", " << pm_point[2] << ")" << endl;
    exit(1);
  }

  vector<pair<double, unsigned int> > tmp_buffer;
  double tmp_dist;    
  double tmp_rsq = pm_r * pm_r;
  size_t i;

  // not all points in the bounding box are within pm_r distance, calculating actual distances and
  // finding those which qualifies
  for (i = 0; i < OutputList.size(); i++){
    tmp_dist  = (pm_point[0] - OutputList[i].first.x()) * (pm_point[0] - OutputList[i].first.x());
    tmp_dist += (pm_point[1] - OutputList[i].first.y()) * (pm_point[1] - OutputList[i].first.y());
    tmp_dist += (pm_point[2] - OutputList[i].first.z()) * (pm_point[2] - OutputList[i].first.z());
    if (tmp_dist <= tmp_rsq){
      pm_out.push_back(OutputList[i].second);
    }
  }
}

// Returns the ID of the nearest point on the ball with respect to point pm_point
unsigned short 
MyBall :: query_nn(const vector<double>&   pm_point) {

  vector<UINT_Key> OutputList;

  double  tmp_local_dist = m_nn_r;
  UINT_Pure_key a = UINT_Pure_key(pm_point[0] - tmp_local_dist, pm_point[1] - tmp_local_dist, pm_point[2] - tmp_local_dist);
  UINT_Pure_key b = UINT_Pure_key(pm_point[0] + tmp_local_dist, pm_point[1] + tmp_local_dist, pm_point[2] + tmp_local_dist);
  UINT_Interval win = UINT_Interval(UINT_Key(a, 0), UINT_Key(b, 0));
  m_index->window_query(win, back_inserter(OutputList));

  if (OutputList.size() == 0){
    cerr << "ERROR: OutputList.size() == 0 in MyBall::query_nn(...)! " << endl;
    cerr << "\twhen query (" << pm_point[0] << ", " << pm_point[1] << ", " << pm_point[2] << ")" << endl;
    cerr << "\tm_nn_r = " << m_nn_r << endl;
    double x = 10;
    x = x / OutputList.size();
    //exit(1);
  }

  vector<pair<double, unsigned int> > tmp_buffer;
  double tmp_dist;    
  size_t i;
  for (i = 0; i < OutputList.size(); i++){
    tmp_dist  = (pm_point[0] - OutputList[i].first.x()) * (pm_point[0] - OutputList[i].first.x());
    tmp_dist += (pm_point[1] - OutputList[i].first.y()) * (pm_point[1] - OutputList[i].first.y());
    tmp_dist += (pm_point[2] - OutputList[i].first.z()) * (pm_point[2] - OutputList[i].first.z());
    tmp_dist = sqrt(tmp_dist);
    tmp_buffer.push_back(pair<double, unsigned int>(tmp_dist, OutputList[i].second));
  }

  sort(tmp_buffer.begin(), tmp_buffer.end(), less_than_univ<double, unsigned int>());

  return tmp_buffer[0].second;
}

void
MyBall_Coarse ::  init() {
  m_nn_r = 0;

  vector<UINT_Key>             InputList;
  vector<double>              tmp_point(3, 0);
  size_t i;

  for (i = 0; i < points_default_coarse; i ++){
    InputList.push_back(UINT_Key(UINT_Pure_key(vertex_default_coarse[i][0], vertex_default_coarse[i][1], vertex_default_coarse[i][2]), i));
    if (vertex_default_coarse[i][3] > m_nn_r) m_nn_r = vertex_default_coarse[i][3];
  }

  m_nn_r *= 1.15;

  m_size += sizeof(double) * 3 * points_default_coarse;

  m_index = new UINT_Range_tree_3_type(InputList.begin(), InputList.end());
}

void
MyBall_Dense:: init() {
  m_nn_r = 0.07;
  vector<UINT_Key>             InputList;
  vector<double>              tmp_point(3, 0);
  size_t i;

  for (i = 0; i < points_default_dense; i ++){
    InputList.push_back(UINT_Key(UINT_Pure_key(vertex_default_dense[i][0], vertex_default_dense[i][1], vertex_default_dense[i][2]), i));
    if (vertex_default_dense[i][3] > m_nn_r) m_nn_r = vertex_default_dense[i][3];
  }

  m_nn_r *= 1.15;

  m_size += sizeof(double) * 3 * points_default_dense;

  m_index = new UINT_Range_tree_3_type(InputList.begin(), InputList.end());
}

void
MyBall_Pose:: init() {
  m_nn_r = 0.036;

  vector<UINT_Key>             InputList;
  vector<double>              tmp_point(3, 0);
  size_t i;

  for (i = 0; i < points_default_pose; i ++){
    InputList.push_back(UINT_Key(UINT_Pure_key(vertex_default_pose[i][0], vertex_default_pose[i][1], vertex_default_pose[i][2]), i));
    if (vertex_default_pose[i][3] > m_nn_r) m_nn_r = vertex_default_pose[i][3];
  }

  m_nn_r *= 1.15;

  m_size += sizeof(double) * 3 * points_default_pose;

  m_index = new UINT_Range_tree_3_type(InputList.begin(), InputList.end());
}

// this routine is used while generating matching table.
// fill out the variable m_arc
void 
MyBall_Pose :: init_pose() {
  size_t i, j;
  vector<double>  curr_point(3, 0);
  vector<double>  rand_point(3, 0);
  double          dist;
  double          theta;
  vector<double>  tmp_Z(3, 0);
  vector<double>  tmp_Y(3, 0);
  vector<double>  tmp_X(3, 0);

  srand(time(NULL));
  vector<vector<double> > tmp_buff;
  vector<vector<double> > tmp_coord;

  m_arcs.clear();
  for (i = 0; i < points_default_pose; i ++){
    curr_point[0] = vertex_default_pose[i][0];
    curr_point[1] = vertex_default_pose[i][1];
    curr_point[2] = vertex_default_pose[i][2];
    j = rand() % points_default_pose;
    rand_point[0] = vertex_default_pose[j][0];
    rand_point[1] = vertex_default_pose[j][1];
    rand_point[2] = vertex_default_pose[j][2];
    dist  = rdl_vector_ssd_sqr<double>(curr_point, rand_point);
    while(dist < 1 || dist > 3){
      j = rand() % points_default_pose;
      rand_point[0] = vertex_default_pose[j][0];
      rand_point[1] = vertex_default_pose[j][1];
      rand_point[2] = vertex_default_pose[j][2];
      dist  = rdl_vector_ssd_sqr<double>(curr_point, rand_point);
    }

    // finding a vector which is perpenticular to the above two vector (curr_point, rand_point)
    /*
    tmp_Y[0] = curr_point[1] * rand_point[2] - curr_point[2] * rand_point[1];
    tmp_Y[1] = curr_point[2] * rand_point[0] - curr_point[0] * rand_point[2];
    tmp_Y[2] = curr_point[0] * rand_point[1] - curr_point[1] * rand_point[0];
    */
    tmp_Y = rdl_cross_3d<double>(curr_point, rand_point);
    normalize_unit_3d(tmp_Y);
    /*
    tmp_X[0] = tmp_Y[1] * curr_point[2] - tmp_Y[2] * curr_point[1];
    tmp_X[1] = tmp_Y[2] * curr_point[0] - tmp_Y[0] * curr_point[2];
    tmp_X[2] = tmp_Y[0] * curr_point[1] - tmp_Y[1] * curr_point[0];
    */
    tmp_X = rdl_cross_3d<double>(tmp_Y, curr_point);
    normalize_unit_3d(tmp_X);
    tmp_buff.clear();
    for (j = 0; j < 360; j ++){
      theta = 2.0*i*M_PI/360;
      // curr_point is rotated so that it is aligned with z-axis, last param is output
      rdl_rotate(curr_point, tmp_X, theta, rand_point);
      tmp_buff.push_back(rand_point);
    }
    random_shuffle(tmp_buff.begin(), tmp_buff.end());
    rand_point = tmp_buff[rand()%tmp_buff.size()];

    tmp_Y[0] = curr_point[1] * rand_point[2] - curr_point[2] * rand_point[1];
    tmp_Y[1] = curr_point[2] * rand_point[0] - curr_point[0] * rand_point[2];
    tmp_Y[2] = curr_point[0] * rand_point[1] - curr_point[1] * rand_point[0];
    normalize_unit_3d(tmp_Y);
    tmp_X[0] = tmp_Y[1] * curr_point[2] - tmp_Y[2] * curr_point[1];
    tmp_X[1] = tmp_Y[2] * curr_point[0] - tmp_Y[0] * curr_point[2];
    tmp_X[2] = tmp_Y[0] * curr_point[1] - tmp_Y[1] * curr_point[0];
    normalize_unit_3d(tmp_X);
    tmp_Z = curr_point;

    tmp_coord.clear();
    tmp_coord.push_back(tmp_X);
    tmp_coord.push_back(tmp_Y);
    tmp_coord.push_back(tmp_Z);
    m_arcs.push_back(tmp_coord);

    tmp_coord.clear();
    tmp_X[0] = 0 - tmp_X[0];
    tmp_X[1] = 0 - tmp_X[1];
    tmp_X[2] = 0 - tmp_X[2];

    tmp_Y[0] = 0 - tmp_Y[0];
    tmp_Y[1] = 0 - tmp_Y[1];
    tmp_Y[2] = 0 - tmp_Y[2];

    tmp_coord.push_back(tmp_X);
    tmp_coord.push_back(tmp_Y);
    tmp_coord.push_back(tmp_Z);
    m_arcs.push_back(tmp_coord);
  }
  m_size += m_arcs.size() * sizeof(double) * 9;
}

unsigned long 
MyUnitBall :: size() {
        unsigned long tmp_size = 0;

        tmp_size += m_unitball_coarse.size() * 3 * sizeof(double);
        tmp_size += m_unitball_dense.size() * 3 * sizeof(double);
        tmp_size += m_unitball_pose.size() * 3 * sizeof(double);
        tmp_size += m_coarse.m_size;
        tmp_size += m_dense.m_size;
        tmp_size += m_pose.m_size;

        tmp_size += m_table.size() * sizeof(unsigned short) * points_default_coarse;
        tmp_size += m_table_coarse2dense.size() * sizeof(unsigned short) * points_default_dense;
        tmp_size += m_table_coarse2dense.size() * sizeof(unsigned short) * points_default_dense;
        tmp_size += m_rmatrix.size() * sizeof(double) * 9;

        return tmp_size;
    }

void 
MyUnitBall :: init() {
        //cout << "load_unitball();" << endl;
        load_unitball();
        compute_farthest_neighbor_dist();
        m_coarse.init();
        m_dense.init();
        m_pose.init();

        stringstream  tmp_fn_table;
        tmp_fn_table << "ssave_matching_table_" << g_matching_table_selector;

        sload_matching_table(tmp_fn_table.str() + ".txt");
        //sload_matching_table_sort(tmp_fn_table.str() + "_sort.txt");

        //cout << "sload_matching_table_dense2coarse(...)" << endl;
        sload_matching_table_dense2coarse(tmp_fn_table.str() + "_dense2coarse.txt");

        //cout << "sload_matching_table_coarse2dense(...)" << endl;
        sload_matching_table_coarse2dense(tmp_fn_table.str() + "_coarse2dense.txt");

        //cout << "sload_matching_table_rmatrix(...)" << endl;
        sload_matching_table_rmatrix(tmp_fn_table.str() + "_rmatrix.txt");

}   
    
void 
MyUnitBall :: sload_matching_table(string pm_filename) {
        cout << "sload_matching_table ... " << pm_filename << endl;
        m_table.clear();
        string in_filename = g_dir_table + pm_filename;
        std::ifstream ifs(in_filename.c_str(), std::ios::binary);
        if (!ifs.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << endl;
            exit(1);
        }

        boost::archive::text_iarchive la(ifs);
        la & m_table;
}

void 
MyUnitBall :: sload_matching_table_coarse2dense(string pm_filename) {
        cout << "sload_matching_table_coarse2dense ... " << pm_filename << endl;
        m_table_coarse2dense.clear();
        string in_filename = g_dir_table + pm_filename;
        std::ifstream ifs(in_filename.c_str(), std::ios::binary);
        if (!ifs.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << endl;
            exit(1);
        }

        boost::archive::text_iarchive la(ifs);
        la & m_table_coarse2dense;
}


void 
MyUnitBall :: sload_matching_table_dense2coarse(string pm_filename) {
        cout << "sload_matching_table_dense2coarse ... " << pm_filename << endl;
        m_table_dense2coarse.clear();
        string in_filename = g_dir_table + pm_filename;
        std::ifstream ifs(in_filename.c_str(), std::ios::binary);
        if (!ifs.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << endl;
            exit(1);
        }

        boost::archive::text_iarchive la(ifs);
        la & m_table_dense2coarse;
}

void 
MyUnitBall :: sload_matching_table_rmatrix(string pm_filename) {
        cout << "sload_matching_table_rmatrix ... " << pm_filename << endl;
        m_rmatrix.clear();
        string in_filename = g_dir_table + pm_filename;
        std::ifstream ifs(in_filename.c_str(), std::ios::binary);
        if (!ifs.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << endl;
            exit(1);
        }

        boost::archive::text_iarchive la(ifs);
        la & m_rmatrix;
}

void 
MyUnitBall :: ssave_matching_table(string pm_filename) {
        std::ofstream ofs(pm_filename.c_str());
        if (!ofs.is_open()){
            cerr << "ERROR: cannot open file " << pm_filename << endl;
            exit(1);
        }

        boost::archive::text_oarchive oa(ofs);
        oa & m_table;
}


void 
MyUnitBall :: ssave_matching_table_coarse2dense(string pm_filename) {
        std::ofstream ofs(pm_filename.c_str());
        if (!ofs.is_open()){
            cerr << "ERROR: cannot open file " << pm_filename << endl;
            exit(1);
        }

        boost::archive::text_oarchive oa(ofs);
        oa & m_table_coarse2dense;
}

void 
MyUnitBall :: ssave_matching_table_dense2coarse(string pm_filename) {
        std::ofstream ofs(pm_filename.c_str());
        if (!ofs.is_open()){
            cerr << "ERROR: cannot open file " << pm_filename << endl;
            exit(1);
        }

		boost::archive::text_oarchive oa(ofs);
        oa & m_table_dense2coarse;
}

void 
MyUnitBall :: ssave_matching_table_rmatrix(string pm_filename) {
        std::ofstream ofs(pm_filename.c_str());
        if (!ofs.is_open()){
            cerr << "ERROR: cannot open file " << pm_filename << endl;
            exit(1);
        }

		boost::archive::text_oarchive oa(ofs);
        oa & m_rmatrix;
}

void 
MyUnitBall :: load_unitball() {
        vector<double>  tmp_point(3, 0);
        size_t i;

        m_unitball_coarse.clear();
        for (i = 0; i < points_default_coarse; i ++){
            tmp_point[0] = vertex_default_coarse[i][0];
            tmp_point[1] = vertex_default_coarse[i][1];
            tmp_point[2] = vertex_default_coarse[i][2];
            m_unitball_coarse.push_back(tmp_point);
        }

        m_unitball_dense.clear();
        for (i = 0; i < points_default_dense; i ++){
            tmp_point[0] = vertex_default_dense[i][0];
            tmp_point[1] = vertex_default_dense[i][1];
            tmp_point[2] = vertex_default_dense[i][2];
            m_unitball_dense.push_back(tmp_point);
        }

        m_unitball_pose.clear();
        for (i = 0; i < points_default_pose; i ++){
            tmp_point[0] = vertex_default_pose[i][0];
            tmp_point[1] = vertex_default_pose[i][1];
            tmp_point[2] = vertex_default_pose[i][2];
            m_unitball_pose.push_back(tmp_point);
        }
}

void 
MyUnitBall :: compute_farthest_neighbor_dist() {
        vector<pair<double, pair<int, int> > >  tmp_edges;
        tmp_edges.clear();
        size_t i, j;
        double tmp_dist;

        for (i = 0; i < points_default_coarse; i ++){
            for (j = i+1; j < points_default_coarse; j ++){
                tmp_dist  = (vertex_default_coarse[i][0] - vertex_default_coarse[j][0])*(vertex_default_coarse[i][0] - vertex_default_coarse[j][0]);
                tmp_dist += (vertex_default_coarse[i][1] - vertex_default_coarse[j][1])*(vertex_default_coarse[i][1] - vertex_default_coarse[j][1]);
                tmp_dist += (vertex_default_coarse[i][2] - vertex_default_coarse[j][2])*(vertex_default_coarse[i][2] - vertex_default_coarse[j][2]);
                tmp_dist = sqrt(tmp_dist);
                tmp_edges.push_back(pair<double, pair<int, int> >(tmp_dist, pair<int, int>(i, j)));
            }
        }
        sort(tmp_edges.begin(), tmp_edges.end(), less_than_univ<double, pair<int, int> >());
        m_farthest_neighbor_dist_coarse = tmp_edges[3*points_default_coarse - 7].first;

        m_coarse.m_nn_r = m_farthest_neighbor_dist_coarse;
        m_dense.m_nn_r = m_farthest_neighbor_dist_coarse;
        m_pose.m_nn_r = m_farthest_neighbor_dist_coarse;

        tmp_edges.clear();
        for (i = 0; i < points_default_dense; i ++){
            for (j = i+1; j < points_default_coarse; j ++){
                tmp_dist  = (vertex_default_dense[i][0] - vertex_default_dense[j][0])*(vertex_default_dense[i][0] - vertex_default_dense[j][0]);
                tmp_dist += (vertex_default_dense[i][1] - vertex_default_dense[j][1])*(vertex_default_dense[i][1] - vertex_default_dense[j][1]);
                tmp_dist += (vertex_default_dense[i][2] - vertex_default_dense[j][2])*(vertex_default_dense[i][2] - vertex_default_dense[j][2]);
                tmp_dist = sqrt(tmp_dist);
                tmp_edges.push_back(pair<double, pair<int, int> >(tmp_dist, pair<int, int>(i, j)));
            }
        }
        sort(tmp_edges.begin(), tmp_edges.end(), less_than_univ<double, pair<int, int> >());
        m_farthest_neighbor_dist_dense = tmp_edges[3*points_default_dense - 7].first;

        tmp_edges.clear();
}

unsigned long 
MyContextBall :: size() {
        unsigned long tmp_sizeof_cs = 0;

        tmp_sizeof_cs += sizeof(m_center_area);
        tmp_sizeof_cs += sizeof(void*); //m_unitball_template
        tmp_sizeof_cs += sizeof(m_ball_id);
        tmp_sizeof_cs += sizeof(m_ball_id_serialno);

        tmp_sizeof_cs += sizeof(double) * 3; //m_ball_center[3];
        tmp_sizeof_cs += sizeof(double) * g_solid_volume_num_dynamicR; //m_solid_volume_dynamicR
        tmp_sizeof_cs += sizeof(double) * 3 * g_solid_volume_num_dynamicR; //m_solid_vector
        tmp_sizeof_cs += sizeof(double) * g_solid_volume_num_dynamicR; //m_solid_angle
        tmp_sizeof_cs += sizeof(double) * 9; //m_rotation_matrix
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_rays.size(); 
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_shell_L2.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_core_L2.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_inner.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_shell_1A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_shell_2A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_shell_3A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_shell_4A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_core_1A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_core_2A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_core_3A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_core_4A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_inner_2A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_inner_core3A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_core1A_shell3A.size();
        tmp_sizeof_cs += sizeof(CONTEXT_RAY_INT_TYPE) * m_context_core1Ashell3A.size();
        tmp_sizeof_cs += sizeof(MyBuriedVert*) * m_buried_verts.size();


        return tmp_sizeof_cs;
}


vector<double>  
MyContextBall :: get_center_point() {
        vector<double>  rtn_point(3, 0);
        rtn_point[0] = m_ball_center[0];
        rtn_point[1] = m_ball_center[1];
        rtn_point[2] = m_ball_center[2];
        return rtn_point;
}

void 
MyContextBall :: sload(boost::archive::text_iarchive &la) {
        la & m_center_area;
        la & m_ball_id;
        la & m_ball_id_serialno;
        la & m_ball_center;
        la & m_solid_volume_dynamicR;
        la & m_solid_vector; //3D vector
        la & m_solid_angle;
        la & m_rotation_matrix;

        la & m_context_rays;
        la & m_context_shell_L2;
        la & m_context_core_L2;
        la & m_context_inner;
	la & m_context_inner_2A;
	la & m_context_inner_core3A;
	la & m_context_core1A_shell3A;
        la & m_context_core1Ashell3A;

        la & m_context_shell_1A;
        la & m_context_shell_2A;
        la & m_context_shell_3A;
        la & m_context_shell_4A;
        la & m_context_core_1A;
        la & m_context_core_2A;
        la & m_context_core_3A;
        la & m_context_core_4A;

	size_t i;
        size_t num_verts;
        la & num_verts;
        m_buried_verts.clear();
        m_buried_verts.reserve(num_verts);
        for (i = 0; i < num_verts; i ++){
            MyBuriedVert* tmp_vert = new MyBuriedVert();
            tmp_vert->sload(la);
            m_buried_verts.push_back(tmp_vert);
        }

}

void 
MyContextBall :: ssave(boost::archive::text_oarchive &la) {
        la & m_center_area;
        la & m_ball_id;
        la & m_ball_id_serialno;
        la & m_ball_center;
        la & m_solid_volume_dynamicR;
        la & m_solid_vector; //3D vector
        la & m_solid_angle;
        la & m_rotation_matrix;

        la & m_context_rays;
        la & m_context_shell_L2;
        la & m_context_core_L2;
        la & m_context_inner;
	la & m_context_inner_2A;
	la & m_context_inner_core3A;
	la & m_context_core1A_shell3A;
        la & m_context_core1Ashell3A;

        la & m_context_shell_1A;
        la & m_context_shell_2A;
        la & m_context_shell_3A;
        la & m_context_shell_4A;
        la & m_context_core_1A;
        la & m_context_core_2A;
        la & m_context_core_3A;
        la & m_context_core_4A;

        size_t i;

        size_t num_verts = m_buried_verts.size();
        la & num_verts;
        for (i = 0; i < num_verts; i ++){
            m_buried_verts[i]->ssave(la);
        }
}

MyContextBall :: MyContextBall(int pm_ball_id, vector<double>& pm_ball_center, 
                               vector<vector<double> >* pm_unitball_template) {
        m_ball_id = pm_ball_id;
        m_ball_id_serialno = 0;
        m_ball_center[0] = pm_ball_center[0];
        m_ball_center[1] = pm_ball_center[1];
        m_ball_center[2] = pm_ball_center[2];
        m_unitball_template = pm_unitball_template;

        //m_context_rays.clear();
        //m_context_shell_L2.clear(); 
        //m_context_core_L2.clear();  
        //m_context_inner.clear();    

        //m_context_shell_1A.clear(); 
        //m_context_shell_2A.clear(); 
        //m_context_shell_3A.clear(); 
        //m_context_shell_4A.clear(); 

        //m_context_core_1A.clear();  
        //m_context_core_2A.clear();  
        //m_context_core_3A.clear();  
        //m_context_core_4A.clear();  


        size_t i, j;
        //set unit matrix
        for (i = 0; i < 3; i ++){
            for (j = 0; j < 3; j ++){
                m_rotation_matrix[i][j] = 0;
            }
            m_rotation_matrix[i][i] = 1;
        }
}

MyContextBall :: MyContextBall() {
        size_t i, j;
        //set unit matrix
        for (i = 0; i < 3; i ++){
            for (j = 0; j < 3; j ++){
                m_rotation_matrix[i][j] = 0;
            }
            m_rotation_matrix[i][i] = 1;
        }
}

void 
MyContextBall :: write_shell_vtk(vector<unsigned int>& pm_rays, vector<vector<double> >& pm_unitball, string out_fn) {
        size_t i, j;
        vector<double>          tmp_unit(3, 0);
        vector<double>          tmp_point_start(3, 0);
        vector<double>          tmp_point_end(3, 0);

        vector<vector<double> > tmp_points;
        vector<pair<int, int> > tmp_lines;

        int id_point = 0;

        vector<double>          tmp_ball_center(3, 0);
        tmp_ball_center[0] = m_ball_center[0];
        tmp_ball_center[1] = m_ball_center[1];
        tmp_ball_center[2] = m_ball_center[2];

        for (i = 0; i < pm_rays.size(); i ++){
            //tmp_unit[0]  = m_rotation_matrix[0][0] * pm_unitball[i][0];
            //tmp_unit[0] += m_rotation_matrix[1][0] * pm_unitball[i][1];
            //tmp_unit[0] += m_rotation_matrix[2][0] * pm_unitball[i][2];

            //tmp_unit[1]  = m_rotation_matrix[0][1] * pm_unitball[i][0];
            //tmp_unit[1] += m_rotation_matrix[1][1] * pm_unitball[i][1];
            //tmp_unit[1] += m_rotation_matrix[2][1] * pm_unitball[i][2];

            //tmp_unit[2]  = m_rotation_matrix[0][2] * pm_unitball[i][0];
            //tmp_unit[2] += m_rotation_matrix[1][2] * pm_unitball[i][1];
            //tmp_unit[2] += m_rotation_matrix[2][2] * pm_unitball[i][2];

            tmp_unit[0]  = m_rotation_matrix[0][0] * pm_unitball[i][0];
            tmp_unit[0] += m_rotation_matrix[0][1] * pm_unitball[i][1];
            tmp_unit[0] += m_rotation_matrix[0][2] * pm_unitball[i][2];

            tmp_unit[1]  = m_rotation_matrix[1][0] * pm_unitball[i][0];
            tmp_unit[1] += m_rotation_matrix[1][1] * pm_unitball[i][1];
            tmp_unit[1] += m_rotation_matrix[1][2] * pm_unitball[i][2];

            tmp_unit[2]  = m_rotation_matrix[2][0] * pm_unitball[i][0];
            tmp_unit[2] += m_rotation_matrix[2][1] * pm_unitball[i][1];
            tmp_unit[2] += m_rotation_matrix[2][2] * pm_unitball[i][2];


            //tmp_unit = pm_unitball[i];

            for (j = 0; j < 32; j ++){
                if (pm_rays[i] & (0x1u << j)){
                    tmp_point_start = tmp_ball_center + (j * g_sample_region_radius / 32.0) * tmp_unit;
                    while(j < 32 && (pm_rays[i] & (0x1u << j))) j ++;
                    tmp_point_end = tmp_ball_center + (j * g_sample_region_radius / 32.0) * tmp_unit;
                    tmp_points.push_back(tmp_point_start);
                    tmp_points.push_back(tmp_point_end);
                    tmp_lines.push_back(pair<int, int>(id_point, id_point + 1));
                    id_point += 2;
                }
            }
        }

        ofstream out_file(out_fn.c_str());
        if (!out_file.is_open()){
            cerr << "ERROR: cannot open file " << out_fn << endl;
            exit(1);
        }

        stringstream out_buffer;
        out_buffer << "# vtk DataFile Version 3.0" << endl;
        out_buffer << "ContextBall_" << m_ball_id << " " << out_fn << endl;
        out_buffer << "ASCII" << endl;
        out_buffer << "DATASET POLYDATA" << endl;

        out_buffer << "POINTS " << tmp_points.size() << " float" << endl;
        for (i = 0; i < tmp_points.size(); i ++){
            out_buffer << tmp_points[i][0] << " "
                << tmp_points[i][1] << " "
                << tmp_points[i][2] << endl;
        }

        out_buffer << "LINES " << tmp_lines.size() << " " << tmp_lines.size() * 3 << endl;
        for (i = 0; i < tmp_lines.size(); i ++){
            out_buffer << "2 " <<  tmp_lines[i].first << " " << tmp_lines[i].second << endl;
        }

        out_file << out_buffer.str();
        out_file.close();

}

void 
MyContextBall :: write_vtk(string out_fn) {
        vector<vector<double> >& pm_unitball = * m_unitball_template;
        size_t i, j;
        vector<double>          tmp_unit(3, 0);
        vector<double>          tmp_point_start(3, 0);
        vector<double>          tmp_point_end(3, 0);
        vector<vector<double> > tmp_points;
        vector<pair<int, int> > tmp_lines;

        int id_point = 0;

        vector<double>          tmp_ball_center(3, 0);
        tmp_ball_center[0] = m_ball_center[0];
        tmp_ball_center[1] = m_ball_center[1];
        tmp_ball_center[2] = m_ball_center[2];

        for (i = 0; i < m_context_rays.size(); i ++){
            //de-normalize the pose
            tmp_unit[0]  = m_rotation_matrix[0][0] * pm_unitball[i][0];
            tmp_unit[0] += m_rotation_matrix[0][1] * pm_unitball[i][1];
            tmp_unit[0] += m_rotation_matrix[0][2] * pm_unitball[i][2];

            tmp_unit[1]  = m_rotation_matrix[1][0] * pm_unitball[i][0];
            tmp_unit[1] += m_rotation_matrix[1][1] * pm_unitball[i][1];
            tmp_unit[1] += m_rotation_matrix[1][2] * pm_unitball[i][2];

            tmp_unit[2]  = m_rotation_matrix[2][0] * pm_unitball[i][0];
            tmp_unit[2] += m_rotation_matrix[2][1] * pm_unitball[i][1];
            tmp_unit[2] += m_rotation_matrix[2][2] * pm_unitball[i][2];

            //the original post, used in dockinfo mode
            //tmp_unit = pm_unitball[i];

            for (j = 0; j < 32; j ++){
                if (m_context_rays[i] & (0x1u << j)){
                    tmp_point_start = tmp_ball_center + (j * g_sample_region_radius / 32.0) * tmp_unit;
                    while(j < 32 && (m_context_rays[i] & (0x1u << j))) j ++;
                    tmp_point_end = tmp_ball_center + (j * g_sample_region_radius / 32.0) * tmp_unit;
                    tmp_points.push_back(tmp_point_start);
                    tmp_points.push_back(tmp_point_end);
                    tmp_lines.push_back(pair<int, int>(id_point, id_point + 1));
                    id_point += 2;
                }
            }
        }

        ofstream out_file(out_fn.c_str());
        if (!out_file.is_open()){
            cerr << "ERROR: cannot open file " << out_fn << endl;
            exit(1);
        }

        stringstream out_buffer;
        out_buffer << "# vtk DataFile Version 3.0" << endl;
        out_buffer << "ContextBall_" << m_ball_id << " " << out_fn << endl;
        out_buffer << "ASCII" << endl;
        out_buffer << "DATASET POLYDATA" << endl;

        out_buffer << "POINTS " << tmp_points.size() << " float" << endl;
        for (i = 0; i < tmp_points.size(); i ++){
            out_buffer << tmp_points[i][0] << " "
                << tmp_points[i][1] << " "
                << tmp_points[i][2] << endl;
        }

        out_buffer << "LINES " << tmp_lines.size() << " " << tmp_lines.size() * 3 << endl;
        for (i = 0; i < tmp_lines.size(); i ++){
            out_buffer << "2 " <<  tmp_lines[i].first << " " << tmp_lines[i].second << endl;
        }

        out_file << out_buffer.str();
        out_file.close();
}

void 
MyContextBall :: write_vtk_dockinfo(vector<vector<double> >& pm_unitball, string out_fn) {
        size_t i, j;
        vector<double>          tmp_unit(3, 0);
        vector<double>          tmp_point_start(3, 0);
        vector<double>          tmp_point_end(3, 0);
        vector<vector<double> > tmp_points;
        vector<pair<int, int> > tmp_lines;

        int id_point = 0;

        vector<double>          tmp_ball_center(3, 0);
        tmp_ball_center[0] = m_ball_center[0];
        tmp_ball_center[1] = m_ball_center[1];
        tmp_ball_center[2] = m_ball_center[2];

        for (i = 0; i < m_context_rays.size(); i ++){

            tmp_unit = pm_unitball[i];

            for (j = 0; j < 32; j ++){
                if (m_context_rays[i] & (0x1u << j)){
                    tmp_point_start = tmp_ball_center + (j * g_sample_region_radius / 32.0) * tmp_unit;
                    while(j < 32 && (m_context_rays[i] & (0x1u << j))) j ++;
                    tmp_point_end = tmp_ball_center + (j * g_sample_region_radius / 32.0) * tmp_unit;
                    tmp_points.push_back(tmp_point_start);
                    tmp_points.push_back(tmp_point_end);
                    tmp_lines.push_back(pair<int, int>(id_point, id_point + 1));
                    id_point += 2;
                }
            }
        }

        ofstream out_file(out_fn.c_str());
        if (!out_file.is_open()){
            cerr << "ERROR: cannot open file " << out_fn << endl;
            exit(1);
        }

        stringstream out_buffer;
        out_buffer << "# vtk DataFile Version 3.0" << endl;
        out_buffer << "ContextBall_" << m_ball_id << " " << out_fn << endl;
        out_buffer << "ASCII" << endl;
        out_buffer << "DATASET POLYDATA" << endl;

        out_buffer << "POINTS " << tmp_points.size() << " float" << endl;
        for (i = 0; i < tmp_points.size(); i ++){
            out_buffer << tmp_points[i][0] << " "
                << tmp_points[i][1] << " "
                << tmp_points[i][2] << endl;
        }

        out_buffer << "LINES " << tmp_lines.size() << " " << tmp_lines.size() * 3 << endl;
        for (i = 0; i < tmp_lines.size(); i ++){
            out_buffer << "2 " <<  tmp_lines[i].first << " " << tmp_lines[i].second << endl;
        }

        out_file << out_buffer.str();
        out_file.close();
}
