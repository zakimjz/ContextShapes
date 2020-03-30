#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include <ext/algorithm>
#include <numeric>
#include <sys/timeb.h>
#include <time.h>
#include <cstdlib>
#include <cstdio>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include "rbt_coarse.h"
#include "rbt_dense.h"
#include "rbt_pose.h"

#include "CGAL/Cartesian.h"
#include "CGAL/Segment_tree_k.h"
#include "CGAL/Range_segment_tree_traits.h"

#include "CGAL/basic.h" 
#include "CGAL/Point_3.h"
#include "CGAL/Range_tree_k.h"

typedef CGAL::Cartesian<double> CTP_Representation;
typedef CGAL::Range_tree_map_traits_3<CTP_Representation, unsigned short> CTP_Traits;
typedef CGAL::Range_tree_3<CTP_Traits> CTP_Range_tree_3_type;
typedef CTP_Traits::Key CTP_Key;
typedef CTP_Traits::Pure_key CTP_Pure_key;
typedef CTP_Traits::Interval CTP_Interval;

using namespace std;


//typedef double RotationMatrix[3][3];
//typedef vector<vector<double> > RotationMatrix(3, vector<double>(3, 0));

template <class FirstType, class SecondType>
struct less_than_univ{
    bool operator()(pair<FirstType, SecondType> pm_A, pair<FirstType, SecondType> pm_B)
    {
        if (pm_A.first < pm_B.first) return true;
        else return false;
    }
};


template <class T>
T rdl_dot_3d(const vector<T>& pm_v1, const vector<T>& pm_v2)
{
    T     r;
    r  = pm_v1[0] * pm_v2[0];
    r += pm_v1[1] * pm_v2[1];
    r += pm_v1[2] * pm_v2[2];

    return r;
}


template <class T>
void rdl_rotate(vector<T> const &pm_axis, \
                vector<T> const &pm_old_point, \
                T  	            pm_angle, \
                vector<T>       &pm_new_point)
{
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

void normalize(vector<double>& v)
{
	double n;
	n = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	if (n < 0.5){
		cerr << "n is very small! = " << n << endl; 
	}
	n = sqrt(n);
	v[0] /= n;
	v[1] /= n;
	v[2] /= n;
}

struct MyBall{
	CTP_Range_tree_3_type*      m_index;

	MyBall() 
	{
		m_index = NULL;
	}

	void query_local(const vector<double>&	pm_point, double pm_r, vector<unsigned short>& pm_out)
	{
		pm_out.clear();

		vector<CTP_Key> OutputList;
		
		CTP_Pure_key a = CTP_Pure_key(pm_point[0] - pm_r, pm_point[1] - pm_r, pm_point[2] - pm_r);
		CTP_Pure_key b = CTP_Pure_key(pm_point[0] + pm_r, pm_point[1] + pm_r, pm_point[2] + pm_r);
		CTP_Interval win = CTP_Interval(CTP_Key(a, 0), CTP_Key(b, 0));
		m_index->window_query(win, back_inserter(OutputList));

		if (OutputList.size() == 0){
			cerr << "ERROR: OutputList.size() == 0 in query_nnb(...)! " << endl;
			cerr << "\twhen query (" << pm_point[0] << ", " << pm_point[1] << ", " << pm_point[2] << ")" << endl;
			exit(1);
		}
		
		vector<pair<double, unsigned short> >	tmp_buffer;
		double tmp_dist;	
		double tmp_rsq = pm_r * pm_r;
		size_t i;
		for (i = 0; i < OutputList.size(); i++){
			tmp_dist  = (pm_point[0] - OutputList[i].first.x()) * (pm_point[0] - OutputList[i].first.x());
			tmp_dist += (pm_point[1] - OutputList[i].first.y()) * (pm_point[1] - OutputList[i].first.y());
			tmp_dist += (pm_point[2] - OutputList[i].first.z()) * (pm_point[2] - OutputList[i].first.z());
			if (tmp_dist <= tmp_rsq){
				pm_out.push_back(OutputList[i].second);
			}
		}
	}

	unsigned short query_nn(const vector<double>&	pm_point)
	{
		vector<CTP_Key> OutputList;
		
		double	tmp_local_dist = 0.09;
		CTP_Pure_key a = CTP_Pure_key(pm_point[0] - tmp_local_dist, pm_point[1] - tmp_local_dist, pm_point[2] - tmp_local_dist);
		CTP_Pure_key b = CTP_Pure_key(pm_point[0] + tmp_local_dist, pm_point[1] + tmp_local_dist, pm_point[2] + tmp_local_dist);
		CTP_Interval win = CTP_Interval(CTP_Key(a, 0), CTP_Key(b, 0));
		m_index->window_query(win, back_inserter(OutputList));

		if (OutputList.size() == 0){
			cerr << "ERROR: OutputList.size() == 0 in query_nnb(...)! " << endl;
			cerr << "\twhen query (" << pm_point[0] << ", " << pm_point[1] << ", " << pm_point[2] << ")" << endl;
			exit(1);
		}
		
		vector<pair<double, unsigned short> >	tmp_buffer;
		double tmp_dist;	
		size_t i;
		for (i = 0; i < OutputList.size(); i++){
			tmp_dist  = (pm_point[0] - OutputList[i].first.x()) * (pm_point[0] - OutputList[i].first.x());
			tmp_dist += (pm_point[1] - OutputList[i].first.y()) * (pm_point[1] - OutputList[i].first.y());
			tmp_dist += (pm_point[2] - OutputList[i].first.z()) * (pm_point[2] - OutputList[i].first.z());
			tmp_buffer.push_back(pair<double, unsigned short>(tmp_dist, OutputList[i].second));
		}

		sort(tmp_buffer.begin(), tmp_buffer.end(), less_than_univ<double, unsigned short>());

		return tmp_buffer[0].second;
	}
};

struct MyBall_Coarse: public MyBall{

	void init()
	{
		vector<CTP_Key>             InputList;
		vector<double>              tmp_point(3, 0);
		size_t i;

		for (i = 0; i < points_default_coarse; i ++){
			InputList.push_back(CTP_Key(CTP_Pure_key(vertex_default_coarse[i][0], vertex_default_coarse[i][1], vertex_default_coarse[i][2]), i));
		}
		
		m_index = new CTP_Range_tree_3_type(InputList.begin(), InputList.end());
	}
};

struct MyBall_Dense: public MyBall{

	void init()
	{
		vector<CTP_Key>             InputList;
		vector<double>              tmp_point(3, 0);
		size_t i;

		for (i = 0; i < points_default_dense; i ++){
			InputList.push_back(CTP_Key(CTP_Pure_key(vertex_default_dense[i][0], vertex_default_dense[i][1], vertex_default_dense[i][2]), i));
		}
		
		m_index = new CTP_Range_tree_3_type(InputList.begin(), InputList.end());
	}
};

struct MyBall_Pose: public MyBall{
	vector<vector<vector<double> > >		m_arcs;

	void init()
	{
		vector<CTP_Key>             InputList;
		vector<double>              tmp_point(3, 0);
		size_t i;

		for (i = 0; i < points_default_pose; i ++){
			InputList.push_back(CTP_Key(CTP_Pure_key(vertex_default_pose[i][0], vertex_default_pose[i][1], vertex_default_pose[i][2]), i));
		}
		
		m_index = new CTP_Range_tree_3_type(InputList.begin(), InputList.end());
	}

	void init_pose()
	{
		size_t i, j;
		vector<double>	curr_point(3, 0);
		vector<double>	rand_point(3, 0);
		double			dist;
		double			theta;
		vector<double>	tmp_Z(3, 0);
		vector<double>	tmp_Y(3, 0);
		vector<double>	tmp_X(3, 0);

		srand(time(NULL));
		vector<vector<double> >	tmp_buff;
		vector<vector<double> >	tmp_coord;

		m_arcs.clear();
		for (i = 0; i < points_default_pose; i ++){
			curr_point[0] = vertex_default_pose[i][0];
			curr_point[1] = vertex_default_pose[i][1];
			curr_point[2] = vertex_default_pose[i][2];
			j = rand() % points_default_pose;
			rand_point[0] = vertex_default_pose[j][0];
			rand_point[1] = vertex_default_pose[j][1];
			rand_point[2] = vertex_default_pose[j][2];
			dist  = (curr_point[0] - rand_point[0]) * (curr_point[0] - rand_point[0]);
			dist += (curr_point[1] - rand_point[1]) * (curr_point[1] - rand_point[1]);
			dist += (curr_point[2] - rand_point[2]) * (curr_point[2] - rand_point[2]);
			while(dist < 1 || dist > 3){
				j = rand() % points_default_pose;
				rand_point[0] = vertex_default_pose[j][0];
				rand_point[1] = vertex_default_pose[j][1];
				rand_point[2] = vertex_default_pose[j][2];
				dist  = (curr_point[0] - rand_point[0]) * (curr_point[0] - rand_point[0]);
				dist += (curr_point[1] - rand_point[1]) * (curr_point[1] - rand_point[1]);
				dist += (curr_point[2] - rand_point[2]) * (curr_point[2] - rand_point[2]);
			}

			tmp_Y[0] = curr_point[1] * rand_point[2] - curr_point[2] * rand_point[1];
			tmp_Y[1] = curr_point[2] * rand_point[0] - curr_point[0] * rand_point[2];
			tmp_Y[2] = curr_point[0] * rand_point[1] - curr_point[1] * rand_point[0];
			normalize(tmp_Y);
			tmp_X[0] = tmp_Y[1] * curr_point[2] - tmp_Y[2] * curr_point[1];
			tmp_X[1] = tmp_Y[2] * curr_point[0] - tmp_Y[0] * curr_point[2];
			tmp_X[2] = tmp_Y[0] * curr_point[1] - tmp_Y[1] * curr_point[0];
			normalize(tmp_X);
			tmp_buff.clear();
			for (j = 0; j < 360; j ++){
				theta = 2.0*i*M_PI/360;
				rdl_rotate(curr_point, tmp_X, theta, rand_point);
				tmp_buff.push_back(rand_point);
			}
			random_shuffle(tmp_buff.begin(), tmp_buff.end());
			rand_point = tmp_buff[rand()%tmp_buff.size()];

			tmp_Y[0] = curr_point[1] * rand_point[2] - curr_point[2] * rand_point[1];
			tmp_Y[1] = curr_point[2] * rand_point[0] - curr_point[0] * rand_point[2];
			tmp_Y[2] = curr_point[0] * rand_point[1] - curr_point[1] * rand_point[0];
			normalize(tmp_Y);
			tmp_X[0] = tmp_Y[1] * curr_point[2] - tmp_Y[2] * curr_point[1];
			tmp_X[1] = tmp_Y[2] * curr_point[0] - tmp_Y[0] * curr_point[2];
			tmp_X[2] = tmp_Y[0] * curr_point[1] - tmp_Y[1] * curr_point[0];
			normalize(tmp_X);
			tmp_Z = curr_point;

			tmp_coord.clear();
			tmp_coord.push_back(tmp_X);
			tmp_coord.push_back(tmp_Y);
			tmp_coord.push_back(tmp_Z);
			m_arcs.push_back(tmp_coord);

            //cout << "X*Y = " << rdl_dot_3d(tmp_X, tmp_Y) << endl;
            //cout << "Y*Z = " << rdl_dot_3d(tmp_Y, tmp_Z) << endl;
            //cout << "Z*X = " << rdl_dot_3d(tmp_Z, tmp_X) << endl;


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
            
            //cout << "X*Y = " << rdl_dot_3d(tmp_X, tmp_Y) << endl;
            //cout << "Y*Z = " << rdl_dot_3d(tmp_Y, tmp_Z) << endl;
            //cout << "Z*X = " << rdl_dot_3d(tmp_Z, tmp_X) << endl << endl;
		}
	}
};


struct MyUnitBall{
	vector<vector<double> >		m_unitball_coarse;
	vector<vector<double> >		m_unitball_dense;

	MyBall_Coarse			m_coarse;
	MyBall_Dense				m_dense;
	MyBall_Pose				m_pose;
	vector<vector<unsigned short> >				m_table;		//[12560*2][1256]
	vector<vector<pair<unsigned short, unsigned short> > >	m_table_sort;		//[12560*2][1256]
	vector<vector<unsigned short> >				m_table_coarse2dense;
	vector<vector<unsigned short> >				m_table_dense2coarse;
	vector<vector<vector<double> > >			m_rmatrix;
    
	void ssave_matching_table(string pm_filename)
	{
		std::ofstream ofs(pm_filename.c_str());
		boost::archive::text_oarchive oa(ofs);
		oa & m_table;
	}
	

	void ssave_matching_table_sort(string pm_filename)
	{
		std::ofstream ofs(pm_filename.c_str());
		boost::archive::text_oarchive oa(ofs);
		oa & m_table_sort;
	}


	void ssave_matching_table_coarse2dense(string pm_filename)
	{
		std::ofstream ofs(pm_filename.c_str());
		boost::archive::text_oarchive oa(ofs);
		oa & m_table_coarse2dense;
	}

	void ssave_matching_table_dense2coarse(string pm_filename)
	{
		std::ofstream ofs(pm_filename.c_str());
		boost::archive::text_oarchive oa(ofs);
		oa & m_table_dense2coarse;
	}

	void ssave_matching_table_rmatrix(string pm_filename)
	{
		std::ofstream ofs(pm_filename.c_str());
		boost::archive::text_oarchive oa(ofs);
		oa & m_rmatrix;
	}
	


	void sload_matching_table(string pm_filename)
	{
		m_table.clear();
		std::ifstream ifs(pm_filename.c_str(), std::ios::binary);
		boost::archive::text_iarchive la(ifs);
		la & m_table;
	}

	void sload_matching_table_sort(string pm_filename)
	{
		m_table_sort.clear();
		std::ifstream ifs(pm_filename.c_str(), std::ios::binary);
		boost::archive::text_iarchive la(ifs);
		la & m_table_sort;
	}

	void sload_matching_table_coarse2dense(string pm_filename)
	{
		m_table_coarse2dense.clear();
		std::ifstream ifs(pm_filename.c_str(), std::ios::binary);
		boost::archive::text_iarchive la(ifs);
		la & m_table_coarse2dense;
	}

	
	void sload_matching_table_dense2coarse(string pm_filename)
	{
		m_table_dense2coarse.clear();
		std::ifstream ifs(pm_filename.c_str(), std::ios::binary);
		boost::archive::text_iarchive la(ifs);
		la & m_table_dense2coarse;
	}

	void sload_matching_table_rmatrix(string pm_filename)
	{
		m_rmatrix.clear();
		std::ifstream ifs(pm_filename.c_str(), std::ios::binary);
		boost::archive::text_iarchive la(ifs);
		la & m_rmatrix;
	}



	void check_table_sort()
	{
		size_t i, j;

		for (i = 0; i < m_table_sort.size(); i ++){
			for (j = 0; j < m_table_sort[i].size(); j ++){
				if (m_table[i][m_table_sort[i][j].first] != m_table_sort[i][j].second){
					cout << "ERROR: " << m_table[i][m_table_sort[i][j].first] << " " << m_table_sort[i][j].second << endl;
 				}
			}
		}
	}

	void write_bin_matching_table_sort(string pm_filename)
	{
        ofstream out_file (pm_filename.c_str(), ios::out|ios::binary);
		unsigned short* out_buffer = new unsigned short [m_table_sort.size() * m_table_sort[0].size() * 2];
		char* head_buffer = reinterpret_cast<char *>(out_buffer);
		size_t i, j;
		
		for (i = 0; i < m_table_sort.size(); i ++){
			for (j = 0; j < m_table_sort[i].size(); j ++){
				*(out_buffer++) = m_table_sort[i][j].first;
				*(out_buffer++) = m_table_sort[i][j].second;
			}
		}


		unsigned int len = m_table_sort.size();
		out_file.write(reinterpret_cast<char *>(&len), sizeof(len));
		len = m_table_sort[0].size();
		out_file.write(reinterpret_cast<char *>(&len), sizeof(len));

		out_file.write(head_buffer, m_table_sort.size() * m_table_sort[0].size() * 2 * sizeof(unsigned short));
		out_file.close();
	}

	void read_bin_matching_table_sort(string pm_filename)
	{
		m_table_dense2coarse.clear();
        ifstream in_file (pm_filename.c_str(), ios::in | ios::binary | ios::ate);
		ifstream::pos_type in_size = in_file.tellg();
		char* head_buffer = new char [in_size];
		in_file.seekg(0, ios::beg);
		in_file.read(head_buffer, in_size);
		in_file.close();


		unsigned int num_row = *(reinterpret_cast<unsigned int *>(head_buffer));
		unsigned int num_column = *(reinterpret_cast<unsigned int *>(head_buffer+sizeof(unsigned int)));
		m_table_sort.insert(m_table_sort.end(), num_row, vector<pair<unsigned short, unsigned short> >(num_column, pair<unsigned short, unsigned short>(0, 0)));
		
		size_t i, j;
		
		for (i = 0; i < num_row; i ++){
			for (j = 0; j < num_column; j ++){
				m_table_sort[i][j].first = (reinterpret_cast<unsigned short *>(head_buffer+2*sizeof(unsigned int)))[i*num_column*2 + j*2];
				m_table_sort[i][j].second = (reinterpret_cast<unsigned short *>(head_buffer+2*sizeof(unsigned int)))[i*num_column*2 + j*2 + 1];
			}
		}
		delete[] head_buffer;
	}



	void write_bin_matching_table_dense2coarse(string pm_filename)
	{
        ofstream out_file (pm_filename.c_str(), ios::out|ios::binary);
		unsigned short* out_buffer = new unsigned short [m_table_dense2coarse.size() * m_table_dense2coarse[0].size()];
		char* head_buffer = reinterpret_cast<char *>(out_buffer);
		size_t i, j;
		
		for (i = 0; i < m_table_dense2coarse.size(); i ++){
			for (j = 0; j < m_table_dense2coarse[i].size(); j ++){
				*(out_buffer) = m_table_dense2coarse[i][j];
				out_buffer ++;
			}
		}


		unsigned int len = m_table_dense2coarse.size();
		out_file.write(reinterpret_cast<char *>(&len), sizeof(len));
		len = m_table_dense2coarse[0].size();
		out_file.write(reinterpret_cast<char *>(&len), sizeof(len));

		out_file.write(head_buffer, m_table_dense2coarse.size() * m_table_dense2coarse[0].size() * sizeof(unsigned short));
		out_file.close();
	}

	void read_bin_matching_table_dense2coarse()
	{
		m_table_dense2coarse.clear();
		ifstream in_file ("matching_table_dense2coarse.bin", ios::in | ios::binary | ios::ate);
		ifstream::pos_type in_size = in_file.tellg();
		char* head_buffer = new char [in_size];
		in_file.seekg(0, ios::beg);
		in_file.read(head_buffer, in_size);
		in_file.close();


		unsigned int num_row = *(reinterpret_cast<unsigned int *>(head_buffer));
		unsigned int num_column = *(reinterpret_cast<unsigned int *>(head_buffer+sizeof(unsigned int)));
		m_table_dense2coarse.insert(m_table_dense2coarse.end(), num_row, vector<unsigned short>(num_column, 0));
		
		size_t i, j;
		
		for (i = 0; i < num_row; i ++){
			for (j = 0; j < num_column; j ++){
				m_table_dense2coarse[i][j] = (reinterpret_cast<unsigned short *>(head_buffer+2*sizeof(unsigned int)))[i*num_column + j];
			}
		}
		delete[] head_buffer;
	}

	void write_bin_matching_table_coarse2dense(string pm_filename)
	{
        ofstream out_file (pm_filename.c_str(), ios::out|ios::binary);
		unsigned short* out_buffer = new unsigned short [m_table_coarse2dense.size() * m_table_coarse2dense[0].size()];
		char* head_buffer = reinterpret_cast<char *>(out_buffer);
		size_t i, j;
		
		for (i = 0; i < m_table_coarse2dense.size(); i ++){
			for (j = 0; j < m_table_coarse2dense[i].size(); j ++){
				*(out_buffer) = m_table_coarse2dense[i][j];
				out_buffer ++;
			}
		}


		unsigned int len = m_table_coarse2dense.size();
		out_file.write(reinterpret_cast<char *>(&len), sizeof(len));
		len = m_table_coarse2dense[0].size();
		out_file.write(reinterpret_cast<char *>(&len), sizeof(len));

		out_file.write(head_buffer, m_table_coarse2dense.size() * m_table_coarse2dense[0].size() * sizeof(unsigned short));
		out_file.close();
	}

	void read_bin_matching_table_coarse2dense()
	{
		m_table_coarse2dense.clear();
		ifstream in_file ("matching_table_coarse2dense.bin", ios::in | ios::binary | ios::ate);
		ifstream::pos_type in_size = in_file.tellg();
		char* head_buffer = new char [in_size];
		in_file.seekg(0, ios::beg);
		in_file.read(head_buffer, in_size);
		in_file.close();
		
		
		unsigned int num_row = *(reinterpret_cast<unsigned int *>(head_buffer));
		unsigned int num_column = *(reinterpret_cast<unsigned int *>(head_buffer+sizeof(unsigned int)));
		m_table_coarse2dense.insert(m_table_coarse2dense.end(), num_row, vector<unsigned short>(num_column, 0));
		
		size_t i, j;
		
		for (i = 0; i < num_row; i ++){
			for (j = 0; j < num_column; j ++){
				m_table_coarse2dense[i][j] = (reinterpret_cast<unsigned short *>(head_buffer+2*sizeof(unsigned int)))[i*num_column + j];
			}
		}
		delete[] head_buffer;
	}

	void write_bin_matching_table(string pm_filename)
	{
        ofstream out_file (pm_filename.c_str(), ios::out|ios::binary);
		unsigned short* out_buffer = new unsigned short [m_table.size() * m_table[0].size()];
		char* head_buffer = reinterpret_cast<char *>(out_buffer);
		size_t i, j;
		
		for (i = 0; i < m_table.size(); i ++){
			for (j = 0; j < m_table[i].size(); j ++){
				*(out_buffer) = m_table[i][j];
				out_buffer ++;
			}
		}


		unsigned int len = m_table.size();
		out_file.write(reinterpret_cast<char *>(&len), sizeof(len));
		len = m_table[0].size();
		out_file.write(reinterpret_cast<char *>(&len), sizeof(len));

		out_file.write(head_buffer, m_table.size() * m_table[0].size() * sizeof(unsigned short));
		out_file.close();
	}

	void read_bin_matching_table()
	{
		m_table.clear();
		ifstream in_file ("matching_table.bin", ios::in | ios::binary | ios::ate);
		ifstream::pos_type in_size = in_file.tellg();
		char* head_buffer = new char [in_size];
		in_file.seekg(0, ios::beg);
		in_file.read(head_buffer, in_size);
		in_file.close();

		unsigned int num_row = *(reinterpret_cast<unsigned int *>(head_buffer));
		unsigned int num_column = *(reinterpret_cast<unsigned int *>(head_buffer+sizeof(unsigned int)));
		m_table.insert(m_table.end(), num_row, vector<unsigned short>(num_column, 0));
		
		size_t i, j;
		
		for (i = 0; i < num_row; i ++){
			for (j = 0; j < num_column; j ++){
				m_table[i][j] = (reinterpret_cast<unsigned short *>(head_buffer+2*sizeof(unsigned int)))[i*num_column + j];
			}
		}

		delete[] head_buffer;
	}


	void write_bin_rotation_matrix(string pm_filename)
	{
        ofstream out_file (pm_filename.c_str(), ios::out|ios::binary);
		double* out_buffer = new double [m_rmatrix.size() * 9];
		char* head_buffer = reinterpret_cast<char *>(out_buffer);
		size_t i, j, k;
		
		for (i = 0; i < m_rmatrix.size(); i ++){
			for (j = 0; j < m_rmatrix[i].size(); j ++){
				for (k = 0; k < m_rmatrix[i][j].size(); k ++){
					*(out_buffer) = m_rmatrix[i][j][k];
					out_buffer ++;
				}
			}
		}
		unsigned int len = m_rmatrix.size();
		out_file.write(reinterpret_cast<char *>(&len), sizeof(len));
		out_file.write(head_buffer, m_rmatrix.size() * 9 * sizeof(double));
		out_file.close();
	}

	void read_bin_rotation_matrix()
	{
		m_rmatrix.clear();
		ifstream in_file ("rotation_matrix.bin", ios::in | ios::binary | ios::ate);
		ifstream::pos_type in_size = in_file.tellg();
		char* head_buffer = new char [in_size];
		in_file.seekg(0, ios::beg);
		in_file.read(head_buffer, in_size);
		in_file.close();

		size_t num_matrix = *(reinterpret_cast<unsigned int *>(head_buffer));
		m_rmatrix.insert(m_rmatrix.end(), num_matrix, vector<vector<double> >(3, vector<double>(3, 0)));
		size_t i, j, k;
		
		for (i = 0; i < num_matrix; i ++){
			for (j = 0; j < 3; j ++){
				for (k = 0; k < 3; k ++){
					m_rmatrix[i][j][k] = (reinterpret_cast<double *>(head_buffer+sizeof(unsigned int)))[i*9 + j*3 + k];
				}
			}
		}
		delete[] head_buffer;
	}

	void print_rotation_matrix()
	{
		size_t i;
		//stringstream out_buffer;
		//ofstream out_file("rotation_matrix.txt");
		//out_buffer << m_rmatrix.size() << endl;
		FILE*	f = fopen ("rotation_matrix.txt","w");
		fprintf(f, "%d\n", m_rmatrix.size());
		for (i = 0; i < m_rmatrix.size(); i ++){
			fprintf(f, "%+1.16e %+1.16e %+1.16e ",m_rmatrix[i][0][0], m_rmatrix[i][0][1], m_rmatrix[i][0][2]);
			fprintf(f, "%+1.16e %+1.16e %+1.16e ",m_rmatrix[i][1][0], m_rmatrix[i][1][1], m_rmatrix[i][1][2]);
			fprintf(f, "%+1.16e %+1.16e %+1.16e\n",m_rmatrix[i][2][0], m_rmatrix[i][2][1], m_rmatrix[i][2][2]);

			//out_buffer << m_rmatrix[i][0][0] << " ";
			//out_buffer << m_rmatrix[i][0][1] << " ";
			//out_buffer << m_rmatrix[i][0][2] << " ";
			//out_buffer << m_rmatrix[i][1][0] << " ";
			//out_buffer << m_rmatrix[i][1][1] << " ";
			//out_buffer << m_rmatrix[i][1][2] << " ";
			//out_buffer << m_rmatrix[i][2][0] << " ";
			//out_buffer << m_rmatrix[i][2][1] << " ";
			//out_buffer << m_rmatrix[i][2][2] << endl;
		}
		//out_file << out_buffer.str();
		//out_file.close();
	}

	void load_rotation_matrix()
	{
		m_rmatrix.clear();

		ifstream in_file("rotation_matrix.txt");
		int num;
		in_file >> num;

		vector<vector<double> > tmp_matrix(3, vector<double>(3, 0));
		int i;

		for (i = 0; i < num; i ++){
			in_file >> tmp_matrix[0][0] 
			        >> tmp_matrix[0][1] 
					>> tmp_matrix[0][2] 
					>> tmp_matrix[1][0] 
			        >> tmp_matrix[1][1] 
					>> tmp_matrix[1][2] 
					>> tmp_matrix[2][0] 
			        >> tmp_matrix[2][1] 
					>> tmp_matrix[2][2]; 
			m_rmatrix.push_back(tmp_matrix);
		}

		in_file.close();
	}

	void load_matching_table()
	{
		m_table.clear();

		ifstream in_file("matching_table.txt");
		int num_pose, num_coarse;
		in_file >> num_pose >> num_coarse;
		
		vector<unsigned short>	tmp_row;
		unsigned short id_ray;
		int i, j;
		for (i = 0; i < num_pose; i ++){
			tmp_row.clear();
			for (j = 0; j < num_coarse; j ++){
				in_file >> id_ray;
				tmp_row.push_back(id_ray);
			}
			m_table.push_back(tmp_row);
		}
		in_file.close();
	}

	void load_matching_table_coarse2dense()
	{
		m_table_coarse2dense.clear();

		ifstream in_file("matching_table_coarse2dense.txt");
		int num_pose, num_coarse;
		in_file >> num_pose >> num_coarse;
		
		vector<unsigned short>	tmp_row;
		unsigned short id_ray;
		int i, j;
		for (i = 0; i < num_pose; i ++){
			tmp_row.clear();
			for (j = 0; j < num_coarse; j ++){
				in_file >> id_ray;
				tmp_row.push_back(id_ray);
			}
			m_table_coarse2dense.push_back(tmp_row);
		}
		in_file.close();
	}

	void load_matching_table_dense2coarse()
	{
		m_table_dense2coarse.clear();

		ifstream in_file("matching_table_dense2coarse.txt");
		int num_pose, num_coarse;
		in_file >> num_pose >> num_coarse;
		
		vector<unsigned short>	tmp_row;
		unsigned short id_ray;
		int i, j;
		for (i = 0; i < num_pose; i ++){
			tmp_row.clear();
			for (j = 0; j < num_coarse; j ++){
				in_file >> id_ray;
				tmp_row.push_back(id_ray);
			}
			m_table_dense2coarse.push_back(tmp_row);
		}
		in_file.close();
	}

	void init()
	{
		load_unitball();
		m_coarse.init();
		m_dense.init();
		m_pose.init();
//cerr << "m_pose.init_pose();" << endl;
		m_pose.init_pose();
	}	

	void load_unitball()
	{
		vector<double>	tmp_point(3, 0);
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

	}

	void count_matching_table()
	{
		vector<int>	tmp_counter(points_default_dense, 0);
		size_t i, j;
		
		for (i = 0; i < m_table.size(); i ++){
			for (j = 0; j < m_table[i].size(); j ++){
				tmp_counter[m_table[i][j]] += 1;
			}
		}

		sort(tmp_counter.begin(), tmp_counter.end());
		copy(tmp_counter.begin(), tmp_counter.end(), ostream_iterator<unsigned short>(cout, "\n"));
	}

	void print_matching_table_coarse2dense()
	{
		size_t i, j;
		stringstream out_buffer;
		out_buffer << m_table_coarse2dense.size() << " " << m_table_coarse2dense[0].size() << endl;
		for (i = 0; i < m_table_coarse2dense.size(); i ++){
			for (j = 0; j < m_table_coarse2dense[i].size(); j ++){
				out_buffer << m_table_coarse2dense[i][j] << " ";
			}
			out_buffer << endl;
		}

		ofstream out_file("matching_table_coarse2dense.txt");
		out_file << out_buffer.str();
		out_file.close();
	}

	void print_matching_table_dense2coarse()
	{
		size_t i, j;
		stringstream out_buffer;
		out_buffer << m_table_dense2coarse.size() << " " << m_table_dense2coarse[0].size() << endl;
		for (i = 0; i < m_table_dense2coarse.size(); i ++){
			for (j = 0; j < m_table_dense2coarse[i].size(); j ++){
				out_buffer << m_table_dense2coarse[i][j] << " ";
			}
			out_buffer << endl;
		}

		ofstream out_file("matching_table_dense2coarse.txt");
		out_file << out_buffer.str();
		out_file.close();
	}

	void print_matching_table()
	{
		size_t i, j;
		stringstream out_buffer;
		out_buffer << m_table.size() << " " << m_table[0].size() << endl;
		for (i = 0; i < m_table.size(); i ++){
			for (j = 0; j < m_table[i].size(); j ++){
				out_buffer << m_table[i][j] << " ";
			}
			out_buffer << endl;
		}

		ofstream out_file("matching_table.txt");
		out_file << out_buffer.str();
		out_file.close();
	}

	void compute_matching_table()
	{
		m_table.clear();
		m_rmatrix.clear();

		vector<unsigned short>	tmp_pose;
		vector<pair<double, pair<unsigned short, unsigned short> > >	tmp_pose_sort;
		vector<pair<unsigned short, unsigned short> >	tmp_pose_sort_row;
		vector<unsigned short>	tmp_coarse2dense;
		vector<unsigned short>	tmp_dense2coarse;
		unsigned short tmp_id_dense;
		size_t i, j;
		//double M[3][3];
		//RotationMatrix M;
		vector<vector<double> > M(3, vector<double>(3, 0));
		vector<double>	tmp_point_old(3, 0);
		vector<double>	tmp_point_new(3, 0);
		double			tmp_local_r_sq = 0.4;
		double			tmp_dist;
		
		for (i = 0; i < m_pose.m_arcs.size(); i ++){
if (i/100*100 == i) cerr << "i = " << i << endl;			
			
			M[0][0] = m_pose.m_arcs[i][0][0];
			M[1][0] = m_pose.m_arcs[i][0][1];
			M[2][0] = m_pose.m_arcs[i][0][2];

			M[0][1] = m_pose.m_arcs[i][1][0];
			M[1][1] = m_pose.m_arcs[i][1][1];
			M[2][1] = m_pose.m_arcs[i][1][2];

			M[0][2] = m_pose.m_arcs[i][2][0];
			M[1][2] = m_pose.m_arcs[i][2][1];
			M[2][2] = m_pose.m_arcs[i][2][2];

			tmp_dist  = M[0][2] * M[0][2];
			tmp_dist += M[1][2] * M[1][2];
			tmp_dist += (M[2][2] - 1) * (M[2][2] - 1);
			if (tmp_dist > tmp_local_r_sq) continue;
		
			m_rmatrix.push_back(M);
			tmp_pose.clear();
			tmp_pose_sort.clear();

			for (j = 0; j < points_default_coarse; j++){
				tmp_point_old[0] = vertex_default_coarse[j][0];
				tmp_point_old[1] = vertex_default_coarse[j][1];
				tmp_point_old[2] = vertex_default_coarse[j][2];
				tmp_point_new[0] = M[0][0] * tmp_point_old[0] + M[0][1] * tmp_point_old[1] + M[0][2] * tmp_point_old[2];
				tmp_point_new[1] = M[1][0] * tmp_point_old[0] + M[1][1] * tmp_point_old[1] + M[1][2] * tmp_point_old[2];
				tmp_point_new[2] = M[2][0] * tmp_point_old[0] + M[2][1] * tmp_point_old[1] + M[2][2] * tmp_point_old[2];
				tmp_id_dense = m_dense.query_nn(tmp_point_new);
				tmp_pose.push_back(tmp_id_dense);
				tmp_pose_sort.push_back(pair<double, pair<unsigned short, unsigned short> >(2 * tmp_point_new[2] - tmp_point_old[2], pair<unsigned short, unsigned short>(j, tmp_id_dense)));
			}
			m_table.push_back(tmp_pose);

			sort(tmp_pose_sort.begin(), tmp_pose_sort.end(), less_than_univ<double, pair<unsigned short, unsigned short> >());
			tmp_pose_sort_row.clear();
			for (j = 0; j < tmp_pose_sort.size(); j ++){
				tmp_pose_sort_row.push_back(tmp_pose_sort[j].second);
			}
			m_table_sort.push_back(tmp_pose_sort_row);	
			
			tmp_coarse2dense.clear();
			for (j = 0; j < points_default_pose; j++){
				tmp_point_old[0] = vertex_default_pose[j][0];
				tmp_point_old[1] = vertex_default_pose[j][1];
				tmp_point_old[2] = vertex_default_pose[j][2];
				tmp_point_new[0] = M[0][0] * tmp_point_old[0] + M[0][1] * tmp_point_old[1] + M[0][2] * tmp_point_old[2];
				tmp_point_new[1] = M[1][0] * tmp_point_old[0] + M[1][1] * tmp_point_old[1] + M[1][2] * tmp_point_old[2];
				tmp_point_new[2] = M[2][0] * tmp_point_old[0] + M[2][1] * tmp_point_old[1] + M[2][2] * tmp_point_old[2];
				tmp_coarse2dense.push_back(m_pose.query_nn(tmp_point_new));
			}
			m_table_coarse2dense.push_back(tmp_coarse2dense);
			tmp_dense2coarse = tmp_coarse2dense;
			for (j = 0; j < tmp_coarse2dense.size(); j ++){
				tmp_dense2coarse[tmp_coarse2dense[j]] = j;
			}
			m_table_dense2coarse.push_back(tmp_dense2coarse);
		}
	}


    //will generate full size table
	void compute_matching_table_v6()
	{
		m_table.clear();
		m_rmatrix.clear();

		vector<unsigned short>	tmp_pose;
		vector<unsigned short>	tmp_coarse2dense;
		vector<unsigned short>	tmp_dense2coarse;
		unsigned short tmp_id_dense;
		size_t i, j;
		vector<vector<double> > M(3, vector<double>(3, 0));
		vector<double>	tmp_point_old(3, 0);
		vector<double>	tmp_point_new(3, 0);
		
		for (i = 0; i < m_pose.m_arcs.size(); i ++){
if (i/100*100 == i) cerr << "i = " << i << endl;			
			
			M[0][0] = m_pose.m_arcs[i][0][0];
			M[1][0] = m_pose.m_arcs[i][0][1];
			M[2][0] = m_pose.m_arcs[i][0][2];

			M[0][1] = m_pose.m_arcs[i][1][0];
			M[1][1] = m_pose.m_arcs[i][1][1];
			M[2][1] = m_pose.m_arcs[i][1][2];

			M[0][2] = m_pose.m_arcs[i][2][0];
			M[1][2] = m_pose.m_arcs[i][2][1];
			M[2][2] = m_pose.m_arcs[i][2][2];

			m_rmatrix.push_back(M);
			tmp_pose.clear();

			for (j = 0; j < points_default_coarse; j++){
				tmp_point_old[0] = vertex_default_coarse[j][0];
				tmp_point_old[1] = vertex_default_coarse[j][1];
				tmp_point_old[2] = vertex_default_coarse[j][2];
				tmp_point_new[0] = M[0][0] * tmp_point_old[0] + M[0][1] * tmp_point_old[1] + M[0][2] * tmp_point_old[2];
				tmp_point_new[1] = M[1][0] * tmp_point_old[0] + M[1][1] * tmp_point_old[1] + M[1][2] * tmp_point_old[2];
				tmp_point_new[2] = M[2][0] * tmp_point_old[0] + M[2][1] * tmp_point_old[1] + M[2][2] * tmp_point_old[2];
				tmp_id_dense = m_dense.query_nn(tmp_point_new);
				tmp_pose.push_back(tmp_id_dense);
			}
			m_table.push_back(tmp_pose);

			tmp_coarse2dense.clear();
			for (j = 0; j < points_default_dense; j++){
				tmp_point_old[0] = vertex_default_dense[j][0];
				tmp_point_old[1] = vertex_default_dense[j][1];
				tmp_point_old[2] = vertex_default_dense[j][2];
				tmp_point_new[0] = M[0][0] * tmp_point_old[0] + M[0][1] * tmp_point_old[1] + M[0][2] * tmp_point_old[2];
				tmp_point_new[1] = M[1][0] * tmp_point_old[0] + M[1][1] * tmp_point_old[1] + M[1][2] * tmp_point_old[2];
				tmp_point_new[2] = M[2][0] * tmp_point_old[0] + M[2][1] * tmp_point_old[1] + M[2][2] * tmp_point_old[2];
				tmp_coarse2dense.push_back(m_dense.query_nn(tmp_point_new));
			}
			m_table_coarse2dense.push_back(tmp_coarse2dense);
			tmp_dense2coarse = tmp_coarse2dense;
			for (j = 0; j < tmp_coarse2dense.size(); j ++){
				tmp_dense2coarse[tmp_coarse2dense[j]] = j;
			}
			m_table_dense2coarse.push_back(tmp_dense2coarse);
		}

        cout << "m_table.size() = " << m_table.size() << endl;
        cout << "m_rmatrix.size() = " << m_rmatrix.size() << endl;
        cout << "m_table_coarse2dense.size() = " << m_table_coarse2dense.size() << endl;
        cout << "m_table_dense2coarse.size() = " << m_table_dense2coarse.size() << endl;
	}


	void compute_matching_table_v6_2(double pm_angle_degree)
	{
		m_table.clear();
		m_rmatrix.clear();

		vector<pair<double, pair<unsigned short, unsigned short> > >	tmp_pose_sort(points_default_coarse, pair<double, pair<unsigned short, unsigned short> >(0, pair<unsigned short, unsigned short>(0, 0)));
		vector<pair<unsigned short, unsigned short> >	tmp_pose_sort_row(points_default_coarse, pair<unsigned short, unsigned short>(0, 0));

		vector<unsigned short>	tmp_pose(points_default_coarse, 0);
		vector<unsigned short>	tmp_coarse2dense(points_default_dense, 0);
		vector<unsigned short>	tmp_dense2coarse(points_default_dense, 0);
		unsigned short tmp_id_dense;
		size_t i, j;
		vector<vector<double> > M(3, vector<double>(3, 0));
		vector<double>	tmp_point_old(3, 0);
		vector<double>	tmp_point_new(3, 0);
		
		double			tmp_local_r_sq = sin(M_PI*pm_angle_degree/360)*sin(M_PI*pm_angle_degree/360)*4;; //0.268 = 30 degree, 0.362 = 35, 0.468 = 40, 0.586 = 45, 
		double			tmp_dist;

		m_table_sort.clear();

        cout << "Please wait ..." << endl;
		for (i = 0; i < m_pose.m_arcs.size(); i ++){
//if (i/100*100 == i) cerr << "i = " << i << endl;			
			
			M[0][0] = m_pose.m_arcs[i][0][0];
			M[1][0] = m_pose.m_arcs[i][0][1];
			M[2][0] = m_pose.m_arcs[i][0][2];

			M[0][1] = m_pose.m_arcs[i][1][0];
			M[1][1] = m_pose.m_arcs[i][1][1];
			M[2][1] = m_pose.m_arcs[i][1][2];

			M[0][2] = m_pose.m_arcs[i][2][0];
			M[1][2] = m_pose.m_arcs[i][2][1];
			M[2][2] = m_pose.m_arcs[i][2][2];

			tmp_dist  = M[0][2] * M[0][2];
			tmp_dist += M[1][2] * M[1][2];
			tmp_dist += (M[2][2] - 1) * (M[2][2] - 1);
			if (tmp_dist > tmp_local_r_sq) continue;

			m_rmatrix.push_back(M);

			for (j = 0; j < points_default_coarse; j++){
				tmp_point_old[0] = vertex_default_coarse[j][0];
				tmp_point_old[1] = vertex_default_coarse[j][1];
				tmp_point_old[2] = vertex_default_coarse[j][2];
				tmp_point_new[0] = M[0][0] * tmp_point_old[0] + M[0][1] * tmp_point_old[1] + M[0][2] * tmp_point_old[2];
				tmp_point_new[1] = M[1][0] * tmp_point_old[0] + M[1][1] * tmp_point_old[1] + M[1][2] * tmp_point_old[2];
				tmp_point_new[2] = M[2][0] * tmp_point_old[0] + M[2][1] * tmp_point_old[1] + M[2][2] * tmp_point_old[2];
				tmp_id_dense = m_dense.query_nn(tmp_point_new);
				tmp_pose[j] = tmp_id_dense;
				tmp_pose_sort[j] = pair<double, pair<unsigned short, unsigned short> >(2 * tmp_point_new[2] - tmp_point_old[2], pair<unsigned short, unsigned short>(j, tmp_id_dense));
			}
			m_table.push_back(tmp_pose);

			sort(tmp_pose_sort.begin(), tmp_pose_sort.end(), less_than_univ<double, pair<unsigned short, unsigned short> >());
			for (j = 0; j < tmp_pose_sort.size(); j ++){
				tmp_pose_sort_row[j] = tmp_pose_sort[j].second;
			}
			m_table_sort.push_back(tmp_pose_sort_row);	

			for (j = 0; j < points_default_dense; j++){
				tmp_point_old[0] = vertex_default_dense[j][0];
				tmp_point_old[1] = vertex_default_dense[j][1];
				tmp_point_old[2] = vertex_default_dense[j][2];
				tmp_point_new[0] = M[0][0] * tmp_point_old[0] + M[0][1] * tmp_point_old[1] + M[0][2] * tmp_point_old[2];
				tmp_point_new[1] = M[1][0] * tmp_point_old[0] + M[1][1] * tmp_point_old[1] + M[1][2] * tmp_point_old[2];
				tmp_point_new[2] = M[2][0] * tmp_point_old[0] + M[2][1] * tmp_point_old[1] + M[2][2] * tmp_point_old[2];
				tmp_coarse2dense[j] = m_dense.query_nn(tmp_point_new);
			}
			m_table_coarse2dense.push_back(tmp_coarse2dense);

			for (j = 0; j < tmp_coarse2dense.size(); j ++){
				tmp_dense2coarse[tmp_coarse2dense[j]] = j;
			}
			m_table_dense2coarse.push_back(tmp_dense2coarse);
		}

        cout << "m_table.size() = " << m_table.size() << endl;
        cout << "m_rmatrix.size() = " << m_rmatrix.size() << endl;
        cout << "m_table_coarse2dense.size() = " << m_table_coarse2dense.size() << endl;
        cout << "m_table_dense2coarse.size() = " << m_table_dense2coarse.size() << endl;
	}

};

void task_map_sphere2circle(int argc, char *argv[])
{
	if (argc < 3){
		cerr << "USAGE: " << argv[0] << " <out_vtk_filename>" << endl;
		return;
	}

    MyUnitBall myBall;

    myBall.init();

	//m_unitball_dense
	size_t i;
	double x, y, z, r, s, u, v;
	vector<vector<double> >	tmp_points;
	vector<double>			tmp_point(2, 0);

	for (i = 0; i < myBall.m_unitball_dense.size(); i ++){
		x = myBall.m_unitball_dense[i][0];
		y = myBall.m_unitball_dense[i][1];
		z = myBall.m_unitball_dense[i][2];
		if (y < 0 || z < 0 || x < 0 || x == -1) continue;
		s = 0.5 * (1 - x);
		r = 2 * sqrt(1 - s);
		u = y / r;
		v = z / r;
		tmp_point[0] = u;
		tmp_point[1] = v;
		tmp_points.push_back(tmp_point);
	}

	//write to vtk file
	string out_fn = argv[2];
	ofstream out_file(out_fn.c_str());
	stringstream out_buffer;

	out_buffer << "# vtk DataFile Version 3.0" << endl;
	out_buffer << argv[2] << endl;
	out_buffer << "ASCII" << endl;
	out_buffer << "DATASET POLYDATA" << endl;
	out_buffer << "POINTS " << tmp_points.size() << " float" << endl;

	size_t j;
	for (j = 0; j < tmp_points.size(); j ++){
		out_buffer << tmp_points[j][0] << " "
				   << tmp_points[j][1] << " "
			       << 0 << endl;
	}

	out_buffer << "VERTICES " << tmp_points.size() << " " << tmp_points.size() * 2 << endl;
	for(j = 0; j < tmp_points.size(); j ++){
		out_buffer << "1 " << j << endl;
	}

	out_file << out_buffer.str();
	out_file.close();


}

void task_mtgen_v6(int argc, char *argv[])
{
    MyUnitBall myBall;

    myBall.init();
cerr << "myBall.compute_matching_table_v6(); " << endl;
	myBall.compute_matching_table_v6();

cerr << "serializing the tables..." << endl;
	myBall.ssave_matching_table("ssave_matching_table_v6.txt");
	myBall.ssave_matching_table_rmatrix("ssave_matching_table_v6_rmatrix.txt");
	//myBall.ssave_matching_table_sort("ssave_matching_table_v6_sort.txt");
	myBall.ssave_matching_table_dense2coarse("ssave_matching_table_v6_dense2coarse.txt");
	myBall.ssave_matching_table_coarse2dense("ssave_matching_table_v6_coarse2dense.txt");


//cerr << "myBall.write_bin_rotation_matrix();" << endl;
//	myBall.write_bin_rotation_matrix("matching_table_rotation_matrix_v6.bin");
//
//cerr << "myBall.write_bin_matching_table();" << endl;
//	myBall.write_bin_matching_table("matching_table_v6.bin");
//
//cerr << "myBall.write_bin_matching_table_coarse2dense();" << endl;
//	myBall.write_bin_matching_table_coarse2dense("matching_table_coarse2dense_v6.bin");
//	
//cerr << "myBall.write_bin_matching_table_dense2coarse();" << endl;	
//	myBall.write_bin_matching_table_dense2coarse("matching_table_dense2coarse_v6.bin");

}


void task_mtgen_v6_2(int argc, char *argv[])
{
    MyUnitBall myBall;

    myBall.init();
//cerr << "myBall.compute_matching_table_v6_2(); " << endl;
	myBall.compute_matching_table_v6_2(atof(argv[1]));

		string fn_table = string("ssave_matching_table_") + string(argv[1]) + string(".txt");
		string fn_table_rmatrix = string("ssave_matching_table_") + string(argv[1]) + string("_rmatrix.txt");
		string fn_table_sort = string("ssave_matching_table_") + string(argv[1]) + string("_sort.txt");
		string fn_table_dense2coarse = string("ssave_matching_table_") + string(argv[1]) + string("_dense2coarse.txt");
		string fn_table_coarse2dense = string("ssave_matching_table_") + string(argv[1]) + string("_coarse2dense.txt");
		
cerr << "serializing the tables..." << endl;
	myBall.ssave_matching_table(fn_table.c_str());
	myBall.ssave_matching_table_rmatrix(fn_table_rmatrix.c_str());
	myBall.ssave_matching_table_sort(fn_table_sort.c_str());
	myBall.ssave_matching_table_dense2coarse(fn_table_dense2coarse.c_str());
	myBall.ssave_matching_table_coarse2dense(fn_table_coarse2dense.c_str());
cerr << "Done!" << endl;

//cerr << "myBall.write_bin_rotation_matrix();" << endl;
//	myBall.write_bin_rotation_matrix("matching_table_rotation_matrix_v6.bin");
//
//cerr << "myBall.write_bin_matching_table();" << endl;
//	myBall.write_bin_matching_table("matching_table_v6.bin");
//
//cerr << "myBall.write_bin_matching_table_coarse2dense();" << endl;
//	myBall.write_bin_matching_table_coarse2dense("matching_table_coarse2dense_v6.bin");
//	
//cerr << "myBall.write_bin_matching_table_dense2coarse();" << endl;	
//	myBall.write_bin_matching_table_dense2coarse("matching_table_dense2coarse_v6.bin");

}


void task_mtgen(int argc, char *argv[])
{
	MyUnitBall myBall;

    myBall.init();
cerr << "myBall.compute_matching_table(); " << endl;
	myBall.compute_matching_table();

cerr << "serializing the tables..." << endl;
	myBall.ssave_matching_table("ssave_matching_table_full.txt");
	myBall.ssave_matching_table_rmatrix("ssave_matching_table_full_rmatrix.txt");
	myBall.ssave_matching_table_sort("ssave_matching_table_full_sort.txt");
	myBall.ssave_matching_table_dense2coarse("ssave_matching_table_full_dense2coarse.txt");
	myBall.ssave_matching_table_coarse2dense("ssave_matching_table_full_coarse2dense.txt");

//cerr << "myBall.write_bin_rotation_matrix();" << endl;
//	myBall.write_bin_rotation_matrix("rotation_matrix.bin");
//
//cerr << "myBall.write_bin_matching_table();" << endl;
//	myBall.write_bin_matching_table("matching_table.bin");
//
//cerr << "myBall.write_bin_matching_table_coarse2dense();" << endl;
//	myBall.write_bin_matching_table_coarse2dense("matching_table_coarse2dense.bin");
//	
//cerr << "myBall.write_bin_matching_table_dense2coarse();" << endl;	
//	myBall.write_bin_matching_table_dense2coarse("matching_table_dense2coarse.bin");
//	
//cerr << "myBall.write_bin_matching_table_sort();" << endl;	
//	myBall.write_bin_matching_table_sort("matching_table_sort.bin");	

}

int main(int argc, char *argv[])
{
    //MyUnitBall myBall;
    //myBall.init();

		if (argc != 2){
				cerr << "Usage: " << argv[0] << " <angle>" << endl;
				return -1; 
		}
		
		task_mtgen_v6_2(argc, argv);

/*		
    if (argc < 2){
        task_mtgen(argc, argv);
    }
    else if (atoi(argv[1]) == 6){
        task_mtgen_v6(argc, argv);
    }
    else if (atoi(argv[1]) == 62){
        task_mtgen_v6_2(argc, argv);
    }
		else if (atoi(argv[1]) == 5){
				task_map_sphere2circle(argc, argv);
		}
*/

/*
	myBall.load_rotation_matrix();
	myBall.load_matching_table();
	myBall.load_matching_table_coarse2dense();
	myBall.load_matching_table_dense2coarse();

	myBall.write_bin_rotation_matrix();
	myBall.write_bin_matching_table();
	myBall.write_bin_matching_table_coarse2dense();
	myBall.write_bin_matching_table_dense2coarse();
*/

/*
	cout << "myBall.read_bin_matching_table();" << endl;
	myBall.read_bin_matching_table();

	//cout << "myBall.print_matching_table();" << endl;
	//myBall.print_matching_table();

	cout << "myBall.read_bin_matching_table();" << endl;
	myBall.read_bin_matching_table_sort();

	cout << "myBall.read_bin_matching_table_coarse2dense();" << endl;
	myBall.read_bin_matching_table_coarse2dense();

	//cout << "myBall.print_matching_table_coarse2dense();" << endl;
	//myBall.print_matching_table_coarse2dense();

	cout << "myBall.read_bin_matching_table_dense2coarse();" << endl;
	myBall.read_bin_matching_table_dense2coarse();

	//cout << "myBall.print_matching_table_dense2coarse();" << endl;
	//myBall.print_matching_table_dense2coarse();

	cout << "myBall.read_bin_rotation_matrix();" << endl;
	myBall.read_bin_rotation_matrix();

	//cout << "myBall.print_rotation_matrix();" << endl;
	//myBall.print_rotation_matrix();

	myBall.check_table_sort();
*/


    return 0;
}
