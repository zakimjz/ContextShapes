#ifndef _MY_COMPARE_SPD_H_
#define _MY_COMPARE_SPD_H_
#include <vector>
#include "global.h"
#include "contextBall.h"
#include "gridSES.h"
#include "util.h"

using namespace std;

extern double weight_in_1_16bits [0x1u << 16] ;
extern double weight_in_2_16bits [0x1u << 16] ;

inline double precomputed32_weight (CONTEXT_RAY_INT_TYPE n) {
    return weight_in_1_16bits[n         & 0x0000ffffu]
    +  weight_in_2_16bits[(n >> 16) & 0x0000ffffu] ;
}

void* mt_subtask_match(void* pm_subtask);

struct MyCompare_SPD;

struct MySubtask{
        unsigned int                                                                            m_subtask_id;
        MyUnitBall*                                                                                     m_unitball;
        MyCompare_SPD*                                                                          m_compare;
    vector<pair<MyContextBall*, MyContextBall*> >               m_cb_pairs;
        vector<MyMatchPair*>                                                            m_rtn_pairs;

};

struct MyCompare_SPD{
    MyGridSES_Dense_SPD*        m_DB_dense; //from DB, densely sampled
    MyGridSES_Coarse_SPD*       m_Query_coarse; //from Query, coarsely sampled, only interaction area used.
    MyUnitBall*					m_unitball;

    //the following for debug purpose only
    int         m_hit_core_3A;
    int         m_hit_core_1A;
    int         m_hit_shell_1A;
    int         m_hit_shell_2A;
    int         m_hit_shell_3A;
    unsigned long long         m_non_hit;

    vector<vector<double> > convert_matrix(double pm_matrix[3][3] )
    {
        vector<vector<double> > M(3, vector<double>(3, 0));
        M[0][0] = pm_matrix[0][0];
        M[0][1] = pm_matrix[0][1];
        M[0][2] = pm_matrix[0][2];

        M[1][0] = pm_matrix[1][0];
        M[1][1] = pm_matrix[1][1];
        M[1][2] = pm_matrix[1][2];

        M[2][0] = pm_matrix[2][0];
        M[2][1] = pm_matrix[2][1];
        M[2][2] = pm_matrix[2][2];

        return M;
    }

	void subtask_match(MySubtask* pm_subtask)
	{
		cout << g_match_pair_str << " SUBTASK_" << pm_subtask->m_subtask_id << ": Starting..." << endl;
		size_t i;
		MyMatchPair*		tmp_pair;
        MyContextBall*      curr_ball_query;
        MyContextBall*      curr_ball_db;

        double time_begin_progress = 0, time_end_progress = 0;
        double time_sum_progress = 0;
        struct timeb tp_progress;
		int	   num_progress = 0;

		ftime(&tp_progress);
        time_begin_progress = tp_progress.time + (tp_progress.millitm * 1.0) / 1000;

		for (i = 0; i < pm_subtask->m_cb_pairs.size(); i ++){
            //cout << "i = " << i << " out of " << pm_subtask->m_cb_pairs.size() << endl;
			curr_ball_db = pm_subtask->m_cb_pairs[i].first; 
			curr_ball_query = pm_subtask->m_cb_pairs[i].second;
			num_progress ++;

			if (num_progress % 2000 == 0){
				cout << g_match_pair_str << " SUBTASK_" << pm_subtask->m_subtask_id << ": progress: " << num_progress << " out of " << pm_subtask->m_cb_pairs.size();
                ftime(&tp_progress);
                time_end_progress = tp_progress.time + (tp_progress.millitm * 1.0) / 1000;
                cout << " \t***ToGO: ";   
                int num_total = pm_subtask->m_cb_pairs.size();
				time_sum_progress = time_end_progress - time_begin_progress;
				cout << num_total * time_sum_progress / (num_progress * 1.0) - time_sum_progress << " seconds." << endl;
            }

            
			tmp_pair = compute_match_pair(curr_ball_db, curr_ball_query);
			if (!tmp_pair) continue;
			pm_subtask->m_rtn_pairs.push_back(tmp_pair);
		}
		cout << g_match_pair_str << " SUBTASK_" << pm_subtask->m_subtask_id << ": Finished!" << endl;
	}


	//
	MyCompare_SPD(MyGridSES_Dense_SPD* A, MyGridSES_Coarse_SPD* B, MyUnitBall* unitball)
    {
        m_DB_dense = A;
        m_Query_coarse = B;
        m_unitball = unitball;

        m_hit_core_3A = 0;
        m_hit_core_1A = 0;
        m_hit_shell_1A = 0;
        m_hit_shell_2A = 0;
        m_hit_shell_3A = 0;
        m_non_hit = 0;

    }

    ~MyCompare_SPD()
    {
    }


	//multithreading version
    void match_SPD_mt()
    {
        double time_begin = 0, time_end = 0;
        struct timeb tp;

  //      double time_begin_progress = 0, time_end_progress = 0;
  //      double time_sum_progress = 0;
  //      //int    time_count_progress = 0;
  //      struct timeb tp_progress;
		//int	   num_progress = 0;

        size_t i, j, k;
        MyContextBall*      curr_ball_query;
        MyContextBall*      curr_ball_db;

        vector<MyContextBall*>      pool_ball_query;
        vector<MyContextBall*>      pool_ball_db;
		vector<pair<MyContextBall*, MyContextBall*> >      pool_cb_pairs;

        MyMatchPair*									tmp_pair;
        vector<pair<double, MyMatchPair*> >          heap_buffer;

        //compute rmsd
        vector<vector<double> >     M_DB(3, vector<double>(3, 0));
        vector<vector<double> >     M_Query(3, vector<double>(3, 0));
        vector<vector<double> >     M_Rotation(3, vector<double>(3, 0));
        vector<vector<double> >     M_DB_R(3, vector<double>(3, 0));
        vector<double>              V_translate_query_to_origin(3, 0);    
        vector<vector<double> >     M_rotate_query(3, vector<double>(3, 0));
        vector<double>              V_translate_query_to_db(3, 0);        


        double          tmp_dist;
        //double          tmp_rmsd_ca = 0;
        //double          tmp_rmsd_local = 0;
        //double          tmp_rmsd_all = 0;
        unsigned int    tmp_pose_id;

        ftime(&tp);
        time_begin = tp.time + (tp.millitm * 1.0) / 1000;


        cout << "begin matching ...." << endl;
        
        
	
		pool_ball_db = m_DB_dense->m_context_balls;
		pool_ball_query = m_Query_coarse->m_context_balls;
		for (i = 0; i < m_Query_coarse->m_context_balls_refine.size(); i ++){
			for (j = 0; j < m_Query_coarse->m_context_balls_refine[i].size(); j ++){
				pool_ball_query.push_back(m_Query_coarse->m_context_balls_refine[i][j]);
			}
		}
		cout << "pool_ball_db.size() = " << pool_ball_db.size() << endl;
		cout << "pool_ball_query.size() = " << pool_ball_query.size() << endl; 

		double tmp_sa = 0;


		//create threads
        for (i = 0; i < pool_ball_query.size(); ++i){
            curr_ball_query = pool_ball_query[i];
            for (j = 0; j < pool_ball_db.size(); ++j){
				curr_ball_db = pool_ball_db[j];

				tmp_sa = 0;
				for (k = 1; k < 5; k ++){
					tmp_sa += curr_ball_db->m_solid_angle[k];
					tmp_sa += curr_ball_query->m_solid_angle[k];
				}
				tmp_sa /= 4.0;
				if (tmp_sa < g_solid_angle_lower_bound || tmp_sa > g_solid_angle_upper_bound) continue;
				
				pool_cb_pairs.push_back(pair<MyContextBall*, MyContextBall*>(curr_ball_db, curr_ball_query));
            }
        }

		vector<pthread_t>	tmp_threads;
		vector<MySubtask*>	tmp_subtasks;
		for (i = 0; i < g_num_threads; i ++){
			tmp_subtasks.push_back(new MySubtask());
			tmp_threads.push_back(pthread_t());
		}

		//partition the task
		unsigned int tmp_job_size = pool_cb_pairs.size() / g_num_threads;
		//cout << "tmp_job_size = " << tmp_job_size << endl;

		for (i = 0; i < g_num_threads; i ++){
			for (j = 0; j < tmp_job_size; j ++){
				tmp_subtasks[i]->m_cb_pairs.push_back(pool_cb_pairs[i*tmp_job_size + j]);
			}
			tmp_subtasks[i]->m_unitball = m_unitball;
			tmp_subtasks[i]->m_subtask_id = i;
			tmp_subtasks[i]->m_compare = this;
		}

		i = 0;
		if (pool_cb_pairs.size() > g_num_threads * (pool_cb_pairs.size() / g_num_threads)){
			for (j = tmp_job_size * g_num_threads; j < pool_cb_pairs.size(); j ++){
				tmp_subtasks[i++]->m_cb_pairs.push_back(pool_cb_pairs[j]);
			}
		}

		for (i = 0; i < g_num_threads; i ++){
			pthread_create( &tmp_threads[i], NULL, mt_subtask_match, tmp_subtasks[i]);
		}

		for (i = 0; i < g_num_threads; i ++){
			pthread_join(tmp_threads[i], NULL);
		}

		cout << "INFO: outputing the rankings..." << endl;
		heap_buffer.clear();
		for (i = 0; i < tmp_subtasks.size(); i ++){
			for (j = 0; j < tmp_subtasks[i]->m_rtn_pairs.size(); j ++)
			heap_buffer.push_back(pair<double, MyMatchPair*>(0 - tmp_subtasks[i]->m_rtn_pairs[j]->m_score, tmp_subtasks[i]->m_rtn_pairs[j]));
		}
        sort(heap_buffer.begin(), heap_buffer.end(), less_than_univ<double, MyMatchPair*>());

		//open the out_file
		string out_filename = g_dir_output
            + m_DB_dense->m_surface->m_filename_pdb.substr(0, m_DB_dense->m_surface->m_filename_pdb.length() - 4)
            + "_vs_"
            + m_Query_coarse->m_surface->m_filename_pdb.substr(0, m_Query_coarse->m_surface->m_filename_pdb.length() - 4)
            + "_SQL_INSERT.sql";
        ofstream out_file(out_filename.c_str());
        if (!out_file.is_open()){
            cerr << "ERROR: cannot open file " << out_filename << endl;
            exit(1);
        }
        print_SQL_CREATE_TABLE(out_file);


		for (i = 0; i < heap_buffer.size(); i ++){
			if (i > g_num_predicates) break;
			tmp_pair = heap_buffer[i].second;
            tmp_dist = tmp_pair->m_dist = rdl_vector_ssd(tmp_pair->m_ball_db->m_ball_center, tmp_pair->m_ball_query->m_ball_center);
            tmp_pose_id = tmp_pair->m_pose_id;

            V_translate_query_to_origin[0] = 0 - tmp_pair->m_ball_query->m_ball_center[0];
            V_translate_query_to_origin[1] = 0 - tmp_pair->m_ball_query->m_ball_center[1];
            V_translate_query_to_origin[2] = 0 - tmp_pair->m_ball_query->m_ball_center[2];

            V_translate_query_to_db[0] = tmp_pair->m_ball_db->m_ball_center[0];
            V_translate_query_to_db[1] = tmp_pair->m_ball_db->m_ball_center[1];
            V_translate_query_to_db[2] = tmp_pair->m_ball_db->m_ball_center[2];

            M_DB[0][0] = tmp_pair->m_ball_db->m_rotation_matrix[0][0];
            M_DB[0][1] = tmp_pair->m_ball_db->m_rotation_matrix[0][1];
            M_DB[0][2] = tmp_pair->m_ball_db->m_rotation_matrix[0][2];
            M_DB[1][0] = tmp_pair->m_ball_db->m_rotation_matrix[1][0];
            M_DB[1][1] = tmp_pair->m_ball_db->m_rotation_matrix[1][1];
            M_DB[1][2] = tmp_pair->m_ball_db->m_rotation_matrix[1][2];
            M_DB[2][0] = tmp_pair->m_ball_db->m_rotation_matrix[2][0];
            M_DB[2][1] = tmp_pair->m_ball_db->m_rotation_matrix[2][1];
            M_DB[2][2] = tmp_pair->m_ball_db->m_rotation_matrix[2][2];

            M_Query[0][0] = tmp_pair->m_ball_query->m_rotation_matrix[0][0];
            M_Query[0][1] = tmp_pair->m_ball_query->m_rotation_matrix[0][1];
            M_Query[0][2] = tmp_pair->m_ball_query->m_rotation_matrix[0][2];
            M_Query[1][0] = tmp_pair->m_ball_query->m_rotation_matrix[1][0];
            M_Query[1][1] = tmp_pair->m_ball_query->m_rotation_matrix[1][1];
            M_Query[1][2] = tmp_pair->m_ball_query->m_rotation_matrix[1][2];
            M_Query[2][0] = tmp_pair->m_ball_query->m_rotation_matrix[2][0];
            M_Query[2][1] = tmp_pair->m_ball_query->m_rotation_matrix[2][1];
            M_Query[2][2] = tmp_pair->m_ball_query->m_rotation_matrix[2][2];

            M_Rotation[0][0] = m_unitball->m_rmatrix[tmp_pose_id][0][0];
            M_Rotation[0][1] = m_unitball->m_rmatrix[tmp_pose_id][0][1];
            M_Rotation[0][2] = m_unitball->m_rmatrix[tmp_pose_id][0][2];
            M_Rotation[1][0] = m_unitball->m_rmatrix[tmp_pose_id][1][0];
            M_Rotation[1][1] = m_unitball->m_rmatrix[tmp_pose_id][1][1];
            M_Rotation[1][2] = m_unitball->m_rmatrix[tmp_pose_id][1][2];
            M_Rotation[2][0] = m_unitball->m_rmatrix[tmp_pose_id][2][0];
            M_Rotation[2][1] = m_unitball->m_rmatrix[tmp_pose_id][2][1];
            M_Rotation[2][2] = m_unitball->m_rmatrix[tmp_pose_id][2][2];
            M_DB_R = rdl_matrix_multiplication_3d(M_DB, M_Rotation);
            M_rotate_query = rdl_matrix_multiplication_3d(M_DB_R, rdl_matrix_transpose_3d(M_Query));
                
            if (tmp_dist <= g_distance_true_pair + 1){
                tmp_pair->m_rmsd_ca = compute_rmsd_ca(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
                tmp_pair->m_rmsd_all = compute_rmsd_all(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
                tmp_pair->m_rmsd_local = compute_rmsd_local(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
            }
            else {
                tmp_pair->m_rmsd_ca = tmp_pair->m_rmsd_all = tmp_pair->m_rmsd_local = 9999;
            }
            tmp_pair->m_translate_query_to_origin = V_translate_query_to_origin;
            tmp_pair->m_rotate_query = M_rotate_query;
            tmp_pair->m_translate_query_to_db = V_translate_query_to_db;

			print_SQL_INSERT_row(tmp_pair, out_file);
		}

		for (i = 0; i < tmp_subtasks.size(); i ++){
			delete tmp_subtasks[i];
		}

        print_SQL_COMMIT(out_file);
        out_file.close();

        ftime(&tp);
        time_end = tp.time + (tp.millitm * 1.0) / 1000;
        cout << "******INFO: time to match " << time_end - time_begin << " seconds." << endl;       

        cout << "end matching ...." << endl;
    }

    void match_SPD_v1_mt()
    {
        double time_begin = 0, time_end = 0;
        struct timeb tp;

        size_t i, j, k;
        MyContextBall*      curr_ball_query;
        MyContextBall*      curr_ball_db;
		vector<pair<MyContextBall*, MyContextBall*> >      pool_cb_pairs;
        
        MyMatchPair*                    tmp_pair;
        vector<MyMatchPair*>            tmp_pairs;
        vector<pair<double, MyMatchPair*> >             heap_buffer;

        //compute rmsd
        //vector<double>              V_Translation(3, 0);
        //vector<double>              V_Touch(3, 0);
        vector<vector<double> >     M_DB(3, vector<double>(3, 0));
        vector<vector<double> >     M_Query(3, vector<double>(3, 0));
        vector<vector<double> >     M_Rotation(3, vector<double>(3, 0));
        vector<vector<double> >     M_DB_R(3, vector<double>(3, 0));
        vector<double>                  V_translate_query_to_origin(3, 0);    
        vector<vector<double> >         M_rotate_query(3, vector<double>(3, 0));
        vector<double>                  V_translate_query_to_db(3, 0);        

        double          tmp_dist;
        unsigned int    tmp_pose_id;
		double			tmp_sa;

        ftime(&tp);
        time_begin = tp.time + (tp.millitm * 1.0) / 1000;

        cout << "begin matching ...." << endl;
        heap_buffer.clear();
        for (i = 0; i < m_Query_coarse->m_context_balls.size(); ++i){
            curr_ball_query = m_Query_coarse->m_context_balls[i];
            for (j = 0; j < m_DB_dense->m_context_balls.size(); ++j){
                curr_ball_db = m_DB_dense->m_context_balls[j];
				
				//filter based on solid angle
				tmp_sa = 0;
				for (k = 1; k < 5; k ++){
					tmp_sa += curr_ball_db->m_solid_angle[k];
					tmp_sa += curr_ball_query->m_solid_angle[k];
				}
				tmp_sa /= 4.0;
				if (tmp_sa < g_solid_angle_lower_bound || tmp_sa > g_solid_angle_upper_bound) continue;
			
				pool_cb_pairs.push_back(pair<MyContextBall*, MyContextBall*>(curr_ball_db, curr_ball_query));
            }
        }

		vector<pthread_t>	tmp_threads;
		vector<MySubtask*>	tmp_subtasks;
		for (i = 0; i < g_num_threads; i ++){
			tmp_subtasks.push_back(new MySubtask());
			tmp_threads.push_back(pthread_t());
		}

		unsigned int tmp_job_size = pool_cb_pairs.size() / g_num_threads;
		for (i = 0; i < g_num_threads; i ++){
			for (j = 0; j < tmp_job_size; j ++){
				tmp_subtasks[i]->m_cb_pairs.push_back(pool_cb_pairs[i*tmp_job_size + j]);
			}
			tmp_subtasks[i]->m_unitball = m_unitball;
			tmp_subtasks[i]->m_subtask_id = i;
			tmp_subtasks[i]->m_compare = this;
		}

		i = 0;
		if (pool_cb_pairs.size() > g_num_threads * (pool_cb_pairs.size() / g_num_threads)){
			for (j = tmp_job_size * g_num_threads; j < pool_cb_pairs.size(); j ++){
				tmp_subtasks[i++]->m_cb_pairs.push_back(pool_cb_pairs[j]);
			}
		}

		for (i = 0; i < g_num_threads; i ++){
			pthread_create( &tmp_threads[i], NULL, mt_subtask_match, tmp_subtasks[i]);
		}

		for (i = 0; i < g_num_threads; i ++){
			pthread_join(tmp_threads[i], NULL);
		}

		cout << "INFO: outputing the rankings..." << endl;
		heap_buffer.clear();
		for (i = 0; i < tmp_subtasks.size(); i ++){
			for (j = 0; j < tmp_subtasks[i]->m_rtn_pairs.size(); j ++)
			heap_buffer.push_back(pair<double, MyMatchPair*>(0 - tmp_subtasks[i]->m_rtn_pairs[j]->m_score, tmp_subtasks[i]->m_rtn_pairs[j]));
		}
        sort(heap_buffer.begin(), heap_buffer.end(), less_than_univ<double, MyMatchPair*>());


        size_t query_ball_id;
        //size_t query_ball_sn;

		pool_cb_pairs.clear();
        for (i = 0; i < heap_buffer.size(); i ++){
            if (i > g_num_predicates) break;
            query_ball_id = heap_buffer[i].second->m_ball_query->m_ball_id;
            curr_ball_db = heap_buffer[i].second->m_ball_db;
            for (j = 0; j < m_Query_coarse->m_context_balls_refine[query_ball_id].size(); j ++){
                curr_ball_query = m_Query_coarse->m_context_balls_refine[query_ball_id][j];
				pool_cb_pairs.push_back(pair<MyContextBall*, MyContextBall*>(curr_ball_db, curr_ball_query));
            }
        }

		tmp_threads.clear();
		for (i = 0; i < g_num_threads; i ++){
			delete tmp_subtasks[i];
		}
		tmp_subtasks.clear();
		for (i = 0; i < g_num_threads; i ++){
			tmp_subtasks.push_back(new MySubtask());
			tmp_threads.push_back(pthread_t());
		}

		tmp_job_size = pool_cb_pairs.size() / g_num_threads;
		for (i = 0; i < g_num_threads; i ++){
			for (j = 0; j < tmp_job_size; j ++){
				tmp_subtasks[i]->m_cb_pairs.push_back(pool_cb_pairs[i*tmp_job_size + j]);
			}
			tmp_subtasks[i]->m_unitball = m_unitball;
			tmp_subtasks[i]->m_subtask_id = i;
			tmp_subtasks[i]->m_compare = this;
		}

		i = 0;
		if (pool_cb_pairs.size() > g_num_threads * (pool_cb_pairs.size() / g_num_threads)){
			for (j = tmp_job_size * g_num_threads; j < pool_cb_pairs.size(); j ++){
				tmp_subtasks[i++]->m_cb_pairs.push_back(pool_cb_pairs[j]);
			}
		}

		for (i = 0; i < g_num_threads; i ++){
			pthread_create( &tmp_threads[i], NULL, mt_subtask_match, tmp_subtasks[i]);
		}

		for (i = 0; i < g_num_threads; i ++){
			pthread_join(tmp_threads[i], NULL);
		}

		for (i = 0; i < tmp_subtasks.size(); i ++){
			for (j = 0; j < tmp_subtasks[i]->m_rtn_pairs.size(); j ++)
			heap_buffer.push_back(pair<double, MyMatchPair*>(0 - tmp_subtasks[i]->m_rtn_pairs[j]->m_score, tmp_subtasks[i]->m_rtn_pairs[j]));
		}
        sort(heap_buffer.begin(), heap_buffer.end(), less_than_univ<double, MyMatchPair*>());

        string out_filename = g_dir_output
            + m_DB_dense->m_surface->m_filename_pdb.substr(0, m_DB_dense->m_surface->m_filename_pdb.length() - 4)
            + "_vs_"
            + m_Query_coarse->m_surface->m_filename_pdb.substr(0, m_Query_coarse->m_surface->m_filename_pdb.length() - 4)
            + "_SQL_INSERT.sql";
        ofstream out_file(out_filename.c_str());
        print_SQL_CREATE_TABLE(out_file);

		for (i = 0; i < heap_buffer.size(); i ++){
			if (i > g_num_predicates) break;
			tmp_pair = heap_buffer[i].second;
            tmp_dist = tmp_pair->m_dist = rdl_vector_ssd(tmp_pair->m_ball_db->m_ball_center, tmp_pair->m_ball_query->m_ball_center);
            tmp_pose_id = tmp_pair->m_pose_id;

            V_translate_query_to_origin[0] = 0 - tmp_pair->m_ball_query->m_ball_center[0];
            V_translate_query_to_origin[1] = 0 - tmp_pair->m_ball_query->m_ball_center[1];
            V_translate_query_to_origin[2] = 0 - tmp_pair->m_ball_query->m_ball_center[2];

            V_translate_query_to_db[0] = tmp_pair->m_ball_db->m_ball_center[0];
            V_translate_query_to_db[1] = tmp_pair->m_ball_db->m_ball_center[1];
            V_translate_query_to_db[2] = tmp_pair->m_ball_db->m_ball_center[2];

            M_DB[0][0] = tmp_pair->m_ball_db->m_rotation_matrix[0][0];
            M_DB[0][1] = tmp_pair->m_ball_db->m_rotation_matrix[0][1];
            M_DB[0][2] = tmp_pair->m_ball_db->m_rotation_matrix[0][2];
            M_DB[1][0] = tmp_pair->m_ball_db->m_rotation_matrix[1][0];
            M_DB[1][1] = tmp_pair->m_ball_db->m_rotation_matrix[1][1];
            M_DB[1][2] = tmp_pair->m_ball_db->m_rotation_matrix[1][2];
            M_DB[2][0] = tmp_pair->m_ball_db->m_rotation_matrix[2][0];
            M_DB[2][1] = tmp_pair->m_ball_db->m_rotation_matrix[2][1];
            M_DB[2][2] = tmp_pair->m_ball_db->m_rotation_matrix[2][2];

            M_Query[0][0] = tmp_pair->m_ball_query->m_rotation_matrix[0][0];
            M_Query[0][1] = tmp_pair->m_ball_query->m_rotation_matrix[0][1];
            M_Query[0][2] = tmp_pair->m_ball_query->m_rotation_matrix[0][2];
            M_Query[1][0] = tmp_pair->m_ball_query->m_rotation_matrix[1][0];
            M_Query[1][1] = tmp_pair->m_ball_query->m_rotation_matrix[1][1];
            M_Query[1][2] = tmp_pair->m_ball_query->m_rotation_matrix[1][2];
            M_Query[2][0] = tmp_pair->m_ball_query->m_rotation_matrix[2][0];
            M_Query[2][1] = tmp_pair->m_ball_query->m_rotation_matrix[2][1];
            M_Query[2][2] = tmp_pair->m_ball_query->m_rotation_matrix[2][2];

            M_Rotation[0][0] = m_unitball->m_rmatrix[tmp_pose_id][0][0];
            M_Rotation[0][1] = m_unitball->m_rmatrix[tmp_pose_id][0][1];
            M_Rotation[0][2] = m_unitball->m_rmatrix[tmp_pose_id][0][2];
            M_Rotation[1][0] = m_unitball->m_rmatrix[tmp_pose_id][1][0];
            M_Rotation[1][1] = m_unitball->m_rmatrix[tmp_pose_id][1][1];
            M_Rotation[1][2] = m_unitball->m_rmatrix[tmp_pose_id][1][2];
            M_Rotation[2][0] = m_unitball->m_rmatrix[tmp_pose_id][2][0];
            M_Rotation[2][1] = m_unitball->m_rmatrix[tmp_pose_id][2][1];
            M_Rotation[2][2] = m_unitball->m_rmatrix[tmp_pose_id][2][2];
            M_DB_R = rdl_matrix_multiplication_3d(M_DB, M_Rotation);
            M_rotate_query = rdl_matrix_multiplication_3d(M_DB_R, rdl_matrix_transpose_3d(M_Query));
                
            //if (tmp_dist <= g_distance_true_pair + 1){
            //    tmp_pair->m_rmsd_ca = compute_rmsd_ca(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
            //    tmp_pair->m_rmsd_all = compute_rmsd_all(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
            //    tmp_pair->m_rmsd_local = compute_rmsd_local(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
            //}
            //else {
            //    tmp_pair->m_rmsd_ca = tmp_pair->m_rmsd_all = tmp_pair->m_rmsd_local = 9999;
            //}
			tmp_pair->m_rmsd_ca = tmp_pair->m_rmsd_all = tmp_pair->m_rmsd_local = 9999;

            tmp_pair->m_translate_query_to_origin = V_translate_query_to_origin;
            tmp_pair->m_rotate_query = M_rotate_query;
            tmp_pair->m_translate_query_to_db = V_translate_query_to_db;

			print_SQL_INSERT_row(tmp_pair, out_file);
		}

		for (i = 0; i < tmp_subtasks.size(); i ++){
			delete tmp_subtasks[i];
		}

        print_SQL_COMMIT(out_file);
        out_file.close();

        ftime(&tp);
        time_end = tp.time + (tp.millitm * 1.0) / 1000;
        cout << "******INFO: time to match " << time_end - time_begin << " seconds." << endl;       

        cout << "end matching ...." << endl;
    }

    void match_SPD()
    {
        double time_begin = 0, time_end = 0;
        struct timeb tp;

        double time_begin_progress = 0, time_end_progress = 0;
        double time_sum_progress = 0;
        struct timeb tp_progress;
		int	   num_progress = 0;

        size_t i, j, k;
        MyContextBall*      curr_ball_query;
        MyContextBall*      curr_ball_db;

        vector<MyContextBall*>      pool_ball_query;
        vector<MyContextBall*>      pool_ball_db;

        MyMatchPair*									tmp_pair;
        vector<pair<double, MyMatchPair*> >          heap_buffer;

        //compute rmsd
        vector<vector<double> >     M_DB(3, vector<double>(3, 0));
        vector<vector<double> >     M_Query(3, vector<double>(3, 0));
        vector<vector<double> >     M_Rotation(3, vector<double>(3, 0));
        vector<vector<double> >     M_DB_R(3, vector<double>(3, 0));
        vector<double>              V_translate_query_to_origin(3, 0);    
        vector<vector<double> >     M_rotate_query(3, vector<double>(3, 0));
        vector<double>              V_translate_query_to_db(3, 0);        


        double          tmp_dist;
        //double          tmp_rmsd_ca = 0;
        //double          tmp_rmsd_local = 0;
        //double          tmp_rmsd_all = 0;
        unsigned int    tmp_pose_id;

        ftime(&tp);
        time_begin = tp.time + (tp.millitm * 1.0) / 1000;

        ftime(&tp_progress);
        time_begin_progress = tp_progress.time + (tp_progress.millitm * 1.0) / 1000;


        string out_filename = g_dir_output
            + m_DB_dense->m_surface->m_filename_pdb.substr(0, m_DB_dense->m_surface->m_filename_pdb.length() - 4)
            + "_vs_"
            + m_Query_coarse->m_surface->m_filename_pdb.substr(0, m_Query_coarse->m_surface->m_filename_pdb.length() - 4)
            + "_SQL_INSERT.sql";

        ofstream out_file(out_filename.c_str());
        if (!out_file.is_open()){
            cerr << "ERROR: cannot open file " << out_filename << endl;
            exit(1);
        }

        print_SQL_CREATE_TABLE(out_file);

        cout << "begin matching ...." << endl;
        
        heap_buffer.clear();
	
		pool_ball_db = m_DB_dense->m_context_balls;
		pool_ball_query = m_Query_coarse->m_context_balls;
		for (i = 0; i < m_Query_coarse->m_context_balls_refine.size(); i ++){
			for (j = 0; j < m_Query_coarse->m_context_balls_refine[i].size(); j ++){
				pool_ball_query.push_back(m_Query_coarse->m_context_balls_refine[i][j]);
			}
		}
		cout << "pool_ball_db.size() = " << pool_ball_db.size() << endl;
		cout << "pool_ball_query.size() = " << pool_ball_query.size() << endl; 

		double tmp_sa = 0;

        for (i = 0; i < pool_ball_query.size(); ++i){
            curr_ball_query = pool_ball_query[i];
            for (j = 0; j < pool_ball_db.size(); ++j){
				curr_ball_db = pool_ball_db[j];
                num_progress ++;

				if (num_progress % 2000 == 0){
                    cout << "progress: " << num_progress << " out of " << pool_ball_query.size() * pool_ball_db.size();
                    ftime(&tp_progress);
                    time_end_progress = tp_progress.time + (tp_progress.millitm * 1.0) / 1000;
                    cout << " \t***ToGO: ";   
                    int num_total = pool_ball_query.size() * pool_ball_db.size();
					time_sum_progress = time_end_progress - time_begin_progress;
					cout << num_total * time_sum_progress / (num_progress * 1.0) - time_sum_progress << " seconds." << endl;
                }

				tmp_sa = 0;
				for (k = 1; k < 5; k ++){
					tmp_sa += curr_ball_db->m_solid_angle[k];
					tmp_sa += curr_ball_query->m_solid_angle[k];
				}
				tmp_sa /= 4.0;
				if (tmp_sa < g_solid_angle_lower_bound || tmp_sa > g_solid_angle_upper_bound) continue;

				tmp_pair = compute_match_pair(curr_ball_db, curr_ball_query);

				if (!tmp_pair) continue;
				heap_buffer.push_back(pair<double, MyMatchPair*>(0 - tmp_pair->m_score, tmp_pair));
            }
        }

        sort(heap_buffer.begin(), heap_buffer.end(), less_than_univ<double, MyMatchPair*>());
		
		for (i = 0; i < heap_buffer.size(); i ++){
			if (i > g_num_predicates) break;
			tmp_pair = heap_buffer[i].second;
            tmp_dist = tmp_pair->m_dist = rdl_vector_ssd(tmp_pair->m_ball_db->m_ball_center, tmp_pair->m_ball_query->m_ball_center);
            tmp_pose_id = tmp_pair->m_pose_id;

            V_translate_query_to_origin[0] = 0 - tmp_pair->m_ball_query->m_ball_center[0];
            V_translate_query_to_origin[1] = 0 - tmp_pair->m_ball_query->m_ball_center[1];
            V_translate_query_to_origin[2] = 0 - tmp_pair->m_ball_query->m_ball_center[2];

            V_translate_query_to_db[0] = tmp_pair->m_ball_db->m_ball_center[0];
            V_translate_query_to_db[1] = tmp_pair->m_ball_db->m_ball_center[1];
            V_translate_query_to_db[2] = tmp_pair->m_ball_db->m_ball_center[2];

            M_DB[0][0] = tmp_pair->m_ball_db->m_rotation_matrix[0][0];
            M_DB[0][1] = tmp_pair->m_ball_db->m_rotation_matrix[0][1];
            M_DB[0][2] = tmp_pair->m_ball_db->m_rotation_matrix[0][2];
            M_DB[1][0] = tmp_pair->m_ball_db->m_rotation_matrix[1][0];
            M_DB[1][1] = tmp_pair->m_ball_db->m_rotation_matrix[1][1];
            M_DB[1][2] = tmp_pair->m_ball_db->m_rotation_matrix[1][2];
            M_DB[2][0] = tmp_pair->m_ball_db->m_rotation_matrix[2][0];
            M_DB[2][1] = tmp_pair->m_ball_db->m_rotation_matrix[2][1];
            M_DB[2][2] = tmp_pair->m_ball_db->m_rotation_matrix[2][2];

            M_Query[0][0] = tmp_pair->m_ball_query->m_rotation_matrix[0][0];
            M_Query[0][1] = tmp_pair->m_ball_query->m_rotation_matrix[0][1];
            M_Query[0][2] = tmp_pair->m_ball_query->m_rotation_matrix[0][2];
            M_Query[1][0] = tmp_pair->m_ball_query->m_rotation_matrix[1][0];
            M_Query[1][1] = tmp_pair->m_ball_query->m_rotation_matrix[1][1];
            M_Query[1][2] = tmp_pair->m_ball_query->m_rotation_matrix[1][2];
            M_Query[2][0] = tmp_pair->m_ball_query->m_rotation_matrix[2][0];
            M_Query[2][1] = tmp_pair->m_ball_query->m_rotation_matrix[2][1];
            M_Query[2][2] = tmp_pair->m_ball_query->m_rotation_matrix[2][2];

            M_Rotation[0][0] = m_unitball->m_rmatrix[tmp_pose_id][0][0];
            M_Rotation[0][1] = m_unitball->m_rmatrix[tmp_pose_id][0][1];
            M_Rotation[0][2] = m_unitball->m_rmatrix[tmp_pose_id][0][2];
            M_Rotation[1][0] = m_unitball->m_rmatrix[tmp_pose_id][1][0];
            M_Rotation[1][1] = m_unitball->m_rmatrix[tmp_pose_id][1][1];
            M_Rotation[1][2] = m_unitball->m_rmatrix[tmp_pose_id][1][2];
            M_Rotation[2][0] = m_unitball->m_rmatrix[tmp_pose_id][2][0];
            M_Rotation[2][1] = m_unitball->m_rmatrix[tmp_pose_id][2][1];
            M_Rotation[2][2] = m_unitball->m_rmatrix[tmp_pose_id][2][2];
            M_DB_R = rdl_matrix_multiplication_3d(M_DB, M_Rotation);
            M_rotate_query = rdl_matrix_multiplication_3d(M_DB_R, rdl_matrix_transpose_3d(M_Query));
                
            if (tmp_dist <= g_distance_true_pair + 1){
                tmp_pair->m_rmsd_ca = compute_rmsd_ca(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
                tmp_pair->m_rmsd_all = compute_rmsd_all(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
                tmp_pair->m_rmsd_local = compute_rmsd_local(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
            }
            else {
                tmp_pair->m_rmsd_ca = tmp_pair->m_rmsd_all = tmp_pair->m_rmsd_local = 9999;
            }
            tmp_pair->m_translate_query_to_origin = V_translate_query_to_origin;
            tmp_pair->m_rotate_query = M_rotate_query;
            tmp_pair->m_translate_query_to_db = V_translate_query_to_db;

			print_SQL_INSERT_row(tmp_pair, out_file);
		}

        print_SQL_COMMIT(out_file);
        out_file.close();

        ftime(&tp);
        time_end = tp.time + (tp.millitm * 1.0) / 1000;
        cout << "******INFO: time to match " << time_end - time_begin << " seconds." << endl;       

        cout << "end matching ...." << endl;
    }

    

    void match_SPD_test()
    {
        double time_begin = 0, time_end = 0;
        struct timeb tp;

        double time_begin_progress = 0, time_end_progress = 0;
        double time_sum_progress = 0;
        struct timeb tp_progress;
		int	   num_progress = 0;

        size_t i, j, k;
        MyContextBall*      curr_ball_query;
        MyContextBall*      curr_ball_db;

        vector<MyContextBall*>      pool_ball_query;
        vector<MyContextBall*>      pool_ball_db;

        MyMatchPair*								tmp_pair;
        vector<pair<double, MyMatchPair*> >         heap_buffer;

        //compute rmsd
        vector<vector<double> >     M_DB(3, vector<double>(3, 0));
        vector<vector<double> >     M_Query(3, vector<double>(3, 0));
        vector<vector<double> >     M_Rotation(3, vector<double>(3, 0));
        vector<vector<double> >     M_DB_R(3, vector<double>(3, 0));
        vector<double>              V_translate_query_to_origin(3, 0);    
        vector<vector<double> >     M_rotate_query(3, vector<double>(3, 0));
        vector<double>              V_translate_query_to_db(3, 0);        


        double          tmp_dist;
        //double          tmp_rmsd_ca = 0;
        //double          tmp_rmsd_local = 0;
        //double          tmp_rmsd_all = 0;
        unsigned int    tmp_pose_id;

        ftime(&tp);
        time_begin = tp.time + (tp.millitm * 1.0) / 1000;

        ftime(&tp_progress);
        time_begin_progress = tp_progress.time + (tp_progress.millitm * 1.0) / 1000;

        cout << "begin matching ...." << endl;
        
        heap_buffer.clear();
	
		pool_ball_db = m_DB_dense->m_context_balls;
		pool_ball_query = m_Query_coarse->m_context_balls;
		for (i = 0; i < m_Query_coarse->m_context_balls_refine.size(); i ++){
			for (j = 0; j < m_Query_coarse->m_context_balls_refine[i].size(); j ++){
				pool_ball_query.push_back(m_Query_coarse->m_context_balls_refine[i][j]);
			}
		}
		cout << "pool_ball_db.size() = " << pool_ball_db.size() << endl;
		cout << "pool_ball_query.size() = " << pool_ball_query.size() << endl; 

		double tmp_sa = 0;


        for (i = 0; i < pool_ball_query.size(); ++i){
            curr_ball_query = pool_ball_query[i];
            for (j = 0; j < pool_ball_db.size(); ++j){
				curr_ball_db = pool_ball_db[j];
                num_progress ++;

				if (num_progress % 2000 == 0){
                    cout << "progress: " << num_progress << " out of " << pool_ball_query.size() * pool_ball_db.size();
                    ftime(&tp_progress);
                    time_end_progress = tp_progress.time + (tp_progress.millitm * 1.0) / 1000;
                    cout << " \t***ToGO: ";   
                    int num_total = pool_ball_query.size() * pool_ball_db.size();
					time_sum_progress = time_end_progress - time_begin_progress;
					cout << num_total * time_sum_progress / (num_progress * 1.0) - time_sum_progress << " seconds." << endl;
                }

				tmp_sa = 0;
				for (k = 1; k < 5; k ++){
					tmp_sa += curr_ball_db->m_solid_angle[k];
					tmp_sa += curr_ball_query->m_solid_angle[k];
				}
				tmp_sa /= 4.0;
				if (tmp_sa < g_solid_angle_lower_bound || tmp_sa > g_solid_angle_upper_bound) continue;

				tmp_pair = compute_match_pair_test(curr_ball_db, curr_ball_query);
            }
        }

		ftime(&tp);
        time_end = tp.time + (tp.millitm * 1.0) / 1000;
        cout << "******INFO: time to match " << time_end - time_begin << " seconds." << endl;       

        cout << "end matching ...." << endl;

    }

    

	void match_SPD_1pair(unsigned int pm_cb_receptor_id, unsigned int pm_cb_ligand_id, unsigned int pm_cb_ligand_idsn)
    {
        double time_begin = 0, time_end = 0;
        struct timeb tp;

        size_t i;
        MyContextBall*      curr_ball_query;
        MyContextBall*      curr_ball_db;

		//compute rmsd
        vector<vector<double> >     M_DB(3, vector<double>(3, 0));
        vector<vector<double> >     M_Query(3, vector<double>(3, 0));
        vector<vector<double> >     M_Rotation(3, vector<double>(3, 0));
        vector<vector<double> >     M_DB_R(3, vector<double>(3, 0));
        vector<double>              V_translate_query_to_origin(3, 0);    
        vector<vector<double> >     M_rotate_query(3, vector<double>(3, 0));
        vector<double>              V_translate_query_to_db(3, 0);        


        double          tmp_dist;
        unsigned int    tmp_pose_id;

        ftime(&tp);
        time_begin = tp.time + (tp.millitm * 1.0) / 1000;

        string out_filename = g_dir_output
            + m_DB_dense->m_surface->m_filename_pdb.substr(0, m_DB_dense->m_surface->m_filename_pdb.length() - 4)
            + "_cb" + rdl_x2str(pm_cb_receptor_id)
			+ "_vs_"
            + m_Query_coarse->m_surface->m_filename_pdb.substr(0, m_Query_coarse->m_surface->m_filename_pdb.length() - 4)
            + "_cb" + rdl_x2str(pm_cb_ligand_id) + "-" + rdl_x2str(pm_cb_ligand_idsn)
			+ "_SQL_INSERT.sql";

        ofstream out_file(out_filename.c_str());
        if (!out_file.is_open()){
            cerr << "ERROR: cannot open file " << out_filename << endl;
            exit(1);
        }

        print_SQL_CREATE_TABLE(out_file);


		curr_ball_db = m_DB_dense->m_context_balls[pm_cb_receptor_id];
		if (pm_cb_ligand_idsn == 0){
			curr_ball_query = m_Query_coarse->m_context_balls[pm_cb_ligand_id];
		}
		else{
			curr_ball_query = m_Query_coarse->m_context_balls_refine[pm_cb_ligand_id][pm_cb_ligand_idsn-1];
		}

		tmp_dist = rdl_vector_ssd(curr_ball_db->m_ball_center, curr_ball_query->m_ball_center);

		vector<MyMatchPair*>	tmp_pairs;
		MyMatchPair*	tmp_pair;

		compute_match_1pair(curr_ball_db, curr_ball_query, tmp_pairs);

		for (i = 0; i < tmp_pairs.size(); i ++){
			tmp_pair = tmp_pairs[i];
			tmp_pose_id = tmp_pair->m_pose_id;
			tmp_pair->m_dist = tmp_dist;

			V_translate_query_to_origin[0] = 0 - tmp_pair->m_ball_query->m_ball_center[0];
			V_translate_query_to_origin[1] = 0 - tmp_pair->m_ball_query->m_ball_center[1];
			V_translate_query_to_origin[2] = 0 - tmp_pair->m_ball_query->m_ball_center[2];

			V_translate_query_to_db[0] = tmp_pair->m_ball_db->m_ball_center[0];
			V_translate_query_to_db[1] = tmp_pair->m_ball_db->m_ball_center[1];
			V_translate_query_to_db[2] = tmp_pair->m_ball_db->m_ball_center[2];

			M_DB[0][0] = tmp_pair->m_ball_db->m_rotation_matrix[0][0];
			M_DB[0][1] = tmp_pair->m_ball_db->m_rotation_matrix[0][1];
			M_DB[0][2] = tmp_pair->m_ball_db->m_rotation_matrix[0][2];
			M_DB[1][0] = tmp_pair->m_ball_db->m_rotation_matrix[1][0];
			M_DB[1][1] = tmp_pair->m_ball_db->m_rotation_matrix[1][1];
			M_DB[1][2] = tmp_pair->m_ball_db->m_rotation_matrix[1][2];
			M_DB[2][0] = tmp_pair->m_ball_db->m_rotation_matrix[2][0];
			M_DB[2][1] = tmp_pair->m_ball_db->m_rotation_matrix[2][1];
			M_DB[2][2] = tmp_pair->m_ball_db->m_rotation_matrix[2][2];

			M_Query[0][0] = tmp_pair->m_ball_query->m_rotation_matrix[0][0];
			M_Query[0][1] = tmp_pair->m_ball_query->m_rotation_matrix[0][1];
			M_Query[0][2] = tmp_pair->m_ball_query->m_rotation_matrix[0][2];
			M_Query[1][0] = tmp_pair->m_ball_query->m_rotation_matrix[1][0];
			M_Query[1][1] = tmp_pair->m_ball_query->m_rotation_matrix[1][1];
			M_Query[1][2] = tmp_pair->m_ball_query->m_rotation_matrix[1][2];
			M_Query[2][0] = tmp_pair->m_ball_query->m_rotation_matrix[2][0];
			M_Query[2][1] = tmp_pair->m_ball_query->m_rotation_matrix[2][1];
			M_Query[2][2] = tmp_pair->m_ball_query->m_rotation_matrix[2][2];

			M_Rotation[0][0] = m_unitball->m_rmatrix[tmp_pose_id][0][0];
			M_Rotation[0][1] = m_unitball->m_rmatrix[tmp_pose_id][0][1];
			M_Rotation[0][2] = m_unitball->m_rmatrix[tmp_pose_id][0][2];
			M_Rotation[1][0] = m_unitball->m_rmatrix[tmp_pose_id][1][0];
			M_Rotation[1][1] = m_unitball->m_rmatrix[tmp_pose_id][1][1];
			M_Rotation[1][2] = m_unitball->m_rmatrix[tmp_pose_id][1][2];
			M_Rotation[2][0] = m_unitball->m_rmatrix[tmp_pose_id][2][0];
			M_Rotation[2][1] = m_unitball->m_rmatrix[tmp_pose_id][2][1];
			M_Rotation[2][2] = m_unitball->m_rmatrix[tmp_pose_id][2][2];
			M_DB_R = rdl_matrix_multiplication_3d(M_DB, M_Rotation);
			M_rotate_query = rdl_matrix_multiplication_3d(M_DB_R, rdl_matrix_transpose_3d(M_Query));
	            
			tmp_pair->m_rmsd_ca = compute_rmsd_ca(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
			tmp_pair->m_rmsd_all = compute_rmsd_all(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
			tmp_pair->m_rmsd_local = compute_rmsd_local(tmp_pair->m_ball_query, V_translate_query_to_origin, M_rotate_query, V_translate_query_to_db);
			tmp_pair->m_translate_query_to_origin = V_translate_query_to_origin;
			tmp_pair->m_rotate_query = M_rotate_query;
			tmp_pair->m_translate_query_to_db = V_translate_query_to_db;

			print_SQL_INSERT_row(tmp_pair, out_file);

		}

        print_SQL_COMMIT(out_file);
        out_file.close();

        ftime(&tp);
        time_end = tp.time + (tp.millitm * 1.0) / 1000;
        cout << "******INFO: time to match " << time_end - time_begin << " seconds." << endl;       

        cout << "end matching ...." << endl;
    }

	void print_SQL_COMMIT(ofstream& out_file)
    {
        out_file << "COMMIT;" << endl;
    }
  
	void print_SQL_CREATE_TABLE(ofstream& out_file)
    {
        out_file << "BEGIN TRANSACTION;" << endl;
        out_file << "CREATE TABLE match_rankings (" << endl;
		out_file << "\tPDB_pair       string," << endl;	
        out_file << "\tball_id_DB       INTEGER," << endl;
        out_file << "\tball_id_Query    INTEGER," << endl;
        out_file << "\tball_id_sn_Query INTEGER," << endl;
        out_file << "\tpose_id          INTEGER," << endl;
        out_file << "\tdist_balls       DOUBLE," << endl;
        out_file << "\tscore            DOUBLE," << endl;

        out_file << "\trmsd_ca          DOUBLE," << endl;
        out_file << "\trmsd_local       DOUBLE," << endl;
        out_file << "\trmsd_all       DOUBLE," << endl;

		out_file << "\t" << "bov_inInner2A_CORE     DOUBLE," << endl;
		out_file << "\t" << "bov_inCore		  DOUBLE," << endl;
		out_file << "\t" << "bov_inInner4A    DOUBLE," << endl;
		out_file << "\t" << "bov_inInner3A    DOUBLE," << endl;
		out_file << "\t" << "bov_inInner2A    DOUBLE," << endl;
		out_file << "\t" << "bov_inInner1A    DOUBLE," << endl;


        out_file << "\tburied_area_inner_inDB      DOUBLE," << endl;
        out_file << "\tburied_area_core_4A_inDB      DOUBLE," << endl;
        out_file << "\tburied_area_core_3A_inDB      DOUBLE," << endl;
        out_file << "\tburied_area_core_2A_inDB      DOUBLE," << endl;
        out_file << "\tburied_area_core_1A_inDB      DOUBLE," << endl;
        out_file << "\tburied_area_shell_1A_inDB      DOUBLE," << endl;
        out_file << "\tburied_area_shell_2A_inDB      DOUBLE," << endl;
        out_file << "\tburied_area_shell_3A_inDB      DOUBLE," << endl;
        out_file << "\tburied_area_shell_4A_inDB      DOUBLE," << endl;
        out_file << "\tburied_area_inner_inQuery      DOUBLE," << endl;
        out_file << "\tburied_area_core_4A_inQuery      DOUBLE," << endl;
        out_file << "\tburied_area_core_3A_inQuery      DOUBLE," << endl;
        out_file << "\tburied_area_core_2A_inQuery      DOUBLE," << endl;
        out_file << "\tburied_area_core_1A_inQuery      DOUBLE," << endl;
        out_file << "\tburied_area_shell_1A_inQuery      DOUBLE," << endl;
        out_file << "\tburied_area_shell_2A_inQuery      DOUBLE," << endl;
        out_file << "\tburied_area_shell_3A_inQuery      DOUBLE," << endl;
        out_file << "\tburied_area_shell_4A_inQuery      DOUBLE," << endl;


        out_file << "\tT0_0        DOUBLE," << endl;
        out_file << "\tT0_1        DOUBLE," << endl;
        out_file << "\tT0_2        DOUBLE," << endl;
        out_file << "\tR_00    DOUBLE," << endl;
        out_file << "\tR_01    DOUBLE," << endl;
        out_file << "\tR_02    DOUBLE," << endl;
        out_file << "\tR_10    DOUBLE," << endl;
        out_file << "\tR_11    DOUBLE," << endl;
        out_file << "\tR_12    DOUBLE," << endl;
        out_file << "\tR_20    DOUBLE," << endl;
        out_file << "\tR_21    DOUBLE," << endl;
        out_file << "\tR_22    DOUBLE," << endl;
        out_file << "\tT1_0        DOUBLE," << endl;
        out_file << "\tT1_1        DOUBLE," << endl;
        out_file << "\tT1_2        DOUBLE," << endl;

        out_file << "\tPRIMARY KEY (ball_id_DB, ball_id_Query, ball_id_sn_Query, pose_id));" << endl;
    }


	void print_SQL_INSERT_row(MyMatchPair* pm_pair, ofstream& out_file)
    {
        size_t i;

        out_file  << "INSERT INTO match_rankings VALUES("
			<< "\t\"" << g_match_pair_str << "\", " << endl
            << "\t" << pm_pair->m_ball_db->m_ball_id << ", " << endl
            << "\t" << pm_pair->m_ball_query->m_ball_id << ", " << endl
            << "\t" << pm_pair->m_ball_query->m_ball_id_serialno << ", " << endl
            << "\t" << pm_pair->m_pose_id << ", " << endl
            << "\t" << pm_pair->m_dist << ", "  << endl
			<< "\t" << pm_pair->m_score << ", "  << endl
            << "\t" << pm_pair->m_rmsd_ca << ", " << endl
            << "\t" << pm_pair->m_rmsd_local << ", " << endl
            << "\t" << pm_pair->m_rmsd_all << ", " << endl
			<< "\t" << pm_pair->m_BOV << ", " << endl;

		for (i = 0; i < pm_pair->m_BOV_v.size(); i ++){
			out_file << "\t" << pm_pair->m_BOV_v[i] << ", " << endl;
		}

		for (i = 0; i < pm_pair->m_buried_area_v.size(); i ++){
            out_file << "\t" << pm_pair->m_buried_area_v[i] << ", " << endl;
        }

        out_file << "\t" << pm_pair->m_translate_query_to_origin[0] << ", " << endl
                 << "\t" << pm_pair->m_translate_query_to_origin[1] << ", " << endl
                 << "\t" << pm_pair->m_translate_query_to_origin[2] << ", " << endl
                 << "\t" << pm_pair->m_rotate_query[0][0] << ", " << endl
                 << "\t" << pm_pair->m_rotate_query[0][1] << ", " << endl
                 << "\t" << pm_pair->m_rotate_query[0][2] << ", " << endl
                 << "\t" << pm_pair->m_rotate_query[1][0] << ", " << endl
                 << "\t" << pm_pair->m_rotate_query[1][1] << ", " << endl
                 << "\t" << pm_pair->m_rotate_query[1][2] << ", " << endl
                 << "\t" << pm_pair->m_rotate_query[2][0] << ", " << endl
                 << "\t" << pm_pair->m_rotate_query[2][1] << ", " << endl
                 << "\t" << pm_pair->m_rotate_query[2][2] << ", " << endl
                 << "\t" << pm_pair->m_translate_query_to_db[0] << ", " << endl
                 << "\t" << pm_pair->m_translate_query_to_db[1] << ", " << endl
                 << "\t" << pm_pair->m_translate_query_to_db[2] << ");" << endl;
    }

    //based on local verts
    double  compute_rmsd_local(MyContextBall* pm_CB_query, 
                            vector<double>& V_translate_query_to_origin, 
                            vector<vector<double> >& M_rotate_query,
                            vector<double>& V_translate_query_to_db)
    {
        vector<double>  tmp_point(3, 0);
        vector<double>  tmp_point_old(3, 0);

        double          tmp_rmsd = 0;
        size_t i;
        int             tmp_num_points = 0;
        for (i = 0; i < pm_CB_query->m_buried_verts.size(); i ++){
            tmp_point_old = pm_CB_query->m_buried_verts[i]->m_xyz;

            tmp_point[0] = tmp_point_old[0] + V_translate_query_to_origin[0];
            tmp_point[1] = tmp_point_old[1] + V_translate_query_to_origin[1];
            tmp_point[2] = tmp_point_old[2] + V_translate_query_to_origin[2];
            tmp_point = rdl_matrix_application_3d(M_rotate_query, tmp_point);  
            tmp_point[0] += V_translate_query_to_db[0];   
            tmp_point[1] += V_translate_query_to_db[1];
            tmp_point[2] += V_translate_query_to_db[2];

            tmp_rmsd += (tmp_point_old[0] - tmp_point[0]) * (tmp_point_old[0] - tmp_point[0]);
            tmp_rmsd += (tmp_point_old[1] - tmp_point[1]) * (tmp_point_old[1] - tmp_point[1]);
            tmp_rmsd += (tmp_point_old[2] - tmp_point[2]) * (tmp_point_old[2] - tmp_point[2]);
            tmp_num_points ++;
        }

        if (tmp_num_points > 0){
            tmp_rmsd = sqrt(tmp_rmsd / tmp_num_points);
        }
        else{
            tmp_rmsd = 9999;
        }

        return tmp_rmsd;
    }

    double  compute_rmsd_ca(MyContextBall* pm_CB_query, 
                            vector<double>& V_translate_query_to_origin, 
                            vector<vector<double> >& M_rotate_query,
                            vector<double>& V_translate_query_to_db)
    {

        vector<double>  tmp_point(3, 0);
        vector<double>  tmp_point_old(3, 0);

        vector<double>  tmp_center_query(3, 0);

        tmp_center_query[0] = pm_CB_query->m_ball_center[0];
        tmp_center_query[1] = pm_CB_query->m_ball_center[1];
        tmp_center_query[2] = pm_CB_query->m_ball_center[2];

        double          tmp_rmsd = 0;
        //double          tmp_dist;
        int             tmp_num_points = 0;
        size_t i;
        for (i = 0; i < m_Query_coarse->m_surface->m_atom_interface_ca.size(); i ++){
            tmp_point_old = m_Query_coarse->m_surface->m_atom_interface_ca[i];
            
            //tmp_dist = rdl_vector_ssd(tmp_center_query, tmp_point_old);
            //if (tmp_dist > g_sample_region_radius) continue;

            tmp_point[0] = tmp_point_old[0] + V_translate_query_to_origin[0];
            tmp_point[1] = tmp_point_old[1] + V_translate_query_to_origin[1];
            tmp_point[2] = tmp_point_old[2] + V_translate_query_to_origin[2];
            tmp_point = rdl_matrix_application_3d(M_rotate_query, tmp_point);  
            tmp_point[0] += V_translate_query_to_db[0];   
            tmp_point[1] += V_translate_query_to_db[1];
            tmp_point[2] += V_translate_query_to_db[2];

            tmp_rmsd += (tmp_point_old[0] - tmp_point[0]) * (tmp_point_old[0] - tmp_point[0]);
            tmp_rmsd += (tmp_point_old[1] - tmp_point[1]) * (tmp_point_old[1] - tmp_point[1]);
            tmp_rmsd += (tmp_point_old[2] - tmp_point[2]) * (tmp_point_old[2] - tmp_point[2]);
            tmp_num_points ++;
        }

        //cout << "tmp_num_points_ca = " << tmp_num_points << endl;
        if (tmp_num_points > 0){
            tmp_rmsd = sqrt(tmp_rmsd / tmp_num_points);
        }
        else{
            tmp_rmsd = 9999;
        }

        return tmp_rmsd;
    }

    double  compute_rmsd_all(MyContextBall* pm_CB_query, 
                            vector<double>& V_translate_query_to_origin, 
                            vector<vector<double> >& M_rotate_query,
                            vector<double>& V_translate_query_to_db)
    {

        

        vector<double>  tmp_point(3, 0);
        vector<double>  tmp_point_old(3, 0);

        vector<double>  tmp_center_query(3, 0);

        tmp_center_query[0] = pm_CB_query->m_ball_center[0];
        tmp_center_query[1] = pm_CB_query->m_ball_center[1];
        tmp_center_query[2] = pm_CB_query->m_ball_center[2];

        double          tmp_rmsd = 0;
        //double          tmp_dist;
        int             tmp_num_points = 0;
        size_t i;
        for (i = 0; i < m_Query_coarse->m_surface->m_atom.size(); i ++){
            tmp_point_old = m_Query_coarse->m_surface->m_atom[i]->point();

            tmp_point[0] = tmp_point_old[0] + V_translate_query_to_origin[0];
            tmp_point[1] = tmp_point_old[1] + V_translate_query_to_origin[1];
            tmp_point[2] = tmp_point_old[2] + V_translate_query_to_origin[2];
            tmp_point = rdl_matrix_application_3d(M_rotate_query, tmp_point);  
            tmp_point[0] += V_translate_query_to_db[0];   
            tmp_point[1] += V_translate_query_to_db[1];
            tmp_point[2] += V_translate_query_to_db[2];

            tmp_rmsd += (tmp_point_old[0] - tmp_point[0]) * (tmp_point_old[0] - tmp_point[0]);
            tmp_rmsd += (tmp_point_old[1] - tmp_point[1]) * (tmp_point_old[1] - tmp_point[1]);
            tmp_rmsd += (tmp_point_old[2] - tmp_point[2]) * (tmp_point_old[2] - tmp_point[2]);
            tmp_num_points ++;
        }

        //cout << "tmp_num_points_ca = " << tmp_num_points << endl;
        if (tmp_num_points > 0){
            tmp_rmsd = sqrt(tmp_rmsd / tmp_num_points);
        }
        else{
            tmp_rmsd = 9999;
        }

        return tmp_rmsd;
    }

    void compute_match_1pair(MyContextBall* pm_db, MyContextBall* pm_query, vector<MyMatchPair*>& pm_out)
    {
		pm_out.clear();

        size_t i;
		double					tmp_score;
		double					tmp_BOV = 0;
		vector<double>			tmp_BOV_v(5, 0);
		vector<double>			tmp_buried_area_v(18, 0);

		if (g_bov_filter_selector == 2){
			for (i = 0; i < m_unitball->m_table.size(); i ++){
				compute_buried_volume_1pair(pm_db, pm_query, i, tmp_BOV_v);
				compute_buried_area(pm_db, pm_query, i, tmp_buried_area_v);
			//rtn_buried_area_v[0] = tmp_buried_area_inner_inDB;
			//rtn_buried_area_v[1] = tmp_buried_area_core_4A_inDB;
			//rtn_buried_area_v[2] = tmp_buried_area_core_3A_inDB;
			//rtn_buried_area_v[3] = tmp_buried_area_core_2A_inDB;
			//rtn_buried_area_v[4] = tmp_buried_area_core_1A_inDB;
			//rtn_buried_area_v[5] = tmp_buried_area_shell_1A_inDB;
			//rtn_buried_area_v[6] = tmp_buried_area_shell_2A_inDB;
			//rtn_buried_area_v[7] = tmp_buried_area_shell_3A_inDB;
			//rtn_buried_area_v[8] = tmp_buried_area_shell_4A_inDB;

			//rtn_buried_area_v[9]  = tmp_buried_area_inner_inQuery;
			//rtn_buried_area_v[10] = tmp_buried_area_core_4A_inQuery;
			//rtn_buried_area_v[11] = tmp_buried_area_core_3A_inQuery;
			//rtn_buried_area_v[12] = tmp_buried_area_core_2A_inQuery;
			//rtn_buried_area_v[13] = tmp_buried_area_core_1A_inQuery;
			//rtn_buried_area_v[14] = tmp_buried_area_shell_1A_inQuery;
			//rtn_buried_area_v[15] = tmp_buried_area_shell_2A_inQuery;
			//rtn_buried_area_v[16] = tmp_buried_area_shell_3A_inQuery;
			//rtn_buried_area_v[17] = tmp_buried_area_shell_4A_inQuery;
				if (tmp_buried_area_v[0] + tmp_buried_area_v[9] > 0) continue;
				if (tmp_buried_area_v[1] + tmp_buried_area_v[10] > 0) continue;

				tmp_score  = 0 - (tmp_buried_area_v[3] + tmp_buried_area_v[12]);
				tmp_score += 4 * (tmp_buried_area_v[4] + tmp_buried_area_v[13]);
				tmp_score += 4 * (tmp_buried_area_v[5] + tmp_buried_area_v[14]);
				tmp_score += (tmp_buried_area_v[6] + tmp_buried_area_v[15]);
				tmp_score += 0.25 * (tmp_buried_area_v[7] + tmp_buried_area_v[16]);
				tmp_score -= 4*(tmp_buried_area_v[2] + tmp_buried_area_v[11]);
				
				MyMatchPair* rtn_pair = new MyMatchPair;
				rtn_pair->m_ball_db = pm_db;
				rtn_pair->m_ball_query = pm_query;
				rtn_pair->m_buried_area_v = tmp_buried_area_v;
				rtn_pair->m_pose_id = i;
				rtn_pair->m_score = tmp_score;
				rtn_pair->m_BOV = tmp_BOV;
				rtn_pair->m_BOV_v = tmp_BOV_v;
				pm_out.push_back(rtn_pair);
			}
		}
		else{ // == 1 or default
			for (i = 0; i < m_unitball->m_table.size(); i ++){
				compute_buried_volume(pm_db, pm_query, i, tmp_BOV);
				compute_buried_area(pm_db, pm_query, i, tmp_buried_area_v);
			//rtn_buried_area_v[0] = tmp_buried_area_inner_inDB;
			//rtn_buried_area_v[1] = tmp_buried_area_core_4A_inDB;
			//rtn_buried_area_v[2] = tmp_buried_area_core_3A_inDB;
			//rtn_buried_area_v[3] = tmp_buried_area_core_2A_inDB;
			//rtn_buried_area_v[4] = tmp_buried_area_core_1A_inDB;
			//rtn_buried_area_v[5] = tmp_buried_area_shell_1A_inDB;
			//rtn_buried_area_v[6] = tmp_buried_area_shell_2A_inDB;
			//rtn_buried_area_v[7] = tmp_buried_area_shell_3A_inDB;
			//rtn_buried_area_v[8] = tmp_buried_area_shell_4A_inDB;

			//rtn_buried_area_v[9]  = tmp_buried_area_inner_inQuery;
			//rtn_buried_area_v[10] = tmp_buried_area_core_4A_inQuery;
			//rtn_buried_area_v[11] = tmp_buried_area_core_3A_inQuery;
			//rtn_buried_area_v[12] = tmp_buried_area_core_2A_inQuery;
			//rtn_buried_area_v[13] = tmp_buried_area_core_1A_inQuery;
			//rtn_buried_area_v[14] = tmp_buried_area_shell_1A_inQuery;
			//rtn_buried_area_v[15] = tmp_buried_area_shell_2A_inQuery;
			//rtn_buried_area_v[16] = tmp_buried_area_shell_3A_inQuery;
			//rtn_buried_area_v[17] = tmp_buried_area_shell_4A_inQuery;

				tmp_score  = (tmp_buried_area_v[3] + tmp_buried_area_v[12]);
				tmp_score += 4 * (tmp_buried_area_v[4] + tmp_buried_area_v[13]);
				tmp_score += 4 * (tmp_buried_area_v[5] + tmp_buried_area_v[14]);
				tmp_score += (tmp_buried_area_v[6] + tmp_buried_area_v[15]);
				tmp_score += 0.25 * (tmp_buried_area_v[7] + tmp_buried_area_v[16]);
				tmp_score -= (tmp_buried_area_v[2] + tmp_buried_area_v[11]);
				tmp_score -= 4* (tmp_buried_area_v[1] + tmp_buried_area_v[10]);

				MyMatchPair* rtn_pair = new MyMatchPair;
				rtn_pair->m_ball_db = pm_db;
				rtn_pair->m_ball_query = pm_query;
				rtn_pair->m_buried_area_v = tmp_buried_area_v;
				rtn_pair->m_pose_id = i;
				rtn_pair->m_score = tmp_score;
				rtn_pair->m_BOV = tmp_BOV;
				rtn_pair->m_BOV_v = tmp_BOV_v;
				pm_out.push_back(rtn_pair);

			}
		}
    }

    MyMatchPair* compute_match_pair(MyContextBall* pm_db, MyContextBall* pm_query)
    {
        size_t i;

		double					tmp_score;
		double					tmp_score2;
		double					tmp_BOV = 0;
		vector<double>			tmp_BOV_v(5, 0);
		vector<double>			tmp_buried_area_v(18, 0);

		double					curr_score = -10;
		unsigned int			curr_pose_id = 0;
		vector<double>			curr_buried_area_v(18, 0);
		double					curr_BOV = 0;
		vector<double>			curr_BOV_v(5, 0);

		if (g_bov_filter_selector == 2){
			for (i = 0; i < m_unitball->m_table.size(); i ++){
				if (compute_buried_volume(pm_db, pm_query, i, tmp_BOV_v) < 0) continue;
				compute_buried_area(pm_db, pm_query, i, tmp_buried_area_v);
			//rtn_buried_area_v[0] = tmp_buried_area_inner_inDB;
			//rtn_buried_area_v[1] = tmp_buried_area_core_4A_inDB;
			//rtn_buried_area_v[2] = tmp_buried_area_core_3A_inDB;
			//rtn_buried_area_v[3] = tmp_buried_area_core_2A_inDB;
			//rtn_buried_area_v[4] = tmp_buried_area_core_1A_inDB;
			//rtn_buried_area_v[5] = tmp_buried_area_shell_1A_inDB;
			//rtn_buried_area_v[6] = tmp_buried_area_shell_2A_inDB;
			//rtn_buried_area_v[7] = tmp_buried_area_shell_3A_inDB;
			//rtn_buried_area_v[8] = tmp_buried_area_shell_4A_inDB;




			//rtn_buried_area_v[9]  = tmp_buried_area_inner_inQuery;
			//rtn_buried_area_v[10] = tmp_buried_area_core_4A_inQuery;
			//rtn_buried_area_v[11] = tmp_buried_area_core_3A_inQuery;
			//rtn_buried_area_v[12] = tmp_buried_area_core_2A_inQuery;
			//rtn_buried_area_v[13] = tmp_buried_area_core_1A_inQuery;
			//rtn_buried_area_v[14] = tmp_buried_area_shell_1A_inQuery;
			//rtn_buried_area_v[15] = tmp_buried_area_shell_2A_inQuery;
			//rtn_buried_area_v[16] = tmp_buried_area_shell_3A_inQuery;
			//rtn_buried_area_v[17] = tmp_buried_area_shell_4A_inQuery;

				//tmp_score  = 0 - 1000 * (tmp_buried_area_v[0] + tmp_buried_area_v[9]);
				//tmp_score += 0 - 4 * (tmp_buried_area_v[1] + tmp_buried_area_v[10]);
				//tmp_score += 0 - (tmp_buried_area_v[2] + tmp_buried_area_v[11]);
				////tmp_score += 0.25 * (tmp_buried_area_v[3] + tmp_buried_area_v[12]);
				//tmp_score += 4*(tmp_buried_area_v[4] + tmp_buried_area_v[13]);
				//tmp_score += 4*(tmp_buried_area_v[5] + tmp_buried_area_v[14]);
				//tmp_score += (tmp_buried_area_v[6] + tmp_buried_area_v[15]);
				//tmp_score += 0.5*(tmp_buried_area_v[7] + tmp_buried_area_v[16]);
				////tmp_score += 0.25 * (tmp_buried_area_v[8] + tmp_buried_area_v[17]);

				//tmp_score2  = 0 - 1000 * (tmp_buried_area_v[0] + tmp_buried_area_v[9]);
				//tmp_score2 += 0 - 4 * (tmp_buried_area_v[1] + tmp_buried_area_v[10]);
				//tmp_score2 += 0.25 * (tmp_buried_area_v[2] + tmp_buried_area_v[11]);
				//tmp_score2 += (tmp_buried_area_v[3] + tmp_buried_area_v[12]);
				//tmp_score2 += 4 * (tmp_buried_area_v[4] + tmp_buried_area_v[13]);
				//tmp_score2 += 4 * (tmp_buried_area_v[5] + tmp_buried_area_v[14]);
				//tmp_score2 += (tmp_buried_area_v[6] + tmp_buried_area_v[15]);
				//tmp_score2 += 0.25 * (tmp_buried_area_v[7] + tmp_buried_area_v[16]);
				////tmp_score2 += 0.25 * (tmp_buried_area_v[8] + tmp_buried_area_v[17]);

				//if (tmp_score2 > tmp_score) tmp_score = tmp_score2;

                tmp_score  = 0 - 100*(tmp_buried_area_v[2] + tmp_buried_area_v[11]);
				tmp_score += 4 * (tmp_buried_area_v[4] + tmp_buried_area_v[13]);
				tmp_score += 4 * (tmp_buried_area_v[5] + tmp_buried_area_v[14]);
				tmp_score += (tmp_buried_area_v[6] + tmp_buried_area_v[15]);
				tmp_score += 0.25 * (tmp_buried_area_v[7] + tmp_buried_area_v[16]);


				if (tmp_score > curr_score){
					curr_score = tmp_score;
					curr_pose_id = i;
					curr_buried_area_v = tmp_buried_area_v;
					curr_BOV_v = tmp_BOV_v;
				}
			}
		}
		else{ // == 1 or default
			for (i = 0; i < m_unitball->m_table.size(); i ++){
				if (compute_buried_volume(pm_db, pm_query, i, tmp_BOV) < 0) continue;
				compute_buried_area(pm_db, pm_query, i, tmp_buried_area_v);
			//rtn_buried_area_v[0] = tmp_buried_area_inner_inDB;
			//rtn_buried_area_v[1] = tmp_buried_area_core_4A_inDB;
			//rtn_buried_area_v[2] = tmp_buried_area_core_3A_inDB;
			//rtn_buried_area_v[3] = tmp_buried_area_core_2A_inDB;
			//rtn_buried_area_v[4] = tmp_buried_area_core_1A_inDB;
			//rtn_buried_area_v[5] = tmp_buried_area_shell_1A_inDB;
			//rtn_buried_area_v[6] = tmp_buried_area_shell_2A_inDB;
			//rtn_buried_area_v[7] = tmp_buried_area_shell_3A_inDB;
			//rtn_buried_area_v[8] = tmp_buried_area_shell_4A_inDB;

			//rtn_buried_area_v[9]  = tmp_buried_area_inner_inQuery;
			//rtn_buried_area_v[10] = tmp_buried_area_core_4A_inQuery;
			//rtn_buried_area_v[11] = tmp_buried_area_core_3A_inQuery;
			//rtn_buried_area_v[12] = tmp_buried_area_core_2A_inQuery;
			//rtn_buried_area_v[13] = tmp_buried_area_core_1A_inQuery;
			//rtn_buried_area_v[14] = tmp_buried_area_shell_1A_inQuery;
			//rtn_buried_area_v[15] = tmp_buried_area_shell_2A_inQuery;
			//rtn_buried_area_v[16] = tmp_buried_area_shell_3A_inQuery;
			//rtn_buried_area_v[17] = tmp_buried_area_shell_4A_inQuery;

				//tmp_score  = 4 * (tmp_buried_area_v[4] + tmp_buried_area_v[13]);
				//tmp_score += 4 * (tmp_buried_area_v[5] + tmp_buried_area_v[14]);
				//tmp_score += (tmp_buried_area_v[6] + tmp_buried_area_v[15]);
				//tmp_score += 0.25 * (tmp_buried_area_v[7] + tmp_buried_area_v[16]);
				//tmp_score -= tmp_BOV;


				//tmp_score  = 0 - 1000 * (tmp_buried_area_v[0] + tmp_buried_area_v[9]);
				//tmp_score += 0 - 4 * (tmp_buried_area_v[1] + tmp_buried_area_v[10]);
				//tmp_score += 0 - (tmp_buried_area_v[2] + tmp_buried_area_v[11]);
				//tmp_score += 0.25 * (tmp_buried_area_v[3] + tmp_buried_area_v[12]);
				//tmp_score +=  (tmp_buried_area_v[4] + tmp_buried_area_v[13]);
				//tmp_score += 4 * (tmp_buried_area_v[5] + tmp_buried_area_v[14]);
				//tmp_score += 4 * (tmp_buried_area_v[6] + tmp_buried_area_v[15]);
				//tmp_score += (tmp_buried_area_v[7] + tmp_buried_area_v[16]);
				//tmp_score += 0.25 * (tmp_buried_area_v[8] + tmp_buried_area_v[17]);


                tmp_score  = 0 - 100*(tmp_buried_area_v[2] + tmp_buried_area_v[11]);
				tmp_score += 4 * (tmp_buried_area_v[4] + tmp_buried_area_v[13]);
				tmp_score += 4 * (tmp_buried_area_v[5] + tmp_buried_area_v[14]);
				tmp_score += (tmp_buried_area_v[6] + tmp_buried_area_v[15]);
				tmp_score += 0.25 * (tmp_buried_area_v[7] + tmp_buried_area_v[16]);


				//tmp_score2  = 0 - 1000 * (tmp_buried_area_v[0] + tmp_buried_area_v[9]);
				//tmp_score2 += 0 - 4 * (tmp_buried_area_v[1] + tmp_buried_area_v[10]);
				//tmp_score2 += 0.25 * (tmp_buried_area_v[2] + tmp_buried_area_v[11]);
				//tmp_score2 += (tmp_buried_area_v[3] + tmp_buried_area_v[12]);
				//tmp_score2 += 4 * (tmp_buried_area_v[4] + tmp_buried_area_v[13]);
				//tmp_score2 += 4 * (tmp_buried_area_v[5] + tmp_buried_area_v[14]);
				//tmp_score2 += (tmp_buried_area_v[6] + tmp_buried_area_v[15]);
				//tmp_score2 += 0.25 * (tmp_buried_area_v[7] + tmp_buried_area_v[16]);
				//tmp_score2 += 0.25 * (tmp_buried_area_v[8] + tmp_buried_area_v[17]);

				//if (tmp_score2 > tmp_score) tmp_score = tmp_score2;

				if (tmp_score > curr_score){
					curr_score = tmp_score;
					curr_pose_id = i;
					curr_buried_area_v = tmp_buried_area_v;
					curr_BOV = tmp_BOV;
				}
			}
		}

        
        if (curr_score < 0) return NULL;

		MyMatchPair* rtn_pair = new MyMatchPair;
		rtn_pair->m_ball_db = pm_db;
		rtn_pair->m_ball_query = pm_query;
		rtn_pair->m_buried_area_v = curr_buried_area_v;
		rtn_pair->m_pose_id = curr_pose_id;
		rtn_pair->m_score = curr_score;
		rtn_pair->m_BOV = curr_BOV;
		rtn_pair->m_BOV_v = curr_BOV_v;

        return rtn_pair;
    }
	


    MyMatchPair* compute_match_pair_test(MyContextBall* pm_db, MyContextBall* pm_query)
    {
        size_t i;

		double					tmp_score;
		double					tmp_BOV = 0;
		vector<double>			tmp_BOV_v(5, 0);
		vector<double>			tmp_buried_area_v(18, 0);

		double					curr_score = -10;
		unsigned int			curr_pose_id = 0;
		vector<double>			curr_buried_area_v(18, 0);
		double					curr_BOV = 0;
		vector<double>			curr_BOV_v(5, 0);

		for (i = 0; i < m_unitball->m_table.size(); i ++){
			//compute_buried_volume(pm_db, pm_query, i, tmp_BOV_v);
			compute_buried_volume(pm_db, pm_query, i, tmp_BOV);
		}

		//for (i = 0; i < m_unitball->m_table.size(); i ++){
		//	compute_buried_area_test(pm_db, pm_query, i, tmp_buried_area_v);
		//}

		return NULL;
    }
	

  int compute_buried_volume(MyContextBall* pm_db, MyContextBall* pm_query, unsigned int pm_post_id, double& rtn_buried_volume)
    {
        CONTEXT_RAY_INT_TYPE curr_solid_mask_db, curr_solid_mask_query;
        rtn_buried_volume = 0;
        size_t i, tmp_jump = 100;

        for (i = 0; i < m_unitball->m_table[pm_post_id].size(); i ++){
            curr_solid_mask_query = pm_query->m_context_rays[i];
            
			curr_solid_mask_db = pm_db->m_context_inner_2A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

            if (i % tmp_jump == 0){
                if (rtn_buried_volume > g_bov_inInner2A_CORE_upper_bound){
                    rtn_buried_volume = rtn_buried_volume + 10;
                    return -1;
                }
            }
        }

        return 1;
    }

	int compute_buried_volume_test(MyContextBall* pm_db, MyContextBall* pm_query, unsigned int pm_post_id, double& rtn_buried_volume)
    {
        CONTEXT_RAY_INT_TYPE curr_solid_mask_db, curr_solid_mask_query;
        rtn_buried_volume = 0;
		double rtn_buried_volume_3A = 0;
        size_t i, tmp_jump = 100;

        for (i = 0; i < m_unitball->m_table[pm_post_id].size(); i ++){
            curr_solid_mask_query = pm_query->m_context_rays[i];
            
			curr_solid_mask_db = pm_db->m_context_core_2A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

			curr_solid_mask_db = pm_db->m_context_inner_core3A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume_3A += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

            if (rtn_buried_volume > 50){
                return -1;
            }

            if (rtn_buried_volume_3A > 15){
                return -1;
            }
		
		}

        return 1;
    }

	//every bit counts... for speed!
	int compute_buried_volume(MyContextBall* pm_db, MyContextBall* pm_query, unsigned int pm_post_id, vector<double>& rtn_buried_volume)
    {
        CONTEXT_RAY_INT_TYPE curr_solid_mask_db, curr_solid_mask_query;
        size_t i, j, tmp_jump = 100;

		rtn_buried_volume[0] = 0; //INNER
		rtn_buried_volume[1] = 0; //CORE_4A		
		rtn_buried_volume[2] = 0; //CORE_3A
		rtn_buried_volume[3] = 0; //CORE_2A
		rtn_buried_volume[4] = 0; //CORE_1A

		for (j = 0; j < m_unitball->m_table[pm_post_id].size() / tmp_jump; j ++){
			for (i = j * tmp_jump; i < (j+1) * tmp_jump; i ++){
				curr_solid_mask_query = pm_query->m_context_rays[i];

				//curr_solid_mask_db = pm_db->m_context_inner[m_unitball->m_table[pm_post_id][i]];
				//rtn_buried_volume[0] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

				//curr_solid_mask_db = pm_db->m_context_core_4A[m_unitball->m_table[pm_post_id][i]];
				//rtn_buried_volume[1] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

				//curr_solid_mask_db = pm_db->m_context_core_3A[m_unitball->m_table[pm_post_id][i]];
				//rtn_buried_volume[2] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

				curr_solid_mask_db = pm_db->m_context_core_2A[m_unitball->m_table[pm_post_id][i]];
				rtn_buried_volume[3] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

				curr_solid_mask_db = pm_db->m_context_core_1A[m_unitball->m_table[pm_post_id][i]];
				rtn_buried_volume[4] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);
			}

			//if (rtn_buried_volume[0] > g_bov_inCore_upper_bound){
   //             //rtn_buried_volume[0] = rtn_buried_volume[0] + 10;
   //             return -1;
   //         }
   //         if (rtn_buried_volume[1] > g_bov_inInner4A_upper_bound){
   //             //rtn_buried_volume[1] = rtn_buried_volume[1] + 10;
   //             return -1;
   //         }
   //         if (rtn_buried_volume[2] > g_bov_inInner3A_upper_bound){
   //             //rtn_buried_volume[2] = rtn_buried_volume[2] + 10;
   //             return -1;
   //         }

            if (rtn_buried_volume[3] > g_bov_inInner2A_upper_bound){
                //rtn_buried_volume[3] = rtn_buried_volume[3] + 10;
                return -1;
            }
            if (rtn_buried_volume[4] > g_bov_inInner1A_upper_bound){
                //rtn_buried_volume[4] = rtn_buried_volume[4] + 10;
                return -1;
            }
		}

		if (m_unitball->m_table[pm_post_id].size() / tmp_jump * tmp_jump < m_unitball->m_table[pm_post_id].size()){
			for (i = m_unitball->m_table[pm_post_id].size() / tmp_jump * tmp_jump; i < m_unitball->m_table[pm_post_id].size(); i ++){
				curr_solid_mask_query = pm_query->m_context_rays[i];

				//curr_solid_mask_db = pm_db->m_context_inner[m_unitball->m_table[pm_post_id][i]];
				//rtn_buried_volume[0] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

				//curr_solid_mask_db = pm_db->m_context_core_4A[m_unitball->m_table[pm_post_id][i]];
				//rtn_buried_volume[1] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

				//curr_solid_mask_db = pm_db->m_context_core_3A[m_unitball->m_table[pm_post_id][i]];
				//rtn_buried_volume[2] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

				curr_solid_mask_db = pm_db->m_context_core_2A[m_unitball->m_table[pm_post_id][i]];
				rtn_buried_volume[3] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

				curr_solid_mask_db = pm_db->m_context_core_1A[m_unitball->m_table[pm_post_id][i]];
				rtn_buried_volume[4] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);
			}

			//if (rtn_buried_volume[0] > g_bov_inCore_upper_bound){
   //             //rtn_buried_volume[0] = rtn_buried_volume[0] + 10;
   //             return -1;
   //         }
   //         if (rtn_buried_volume[1] > g_bov_inInner4A_upper_bound){
   //             //rtn_buried_volume[1] = rtn_buried_volume[1] + 10;
   //             return -1;
   //         }
   //         if (rtn_buried_volume[2] > g_bov_inInner3A_upper_bound){
   //             //rtn_buried_volume[2] = rtn_buried_volume[2] + 10;
   //             return -1;
   //         }

            if (rtn_buried_volume[3] > g_bov_inInner2A_upper_bound){
                //rtn_buried_volume[3] = rtn_buried_volume[3] + 10;
                return -1;
            }
            if (rtn_buried_volume[4] > g_bov_inInner1A_upper_bound){
                //rtn_buried_volume[4] = rtn_buried_volume[4] + 10;
                return -1;
            }
		}

        return 1;
    }

	int compute_buried_volume_1pair(MyContextBall* pm_db, MyContextBall* pm_query, unsigned int pm_post_id, vector<double>& rtn_buried_volume)
    {
        CONTEXT_RAY_INT_TYPE curr_solid_mask_db, curr_solid_mask_query;
        size_t i, tmp_jump = 100;

		rtn_buried_volume[0] = 0; //INNER
		rtn_buried_volume[1] = 0; //CORE_4A		
		rtn_buried_volume[2] = 0; //CORE_3A
		rtn_buried_volume[3] = 0; //CORE_2A
		rtn_buried_volume[4] = 0; //CORE_1A

        for (i = 0; i < m_unitball->m_table[pm_post_id].size(); i ++){
            curr_solid_mask_query = pm_query->m_context_rays[i];

            curr_solid_mask_db = pm_db->m_context_inner[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume[0] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

            curr_solid_mask_db = pm_db->m_context_core_4A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume[1] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

            curr_solid_mask_db = pm_db->m_context_core_3A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume[2] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

            curr_solid_mask_db = pm_db->m_context_core_2A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume[3] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

            curr_solid_mask_db = pm_db->m_context_core_1A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume[4] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

   //         if (i % tmp_jump == 0){
   //             //if (rtn_buried_volume[0] > g_bov_inCore_upper_bound){
   //             //    rtn_buried_volume[0] = rtn_buried_volume[0] + 10;
   //             //    return -1;
   //             //}
   //             //if (rtn_buried_volume[1] > g_bov_inInner4A_upper_bound){
   //             //    rtn_buried_volume[1] = rtn_buried_volume[1] + 10;
   //             //    return -1;
   //             //}
   //             if (rtn_buried_volume[2] > g_bov_inInner3A_upper_bound){
   //                 //rtn_buried_volume[2] = rtn_buried_volume[2] + 10;
   //                 return -1;
   //             }
   //             if (rtn_buried_volume[3] > g_bov_inInner2A_upper_bound){
   //                 //rtn_buried_volume[3] = rtn_buried_volume[3] + 10;
   //                 return -1;
   //             }
   //             if (rtn_buried_volume[4] > g_bov_inInner1A_upper_bound){
   //                 //rtn_buried_volume[4] = rtn_buried_volume[4] + 10;
   //                 return -1;
   //             }
			//}
        }

        return 1;
    }

	int compute_buried_volume_test(MyContextBall* pm_db, MyContextBall* pm_query, unsigned int pm_post_id, vector<double>& rtn_buried_volume)
    {
        CONTEXT_RAY_INT_TYPE curr_solid_mask_db, curr_solid_mask_query;
        size_t i, tmp_jump = 100;

		rtn_buried_volume[0] = 0; //INNER
		rtn_buried_volume[1] = 0; //CORE_4A		
		rtn_buried_volume[2] = 0; //CORE_3A
		rtn_buried_volume[3] = 0; //CORE_2A
		rtn_buried_volume[4] = 0; //CORE_1A

        for (i = 0; i < m_unitball->m_table[pm_post_id].size(); i ++){
            curr_solid_mask_query = pm_query->m_context_rays[i];

            curr_solid_mask_db = pm_db->m_context_inner[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume[0] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

            curr_solid_mask_db = pm_db->m_context_core_4A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume[1] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

            curr_solid_mask_db = pm_db->m_context_core_3A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume[2] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

            curr_solid_mask_db = pm_db->m_context_core_2A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume[3] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);

            curr_solid_mask_db = pm_db->m_context_core_1A[m_unitball->m_table[pm_post_id][i]];
            rtn_buried_volume[4] += GET_WEIGHT(curr_solid_mask_query & curr_solid_mask_db);
        }

        return 1;
    }


	void compute_buried_area_test(MyContextBall* pm_db, MyContextBall* pm_query, unsigned int pm_post_id, vector<double>& rtn_buried_area_v)
    {
        CONTEXT_RAY_INT_TYPE curr_solid_mask_db, curr_solid_mask_query;
        size_t i;

        double  tmp_bsa_core1A_shell3A_inQuery = 0;
		double  tmp_bsa_core2A_inQuery = 0;
		double  tmp_bsa_inner_core3A_inQuery = 0;

        double  tmp_bsa_core1A_shell3A_inDB = 0;
		double  tmp_bsa_core2A_inDB = 0;
		double  tmp_bsa_inner_core3A_inDB = 0;

        double  tmp_bsa_core1A_shell3A = 0;
		double  tmp_bsa_core2A = 0;
		double  tmp_bsa_inner_core3A = 0;

		for (i = 0; i < pm_db->m_buried_verts.size(); i ++){
			if (i % 200 == 0){
				if (tmp_bsa_inner_core3A > 25) return;
				if (tmp_bsa_core2A > 150) return;
			}
			curr_solid_mask_db = pm_db->m_buried_verts[i]->m_vert_mask;
            curr_solid_mask_query = pm_query->m_context_inner_core3A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
            if (curr_solid_mask_db & curr_solid_mask_query){
                tmp_bsa_inner_core3A += pm_db->m_buried_verts[i]->m_area;
                continue;
            }

			curr_solid_mask_query = pm_query->m_context_core_2A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
            if (curr_solid_mask_db & curr_solid_mask_query){
                tmp_bsa_core2A += pm_db->m_buried_verts[i]->m_area;
            }

		}

        for (i = 0; i < pm_query->m_buried_verts.size(); i ++){
			if (i % 200 == 0){
				if (tmp_bsa_inner_core3A > 25) return;
				if (tmp_bsa_core2A > 150) return;
			}

            curr_solid_mask_query = pm_query->m_buried_verts[i]->m_vert_mask;
            curr_solid_mask_db = pm_db->m_context_inner_core3A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
            if (curr_solid_mask_query & curr_solid_mask_db){
                tmp_bsa_inner_core3A += pm_query->m_buried_verts[i]->m_area;
                continue;
            }
            curr_solid_mask_db = pm_db->m_context_core_2A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
            if (curr_solid_mask_query & curr_solid_mask_db){
                tmp_bsa_core2A += pm_query->m_buried_verts[i]->m_area;
            }
		}

		for (i = 0; i < pm_db->m_buried_verts.size(); i ++){
            curr_solid_mask_db = pm_db->m_buried_verts[i]->m_vert_mask;

            curr_solid_mask_query = pm_query->m_context_core1A_shell3A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
            if (curr_solid_mask_db & curr_solid_mask_query){
                tmp_bsa_core1A_shell3A += pm_db->m_buried_verts[i]->m_area;
            }
		}

  //      for (i = 0; i < pm_query->m_buried_verts.size(); i ++){
  //          curr_solid_mask_query = pm_query->m_buried_verts[i]->m_vert_mask;
  //          
  //          curr_solid_mask_db = pm_db->m_context_core1A_shell3A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
  //          if (curr_solid_mask_query & curr_solid_mask_db){
  //              tmp_bsa_core1A_shell3A += pm_query->m_buried_verts[i]->m_area;
  //          }
		//}


        rtn_buried_area_v[0] = tmp_bsa_inner_core3A_inDB;
        rtn_buried_area_v[1] = tmp_bsa_core2A_inDB;
        rtn_buried_area_v[2] = tmp_bsa_core1A_shell3A_inDB;
        rtn_buried_area_v[3] = tmp_bsa_inner_core3A_inQuery;
        rtn_buried_area_v[4] = tmp_bsa_core2A_inQuery;
        rtn_buried_area_v[5] = tmp_bsa_core1A_shell3A_inQuery;

		return;
    }


	//void compute_buried_area(MyContextBall* pm_db, MyContextBall* pm_query, unsigned int pm_post_id, vector<double>& rtn_buried_area_v)
 //   {
 //       CONTEXT_RAY_INT_TYPE curr_solid_mask_db, curr_solid_mask_query;
 //       size_t i;

 //       double  tmp_buried_area_inner_inDB = 0;
 //       double  tmp_buried_area_core_4A_inDB = 0;
 //       double  tmp_buried_area_core_3A_inDB = 0;
 //       double  tmp_buried_area_core_2A_inDB = 0;
 //       double  tmp_buried_area_core_1A_inDB = 0;
 //       double  tmp_buried_area_shell_1A_inDB = 0;
 //       double  tmp_buried_area_shell_2A_inDB = 0;
 //       double  tmp_buried_area_shell_3A_inDB = 0;
 //       double  tmp_buried_area_shell_4A_inDB = 0;

 //       double  tmp_buried_area_inner_inQuery = 0;
 //       double  tmp_buried_area_core_4A_inQuery = 0;
 //       double  tmp_buried_area_core_3A_inQuery = 0;
 //       double  tmp_buried_area_core_2A_inQuery = 0;
 //       double  tmp_buried_area_core_1A_inQuery = 0;
 //       double  tmp_buried_area_shell_1A_inQuery = 0;
 //       double  tmp_buried_area_shell_2A_inQuery = 0;
 //       double  tmp_buried_area_shell_3A_inQuery = 0;
 //       double  tmp_buried_area_shell_4A_inQuery = 0;

 //       for (i = 0; i < pm_db->m_buried_verts.size(); i ++){
 //           curr_solid_mask_db = pm_db->m_buried_verts[i]->m_vert_mask;
 //           curr_solid_mask_query = pm_query->m_context_shell_L2[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //           if (curr_solid_mask_db & curr_solid_mask_query){
 //               curr_solid_mask_query = pm_query->m_context_shell_4A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //               if (curr_solid_mask_db & curr_solid_mask_query){
 //                   tmp_buried_area_shell_4A_inQuery += pm_db->m_buried_verts[i]->m_area;
 //                   continue;
 //               }

 //               curr_solid_mask_query = pm_query->m_context_shell_3A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //               if (curr_solid_mask_db & curr_solid_mask_query){
 //                   tmp_buried_area_shell_3A_inQuery += pm_db->m_buried_verts[i]->m_area;
 //                   continue;
 //               }
	//			
 //               curr_solid_mask_query = pm_query->m_context_shell_2A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //               if (curr_solid_mask_db & curr_solid_mask_query){
 //                   tmp_buried_area_shell_2A_inQuery += pm_db->m_buried_verts[i]->m_area;
 //                   continue;
 //               }

	//			curr_solid_mask_query = pm_query->m_context_shell_1A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //               if (curr_solid_mask_db & curr_solid_mask_query){
 //                   tmp_buried_area_shell_1A_inQuery += pm_db->m_buried_verts[i]->m_area;
 //                   continue;
 //               }



 //           }
 //           else{
 //               curr_solid_mask_query = pm_query->m_context_core_L2[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //               if (curr_solid_mask_db & curr_solid_mask_query){
 //                   curr_solid_mask_query = pm_query->m_context_core_1A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //                   if (curr_solid_mask_db & curr_solid_mask_query){
 //                       tmp_buried_area_core_1A_inQuery += pm_db->m_buried_verts[i]->m_area;
 //                       continue;
 //                   }

 //                   curr_solid_mask_query = pm_query->m_context_core_2A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //                   if (curr_solid_mask_db & curr_solid_mask_query){
 //                       tmp_buried_area_core_2A_inQuery += pm_db->m_buried_verts[i]->m_area;
 //                       continue;
 //                   }

 //                   curr_solid_mask_query = pm_query->m_context_core_3A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //                   if (curr_solid_mask_db & curr_solid_mask_query){
 //                       tmp_buried_area_core_3A_inQuery += pm_db->m_buried_verts[i]->m_area;
 //                       continue;
 //                   }

 //                   curr_solid_mask_query = pm_query->m_context_core_4A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //                   if (curr_solid_mask_db & curr_solid_mask_query){
 //                       tmp_buried_area_core_4A_inQuery += pm_db->m_buried_verts[i]->m_area;
 //                       continue;
 //                   }

 //                   curr_solid_mask_query = pm_query->m_context_inner[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
 //                   if (curr_solid_mask_db & curr_solid_mask_query){
 //                       tmp_buried_area_inner_inQuery += pm_db->m_buried_verts[i]->m_area;
 //                       continue;
 //                   }
 //               }
 //           }
 //       }

 //       for (i = 0; i < pm_query->m_buried_verts.size(); i ++){
 //           curr_solid_mask_query = pm_query->m_buried_verts[i]->m_vert_mask;
 //           curr_solid_mask_db = pm_db->m_context_shell_L2[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //           if (curr_solid_mask_query & curr_solid_mask_db){
 //               curr_solid_mask_db = pm_db->m_context_shell_4A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //               if (curr_solid_mask_query & curr_solid_mask_db){
 //                   tmp_buried_area_shell_4A_inDB += pm_query->m_buried_verts[i]->m_area;
 //                   continue;
 //               }

 //               curr_solid_mask_db = pm_db->m_context_shell_3A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //               if (curr_solid_mask_query & curr_solid_mask_db){
 //                   tmp_buried_area_shell_3A_inDB += pm_query->m_buried_verts[i]->m_area;
 //                   continue;
 //               }

 //               curr_solid_mask_db = pm_db->m_context_shell_2A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //               if (curr_solid_mask_query & curr_solid_mask_db){
 //                   tmp_buried_area_shell_2A_inDB += pm_query->m_buried_verts[i]->m_area;
 //                   continue;
 //               }

	//			curr_solid_mask_db = pm_db->m_context_shell_1A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //               if (curr_solid_mask_query & curr_solid_mask_db){
 //                   tmp_buried_area_shell_1A_inDB += pm_query->m_buried_verts[i]->m_area;
 //                   continue;
 //               }
 //           }
 //           else{
 //               curr_solid_mask_db = pm_db->m_context_core_L2[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //               if (curr_solid_mask_query & curr_solid_mask_db){
 //                   curr_solid_mask_db = pm_db->m_context_core_1A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //                   if (curr_solid_mask_query & curr_solid_mask_db){
 //                       tmp_buried_area_core_1A_inDB += pm_query->m_buried_verts[i]->m_area;
 //                       continue;
 //                   }

 //                   curr_solid_mask_db = pm_db->m_context_core_2A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //                   if (curr_solid_mask_query & curr_solid_mask_db){
 //                       tmp_buried_area_core_2A_inDB += pm_query->m_buried_verts[i]->m_area;
 //                       continue;
 //                   }

 //                   curr_solid_mask_db = pm_db->m_context_core_3A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //                   if (curr_solid_mask_query & curr_solid_mask_db){
 //                       tmp_buried_area_core_3A_inDB += pm_query->m_buried_verts[i]->m_area;
 //                       continue;
 //                   }

 //                   curr_solid_mask_db = pm_db->m_context_core_4A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //                   if (curr_solid_mask_query & curr_solid_mask_db){
 //                       tmp_buried_area_core_4A_inDB += pm_query->m_buried_verts[i]->m_area;
 //                       continue;
 //                   }

 //                   curr_solid_mask_db = pm_db->m_context_inner[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
 //                   if (curr_solid_mask_query & curr_solid_mask_db){
 //                       tmp_buried_area_inner_inDB += pm_query->m_buried_verts[i]->m_area;
 //                       continue;
 //                   }

 //               }
 //           }
 //       }

 //       rtn_buried_area_v[0] = tmp_buried_area_inner_inDB;
 //       rtn_buried_area_v[1] = tmp_buried_area_core_4A_inDB;
 //       rtn_buried_area_v[2] = tmp_buried_area_core_3A_inDB;
 //       rtn_buried_area_v[3] = tmp_buried_area_core_2A_inDB;
 //       rtn_buried_area_v[4] = tmp_buried_area_core_1A_inDB;
 //       rtn_buried_area_v[5] = tmp_buried_area_shell_1A_inDB;
 //       rtn_buried_area_v[6] = tmp_buried_area_shell_2A_inDB;
 //       rtn_buried_area_v[7] = tmp_buried_area_shell_3A_inDB;
 //       rtn_buried_area_v[8] = tmp_buried_area_shell_4A_inDB;

 //       rtn_buried_area_v[9]  = tmp_buried_area_inner_inQuery;
 //       rtn_buried_area_v[10] = tmp_buried_area_core_4A_inQuery;
 //       rtn_buried_area_v[11] = tmp_buried_area_core_3A_inQuery;
 //       rtn_buried_area_v[12] = tmp_buried_area_core_2A_inQuery;
 //       rtn_buried_area_v[13] = tmp_buried_area_core_1A_inQuery;
 //       rtn_buried_area_v[14] = tmp_buried_area_shell_1A_inQuery;
 //       rtn_buried_area_v[15] = tmp_buried_area_shell_2A_inQuery;
 //       rtn_buried_area_v[16] = tmp_buried_area_shell_3A_inQuery;
 //       rtn_buried_area_v[17] = tmp_buried_area_shell_4A_inQuery;

 //       return;
 //   }



	void compute_buried_area(MyContextBall* pm_db, MyContextBall* pm_query, unsigned int pm_post_id, vector<double>& rtn_buried_area_v)
    {
        CONTEXT_RAY_INT_TYPE curr_solid_mask_db, curr_solid_mask_query;
        size_t i;

        //double  tmp_buried_area_inner_inDB = 0;
        //double  tmp_buried_area_core_4A_inDB = 0;
        double  tmp_buried_area_core_3A_inDB = 0;
        //double  tmp_buried_area_core_2A_inDB = 0;

        double  tmp_buried_area_core_1A_inDB = 0;
        double  tmp_buried_area_shell_1A_inDB = 0;
        double  tmp_buried_area_shell_2A_inDB = 0;
        double  tmp_buried_area_shell_3A_inDB = 0;
        
        //double  tmp_buried_area_shell_4A_inDB = 0;

        //double  tmp_buried_area_inner_inQuery = 0;
        //double  tmp_buried_area_core_4A_inQuery = 0;
        double  tmp_buried_area_core_3A_inQuery = 0;
        //double  tmp_buried_area_core_2A_inQuery = 0;
        
        double  tmp_buried_area_core_1A_inQuery = 0;
        double  tmp_buried_area_shell_1A_inQuery = 0;
        double  tmp_buried_area_shell_2A_inQuery = 0;
        double  tmp_buried_area_shell_3A_inQuery = 0;
        
        //double  tmp_buried_area_shell_4A_inQuery = 0;

        for (i = 0; i < pm_db->m_buried_verts.size(); i ++){
            curr_solid_mask_db = pm_db->m_buried_verts[i]->m_vert_mask;
            curr_solid_mask_query = pm_query->m_context_core1A_shell3A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];

            //cout << "pm_query->m_context_core1Ashell3A.size() = " << pm_query->m_context_core1Ashell3A.size() << endl;

            if (curr_solid_mask_db & curr_solid_mask_query){
                curr_solid_mask_query = pm_query->m_context_core1Ashell3A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
                if (curr_solid_mask_db & curr_solid_mask_query){
                    curr_solid_mask_query = pm_query->m_context_shell_3A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
                    if (curr_solid_mask_db & curr_solid_mask_query){
                        tmp_buried_area_shell_3A_inQuery += pm_db->m_buried_verts[i]->m_area;
                        //m_hit_shell_3A ++;
                        continue;
                    }

                    tmp_buried_area_core_1A_inQuery += pm_db->m_buried_verts[i]->m_area;
                    //m_hit_core_1A ++;

                }
                else{

                    curr_solid_mask_query = pm_query->m_context_shell_2A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
                    if (curr_solid_mask_db & curr_solid_mask_query){
                        tmp_buried_area_shell_2A_inQuery += pm_db->m_buried_verts[i]->m_area;
                        //m_hit_shell_2A ++;
                        continue;
                    }

				    curr_solid_mask_query = pm_query->m_context_shell_1A[m_unitball->m_table_dense2coarse[pm_post_id][pm_db->m_buried_verts[i]->m_ray_id]];
                    if (curr_solid_mask_db & curr_solid_mask_query){
                        tmp_buried_area_shell_1A_inQuery += pm_db->m_buried_verts[i]->m_area;
                        //m_hit_shell_1A ++;
                        continue;
                    }


                    tmp_buried_area_core_3A_inQuery += pm_db->m_buried_verts[i]->m_area;
                    //m_hit_core_3A ++;
                }
            }
        }

        for (i = 0; i < pm_query->m_buried_verts.size(); i ++){
            curr_solid_mask_query = pm_query->m_buried_verts[i]->m_vert_mask;
            curr_solid_mask_db = pm_db->m_context_core1A_shell3A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
            
            //cout << "pm_db->m_context_core1Ashell3A.size() = " << pm_db->m_context_core1Ashell3A.size() << endl;
            if (curr_solid_mask_query & curr_solid_mask_db){
                curr_solid_mask_db = pm_db->m_context_core1Ashell3A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
                if (curr_solid_mask_query & curr_solid_mask_db){
                    curr_solid_mask_db = pm_db->m_context_shell_3A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
                    if (curr_solid_mask_query & curr_solid_mask_db){
                        tmp_buried_area_shell_3A_inDB += pm_query->m_buried_verts[i]->m_area;
                        //m_hit_shell_3A ++;
                        continue;
                    }

                    tmp_buried_area_core_1A_inDB += pm_query->m_buried_verts[i]->m_area;
                    //m_hit_core_1A ++;
                }
                else{
                    curr_solid_mask_db = pm_db->m_context_shell_2A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
                    if (curr_solid_mask_query & curr_solid_mask_db){
                        tmp_buried_area_shell_2A_inDB += pm_query->m_buried_verts[i]->m_area;
                        //m_hit_shell_2A ++;
                        continue;
                    }

				    curr_solid_mask_db = pm_db->m_context_shell_1A[m_unitball->m_table_coarse2dense[pm_post_id][pm_query->m_buried_verts[i]->m_ray_id]];
                    if (curr_solid_mask_query & curr_solid_mask_db){
                        tmp_buried_area_shell_1A_inDB += pm_query->m_buried_verts[i]->m_area;
                        //m_hit_shell_1A ++;
                        continue;
                    }

                    tmp_buried_area_core_3A_inDB += pm_query->m_buried_verts[i]->m_area;
                    //m_hit_core_3A ++;
                }
            }
        }

        //rtn_buried_area_v[0] = tmp_buried_area_inner_inDB;
        //rtn_buried_area_v[1] = tmp_buried_area_core_4A_inDB;
        rtn_buried_area_v[2] = tmp_buried_area_core_3A_inDB;
        //rtn_buried_area_v[3] = tmp_buried_area_core_2A_inDB;

        rtn_buried_area_v[4] = tmp_buried_area_core_1A_inDB;
        rtn_buried_area_v[5] = tmp_buried_area_shell_1A_inDB;
        rtn_buried_area_v[6] = tmp_buried_area_shell_2A_inDB;
        rtn_buried_area_v[7] = tmp_buried_area_shell_3A_inDB;
        
        //rtn_buried_area_v[8] = tmp_buried_area_shell_4A_inDB;

        //rtn_buried_area_v[9]  = tmp_buried_area_inner_inQuery;
        //rtn_buried_area_v[10] = tmp_buried_area_core_4A_inQuery;
        rtn_buried_area_v[11] = tmp_buried_area_core_3A_inQuery;
        //rtn_buried_area_v[12] = tmp_buried_area_core_2A_inQuery;

        rtn_buried_area_v[13] = tmp_buried_area_core_1A_inQuery;
        rtn_buried_area_v[14] = tmp_buried_area_shell_1A_inQuery;
        rtn_buried_area_v[15] = tmp_buried_area_shell_2A_inQuery;
        rtn_buried_area_v[16] = tmp_buried_area_shell_3A_inQuery;
        
        //rtn_buried_area_v[17] = tmp_buried_area_shell_4A_inQuery;

        return;
    }



};

#endif
