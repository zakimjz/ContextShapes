#include "MyDockinfo.h"
#include "util.h"

using namespace std;

#define GET_WEIGHT precomputed32_weight
extern string g_match_pair_str;

double iterated_weight_16bits (CONTEXT_RAY_INT_TYPE n, unsigned int b, unsigned int e) {
    unsigned int count=0;
    double  tmp_weight = 0;
    double  tmp_r1;
    double  tmp_r2;
    double  tmp_c = 4.0 * M_PI / (3.0 * points_default_coarse);
    for (count = b; count < e; count ++){
        if ((n >> count) & 0x1u){
            tmp_r2 = (count + 1) * g_sample_region_radius / CONTEXT_RAY_INT_SIZE; 
            tmp_r1 = count * g_sample_region_radius / CONTEXT_RAY_INT_SIZE; 
            tmp_weight += tmp_c * (tmp_r2 * tmp_r2 * tmp_r2 - tmp_r1 * tmp_r1 * tmp_r1);
        }
    }

    return tmp_weight;
}

void compute_weight_in_16bits() {
    unsigned int i ;    
    CONTEXT_RAY_INT_TYPE n;

    if (CONTEXT_RAY_INT_SIZE == 32){
        for (i = 0; i < (0x1u<<16); i++){
            n = i;
            weight_in_1_16bits[i] = iterated_weight_16bits (n, 0, 16) ;
            weight_in_2_16bits[i] = iterated_weight_16bits (n << 16, 16, 32) ;
        }
    }
    else{
        cerr << "ERROR: wrong CONTEXT_RAY_INT_SIZE in void compute_weight_in_16bits()!" << endl;
    }

    return ;

}

inline double precomputed32_weight (CONTEXT_RAY_INT_TYPE n) {
    return weight_in_1_16bits[n         & 0x0000ffffu]
    +  weight_in_2_16bits[(n >> 16) & 0x0000ffffu] ;
}

void test_weight_table_32() {
    size_t i;
    double w1, w2;
    for (i = 0; i <= 0xffffffffu; i ++){
        if (i/10000*10000 == i) cout << "i = " << i << endl;

        w1 = iterated_weight_16bits(i, 0, 32);
        w2 = precomputed32_weight(i);
        if (fabs(w1 - w2) > SMALL_NUM){
            cerr << "ERROR: inconsistency in weight table! i = " << i << "; w1 = " << w1 << " w2 = " << w2 << endl;
        }
    }
}

double
MyDockinfo :: compute_rmsd_local(MyContextBall* pm_CB_ligand, 
                            vector<double>& V_translate_ligand_to_origin, 
                            vector<vector<double> >& M_rotate_ligand,
                            vector<double>& V_translate_ligand_to_receptor)
{
        vector<double>  tmp_point(3, 0);
        vector<double>  tmp_point_old(3, 0);

        double          tmp_rmsd = 0;
        size_t i;
        int             tmp_num_points = 0;
        for (i = 0; i < pm_CB_ligand->m_buried_verts.size(); i ++){
            tmp_point_old = pm_CB_ligand->m_buried_verts[i]->m_xyz;

            tmp_point[0] = tmp_point_old[0] + V_translate_ligand_to_origin[0];
            tmp_point[1] = tmp_point_old[1] + V_translate_ligand_to_origin[1];
            tmp_point[2] = tmp_point_old[2] + V_translate_ligand_to_origin[2];
            tmp_point = rdl_matrix_application_3d(M_rotate_ligand, tmp_point);  
            tmp_point[0] += V_translate_ligand_to_receptor[0];   
            tmp_point[1] += V_translate_ligand_to_receptor[1];
            tmp_point[2] += V_translate_ligand_to_receptor[2];

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
void 
MyDockinfo :: compute_BSA_Grid(MyGridSES_SPD* pm_GridSES_Cells, MyContextBall* pm_CB_Cells, 
		     MyGridSES_SPD* pm_GridSES_Grid, MyContextBall* pm_CB_Grid, 
			 vector<double>& pm_BSA_inGrid)
{
		pm_BSA_inGrid[0] = 0; //CORE	INNER
		pm_BSA_inGrid[1] = 0; //InL3	CORE_4A
		pm_BSA_inGrid[2] = 0; //InL2	CORE_3A
		pm_BSA_inGrid[3] = 0; //InL1	CORE_2A
		pm_BSA_inGrid[4] = 0; //OutL1	CORE_1A

		pm_BSA_inGrid[5] = 0; //OutL2	SHELL_1A
		pm_BSA_inGrid[6] = 0; //OutL3	SHELL_2A
		pm_BSA_inGrid[7] = 0; //SIS		SHELL_3A
		pm_BSA_inGrid[8] = 0; //SIS		SHELL_4A

		unsigned int			tmp_cell_type;
		vector<double>			tmp_point(3, 0);
		size_t k;

		for (k = 0; k < pm_CB_Cells->m_buried_verts.size(); k ++){
			tmp_point = pm_CB_Cells->m_buried_verts[k]->m_xyz;
			tmp_cell_type = pm_GridSES_Grid->check_cell_type(tmp_point);
			
			if ((tmp_cell_type & RDL_MASK) == RDL_INNER){
				pm_BSA_inGrid[0] += pm_CB_Cells->m_buried_verts[k]->m_area;	
			}
			else if ((tmp_cell_type & RDL_MASK) == RDL_CORE_L2){
				double tmp_mc_dist = sqrt(pm_GridSES_Grid->mc_get_x(tmp_cell_type)*pm_GridSES_Grid->mc_get_x(tmp_cell_type)*1.0 + 
					pm_GridSES_Grid->mc_get_y(tmp_cell_type)*pm_GridSES_Grid->mc_get_y(tmp_cell_type)*1.0 + 
					pm_GridSES_Grid->mc_get_z(tmp_cell_type)*pm_GridSES_Grid->mc_get_z(tmp_cell_type)*1.0) * g_grid_cell_len;
				if (tmp_mc_dist <= 1){
					pm_BSA_inGrid[4] += pm_CB_Cells->m_buried_verts[k]->m_area;	
				}
				else if (tmp_mc_dist <= 2){
					pm_BSA_inGrid[3] += pm_CB_Cells->m_buried_verts[k]->m_area;	
				}
				else if (tmp_mc_dist <= 3){
					pm_BSA_inGrid[2] += pm_CB_Cells->m_buried_verts[k]->m_area;	
					//pm_BSA_inGrid[8] += pm_CB_Cells->m_buried_verts[k]->m_area;	
				}
				else{
					pm_BSA_inGrid[1] += pm_CB_Cells->m_buried_verts[k]->m_area;	
				}
			}
			else if ((tmp_cell_type & RDL_MASK) == RDL_SHELL_L2){
				double tmp_mc_dist = sqrt(pm_GridSES_Grid->mc_get_x(tmp_cell_type)*pm_GridSES_Grid->mc_get_x(tmp_cell_type)*1.0 + 
					pm_GridSES_Grid->mc_get_y(tmp_cell_type)*pm_GridSES_Grid->mc_get_y(tmp_cell_type)*1.0 + 
					pm_GridSES_Grid->mc_get_z(tmp_cell_type)*pm_GridSES_Grid->mc_get_z(tmp_cell_type)*1.0) * g_grid_cell_len;
				if (tmp_mc_dist < 1){
					pm_BSA_inGrid[5] += pm_CB_Cells->m_buried_verts[k]->m_area;	
				}
				else if (tmp_mc_dist < 2){
					pm_BSA_inGrid[6] += pm_CB_Cells->m_buried_verts[k]->m_area;	
				}
				else if (tmp_mc_dist < 3){
					pm_BSA_inGrid[7] += pm_CB_Cells->m_buried_verts[k]->m_area;	
				}
				else{
					pm_BSA_inGrid[8] += pm_CB_Cells->m_buried_verts[k]->m_area;	
				}
			}
			else if((tmp_cell_type & RDL_MASK) == RDL_SES){
				pm_BSA_inGrid[5] += pm_CB_Cells->m_buried_verts[k]->m_area;	
			}
		}
}

void 
MyDockinfo ::	compute_BSA_CB(MyContextBall* pm_CB_receptor, 
					    MyContextBall* pm_CB_ligand, 
						vector<double>& pm_BSA_CB)
	{
		pm_BSA_CB[0] = 0; //inR_INNER
		pm_BSA_CB[1] = 0; //inR_CORE_4A
		pm_BSA_CB[2] = 0; //inR_CORE_3A
		pm_BSA_CB[3] = 0; //inR_CORE_2A
		pm_BSA_CB[4] = 0; //inR_CORE_1A
		pm_BSA_CB[5] = 0; //inR_SHELL_1A
		pm_BSA_CB[6] = 0; //inR_SHELL_2A
		pm_BSA_CB[7] = 0; //inR_SHELL_3A
		pm_BSA_CB[8] = 0; //inR_SHELL_4A

		pm_BSA_CB[9] = 0;  //inL_INNER
		pm_BSA_CB[10] = 0; //inL_CORE_4A
		pm_BSA_CB[11] = 0; //inL_CORE_3A
		pm_BSA_CB[12] = 0; //inL_CORE_2A
		pm_BSA_CB[13] = 0; //inL_CORE_1A
		pm_BSA_CB[14] = 0; //inL_SHELL_1A
		pm_BSA_CB[15] = 0; //inL_SHELL_2A
		pm_BSA_CB[16] = 0; //inL_SHELL_3A
		pm_BSA_CB[17] = 0; //inL_SHELL_4A


		CONTEXT_RAY_INT_TYPE	tmp_ray_ligand, tmp_ray_receptor;
		unsigned int			tmp_ray_id_ligand, tmp_ray_id_receptor;
		double					tmp_area;
		size_t k;

		//inR_ 0_
		for (k = 0; k < pm_CB_ligand->m_buried_verts.size(); k ++){
			tmp_ray_ligand    = pm_CB_ligand->m_buried_verts[k]->m_vert_mask;
			tmp_ray_id_ligand = pm_CB_ligand->m_buried_verts[k]->m_ray_id;
			tmp_area          = pm_CB_ligand->m_buried_verts[k]->m_area;
			
			//0_
			tmp_ray_receptor = pm_CB_receptor->m_context_inner[tmp_ray_id_ligand];
			if (tmp_ray_ligand & tmp_ray_receptor){
				pm_BSA_CB[0] += tmp_area;
				continue;
			}

			tmp_ray_receptor = pm_CB_receptor->m_context_core_4A[tmp_ray_id_ligand];
			if (tmp_ray_ligand & tmp_ray_receptor){
				pm_BSA_CB[1] += tmp_area;
				continue;
			}

			tmp_ray_receptor = pm_CB_receptor->m_context_core_3A[tmp_ray_id_ligand];
			if (tmp_ray_ligand & tmp_ray_receptor){
				pm_BSA_CB[2] += tmp_area;
				continue;
			}

			tmp_ray_receptor = pm_CB_receptor->m_context_core_2A[tmp_ray_id_ligand];
			if (tmp_ray_ligand & tmp_ray_receptor){
				pm_BSA_CB[3] += tmp_area;
				continue;
			}

			tmp_ray_receptor = pm_CB_receptor->m_context_core_1A[tmp_ray_id_ligand];
			if (tmp_ray_ligand & tmp_ray_receptor){
				pm_BSA_CB[4] += tmp_area;
				continue;
			}

			tmp_ray_receptor = pm_CB_receptor->m_context_shell_1A[tmp_ray_id_ligand];
			if (tmp_ray_ligand & tmp_ray_receptor){
				pm_BSA_CB[5] += tmp_area;
				continue;
			}

			tmp_ray_receptor = pm_CB_receptor->m_context_shell_2A[tmp_ray_id_ligand];
			if (tmp_ray_ligand & tmp_ray_receptor){
				pm_BSA_CB[6] += tmp_area;
				continue;
			}

			tmp_ray_receptor = pm_CB_receptor->m_context_shell_3A[tmp_ray_id_ligand];
			if (tmp_ray_ligand & tmp_ray_receptor){
				pm_BSA_CB[7] += tmp_area;
				continue;
			}

			tmp_ray_receptor = pm_CB_receptor->m_context_shell_4A[tmp_ray_id_ligand];
			if (tmp_ray_ligand & tmp_ray_receptor){
				pm_BSA_CB[8] += tmp_area;
				continue;
			}

		}

		//inL_ 0_
		for (k = 0; k < pm_CB_receptor->m_buried_verts.size(); k ++){
			tmp_ray_receptor	= pm_CB_receptor->m_buried_verts[k]->m_vert_mask;
			tmp_ray_id_receptor = pm_CB_receptor->m_buried_verts[k]->m_ray_id;
			tmp_area			= pm_CB_receptor->m_buried_verts[k]->m_area;
			
			//0_
			tmp_ray_ligand = pm_CB_ligand->m_context_inner[tmp_ray_id_receptor];
			if (tmp_ray_receptor & tmp_ray_ligand){
				pm_BSA_CB[9] += tmp_area;
			}

			tmp_ray_ligand = pm_CB_ligand->m_context_core_4A[tmp_ray_id_receptor];
			if (tmp_ray_receptor & tmp_ray_ligand){
				pm_BSA_CB[10] += tmp_area;
				continue;
			}

			tmp_ray_ligand = pm_CB_ligand->m_context_core_3A[tmp_ray_id_receptor];
			if (tmp_ray_receptor & tmp_ray_ligand){
				pm_BSA_CB[11] += tmp_area;
				continue;
			}

			tmp_ray_ligand = pm_CB_ligand->m_context_core_2A[tmp_ray_id_receptor];
			if (tmp_ray_receptor & tmp_ray_ligand){
				pm_BSA_CB[12] += tmp_area;
				continue;
			}

			tmp_ray_ligand = pm_CB_ligand->m_context_core_1A[tmp_ray_id_receptor];
			if (tmp_ray_receptor & tmp_ray_ligand){
				pm_BSA_CB[13] += tmp_area;
				continue;
			}

			tmp_ray_ligand = pm_CB_ligand->m_context_shell_1A[tmp_ray_id_receptor];
			if (tmp_ray_receptor & tmp_ray_ligand){
				pm_BSA_CB[14] += tmp_area;
				continue;
			}

			tmp_ray_ligand = pm_CB_ligand->m_context_shell_2A[tmp_ray_id_receptor];
			if (tmp_ray_receptor & tmp_ray_ligand){
				pm_BSA_CB[15] += tmp_area;
				continue;
			}

			tmp_ray_ligand = pm_CB_ligand->m_context_shell_3A[tmp_ray_id_receptor];
			if (tmp_ray_receptor & tmp_ray_ligand){
				pm_BSA_CB[16] += tmp_area;
				continue;
			}

			tmp_ray_ligand = pm_CB_ligand->m_context_shell_4A[tmp_ray_id_receptor];
			if (tmp_ray_receptor & tmp_ray_ligand){
				pm_BSA_CB[17] += tmp_area;
				continue;
			}

		}

}


void 
MyDockinfo :: compute_BOV_CB(MyContextBall* pm_CB_receptor,
						MyContextBall* pm_CB_ligand, 
						vector<double>& pm_BOV_CB)
{
		pm_BOV_CB[0] = 0; //inR_INNER
		pm_BOV_CB[1] = 0; //inR_CORE_4A
		pm_BOV_CB[2] = 0; //inR_CORE_3A
		pm_BOV_CB[3] = 0; //inR_CORE_2A
		pm_BOV_CB[4] = 0; //inL_CORE_1A

		pm_BOV_CB[5] = 0; //inL_INNER
		pm_BOV_CB[6] = 0; //inL_CORE_4A
		pm_BOV_CB[7] = 0; //inL_CORE_3A
		pm_BOV_CB[8] = 0; //inL_CORE_2A
		pm_BOV_CB[9] = 0; //inL_CORE_1A


		size_t i;
		//unsigned int			tmp_ray_id_cells, tmp_ray_id_grid;
		CONTEXT_RAY_INT_TYPE	tmp_ray_cells, tmp_ray_grid;

		//inR_ 0_
		for (i = 0; i < pm_CB_ligand->m_context_rays.size(); i ++){
			tmp_ray_cells = pm_CB_ligand->m_context_rays[i];
			
			tmp_ray_grid = pm_CB_receptor->m_context_inner[i];
			if (tmp_ray_cells & tmp_ray_grid){
				pm_BOV_CB[0] += GET_WEIGHT(tmp_ray_cells & tmp_ray_grid);
			}

			tmp_ray_grid = pm_CB_receptor->m_context_core_4A[i];
			if (tmp_ray_cells & tmp_ray_grid){
				pm_BOV_CB[1] += GET_WEIGHT(tmp_ray_cells & tmp_ray_grid);
			}

			tmp_ray_grid = pm_CB_receptor->m_context_core_3A[i];
			if (tmp_ray_cells & tmp_ray_grid){
				pm_BOV_CB[2] += GET_WEIGHT(tmp_ray_cells & tmp_ray_grid);
			}

			tmp_ray_grid = pm_CB_receptor->m_context_core_2A[i];
			if (tmp_ray_cells & tmp_ray_grid){
				pm_BOV_CB[3] += GET_WEIGHT(tmp_ray_cells & tmp_ray_grid);
			}

			tmp_ray_grid = pm_CB_receptor->m_context_core_1A[i];
			if (tmp_ray_cells & tmp_ray_grid){
				pm_BOV_CB[4] += GET_WEIGHT(tmp_ray_cells & tmp_ray_grid);
			}

		}


		//inL_  0_
		for (i = 0; i < pm_CB_receptor->m_context_rays.size(); i ++){
			tmp_ray_cells = pm_CB_receptor->m_context_rays[i];
			
			tmp_ray_grid = pm_CB_ligand->m_context_inner[i];
			if (tmp_ray_cells & tmp_ray_grid){
				pm_BOV_CB[5] += GET_WEIGHT(tmp_ray_cells & tmp_ray_grid);
			}

			tmp_ray_grid = pm_CB_ligand->m_context_core_4A[i];
			if (tmp_ray_cells & tmp_ray_grid){
				pm_BOV_CB[6] += GET_WEIGHT(tmp_ray_cells & tmp_ray_grid);
			}

			tmp_ray_grid = pm_CB_ligand->m_context_core_3A[i];
			if (tmp_ray_cells & tmp_ray_grid){
				pm_BOV_CB[7] += GET_WEIGHT(tmp_ray_cells & tmp_ray_grid);
			}

			tmp_ray_grid = pm_CB_ligand->m_context_core_2A[i];
			if (tmp_ray_cells & tmp_ray_grid){
				pm_BOV_CB[8] += GET_WEIGHT(tmp_ray_cells & tmp_ray_grid);
			}

			tmp_ray_grid = pm_CB_ligand->m_context_core_1A[i];
			if (tmp_ray_cells & tmp_ray_grid){
				pm_BOV_CB[9] += GET_WEIGHT(tmp_ray_cells & tmp_ray_grid);
			}

		}
}

void 
MyDockinfo :: compute_BOV_Grid(MyGridSES_SPD* pm_GridSES_Cells, MyContextBall* pm_CB_Cells, 
		     MyGridSES_SPD* pm_GridSES_Grid, MyContextBall* pm_CB_Grid, 
			 vector<double>& pm_BOV_inGrid)
{
		pm_BOV_inGrid[0] = 0; //INNER
		pm_BOV_inGrid[1] = 0; //CORE_4A
		pm_BOV_inGrid[2] = 0; //CORE_3A
		pm_BOV_inGrid[3] = 0; //CORE_2A
		pm_BOV_inGrid[4] = 0; //CORE_1A

		unsigned int			tmp_cell_type;
		vector<double>			tmp_point(3, 0);
		//size_t k;

        int ix_center, iy_center, iz_center;
        vector<double>  tmp_center(3, 0);
        tmp_center[0] = pm_CB_Cells->m_ball_center[0];
        tmp_center[1] = pm_CB_Cells->m_ball_center[1];
        tmp_center[2] = pm_CB_Cells->m_ball_center[2];
        
		ix_center = static_cast<int>((tmp_center[0] - pm_GridSES_Cells->m_min_x) * pm_GridSES_Cells->m_step_multiplier);
        //if (tmp_center[0] - (ix_center * pm_GridSES_Cells->m_step_len + pm_GridSES_Cells->m_min_x) > pm_GridSES_Cells->m_step_len * 0.5) ix_center ++;
        iy_center = static_cast<int>((tmp_center[1] - pm_GridSES_Cells->m_min_y) * pm_GridSES_Cells->m_step_multiplier);
        //if (tmp_center[1] - (iy_center * pm_GridSES_Cells->m_step_len + pm_GridSES_Cells->m_min_y) > pm_GridSES_Cells->m_step_len * 0.5) iy_center ++;
        iz_center = static_cast<int>((tmp_center[2] - pm_GridSES_Cells->m_min_z) * pm_GridSES_Cells->m_step_multiplier);
        //if (tmp_center[2] - (iz_center * pm_GridSES_Cells->m_step_len + pm_GridSES_Cells->m_min_z) > pm_GridSES_Cells->m_step_len * 0.5) iz_center ++;


		int iix, iiy, iiz;
		int	iilen = static_cast<int>(pm_GridSES_Cells->m_step_multiplier * g_sample_region_radius);
		double	ii_dist;
		double  tmp_cell_volume = pm_GridSES_Cells->m_step_len * pm_GridSES_Cells->m_step_len * pm_GridSES_Cells->m_step_len;
		for (iix = (-1) * iilen; iix <= iilen; iix ++){
			if (iix + ix_center < 0 ) continue;
			if (iix + ix_center >= pm_GridSES_Cells->m_size_x) continue;

			for (iiy = (-1) * iilen; iiy <= iilen; iiy ++){
				if (iiy + iy_center < 0 ) continue;
				if (iiy + iy_center >= pm_GridSES_Cells->m_size_y) continue;

				for (iiz = (-1) * iilen; iiz <= iilen; iiz ++){
					if (iiz + iz_center < 0 ) continue;
					if (iiz + iz_center >= pm_GridSES_Cells->m_size_z) continue;
					
					tmp_cell_type = pm_GridSES_Cells->m_grid[iix + ix_center][iiy + iy_center][iiz + iz_center];
					if ((tmp_cell_type & RDL_MASK) == RDL_OUTER
						|| (tmp_cell_type & RDL_MASK) == RDL_SHELL_L2){
						continue;
					}

					ii_dist = sqrt((iix*iix + iiy*iiy + iiz*iiz) * (pm_GridSES_Cells->m_step_len * pm_GridSES_Cells->m_step_len));
					if (ii_dist > g_sample_region_radius) continue;
					tmp_point[0] = iix * pm_GridSES_Cells->m_step_len + pm_CB_Cells->m_ball_center[0];
					tmp_point[1] = iiy * pm_GridSES_Cells->m_step_len + pm_CB_Cells->m_ball_center[1];
					tmp_point[2] = iiz * pm_GridSES_Cells->m_step_len + pm_CB_Cells->m_ball_center[2];
					
					tmp_cell_type = pm_GridSES_Grid->check_cell_type(tmp_point);

					if ((tmp_cell_type & RDL_MASK) == RDL_OUTER){
						continue;
					}
					else if ((tmp_cell_type & RDL_MASK) == RDL_INNER){
						pm_BOV_inGrid[0] += tmp_cell_volume;
					}
					else if((tmp_cell_type & RDL_MASK) == RDL_CORE_L2){
						double tmp_mc_dist = sqrt(pm_GridSES_Grid->mc_get_x(tmp_cell_type)*pm_GridSES_Grid->mc_get_x(tmp_cell_type)*1.0 + 
							pm_GridSES_Grid->mc_get_y(tmp_cell_type)*pm_GridSES_Grid->mc_get_y(tmp_cell_type)*1.0 + 
							pm_GridSES_Grid->mc_get_z(tmp_cell_type)*pm_GridSES_Grid->mc_get_z(tmp_cell_type)*1.0) * g_grid_cell_len;
						if (tmp_mc_dist <= 1){
							pm_BOV_inGrid[4] += tmp_cell_volume;
						}
						else if (tmp_mc_dist <= 2){
							pm_BOV_inGrid[3] += tmp_cell_volume;
						}
						else if (tmp_mc_dist <= 3){
							pm_BOV_inGrid[2] += tmp_cell_volume;
						}
						else{
							pm_BOV_inGrid[1] += tmp_cell_volume;
						}
					}
					else if ((tmp_cell_type & RDL_MASK) == RDL_SES){
						pm_BOV_inGrid[4] += tmp_cell_volume;	
					}
				}
			}
		}
}


void 
MyDockinfo :: match_dockinfo()
{
		size_t i, j, k;

		MyContextBall* cb_receptor;
		MyContextBall* cb_ligand;
		vector<MyContextBall*> pool_cb_receptors;
		vector<MyContextBall*> pool_cb_ligands;

		pool_cb_receptors = m_receptor->m_context_balls;
		pool_cb_ligands = m_ligand->m_context_balls;

		for (i = 0; i < m_ligand->m_context_balls_refine.size(); i ++){
			for (j = 0; j < m_ligand->m_context_balls_refine[i].size(); j ++){
				pool_cb_ligands.push_back(m_ligand->m_context_balls_refine[i][j]);
			}
		}
		

		vector<double>	tmp_solid_angles(4, 0);

		string tmp_fn_dockinfo_sql_INSERT = g_dir_output + g_match_pair_str + "_dockinfo_INSERT.sql";
		string tmp_fn_dockinfo_sql_CREATE_TABLE = g_dir_output + g_match_pair_str + "_dockinfo_CREATE_TABLE.sql";
		
		string tmp_fn_dockinfo_svsa_sql_INSERT = g_dir_output + g_match_pair_str + "_dockinfo_svsa_INSERT.sql";
		string tmp_fn_dockinfo_svsa_sql_CREATE_TABLE = g_dir_output + g_match_pair_str + "_dockinfo_svsa_CREATE_TABLE.sql";


		ofstream out_dockinfo_table(tmp_fn_dockinfo_sql_CREATE_TABLE.c_str());
		if (!out_dockinfo_table.is_open()){
			cerr << "ERROR: can not open file " << tmp_fn_dockinfo_sql_CREATE_TABLE << endl;
			exit(1);
		}
		out_dockinfo_table << "CREATE TABLE dockinfo(" << endl
						 << "\t" << "PDB_PAIR STRING," << endl
						 <<	"\t" << "CB_ID_Receptor INTEGER, " << endl
						 << "\t" << "CB_ID_Ligand INTEGER, " << endl
						 << "\t" << "CB_ID_Ligand_sn INTEGER, " << endl
						 << "\t" << "DIST DOUBLE, " << endl
						 << "\t" << "SVA_0 DOUBLE, " << endl
						 << "\t" << "SVA_1 DOUBLE, " << endl
						 << "\t" << "SVA_2 DOUBLE, " << endl
						 << "\t" << "SVA_3 DOUBLE, " << endl
						 << "\t" << "SVA_4 DOUBLE, " << endl
						 << "\t" << "SVA_5 DOUBLE, " << endl
						 << "\t" << "SVA_6 DOUBLE, " << endl

						 << "\t" << "SA_0 DOUBLE, " << endl
						 << "\t" << "SA_1 DOUBLE, " << endl
						 << "\t" << "SA_2 DOUBLE, " << endl
						 << "\t" << "SA_3 DOUBLE, " << endl
						 << "\t" << "SA_4 DOUBLE, " << endl
						 << "\t" << "SA_5 DOUBLE, " << endl
						 << "\t" << "SA_6 DOUBLE, " << endl

						 << "\t" << "PRJ_DIST DOUBLE, " << endl

						 << "\t" << "BSA_inL_Grid_INNER		DOUBLE, " << endl
						 << "\t" << "BSA_inL_Grid_CORE_4A	DOUBLE, " << endl
						 << "\t" << "BSA_inL_Grid_CORE_3A	DOUBLE, " << endl
						 << "\t" << "BSA_inL_Grid_CORE_2A	DOUBLE, " << endl
						 << "\t" << "BSA_inL_Grid_CORE_1A	DOUBLE, " << endl
						 << "\t" << "BSA_inL_Grid_SHELL_1A	DOUBLE, " << endl
						 << "\t" << "BSA_inL_Grid_SHELL_2A	DOUBLE, " << endl
						 << "\t" << "BSA_inL_Grid_SHELL_3A	DOUBLE, " << endl
						 << "\t" << "BSA_inL_Grid_SHELL_4A	DOUBLE, " << endl

						 << "\t" << "BSA_inR_Grid_INNER		DOUBLE, " << endl
						 << "\t" << "BSA_inR_Grid_CORE_4A	DOUBLE, " << endl
						 << "\t" << "BSA_inR_Grid_CORE_3A	DOUBLE, " << endl
						 << "\t" << "BSA_inR_Grid_CORE_2A	DOUBLE, " << endl
						 << "\t" << "BSA_inR_Grid_CORE_1A	DOUBLE, " << endl
						 << "\t" << "BSA_inR_Grid_SHELL_1A	DOUBLE, " << endl
						 << "\t" << "BSA_inR_Grid_SHELL_2A	DOUBLE, " << endl
						 << "\t" << "BSA_inR_Grid_SHELL_3A	DOUBLE, " << endl
						 << "\t" << "BSA_inR_Grid_SHELL_4A	DOUBLE, " << endl

						 << "\t" << "BOV_inL_Grid_INNER		DOUBLE, " << endl
						 << "\t" << "BOV_inL_Grid_CORE_4A	DOUBLE, " << endl
						 << "\t" << "BOV_inL_Grid_CORE_3A	DOUBLE, " << endl
						 << "\t" << "BOV_inL_Grid_CORE_2A	DOUBLE, " << endl
						 << "\t" << "BOV_inL_Grid_CORE_1A	DOUBLE, " << endl

						 << "\t" << "BOV_inR_Grid_INNER		DOUBLE, " << endl
						 << "\t" << "BOV_inR_Grid_CORE_4A	DOUBLE, " << endl
						 << "\t" << "BOV_inR_Grid_CORE_3A	DOUBLE, " << endl
						 << "\t" << "BOV_inR_Grid_CORE_2A	DOUBLE, " << endl
						 << "\t" << "BOV_inR_Grid_CORE_1A	DOUBLE, " << endl


						 << "\t" << "BSA_inL_CB_INNER		DOUBLE, " << endl
						 << "\t" << "BSA_inL_CB_CORE_4A		DOUBLE, " << endl
						 << "\t" << "BSA_inL_CB_CORE_3A		DOUBLE, " << endl
						 << "\t" << "BSA_inL_CB_CORE_2A		DOUBLE, " << endl
						 << "\t" << "BSA_inL_CB_CORE_1A		DOUBLE, " << endl
						 << "\t" << "BSA_inL_CB_SHELL_1A	DOUBLE, " << endl
						 << "\t" << "BSA_inL_CB_SHELL_2A	DOUBLE, " << endl
						 << "\t" << "BSA_inL_CB_SHELL_3A	DOUBLE, " << endl
						 << "\t" << "BSA_inL_CB_SHELL_4A	DOUBLE, " << endl

						 << "\t" << "BSA_inR_CB_INNER		DOUBLE, " << endl
						 << "\t" << "BSA_inR_CB_CORE_4A		DOUBLE, " << endl
						 << "\t" << "BSA_inR_CB_CORE_3A		DOUBLE, " << endl
						 << "\t" << "BSA_inR_CB_CORE_2A		DOUBLE, " << endl
						 << "\t" << "BSA_inR_CB_CORE_1A		DOUBLE, " << endl
						 << "\t" << "BSA_inR_CB_SHELL_1A	DOUBLE, " << endl
						 << "\t" << "BSA_inR_CB_SHELL_2A	DOUBLE, " << endl
						 << "\t" << "BSA_inR_CB_SHELL_3A	DOUBLE, " << endl
						 << "\t" << "BSA_inR_CB_SHELL_4A	DOUBLE, " << endl

						 << "\t" << "BOV_inL_CB_INNER		DOUBLE, " << endl
						 << "\t" << "BOV_inL_CB_CORE_4A	DOUBLE, " << endl
						 << "\t" << "BOV_inL_CB_CORE_3A	DOUBLE, " << endl
						 << "\t" << "BOV_inL_CB_CORE_2A	DOUBLE, " << endl
						 << "\t" << "BOV_inL_CB_CORE_1A	DOUBLE, " << endl

						 << "\t" << "BOV_inR_CB_INNER		DOUBLE, " << endl
						 << "\t" << "BOV_inR_CB_CORE_4A	DOUBLE, " << endl
						 << "\t" << "BOV_inR_CB_CORE_3A	DOUBLE, " << endl
						 << "\t" << "BOV_inR_CB_CORE_2A	DOUBLE, " << endl
						 << "\t" << "BOV_inR_CB_CORE_1A	DOUBLE, " << endl



						 << "\t" << "R_x DOUBLE, " << endl
						 << "\t" << "R_y DOUBLE, " << endl
						 << "\t" << "R_z DOUBLE, " << endl
						 << "\t" << "L_x DOUBLE, " << endl
						 << "\t" << "L_y DOUBLE, " << endl
						 << "\t" << "L_z DOUBLE, " << endl

						 << "\t" << "PRIMARY KEY(PDB_PAIR, CB_ID_Receptor, CB_ID_Ligand, CB_ID_Ligand_sn));" << endl;
		out_dockinfo_table.close();

		ofstream out_dockinfo_insert(tmp_fn_dockinfo_sql_INSERT.c_str());
		if (!out_dockinfo_insert.is_open()){
			cerr << "ERROR: can not open file " << tmp_fn_dockinfo_sql_INSERT << endl;
			exit(1);
		}
		out_dockinfo_insert << "BEGIN TRANSACTION;" << endl;


		vector<double>	tmp_point_receptor(3, 0);
		vector<double>	tmp_point_ligand(3, 0);
		vector<double>	tmp_point(3, 0);
		vector<double>	tmp_point_new(3, 0);
		vector<double>	tmp_point_diff(3, 0);
		double			tmp_dist;

		vector<double>	tmp_BSA_inR_Grid(9, 0);
		vector<double>	tmp_BSA_inL_Grid(9, 0);
		vector<double>  tmp_BOV_inR_Grid(5, 0);
		vector<double>  tmp_BOV_inL_Grid(5, 0);

		vector<double>	tmp_BSA_CB(18, 0); 
		vector<double>	tmp_BOV_CB(10, 0); 
		
		
		vector<unsigned int>	tmp_local_verts;
		//unsigned int  tmp_cell_type;
		
        vector<vector<double> >			M_Receptor(3, vector<double>(3, 0));
        vector<vector<double> >			M_Ligand(3, vector<double>(3, 0));
        vector<vector<double> >			M_Rotation(3, vector<double>(3, 0));
        vector<vector<double> >			M_Receptor_R(3, vector<double>(3, 0));
        vector<double>                  V_translate_ligand_to_origin(3, 0);    
        vector<vector<double> >         M_rotate_ligand(3, vector<double>(3, 0));
        vector<double>                  V_translate_ligand_to_receptor(3, 0);        

		vector<double>	tmp_sv_ligand(3,0);
		vector<double>	tmp_sv_receptor(3,0);
		double	tmp_SV_angle; //solid vector angle

		vector<double>	tmp_sn_ligand(3, 0);
		vector<double>	tmp_sn_receptor(3, 0);
		//double	tmp_SN_angle; //surface normal angle


		for (i = 0; i < pool_cb_receptors.size(); i ++){
			cb_receptor = pool_cb_receptors[i];
			tmp_point_receptor[0] = cb_receptor->m_ball_center[0];
			tmp_point_receptor[1] = cb_receptor->m_ball_center[1];
			tmp_point_receptor[2] = cb_receptor->m_ball_center[2];
			for (j = 0; j < pool_cb_ligands.size(); j ++){
				cb_ligand = pool_cb_ligands[j];
				tmp_point_ligand[0] = cb_ligand->m_ball_center[0];
				tmp_point_ligand[1] = cb_ligand->m_ball_center[1];
				tmp_point_ligand[2] = cb_ligand->m_ball_center[2];
				tmp_dist = rdl_vector_ssd(tmp_point_receptor, tmp_point_ligand);
				if ( tmp_dist > 3) continue;
				
				//check each surface point of receptor, to see in which layer of ligand it is located
				//cout << "compute_BSA_Grid(m_receptor, cb_receptor, m_ligand, cb_ligand, tmp_BSA_inL_Grid);" << endl << flush;
				compute_BSA_Grid(m_receptor, cb_receptor, m_ligand, cb_ligand, tmp_BSA_inL_Grid);
				
				//check each surface point of ligand, to see in which layer of receptor it is located
				//cout << "compute_BSA_Grid(m_ligand, cb_ligand, m_receptor, cb_receptor, tmp_BSA_inR_Grid);" << endl << flush;
				compute_BSA_Grid(m_ligand, cb_ligand, m_receptor, cb_receptor, tmp_BSA_inR_Grid);
				
				//check each local cell of a receptor CB, to see in which layer of ligand it is located
				//cout << "compute_BOV_Grid(m_receptor, cb_receptor, m_ligand, cb_ligand, tmp_BOV_inL_Grid);" << endl << flush;
				compute_BOV_Grid(m_receptor, cb_receptor, m_ligand, cb_ligand, tmp_BOV_inL_Grid);
				
				//check each local cell of a ligand CB, to see in which layer of receptor it is located
				//cout << "compute_BOV_Grid(m_ligand, cb_ligand, m_receptor, cb_receptor, tmp_BOV_inR_Grid);" << endl << flush;
				compute_BOV_Grid(m_ligand, cb_ligand, m_receptor, cb_receptor, tmp_BOV_inR_Grid);

				//cout << "compute_BSA_CB(cb_receptor, cb_ligand, tmp_BSA_CB);" << endl << flush;
				compute_BSA_CB(cb_receptor, cb_ligand, tmp_BSA_CB);
				
				//cout << "compute_BOV_CB(cb_receptor, cb_ligand, tmp_BOV_CB);" << endl << flush;
				compute_BOV_CB(cb_receptor, cb_ligand, tmp_BOV_CB);

				//best pose in rotational space, use ov?
				//quality to use critical points to calculate BSA

				out_dockinfo_insert << "INSERT INTO dockinfo VALUES(" << endl
					                << "\t\"" << g_match_pair_str << "\", " << endl
								    << "\t" << cb_receptor->m_ball_id << ", " << endl
								    << "\t" << cb_ligand->m_ball_id << ", " << endl
									<< "\t" << cb_ligand->m_ball_id_serialno << ", " << endl
									<< "\t" << tmp_dist << ", " << endl;
				
				for (k = 0; k < g_dynamic_radius.size(); k ++){
					tmp_sv_receptor[0] = cb_receptor->m_solid_vector[k][0];
					tmp_sv_receptor[1] = cb_receptor->m_solid_vector[k][1];
					tmp_sv_receptor[2] = cb_receptor->m_solid_vector[k][2];

					tmp_sv_ligand[0] = cb_ligand->m_solid_vector[k][0];
					tmp_sv_ligand[1] = cb_ligand->m_solid_vector[k][1];
					tmp_sv_ligand[2] = cb_ligand->m_solid_vector[k][2];

					tmp_SV_angle = acos(rdl_dot_3d(tmp_sv_ligand, tmp_sv_receptor)/(norm_3d(tmp_sv_ligand) * norm_3d(tmp_sv_receptor)));
					tmp_SV_angle = tmp_SV_angle * 180 / M_PI;

					out_dockinfo_insert << "\t" << tmp_SV_angle << ", " << endl;
				}

				for (k = 0; k < g_dynamic_radius.size(); k ++){
					out_dockinfo_insert << "\t" << cb_receptor->m_solid_angle[k] + cb_ligand->m_solid_angle[k] << ", " << endl;
				}

				tmp_sn_receptor[0] =  cb_receptor->m_solid_vector[g_solid_volume_radius_selector][0];
				tmp_sn_receptor[1] =  cb_receptor->m_solid_vector[g_solid_volume_radius_selector][1];
				tmp_sn_receptor[2] =  cb_receptor->m_solid_vector[g_solid_volume_radius_selector][2];
				vector<double>	tmp_vec_centers(3, 0);
				tmp_vec_centers[0] = cb_ligand->m_ball_center[0] - cb_receptor->m_ball_center[0];
				tmp_vec_centers[1] = cb_ligand->m_ball_center[1] - cb_receptor->m_ball_center[1];
				tmp_vec_centers[2] = cb_ligand->m_ball_center[2] - cb_receptor->m_ball_center[2];
				if (norm_3d(tmp_vec_centers) == 0){
					cerr << "ERROR: norm_3d(tmp_vec_centers) == 0" << endl;
					exit(1);
				}
				else if (norm_3d(tmp_sn_receptor) == 0){
					cerr << "ERROR: norm_3d(tmp_sn_receptor) == 0" << endl;
					exit(1);
				}
				else{
					double tmp_prj_dist = norm_3d(tmp_vec_centers) * (rdl_dot_3d(tmp_vec_centers, tmp_sn_receptor)/(norm_3d(tmp_vec_centers) * norm_3d(tmp_sn_receptor)));
					out_dockinfo_insert << "\t" << tmp_prj_dist << ", " << endl;
				}

				for (k = 0; k < tmp_BSA_inL_Grid.size(); k ++){
					out_dockinfo_insert << "\t" << tmp_BSA_inL_Grid[k] << ", " << endl;
				}
				for (k = 0; k < tmp_BSA_inR_Grid.size(); k ++){
					out_dockinfo_insert << "\t" << tmp_BSA_inR_Grid[k] << ", " << endl;
				}
				out_dockinfo_insert << "\t" << tmp_BOV_inL_Grid[0]/3 << ", " << endl
									<< "\t" << tmp_BOV_inL_Grid[1]/3 << ", " << endl
									<< "\t" << tmp_BOV_inL_Grid[2]/3 << ", " << endl
									<< "\t" << tmp_BOV_inL_Grid[3]/3 << ", " << endl
									<< "\t" << tmp_BOV_inL_Grid[4]/3 << ", " << endl;

				out_dockinfo_insert << "\t" << tmp_BOV_inR_Grid[0]/3 << ", " << endl
									<< "\t" << tmp_BOV_inR_Grid[1]/3 << ", " << endl
									<< "\t" << tmp_BOV_inR_Grid[2]/3 << ", " << endl
									<< "\t" << tmp_BOV_inR_Grid[3]/3 << ", " << endl
									<< "\t" << tmp_BOV_inR_Grid[4]/3 << ", " << endl;

				for (k = 0; k < tmp_BSA_CB.size(); k ++){
					out_dockinfo_insert << "\t" << tmp_BSA_CB[k] << ", " << endl;
				}

				for (k = 0; k < tmp_BOV_CB.size(); k ++){
					out_dockinfo_insert << "\t" << tmp_BOV_CB[k] << ", " << endl;
				}

				out_dockinfo_insert << "\t" << cb_receptor->m_ball_center[0] << ", " << endl
									<< "\t" << cb_receptor->m_ball_center[1] << ", " << endl
									<< "\t" << cb_receptor->m_ball_center[2] << ", " << endl
									<< "\t" << cb_ligand->m_ball_center[0] << ", " << endl
									<< "\t" << cb_ligand->m_ball_center[1] << ", " << endl
									<< "\t" << cb_ligand->m_ball_center[2] << ");" << endl;


			}
		}
		
		out_dockinfo_insert << "COMMIT;" << endl;
		out_dockinfo_insert.close();


}

	
MyDockinfo :: ~MyDockinfo()
{
		//delete m_ligand->m_unitball;
		//delete m_ligand->m_surface;
		//delete m_receptor->m_surface;
		//delete m_ligand;
		//delete m_receptor;
}

MyDockinfo :: MyDockinfo(MyGridSES_Dense_SPD* A, MyGridSES_Coarse_SPD* B, MyUnitBall* unitball)
{
        m_receptor = A;
        m_ligand = B;
        m_unitball = unitball;
}
