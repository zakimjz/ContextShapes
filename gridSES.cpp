#include "gridSES.h"
#include "memusage.h"

using namespace std;


typedef unsigned int    TID_VERT;
typedef unsigned int    TID_VERT_TYPE;
typedef unsigned int    TID_ATOM;


unsigned int
MyGridSES_SPD :: check_cell_type(vector<double>& pm_point) {

  vector<double> tmp_point = pm_point;
  int ix, iy, iz;

  if (tmp_point[0] - m_min_x < 0){
    return RDL_OUTER;
  }
  ix = static_cast<int>((tmp_point[0] - m_min_x) * m_step_multiplier);
    if (ix > m_size_x - 1){
    return RDL_OUTER;
  }

  if (tmp_point[1] - m_min_y < 0){
    return RDL_OUTER;
  }
  iy = static_cast<int>((tmp_point[1] - m_min_y) * m_step_multiplier);
    if (iy > m_size_y - 1){
    return RDL_OUTER;
  }

  if (tmp_point[2] - m_min_z < 0){
    return RDL_OUTER;
  }
  iz = static_cast<int>((tmp_point[2] - m_min_z) * m_step_multiplier);
    if (iz > m_size_z - 1){
    return RDL_OUTER;
  }

  return m_grid[ix][iy][iz];
}
    
//inline size_t 
size_t 
MyGridSES_SPD :: mc_get_x(unsigned int pm_label) {
  return (pm_label & 0x00FF0000) >> 16;
} 

//inline size_t 
size_t 
MyGridSES_SPD :: mc_get_y(unsigned int pm_label) {
  return (pm_label & 0x0000FF00) >> 8;
}

//inline size_t 
size_t 
MyGridSES_SPD :: mc_get_z(unsigned int pm_label) {
  return pm_label & 0x000000FF;
}

void 
MyGridSES_SPD :: mc_init() {
        size_t  mc_size;
        size_t i, j, k;

        //mc_table_SHELL_L2
        m_mc_localR = 4;
        mc_size = m_mc_size = static_cast<size_t>(m_step_multiplier * 4.0) + 1;
        m_mc_table.clear();
        m_mc_table.insert(m_mc_table.begin(), mc_size, vector<vector<double> >(mc_size, vector<double>(mc_size, 0.0)));
        for (i = 0; i < mc_size; i ++){
            for (j = 0; j < mc_size; j ++){
                for (k = 0; k < mc_size; k ++){
                    m_mc_table[i][j][k] = sqrt((i*i + j*j + k*k) * (m_step_len*m_step_len));
                }
            }
        }

}

void 
MyGridSES_SPD :: sload_context_balls(boost::archive::text_iarchive &la) {
        la & m_step_len;
        la & m_step_multiplier;
        la & m_min_x;
        la & m_min_y;
        la & m_min_z;
        la & m_max_x;
        la & m_max_y;
        la & m_max_z;
        la & m_size_x;
        la & m_size_y;
        la & m_size_z;

        size_t num_cb, num_cb_refine;
        size_t i, j;

        la & num_cb;
        m_context_balls.clear();
        m_context_balls.reserve(num_cb);
        for (i = 0; i < num_cb; i ++){
            MyContextBall* tmp_ball = new MyContextBall;
            tmp_ball->sload(la);
            m_context_balls.push_back(tmp_ball);
        }

        la & num_cb;
        m_context_balls_refine.clear();
        m_context_balls_refine.insert(m_context_balls_refine.end(), num_cb, vector<MyContextBall*>());
        for (i = 0; i < num_cb;  i ++){
            la & num_cb_refine;
            for (j = 0; j < num_cb_refine; j ++){
                MyContextBall* tmp_ball = new MyContextBall;
                tmp_ball->sload(la);
                m_context_balls_refine[i].push_back(tmp_ball);
            }
        }
}

void 
MyGridSES_SPD :: ssave_context_balls(boost::archive::text_oarchive &la) {
        la & m_step_len;
        la & m_step_multiplier;
        la & m_min_x;
        la & m_min_y;
        la & m_min_z;
        la & m_max_x;
        la & m_max_y;
        la & m_max_z;
        la & m_size_x;
        la & m_size_y;
        la & m_size_z;

        size_t num_cb, num_cb_refine;
        size_t i, j;

        num_cb = m_context_balls.size();
        la & num_cb;

        for (i = 0; i < num_cb; i ++){
            m_context_balls[i]->ssave(la);
        }

        num_cb = m_context_balls_refine.size();
        la & num_cb;
        for (i = 0; i < num_cb; i ++){
            num_cb_refine = m_context_balls_refine[i].size();
            la & num_cb_refine;
            for (j = 0; j < num_cb_refine; j ++){
                m_context_balls_refine[i][j]->ssave(la);
            }
        }
}

MyGridSES_SPD :: ~MyGridSES_SPD() {
        size_t i;
        for (i = 0; i < m_context_balls.size(); i ++){
            delete m_context_balls[i];
        }
}

MyGridSES_SPD :: MyGridSES_SPD(MySES* pm_ses, MyUnitBall* pm_ball) {
        m_step_len = g_grid_cell_len;
        m_step_multiplier = 1.0 / g_grid_cell_len;

        m_surface = pm_ses;
        m_unitball = pm_ball;
}
	
void 
MyGridSES_SPD :: write_context_ball_vtk_criticalPointBased(string pm_protein_name) {
        size_t i;
        vector<double>  tmp_center(3, 0);

        for (i = 0; i < m_context_balls.size(); i ++){
            stringstream out_fn;
            out_fn << g_dir_output << pm_protein_name << "_CB_" << m_context_balls[i]->m_ball_id << ".vtk";
            m_context_balls[i]->write_vtk(out_fn.str());
        }
}

void 
MyGridSES_SPD :: write_context_ball_vtk_localverts(string pm_protein_name) {
        size_t i, j;
        vector<double>  tmp_center(3, 0);

        for (i = 0; i < m_context_balls.size(); i ++){
            stringstream out_fn;
            out_fn << g_dir_output << pm_protein_name << "_CBverts_" << m_context_balls[i]->m_ball_id << ".vtk";

            string out_filename = out_fn.str();

            stringstream out_buffer;
            out_buffer << "# vtk DataFile Version 3.0" << endl;
            out_buffer << out_filename << endl;
            out_buffer << "ASCII" << endl;
            out_buffer << "DATASET POLYDATA" << endl;

            vector<vector<double> > tmp_points;
            tmp_points.clear();
            tmp_points.reserve(m_context_balls[i]->m_buried_verts.size());

            for (j = 0; j < m_context_balls[i]->m_buried_verts.size(); j ++){
                tmp_points.push_back(m_surface->m_vert_xyz[m_context_balls[i]->m_buried_verts[j]->m_vert_id]);
            }

            out_buffer << "POINTS " << tmp_points.size() << " float" << endl;
            size_t m;
            for(m = 0; m < tmp_points.size(); m ++){
                out_buffer << tmp_points[m][0] << " "
                    << tmp_points[m][1] << " "
                    << tmp_points[m][2] << endl;
            }

            out_buffer << "VERTICES " << tmp_points.size() << " " << tmp_points.size() * 2 << endl;
            for(m = 0; m < tmp_points.size(); m ++){
                out_buffer << "1 " << m << endl;
            }

            ofstream out_file(out_filename.c_str());
			if (!out_file.is_open()){
				cerr << "ERROR: cannot open file " << out_filename << endl;
				exit(1);
			}

            out_file << out_buffer.str();
            out_file.close();
        }
}

void 
MyGridSES_SPD :: write_grid_vtk() {
        string tmp_fn_ses = g_dir_output + m_surface->m_filename_pdb + "_ses.vtk";
        cout << "writing SES, saving into: " << tmp_fn_ses << endl;
        write_grid_vtk(RDL_SES, tmp_fn_ses);

        //string tmp_fn_outer = g_dir_output + m_surface->m_filename_pdb + "_outer.vtk";
        //cout << "writing outer, saving into: " << tmp_fn_outer << endl;
        //write_grid_vtk(RDL_OUTER, tmp_fn_outer);


        string tmp_fn_inner = g_dir_output + m_surface->m_filename_pdb + "_inner.vtk";
        cout << "writing inner, saving into: " << tmp_fn_inner << endl;
        write_grid_vtk(RDL_INNER, tmp_fn_inner);


        string tmp_fn_core_L2 = g_dir_output + m_surface->m_filename_pdb + "_core_L2.vtk";
        cout << "writing core_L2, saving into: " << tmp_fn_core_L2 << endl;
        write_grid_vtk(RDL_CORE_L2, tmp_fn_core_L2);


        string tmp_fn_shell_L2 = g_dir_output + m_surface->m_filename_pdb + "_shell_L2.vtk";
        cout << "writing shell_L2, saving into: " << tmp_fn_shell_L2 << endl;
        write_grid_vtk(RDL_SHELL_L2, tmp_fn_shell_L2);

}

void 
MyGridSES_SPD :: write_grid_vtk(unsigned int pm_type, string pm_fn) {
        stringstream out_buffer;
        out_buffer << "# vtk DataFile Version 3.0" << endl;
        out_buffer << m_surface->m_filename_pdb.substr(0, m_surface->m_filename_pdb.size() - 4) << endl;
        out_buffer << "ASCII" << endl;
        out_buffer << "DATASET POLYDATA" << endl;
        //out_buffer << "DIMENSIONS " << m_size_x << " " << m_size_y << " " << m_size_z << endl;

        int i, j, k;
        vector<vector<double> > tmp_points;
        vector<double>          tmp_point(3, 0);
        cout << "\tget the points..." << endl;
        for (i = 0; i < m_size_x; i ++){
            for (j = 0; j < m_size_y; j ++){
                for (k = 0; k < m_size_z; k ++){
                    if ((m_grid[i][j][k] & RDL_MASK) == pm_type){
                        tmp_point[0] = i * m_step_len + m_min_x + m_step_len / 2;
                        tmp_point[1] = j * m_step_len + m_min_y + m_step_len / 2;
                        tmp_point[2] = k * m_step_len + m_min_z + m_step_len / 2;
                        tmp_points.push_back(tmp_point);
                    }
                }
            }
        }

        cout << "\toutput the points..." << endl;
        out_buffer << "POINTS " << tmp_points.size() << " float" << endl;
        size_t m;
        for(m = 0; m < tmp_points.size(); m ++){
            out_buffer << tmp_points[m][0] << " "
                << tmp_points[m][1] << " "
                << tmp_points[m][2] << endl;
        }

        cout << "\tsetup the vertices..." << endl;
        out_buffer << "VERTICES " << tmp_points.size() << " " << tmp_points.size() * 2 << endl;
        for(m = 0; m < tmp_points.size(); m ++){
            out_buffer << "1 " << m << endl;
        }

        //string out_fn = m_surface->m_filename_pdb.substr(0, m_surface->m_filename_pdb.size() - 4) + "_" + pm_fn + ".vtk";
        cout << "\tclosing the file..." << endl;
        ofstream out_file(pm_fn.c_str());
        if (!out_file.is_open()){
            cerr << "ERROR: cannot open file " << pm_fn << endl;
            exit(1);
        }

        out_file << out_buffer.str();
        out_file.close();

}


void 
MyGridSES_SPD :: write_critical_points(string pm_fn) {
        //write_critical_points(m_critical_belt, pm_fn + "_belt.vtk");
        write_critical_points(m_critical_concave, pm_fn + "_concave.vtk");
        write_critical_points(m_critical_convex, pm_fn + "_convex.vtk");
}

void 
MyGridSES_SPD :: write_critical_points(vector<vector<double> >& pm_points, string pm_fn) {
        stringstream out_buffer;
        out_buffer << "# vtk DataFile Version 3.0" << endl;
        out_buffer << pm_fn << endl;
        out_buffer << "ASCII" << endl;
        out_buffer << "DATASET POLYDATA" << endl;

        out_buffer << "POINTS " << pm_points.size() << " float" << endl;
        size_t m;
        for(m = 0; m < pm_points.size(); m ++){
            out_buffer << pm_points[m][0] << " "
                << pm_points[m][1] << " "
                << pm_points[m][2] << endl;
        }

        out_buffer << "VERTICES " << pm_points.size() << " " << pm_points.size() * 2 << endl;
        for(m = 0; m < pm_points.size(); m ++){
            out_buffer << "1 " << m << endl;
        }

        ofstream out_file(pm_fn.c_str());
        if (!out_file.is_open()){
            cerr << "ERROR: cannot open file " << pm_fn << endl;
            exit(1);
        }

        out_file << out_buffer.str();
        out_file.close();
}

void 
MyGridSES_SPD :: compute_context_balls_criticalPointBased() {
        size_t i, j;
        size_t num_balls = 0;
        size_t num_balls_total = m_context_balls.size();

        //get the total number of balls to build, for the use of progress bar.
        for (i = 0; i < m_context_balls_refine.size(); i ++){
            num_balls_total += m_context_balls_refine[i].size();
        }

        for (i = 0; i < m_context_balls.size(); i ++){
            if (num_balls % 100 == 0)
                cout << "ball# = " << num_balls << " out of " << num_balls_total << "(" << m_context_balls.size() << ")" << endl;
            compute_CB_solid_volume(m_context_balls[i]);
            compute_CB_solid_angle(m_context_balls[i]);
            compute_CB_rotation_matrix(m_context_balls[i], g_solid_volume_radius_selector); 
            compute_CB_context_rays(m_context_balls[i]);
            compute_CB_context_shell_L2(m_context_balls[i]);
            compute_CB_buried_verts(m_context_balls[i]);
            num_balls ++;
        }

        for (i = 0; i < m_context_balls_refine.size(); i ++){
            for (j = 0; j < m_context_balls_refine[i].size(); j ++){
                if (num_balls % 100 == 0)
                    cout << "ball# = " << num_balls << " out of " << num_balls_total << "(" << m_context_balls.size() << ")" << endl;

                compute_CB_solid_volume(m_context_balls_refine[i][j]);
                compute_CB_solid_angle(m_context_balls_refine[i][j]);
                compute_CB_rotation_matrix(m_context_balls_refine[i][j], g_solid_volume_radius_selector); 
                compute_CB_context_rays(m_context_balls_refine[i][j]);
                compute_CB_context_shell_L2(m_context_balls_refine[i][j]);
                compute_CB_buried_verts(m_context_balls_refine[i][j]);
                num_balls ++;
            }
        }
}


void 
MyGridSES_SPD :: compute_context_balls_criticalPointBased_dockinfo() {
        size_t i, j;
        size_t num_balls = 0;
        size_t num_balls_total = m_context_balls.size();

        //get the total number of balls to build, for the use of progress bar.
        for (i = 0; i < m_context_balls_refine.size(); i ++){
            num_balls_total += m_context_balls_refine[i].size();
        }

        for (i = 0; i < m_context_balls.size(); i ++){
            if (num_balls % 100 == 0)
                cout << "ball# = " << num_balls << " out of " << num_balls_total << "(" << m_context_balls.size() << ")" << endl;
            compute_CB_solid_volume(m_context_balls[i]);
            compute_CB_solid_angle(m_context_balls[i]);
            compute_CB_rotation_matrix(m_context_balls[i], g_solid_volume_radius_selector); 
            compute_CB_context_rays_dockinfo(m_context_balls[i]);
            compute_CB_context_shell_L2_dockinfo(m_context_balls[i]);
            compute_CB_buried_verts_dockinfo(m_context_balls[i]);
            num_balls ++;
        }

        for (i = 0; i < m_context_balls_refine.size(); i ++){
            for (j = 0; j < m_context_balls_refine[i].size(); j ++){
                if (num_balls % 100 == 0)
                    cout << "ball# = " << num_balls << " out of " << num_balls_total << "(" << m_context_balls.size() << ")" << endl;

                compute_CB_solid_volume(m_context_balls_refine[i][j]);
                compute_CB_solid_angle(m_context_balls_refine[i][j]);
                compute_CB_rotation_matrix(m_context_balls_refine[i][j], g_solid_volume_radius_selector); 
                compute_CB_context_rays_dockinfo(m_context_balls_refine[i][j]);
                compute_CB_context_shell_L2_dockinfo(m_context_balls_refine[i][j]);
                compute_CB_buried_verts_dockinfo(m_context_balls_refine[i][j]);
                num_balls ++;
            }
        }
}
void 
MyGridSES_SPD :: compute_critical_points() {
        m_critical_points.clear();
        m_critical_points.insert(m_critical_points.end(), m_critical_concave.begin(), m_critical_concave.end());
        m_critical_points.insert(m_critical_points.end(), m_critical_convex.begin(), m_critical_convex.end());

        m_critical_points_area.clear();
        m_critical_points_area.insert(m_critical_points_area.end(), m_critical_concave_area.begin(), m_critical_concave_area.end());
        m_critical_points_area.insert(m_critical_points_area.end(), m_critical_convex_area.begin(), m_critical_convex_area.end());
}

void 
MyGridSES_SPD :: compute_critical_points_convex() {
        size_t i, j, ball_id = 0;
        vector<double>  v1_xyz(3, 0), v2_xyz(3, 0), v3_xyz(3, 0);
        TID_VERT        v1_id, v2_id, v3_id;
        int             atom_id;

        vector<vector<double> > tmp_rmatrix(3, vector<double>(3, 0));
        vector<double>  tmp_x_axis(3, 0);
        vector<double>  tmp_y_axis(3, 0);
        vector<double>  tmp_z_axis(3, 0);
        vector<double>  tmp_vector(3, 0);

        m_critical_convex.clear();

        for (i = 0; i < m_surface->m_patch_all.size(); i++){
            if (m_surface->m_patch_all[i]->m_patch_type != 3) continue;
            if (m_surface->m_patch_all[i]->m_area < 1) continue;
            ball_id ++;
            vector<double>  tmp_xyz(3, 0);
            vector<double>  tmp_center_xyz(3, 0);
            double          tmp_weight = 0;
            double          tmp_triangle_area = 0;

            for (j = 0; j < m_surface->m_patch_all[i]->m_triangles.size(); j ++){
                v1_id = m_surface->m_patch_all[i]->m_triangles[j][0];
                v2_id = m_surface->m_patch_all[i]->m_triangles[j][1];
                v3_id = m_surface->m_patch_all[i]->m_triangles[j][2];
                v1_xyz = m_surface->m_vert_xyz[v1_id];
                v2_xyz = m_surface->m_vert_xyz[v2_id];
                v3_xyz = m_surface->m_vert_xyz[v3_id];

                tmp_triangle_area = compute_triangle_area(v1_xyz, v2_xyz, v3_xyz);
                tmp_xyz[0] = (v1_xyz[0] + v2_xyz[0] + v3_xyz[0])/3;
                tmp_xyz[1] = (v1_xyz[1] + v2_xyz[1] + v3_xyz[1])/3;
                tmp_xyz[2] = (v1_xyz[2] + v2_xyz[2] + v3_xyz[2])/3;

                tmp_weight += tmp_triangle_area;
                tmp_center_xyz[0] += tmp_triangle_area * tmp_xyz[0];
                tmp_center_xyz[1] += tmp_triangle_area * tmp_xyz[1];
                tmp_center_xyz[2] += tmp_triangle_area * tmp_xyz[2];
            }
            tmp_center_xyz[0] /= tmp_weight;
            tmp_center_xyz[1] /= tmp_weight;
            tmp_center_xyz[2] /= tmp_weight;


            // normalize to surface
            v1_id = m_surface->m_patch_all[i]->m_triangles[0][0];
            atom_id = m_surface->m_vert_atom_id[v1_id];

            vector<double>  tmp_atom_center(3, 0);
            tmp_atom_center = m_surface->m_atom_xyz[atom_id];
            double atom_radius = m_surface->m_atom_radius[atom_id];

            tmp_xyz = tmp_center_xyz - tmp_atom_center;
            normalize_unit_3d(tmp_xyz);

            if (m_surface->m_patch_all[i]->m_area < 4){
                tmp_xyz = atom_radius * tmp_xyz;
                tmp_xyz = tmp_xyz + tmp_atom_center;

                m_critical_convex.push_back(tmp_xyz);
                m_critical_convex_area.push_back(m_surface->m_patch_all[i]->m_area);
                m_critical_convex_ids.push_back(pair<size_t, size_t>(ball_id - 1, 0));

                continue;
            }

            ///compute the rotation matrix
            tmp_z_axis = tmp_xyz;
            ////choose the best axis to generate the y axis
            double  curr_min = 0, tmp_dot;
            vector<double>  sel_vector(3, 0);
            vector<double>  tmp_vector(3, 0);
            tmp_vector[0] = 1;
            tmp_vector[1] = 0;
            tmp_vector[2] = 0;
            sel_vector = tmp_vector;
            curr_min = tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
            if (curr_min < 0){
                tmp_vector[0] = -1;
                tmp_vector[1] = 0;
                tmp_vector[2] = 0;
                sel_vector = tmp_vector;
                curr_min = tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
            }
            tmp_vector[0] = 0;
            tmp_vector[1] = 1;
            tmp_vector[2] = 0;
            tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
            if (tmp_dot >= 0 && tmp_dot < curr_min){
                sel_vector = tmp_vector;
                curr_min = tmp_dot;
            }
            else{
                tmp_vector[0] = 0;
                tmp_vector[1] = -1;
                tmp_vector[2] = 0;
                tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
                if (tmp_dot >= 0 && tmp_dot < curr_min){
                    sel_vector = tmp_vector;
                    curr_min = tmp_dot;
                }
            }
            tmp_vector[0] = 0;
            tmp_vector[1] = 0;
            tmp_vector[2] = 1;
            tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
            if (tmp_dot >= 0 && tmp_dot < curr_min){
                sel_vector = tmp_vector;
                curr_min = tmp_dot;
            }
            else{
                tmp_vector[0] = 0;
                tmp_vector[1] = 0;
                tmp_vector[2] = -1;
                tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
                if (tmp_dot >= 0 && tmp_dot < curr_min){
                    sel_vector = tmp_vector;
                    curr_min = tmp_dot;
                }
            }
            ////decide the y axis
            tmp_y_axis = rdl_cross_3d(tmp_z_axis, sel_vector);
            normalize_unit_3d(tmp_y_axis);
            ////decide the x axis
            tmp_x_axis = rdl_cross_3d(tmp_y_axis, tmp_z_axis);
            normalize_unit_3d(tmp_x_axis);

            tmp_rmatrix[0][0] = tmp_x_axis[0];
            tmp_rmatrix[1][0] = tmp_x_axis[1];
            tmp_rmatrix[2][0] = tmp_x_axis[2];
            tmp_rmatrix[0][1] = tmp_y_axis[0];
            tmp_rmatrix[1][1] = tmp_y_axis[1];
            tmp_rmatrix[2][1] = tmp_y_axis[2];
            tmp_rmatrix[0][2] = tmp_z_axis[0];
            tmp_rmatrix[1][2] = tmp_z_axis[1];
            tmp_rmatrix[2][2] = tmp_z_axis[2];

            size_t num_points = static_cast<size_t>(m_surface->m_patch_all[i]->m_area * g_uniform_CP_template_size / (4 * M_PI * atom_radius * atom_radius));

            tmp_xyz = atom_radius * tmp_xyz;
            tmp_xyz = tmp_xyz + tmp_atom_center;
            m_critical_convex.push_back(tmp_xyz);
            m_critical_convex_area.push_back(m_surface->m_patch_all[i]->m_area);
            m_critical_convex_ids.push_back(pair<size_t, size_t>(ball_id - 1, 0));
            for (j = 1; j < num_points; j ++){
                tmp_xyz[0]  = tmp_rmatrix[0][0] * g_uniform_CP_template[j][0];
                tmp_xyz[0] += tmp_rmatrix[0][1] * g_uniform_CP_template[j][1];
                tmp_xyz[0] += tmp_rmatrix[0][2] * g_uniform_CP_template[j][2];
                tmp_xyz[1]  = tmp_rmatrix[1][0] * g_uniform_CP_template[j][0];
                tmp_xyz[1] += tmp_rmatrix[1][1] * g_uniform_CP_template[j][1];
                tmp_xyz[1] += tmp_rmatrix[1][2] * g_uniform_CP_template[j][2];
                tmp_xyz[2]  = tmp_rmatrix[2][0] * g_uniform_CP_template[j][0];
                tmp_xyz[2] += tmp_rmatrix[2][1] * g_uniform_CP_template[j][1];
                tmp_xyz[2] += tmp_rmatrix[2][2] * g_uniform_CP_template[j][2];

                tmp_xyz = atom_radius * tmp_xyz;
                tmp_xyz = tmp_xyz + tmp_atom_center;
                m_critical_convex.push_back(tmp_xyz);
                m_critical_convex_area.push_back(m_surface->m_patch_all[i]->m_area / num_points);
                m_critical_convex_ids.push_back(pair<size_t, size_t>(ball_id - 1, j));
            }
        }
}

void 
MyGridSES_SPD :: compute_critical_points_convex_as_concave() {
        size_t i, j, ball_id = 0;
        vector<double>  v1_xyz(3, 0), v2_xyz(3, 0), v3_xyz(3, 0);
        TID_VERT        v1_id, v2_id, v3_id;
        int             atom_id;

        vector<vector<double> > tmp_rmatrix(3, vector<double>(3, 0));
        vector<double>  tmp_x_axis(3, 0);
        vector<double>  tmp_y_axis(3, 0);
        vector<double>  tmp_z_axis(3, 0);
        vector<double>  tmp_vector(3, 0);

        m_critical_concave.clear();

        for (i = 0; i < m_surface->m_patch_all.size(); i++){
            if (m_surface->m_patch_all[i]->m_patch_type != 3) continue;
            if (m_surface->m_patch_all[i]->m_area < 1) continue;
            ball_id ++;
            vector<double>  tmp_xyz(3, 0);
            vector<double>  tmp_center_xyz(3, 0);
            double          tmp_weight = 0;
            double          tmp_triangle_area = 0;

            for (j = 0; j < m_surface->m_patch_all[i]->m_triangles.size(); j ++){
                v1_id = m_surface->m_patch_all[i]->m_triangles[j][0];
                v2_id = m_surface->m_patch_all[i]->m_triangles[j][1];
                v3_id = m_surface->m_patch_all[i]->m_triangles[j][2];
                v1_xyz = m_surface->m_vert_xyz[v1_id];
                v2_xyz = m_surface->m_vert_xyz[v2_id];
                v3_xyz = m_surface->m_vert_xyz[v3_id];

                tmp_triangle_area = compute_triangle_area(v1_xyz, v2_xyz, v3_xyz);
                tmp_xyz[0] = (v1_xyz[0] + v2_xyz[0] + v3_xyz[0])/3;
                tmp_xyz[1] = (v1_xyz[1] + v2_xyz[1] + v3_xyz[1])/3;
                tmp_xyz[2] = (v1_xyz[2] + v2_xyz[2] + v3_xyz[2])/3;

                tmp_weight += tmp_triangle_area;
                tmp_center_xyz[0] += tmp_triangle_area * tmp_xyz[0];
                tmp_center_xyz[1] += tmp_triangle_area * tmp_xyz[1];
                tmp_center_xyz[2] += tmp_triangle_area * tmp_xyz[2];
            }
            tmp_center_xyz[0] /= tmp_weight;
            tmp_center_xyz[1] /= tmp_weight;
            tmp_center_xyz[2] /= tmp_weight;

            // normalize to surface
            v1_id = m_surface->m_patch_all[i]->m_triangles[0][0];
            atom_id = m_surface->m_vert_atom_id[v1_id];

            vector<double>  tmp_atom_center(3, 0);
            tmp_atom_center = m_surface->m_atom_xyz[atom_id];
            double atom_radius = m_surface->m_atom_radius[atom_id];

            tmp_xyz = tmp_center_xyz - tmp_atom_center;
            normalize_unit_3d(tmp_xyz);

            tmp_xyz = atom_radius * tmp_xyz;
            tmp_xyz = tmp_xyz + tmp_atom_center;

            m_critical_concave.push_back(tmp_xyz);
            m_critical_concave_area.push_back(m_surface->m_patch_all[i]->m_area);
        }
}

void 
MyGridSES_SPD :: compute_critical_points_concave() {
        size_t i, j;
        vector<double>  v1_xyz(3, 0), v2_xyz(3, 0), v3_xyz(3, 0);
        TID_VERT        v1_id, v2_id, v3_id;

        m_critical_concave.clear();

        for (i = 0; i < m_surface->m_patch_all.size(); i++){
            if (m_surface->m_patch_all[i]->m_patch_type != 2) continue;
            if (m_surface->m_patch_all[i]->m_area < 1) continue;
            //if (m_surface->m_patch_all[i]->m_area > 3) continue;

            vector<double>  tmp_xyz(3, 0);
            vector<double>  tmp_center_xyz(3, 0);
            double          tmp_weight = 0;
            double          tmp_triangle_area = 0;

            for (j = 0; j < m_surface->m_patch_all[i]->m_triangles.size(); j ++){
                v1_id = m_surface->m_patch_all[i]->m_triangles[j][0];
                v2_id = m_surface->m_patch_all[i]->m_triangles[j][1];
                v3_id = m_surface->m_patch_all[i]->m_triangles[j][2];
                v1_xyz = m_surface->m_vert_xyz[v1_id];
                v2_xyz = m_surface->m_vert_xyz[v2_id];
                v3_xyz = m_surface->m_vert_xyz[v3_id];

                tmp_triangle_area = compute_triangle_area(v1_xyz, v2_xyz, v3_xyz);
                tmp_xyz[0] = (v1_xyz[0] + v2_xyz[0] + v3_xyz[0])/3;
                tmp_xyz[1] = (v1_xyz[1] + v2_xyz[1] + v3_xyz[1])/3;
                tmp_xyz[2] = (v1_xyz[2] + v2_xyz[2] + v3_xyz[2])/3;

                tmp_weight += tmp_triangle_area;
                tmp_center_xyz[0] += tmp_triangle_area * tmp_xyz[0];
                tmp_center_xyz[1] += tmp_triangle_area * tmp_xyz[1];
                tmp_center_xyz[2] += tmp_triangle_area * tmp_xyz[2];
            }
            tmp_center_xyz[0] /= tmp_weight;
            tmp_center_xyz[1] /= tmp_weight;
            tmp_center_xyz[2] /= tmp_weight;

            //normalization is not necessary

            m_critical_concave.push_back(tmp_center_xyz);
            m_critical_concave_area.push_back(m_surface->m_patch_all[i]->m_area);
        }
}

void 
MyGridSES_SPD :: print_critical_concave_vs_ses_vertex() {
        int num_vertex = m_surface->m_vert.size();
        int num_critical = m_critical_concave.size();
        cout << "#CP=" 
             << num_critical 
             << " out of #VERT=" 
             << num_vertex 
             << " PRUNE(percent)= "
             << 1 - num_critical * 1.0 / num_vertex << endl;
}

void 
MyGridSES_SPD :: print_critical_convex_vs_ses_vertex() {
        int num_vertex = m_surface->m_vert.size();
        int num_critical = m_critical_convex.size();
        cout << "#CP=" 
             << num_critical 
             << " out of #VERT=" 
             << num_vertex 
             << " PRUNE(percent)= "
             << 1 - num_critical * 1.0 / num_vertex << endl;
}

void 
MyGridSES_SPD :: create_context_balls(vector<vector<double> >* pm_unitball_template) {
        size_t i, ball_id = 0;
        MyContextBall*  tmp_context_ball;

        //cout << "m_critical_points.size() = " << m_critical_points.size() << endl;
        //cout << "m_critical_points_area.size() = " << m_critical_points_area.size() << endl;

        m_context_balls.clear();
        m_context_balls.reserve(m_critical_points.size());

        if (m_critical_convex.size() > 0 && m_critical_concave.size() > 0){
            cerr << "ERROR: both convex and concave are counted in???" << endl;
        }

        if (m_critical_convex.size() > 0){
            for (i = 0; i < m_critical_convex_ids.size(); i ++){
                if (m_critical_convex_ids[i].second != 0) continue;
                if (m_critical_convex_ids[i].first != ball_id){
                    cerr << "ERROR: ball_id disordered!" << endl;
                }
                tmp_context_ball = new MyContextBall(m_critical_convex_ids[i].first, m_critical_convex[i], pm_unitball_template);
                tmp_context_ball->m_center_area = m_critical_convex_area[i];
                tmp_context_ball->m_ball_id_serialno = m_critical_convex_ids[i].second;
                m_context_balls.push_back(tmp_context_ball);
                m_context_balls_refine.push_back(vector<MyContextBall*>());
                ball_id ++;
            }

            for (i = 0; i < m_critical_convex_ids.size(); i ++){
                if (m_critical_convex_ids[i].second == 0) continue;
                tmp_context_ball = new MyContextBall(m_critical_convex_ids[i].first, m_critical_convex[i], pm_unitball_template);
                tmp_context_ball->m_center_area = m_critical_convex_area[i];
                tmp_context_ball->m_ball_id_serialno = m_critical_convex_ids[i].second;
                m_context_balls_refine[m_critical_convex_ids[i].first].push_back(tmp_context_ball);
                ball_id ++;
            }

            //for (i = 0; i < m_critical_convex.size(); i ++){
            //  tmp_context_ball = new MyContextBall(m_critical_convex_ids[i].first, m_critical_convex[i], pm_unitball_template);
            //  tmp_context_ball->m_center_area = m_critical_convex_area[i];
            //  tmp_context_ball->m_ball_id_serialno = m_critical_convex_ids[i].second;
            //  m_context_balls.push_back(tmp_context_ball);
            //  ball_id ++;
            //}

        }

        if (m_critical_concave.size() > 0){
            for (i = 0; i < m_critical_concave.size(); i ++){
                tmp_context_ball = new MyContextBall(ball_id, m_critical_concave[i], pm_unitball_template);
                tmp_context_ball->m_center_area = m_critical_concave_area[i];
                tmp_context_ball->m_ball_id_serialno = 0;
                m_context_balls.push_back(tmp_context_ball);
                ball_id ++;
            }

        }
}

//compute solid volume and solid vector
void 
MyGridSES_SPD :: compute_CB_solid_volume(MyContextBall* pm_cb) {
        int ix, iy, iz;
        int ix_center, iy_center, iz_center;

        vector<double>  tmp_center(3, 0);
        tmp_center[0] = pm_cb->m_ball_center[0];
        tmp_center[1] = pm_cb->m_ball_center[1];
        tmp_center[2] = pm_cb->m_ball_center[2];

        ix_center = static_cast<int>((tmp_center[0] - m_min_x) * m_step_multiplier);
        if (tmp_center[0] - (ix_center * m_step_len + m_min_x) > m_step_len * 0.5) ix_center ++;

        iy_center = static_cast<int>((tmp_center[1] - m_min_y) * m_step_multiplier);
        if (tmp_center[1] - (iy_center * m_step_len + m_min_y) > m_step_len * 0.5) iy_center ++;

        iz_center = static_cast<int>((tmp_center[2] - m_min_z) * m_step_multiplier);
        if (tmp_center[2] - (iz_center * m_step_len + m_min_z) > m_step_len * 0.5) iz_center ++;

        size_t i, j;

        double curr_volume = 0;
        double curr_r = 0;
        double cell_volume = m_step_len * m_step_len * m_step_len;
        double cell_volume_half = 0.5 * cell_volume;
        vector<double>  curr_solid_volume(3, 0);
        double          curr_solid_volume_count = 0;

        for (i = 0; i < g_dynamic_radius.size(); i ++){
            curr_r = g_dynamic_radius[i].first;

            //the inner part
            for (j = 0; j < g_dynamic_radius[i].second.first.size(); j ++){
                ix = ix_center + g_dynamic_radius[i].second.first[j][0];
                if (ix < 0 || ix >= m_size_x) continue;

                iy = iy_center + g_dynamic_radius[i].second.first[j][1];
                if (iy < 0 || iy >= m_size_y) continue;

                iz = iz_center + g_dynamic_radius[i].second.first[j][2];
                if (iz < 0 || iz >= m_size_z) continue;

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_INNER || (m_grid[ix][iy][iz] & RDL_MASK) == RDL_CORE_L2){
                    curr_volume += cell_volume;
                    curr_solid_volume[0] += ix * m_step_len + m_min_x;
                    curr_solid_volume[1] += iy * m_step_len + m_min_y;
                    curr_solid_volume[2] += iz * m_step_len + m_min_z;
                    curr_solid_volume_count += 1;
                }
                else if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_SES){
                    curr_volume += cell_volume_half;
                    curr_solid_volume[0] += (ix * m_step_len + m_min_x) * 0.5;
                    curr_solid_volume[1] += (iy * m_step_len + m_min_y) * 0.5;
                    curr_solid_volume[2] += (iz * m_step_len + m_min_z) * 0.5;
                    curr_solid_volume_count += 0.5;
                }
            }


            //the surface (on the template ball) part
            for (j = 0; j < g_dynamic_radius[i].second.second.size(); j ++){
                ix = ix_center + g_dynamic_radius[i].second.second[j][0];
                if (ix < 0 || ix >= m_size_x) continue;

                iy = iy_center + g_dynamic_radius[i].second.second[j][1];
                if (iy < 0 || iy >= m_size_y) continue;

                iz = iz_center + g_dynamic_radius[i].second.second[j][2];
                if (iz < 0 || iz >= m_size_z) continue;

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_INNER || (m_grid[ix][iy][iz] & RDL_MASK) == RDL_CORE_L2 || (m_grid[ix][iy][iz] & RDL_MASK) == RDL_SES){
                    curr_volume += cell_volume_half;
                    curr_solid_volume[0] += (ix * m_step_len + m_min_x) * 0.5;
                    curr_solid_volume[1] += (iy * m_step_len + m_min_y) * 0.5;
                    curr_solid_volume[2] += (iz * m_step_len + m_min_z) * 0.5;
                    curr_solid_volume_count += 0.5;
                }
            }

            pm_cb->m_solid_volume_dynamicR[i] = curr_volume * 3.0 / (4.0 * M_PI * curr_r * curr_r * curr_r);
            pm_cb->m_solid_vector[i][0] = curr_solid_volume[0] / curr_solid_volume_count - pm_cb->m_ball_center[0];
            pm_cb->m_solid_vector[i][1] = curr_solid_volume[1] / curr_solid_volume_count - pm_cb->m_ball_center[1];
            pm_cb->m_solid_vector[i][2] = curr_solid_volume[2] / curr_solid_volume_count - pm_cb->m_ball_center[2];
        }

        return;
}

//compute solid angle
void 
MyGridSES_SPD :: compute_CB_solid_angle(MyContextBall* pm_cb) {
        const vector<vector<double> >& tmp_unitball = m_unitball->m_unitball_dense;
        double r;
        double tmp_solid_angle;
        size_t i, j;
        vector<double>  tmp_point(3, 0);
        vector<double>  tmp_unit(3, 0);
        int ix, iy, iz;

        vector<double>  tmp_center(3, 0);
        tmp_center[0] = pm_cb->m_ball_center[0];
        tmp_center[1] = pm_cb->m_ball_center[1];
        tmp_center[2] = pm_cb->m_ball_center[2];

        for (j = 0; j < g_dynamic_radius.size(); j ++){
            r = g_dynamic_radius[j].first;
            tmp_solid_angle = 0;
            for (i = 0; i < tmp_unitball.size(); i ++){
                tmp_unit = tmp_unitball[i];
                tmp_point = tmp_center + r * tmp_unit;
                if (tmp_point[0] - m_min_x < 0) continue;
                ix = static_cast<int>((tmp_point[0] - m_min_x) * m_step_multiplier);
                if (ix > m_size_x - 1) continue;

                if (tmp_point[1] - m_min_y < 0) continue;
                iy = static_cast<int>((tmp_point[1] - m_min_y) * m_step_multiplier);
                if (iy > m_size_y - 1) continue;

                if (tmp_point[2] - m_min_z < 0) continue;
                iz = static_cast<int>((tmp_point[2] - m_min_z) * m_step_multiplier);
                if (iz > m_size_z - 1) continue;

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_INNER || (m_grid[ix][iy][iz] & RDL_MASK) == RDL_CORE_L2) tmp_solid_angle += 1;
                //else if (m_grid[ix][iy][iz] == RDL_SES) tmp_solid_angle += 0.5;
            }
            pm_cb->m_solid_angle[j] = tmp_solid_angle / (tmp_unitball.size() * 1.0);
        }
        return ;
}

//compute Context Ball's m_rotation_matrix
void 
MyGridSES_SPD :: compute_CB_rotation_matrix(MyContextBall* pm_cb, int pm_sv) {
        bool is_coarse = false;
        if (pm_cb->m_unitball_template->size() == points_default_coarse){
            is_coarse = true;
        }

        vector<double>  tmp_x_axis(3, 0);
        vector<double>  tmp_y_axis(3, 0);
        vector<double>  tmp_z_axis(3, 0);
        vector<double>  tmp_vector(3, 0);
        if (is_coarse){
            tmp_z_axis[0] = 0 - pm_cb->m_solid_vector[pm_sv][0];
            tmp_z_axis[1] = 0 - pm_cb->m_solid_vector[pm_sv][1];
            tmp_z_axis[2] = 0 - pm_cb->m_solid_vector[pm_sv][2];
        }
        else {
            tmp_z_axis[0] = pm_cb->m_solid_vector[pm_sv][0];
            tmp_z_axis[1] = pm_cb->m_solid_vector[pm_sv][1];
            tmp_z_axis[2] = pm_cb->m_solid_vector[pm_sv][2];
        }

        vector<double> tmp_z_old = tmp_z_axis;
        normalize_unit_3d(tmp_z_axis);
        if (rdl_dot_3d(tmp_z_old, tmp_z_axis) < 0.95 * norm_3d(tmp_z_old)){
            cerr << "ERROR: in normalize_unit_3d(tmp_z_axis);" << endl;
        }

        //choose the best axis to generate the y axis
        double  curr_min = 0, tmp_dot;
        vector<double>  sel_vector(3, 0);
        tmp_vector[0] = 1;
        tmp_vector[1] = 0;
        tmp_vector[2] = 0;
        sel_vector = tmp_vector;
        curr_min = tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
        if (curr_min < 0){
            tmp_vector[0] = -1;
            tmp_vector[1] = 0;
            tmp_vector[2] = 0;
            sel_vector = tmp_vector;
            curr_min = tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
        }
        tmp_vector[0] = 0;
        tmp_vector[1] = 1;
        tmp_vector[2] = 0;
        tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
        if (tmp_dot >= 0 && tmp_dot < curr_min){
            sel_vector = tmp_vector;
            curr_min = tmp_dot;
        }
        else{
            tmp_vector[0] = 0;
            tmp_vector[1] = -1;
            tmp_vector[2] = 0;
            tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
            if (tmp_dot >= 0 && tmp_dot < curr_min){
                sel_vector = tmp_vector;
                curr_min = tmp_dot;
            }
        }
        tmp_vector[0] = 0;
        tmp_vector[1] = 0;
        tmp_vector[2] = 1;
        tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
        if (tmp_dot >= 0 && tmp_dot < curr_min){
            sel_vector = tmp_vector;
            curr_min = tmp_dot;
        }
        else{
            tmp_vector[0] = 0;
            tmp_vector[1] = 0;
            tmp_vector[2] = -1;
            tmp_dot = rdl_dot_3d(tmp_z_axis, tmp_vector);
            if (tmp_dot >= 0 && tmp_dot < curr_min){
                sel_vector = tmp_vector;
                curr_min = tmp_dot;
            }
        }

        tmp_y_axis = rdl_cross_3d(tmp_z_axis, sel_vector);
        normalize_unit_3d(tmp_y_axis);
        tmp_x_axis = rdl_cross_3d(tmp_y_axis, tmp_z_axis);
        normalize_unit_3d(tmp_x_axis);

        pm_cb->m_rotation_matrix[0][0] = tmp_x_axis[0];
        pm_cb->m_rotation_matrix[1][0] = tmp_x_axis[1];
        pm_cb->m_rotation_matrix[2][0] = tmp_x_axis[2];

        pm_cb->m_rotation_matrix[0][1] = tmp_y_axis[0];
        pm_cb->m_rotation_matrix[1][1] = tmp_y_axis[1];
        pm_cb->m_rotation_matrix[2][1] = tmp_y_axis[2];

        pm_cb->m_rotation_matrix[0][2] = tmp_z_axis[0];
        pm_cb->m_rotation_matrix[1][2] = tmp_z_axis[1];
        pm_cb->m_rotation_matrix[2][2] = tmp_z_axis[2];
}

// compute Context Ball's m_context_rays 
//  and m_context_shell_L2 based on the same unitball_template
void 
MyGridSES_SPD :: compute_CB_context_rays(MyContextBall* pm_cb) {
        pm_cb->m_context_rays.clear();
        pm_cb->m_context_rays.reserve(pm_cb->m_unitball_template->size());

        vector<double>  tmp_ray_start(3, 0); 
        vector<double>  tmp_unit(3, 0); 
        CONTEXT_RAY_INT_TYPE    tmp_context_ray;
        size_t          i, j;
        int             ix, iy, iz;
        vector<double>  tmp_point(3, 0);

        tmp_ray_start[0] = pm_cb->m_ball_center[0];
        tmp_ray_start[1] = pm_cb->m_ball_center[1];
        tmp_ray_start[2] = pm_cb->m_ball_center[2];

        for (i = 0; i < pm_cb->m_unitball_template->size(); i ++){
            tmp_context_ray = 0;

            tmp_unit[0]  = pm_cb->m_rotation_matrix[0][0] * (*(pm_cb->m_unitball_template))[i][0];
            tmp_unit[0] += pm_cb->m_rotation_matrix[0][1] * (*(pm_cb->m_unitball_template))[i][1];
            tmp_unit[0] += pm_cb->m_rotation_matrix[0][2] * (*(pm_cb->m_unitball_template))[i][2];
            tmp_unit[1]  = pm_cb->m_rotation_matrix[1][0] * (*(pm_cb->m_unitball_template))[i][0];
            tmp_unit[1] += pm_cb->m_rotation_matrix[1][1] * (*(pm_cb->m_unitball_template))[i][1];
            tmp_unit[1] += pm_cb->m_rotation_matrix[1][2] * (*(pm_cb->m_unitball_template))[i][2];
            tmp_unit[2]  = pm_cb->m_rotation_matrix[2][0] * (*(pm_cb->m_unitball_template))[i][0];
            tmp_unit[2] += pm_cb->m_rotation_matrix[2][1] * (*(pm_cb->m_unitball_template))[i][1];
            tmp_unit[2] += pm_cb->m_rotation_matrix[2][2] * (*(pm_cb->m_unitball_template))[i][2];

            for (j = 0; j < CONTEXT_RAY_INT_SIZE; j ++){
                tmp_point = tmp_ray_start + (((j + 0.5)* g_sample_region_radius) / (CONTEXT_RAY_INT_SIZE * 1.0)) * tmp_unit;

                if (tmp_point[0] - m_min_x < 0){
                    continue;
                }
                ix = static_cast<int>((tmp_point[0] - m_min_x) * m_step_multiplier);
                if (ix > m_size_x - 1){
                    continue;
                }

                if (tmp_point[1] - m_min_y < 0){
                    continue;
                }
                iy = static_cast<int>((tmp_point[1] - m_min_y) * m_step_multiplier);
                if (iy > m_size_y - 1){
                    continue;
                }

                if (tmp_point[2] - m_min_z < 0){
                    continue;
                }
                iz = static_cast<int>((tmp_point[2] - m_min_z) * m_step_multiplier);
                if (iz > m_size_z - 1){
                    continue;
                }

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_INNER || (m_grid[ix][iy][iz] & RDL_MASK) == RDL_CORE_L2){
                    tmp_context_ray |= (0x1u << j);
                }
            }
            pm_cb->m_context_rays.push_back(tmp_context_ray);
        }
}

void 
MyGridSES_SPD :: compute_CB_context_rays_dockinfo(MyContextBall* pm_cb) {
        pm_cb->m_context_rays.clear();
        pm_cb->m_context_rays.reserve(pm_cb->m_unitball_template->size());

        vector<double>  tmp_ray_start(3, 0); 
        vector<double>  tmp_unit(3, 0); 
        CONTEXT_RAY_INT_TYPE    tmp_context_ray;
        size_t          i, j;
        int             ix, iy, iz;
        vector<double>  tmp_point(3, 0);

        tmp_ray_start[0] = pm_cb->m_ball_center[0];
        tmp_ray_start[1] = pm_cb->m_ball_center[1];
        tmp_ray_start[2] = pm_cb->m_ball_center[2];

        for (i = 0; i < pm_cb->m_unitball_template->size(); i ++){
            tmp_context_ray = 0;

			tmp_unit[0] = (*(pm_cb->m_unitball_template))[i][0];
			tmp_unit[1] = (*(pm_cb->m_unitball_template))[i][1];
			tmp_unit[2] = (*(pm_cb->m_unitball_template))[i][2];


            //tmp_unit[0]  = pm_cb->m_rotation_matrix[0][0] * (*(pm_cb->m_unitball_template))[i][0];
            //tmp_unit[0] += pm_cb->m_rotation_matrix[0][1] * (*(pm_cb->m_unitball_template))[i][1];
            //tmp_unit[0] += pm_cb->m_rotation_matrix[0][2] * (*(pm_cb->m_unitball_template))[i][2];
            //tmp_unit[1]  = pm_cb->m_rotation_matrix[1][0] * (*(pm_cb->m_unitball_template))[i][0];
            //tmp_unit[1] += pm_cb->m_rotation_matrix[1][1] * (*(pm_cb->m_unitball_template))[i][1];
            //tmp_unit[1] += pm_cb->m_rotation_matrix[1][2] * (*(pm_cb->m_unitball_template))[i][2];
            //tmp_unit[2]  = pm_cb->m_rotation_matrix[2][0] * (*(pm_cb->m_unitball_template))[i][0];
            //tmp_unit[2] += pm_cb->m_rotation_matrix[2][1] * (*(pm_cb->m_unitball_template))[i][1];
            //tmp_unit[2] += pm_cb->m_rotation_matrix[2][2] * (*(pm_cb->m_unitball_template))[i][2];

            for (j = 0; j < CONTEXT_RAY_INT_SIZE; j ++){
                tmp_point = tmp_ray_start + (((j + 0.5)* g_sample_region_radius) / (CONTEXT_RAY_INT_SIZE * 1.0)) * tmp_unit;

                if (tmp_point[0] - m_min_x < 0){
                    continue;
                }
                ix = static_cast<int>((tmp_point[0] - m_min_x) * m_step_multiplier);
                if (ix > m_size_x - 1){
                    continue;
                }

                if (tmp_point[1] - m_min_y < 0){
                    continue;
                }
                iy = static_cast<int>((tmp_point[1] - m_min_y) * m_step_multiplier);
                if (iy > m_size_y - 1){
                    continue;
                }

                if (tmp_point[2] - m_min_z < 0){
                    continue;
                }
                iz = static_cast<int>((tmp_point[2] - m_min_z) * m_step_multiplier);
                if (iz > m_size_z - 1){
                    continue;
                }

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_INNER || (m_grid[ix][iy][iz] & RDL_MASK) == RDL_CORE_L2){
                    tmp_context_ray |= (0x1u << j);
                }
            }
            pm_cb->m_context_rays.push_back(tmp_context_ray);
        }
}


//like buried verts, shell_L2 should be always densely sampled
void 
MyGridSES_SPD :: compute_CB_context_shell_L2(MyContextBall* pm_cb) {
        pm_cb->m_context_shell_L2.clear();
        pm_cb->m_context_core_L2.clear();
        pm_cb->m_context_inner.clear();
		pm_cb->m_context_inner_2A.clear();
		pm_cb->m_context_inner_core3A.clear();
		pm_cb->m_context_core1A_shell3A.clear();
        pm_cb->m_context_core1Ashell3A.clear();

        pm_cb->m_context_shell_1A.clear(); 
        pm_cb->m_context_shell_2A.clear(); 
        pm_cb->m_context_shell_3A.clear(); 
        pm_cb->m_context_shell_4A.clear(); 
        pm_cb->m_context_core_1A.clear();  
        pm_cb->m_context_core_2A.clear();  
        pm_cb->m_context_core_3A.clear();  
        pm_cb->m_context_core_4A.clear();  

        pm_cb->m_context_shell_L2.reserve(m_unitball->m_unitball_dense.size());
        pm_cb->m_context_core_L2.reserve(m_unitball->m_unitball_dense.size());
        pm_cb->m_context_inner.reserve(m_unitball->m_unitball_dense.size());
		pm_cb->m_context_inner_2A.reserve(m_unitball->m_unitball_dense.size());
		pm_cb->m_context_inner_core3A.reserve(m_unitball->m_unitball_dense.size());
		pm_cb->m_context_core1A_shell3A.reserve(m_unitball->m_unitball_dense.size());
        pm_cb->m_context_core1Ashell3A.reserve(m_unitball->m_unitball_dense.size());

        pm_cb->m_context_shell_1A.reserve(m_unitball->m_unitball_dense.size()); 
        pm_cb->m_context_shell_2A.reserve(m_unitball->m_unitball_dense.size()); 
        pm_cb->m_context_shell_3A.reserve(m_unitball->m_unitball_dense.size()); 
        pm_cb->m_context_shell_4A.reserve(m_unitball->m_unitball_dense.size()); 
        pm_cb->m_context_core_1A.reserve(m_unitball->m_unitball_dense.size());  
        pm_cb->m_context_core_2A.reserve(m_unitball->m_unitball_dense.size());  
        pm_cb->m_context_core_3A.reserve(m_unitball->m_unitball_dense.size());  
        pm_cb->m_context_core_4A.reserve(m_unitball->m_unitball_dense.size());  

        vector<double>  tmp_ray_start(3, 0); 
        vector<double>  tmp_unit(3, 0); 
        double          tmp_dist_sq;
        
        CONTEXT_RAY_INT_TYPE    tmp_context_shell_L2;
        CONTEXT_RAY_INT_TYPE    tmp_context_core_L2;
        CONTEXT_RAY_INT_TYPE    tmp_context_inner;
		CONTEXT_RAY_INT_TYPE    tmp_context_inner_2A;
		CONTEXT_RAY_INT_TYPE    tmp_context_inner_core3A;
		CONTEXT_RAY_INT_TYPE    tmp_context_core1A_shell3A;
        CONTEXT_RAY_INT_TYPE    tmp_context_core1Ashell3A;

        CONTEXT_RAY_INT_TYPE    tmp_context_shell_1A;
        CONTEXT_RAY_INT_TYPE    tmp_context_shell_2A;
        CONTEXT_RAY_INT_TYPE    tmp_context_shell_3A;
        CONTEXT_RAY_INT_TYPE    tmp_context_shell_4A;
        CONTEXT_RAY_INT_TYPE    tmp_context_core_1A;
        CONTEXT_RAY_INT_TYPE    tmp_context_core_2A;
        CONTEXT_RAY_INT_TYPE    tmp_context_core_3A;
        CONTEXT_RAY_INT_TYPE    tmp_context_core_4A;

        size_t          i, j;
        size_t          ii, jj, kk; //index into m_mc_table
        int             ix, iy, iz;
        vector<double>  tmp_point(3, 0);

        tmp_ray_start[0] = pm_cb->m_ball_center[0];
        tmp_ray_start[1] = pm_cb->m_ball_center[1];
        tmp_ray_start[2] = pm_cb->m_ball_center[2];


        for (i = 0; i < m_unitball->m_unitball_dense.size(); i ++){
            tmp_context_shell_L2 = 0;
            tmp_context_core_L2 = 0;
            tmp_context_inner = 0;
            tmp_context_shell_1A = 0;
            tmp_context_shell_2A = 0;
            tmp_context_shell_3A = 0;
            tmp_context_shell_4A = 0;
            tmp_context_core_1A = 0;
            tmp_context_core_2A = 0;
            tmp_context_core_3A = 0;
            tmp_context_core_4A = 0;

            tmp_unit[0]  = pm_cb->m_rotation_matrix[0][0] * m_unitball->m_unitball_dense[i][0];
            tmp_unit[0] += pm_cb->m_rotation_matrix[0][1] * m_unitball->m_unitball_dense[i][1];
            tmp_unit[0] += pm_cb->m_rotation_matrix[0][2] * m_unitball->m_unitball_dense[i][2];
            tmp_unit[1]  = pm_cb->m_rotation_matrix[1][0] * m_unitball->m_unitball_dense[i][0];
            tmp_unit[1] += pm_cb->m_rotation_matrix[1][1] * m_unitball->m_unitball_dense[i][1];
            tmp_unit[1] += pm_cb->m_rotation_matrix[1][2] * m_unitball->m_unitball_dense[i][2];
            tmp_unit[2]  = pm_cb->m_rotation_matrix[2][0] * m_unitball->m_unitball_dense[i][0];
            tmp_unit[2] += pm_cb->m_rotation_matrix[2][1] * m_unitball->m_unitball_dense[i][1];
            tmp_unit[2] += pm_cb->m_rotation_matrix[2][2] * m_unitball->m_unitball_dense[i][2];

            for (j = 0; j < CONTEXT_RAY_INT_SIZE; j ++){
                tmp_point = tmp_ray_start + (((j + 0.5) * g_sample_region_radius) / (CONTEXT_RAY_INT_SIZE * 1.0)) * tmp_unit;

                if (tmp_point[0] - m_min_x < 0){
                    continue;
                }
                ix = static_cast<int>((tmp_point[0] - m_min_x) * m_step_multiplier);
                if (ix > m_size_x - 1){
                    continue;
                }

                if (tmp_point[1] - m_min_y < 0){
                    continue;
                }
                iy = static_cast<int>((tmp_point[1] - m_min_y) * m_step_multiplier);
                if (iy > m_size_y - 1){
                    continue;
                }

                if (tmp_point[2] - m_min_z < 0){
                    continue;
                }
                iz = static_cast<int>((tmp_point[2] - m_min_z) * m_step_multiplier);
                if (iz > m_size_z - 1){
                    continue;
                }

                ii = mc_get_x(m_grid[ix][iy][iz]);
                jj = mc_get_y(m_grid[ix][iy][iz]);
                kk = mc_get_z(m_grid[ix][iy][iz]);

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_SES || 
                    (m_grid[ix][iy][iz] & RDL_MASK) == RDL_SHELL_L2 || 
                    (m_grid[ix][iy][iz] & RDL_MASK) == RDL_SHELL_L2_SURF){
                    
                    tmp_context_shell_L2 |= (0x1u << j);

                    ii = mc_get_x(m_grid[ix][iy][iz]);
                    jj = mc_get_y(m_grid[ix][iy][iz]);
                    kk = mc_get_z(m_grid[ix][iy][iz]);

                    tmp_dist_sq = (ii*ii + jj*jj + kk*kk) * (m_step_len*m_step_len);
                    if (tmp_dist_sq >= 0 && tmp_dist_sq < 1){
                        tmp_context_shell_1A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 1 && tmp_dist_sq < 4){
                        tmp_context_shell_2A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 4 && tmp_dist_sq < 9){
                        tmp_context_shell_3A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 9){
                        tmp_context_shell_4A |= (0x1u << j);
                    }
                }

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_CORE_L2){
                    tmp_context_core_L2 |= (0x1u << j);

                    ii = mc_get_x(m_grid[ix][iy][iz]);
                    jj = mc_get_y(m_grid[ix][iy][iz]);
                    kk = mc_get_z(m_grid[ix][iy][iz]);

                    tmp_dist_sq = (ii*ii + jj*jj + kk*kk) * (m_step_len*m_step_len);
                    if (tmp_dist_sq >= 0 && tmp_dist_sq < 1){
                        tmp_context_core_1A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 1 && tmp_dist_sq < 4){
                        tmp_context_core_2A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 4 && tmp_dist_sq < 9){
                        tmp_context_core_3A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 9){
                        tmp_context_core_4A |= (0x1u << j);
                    }
                }

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_INNER){
                    tmp_context_inner |= (0x1u << j);
                }
            }
			
			tmp_context_inner_2A = tmp_context_inner | tmp_context_core_4A | tmp_context_core_3A | tmp_context_core_2A;
			tmp_context_inner_core3A = tmp_context_inner | tmp_context_core_4A | tmp_context_core_3A;
            tmp_context_core1A_shell3A = tmp_context_core_3A | tmp_context_core_1A | tmp_context_shell_1A | tmp_context_shell_2A | tmp_context_shell_3A;
            tmp_context_core1Ashell3A = tmp_context_core_1A | tmp_context_shell_3A;

			pm_cb->m_context_shell_L2.push_back(tmp_context_shell_L2);
            pm_cb->m_context_core_L2.push_back(tmp_context_core_L2);
            pm_cb->m_context_inner.push_back(tmp_context_inner);
			pm_cb->m_context_inner_2A.push_back(tmp_context_inner_2A);
			pm_cb->m_context_inner_core3A.push_back(tmp_context_inner_core3A);
			pm_cb->m_context_core1A_shell3A.push_back(tmp_context_core1A_shell3A);
            pm_cb->m_context_core1Ashell3A.push_back(tmp_context_core1Ashell3A);

            pm_cb->m_context_shell_1A.push_back(tmp_context_shell_1A); 
            pm_cb->m_context_shell_2A.push_back(tmp_context_shell_2A); 
            pm_cb->m_context_shell_3A.push_back(tmp_context_shell_3A); 
            pm_cb->m_context_shell_4A.push_back(tmp_context_shell_4A); 
            pm_cb->m_context_core_1A.push_back(tmp_context_core_1A); 
            pm_cb->m_context_core_2A.push_back(tmp_context_core_2A);
            pm_cb->m_context_core_3A.push_back(tmp_context_core_3A);
            pm_cb->m_context_core_4A.push_back(tmp_context_core_4A);
        }
}

void 
MyGridSES_SPD :: compute_CB_context_shell_L2_dockinfo(MyContextBall* pm_cb) {
        pm_cb->m_context_shell_L2.clear();
        pm_cb->m_context_core_L2.clear();
        pm_cb->m_context_inner.clear();

        pm_cb->m_context_shell_1A.clear(); 
        pm_cb->m_context_shell_2A.clear(); 
        pm_cb->m_context_shell_3A.clear(); 
        pm_cb->m_context_shell_4A.clear(); 
        pm_cb->m_context_core_1A.clear();  
        pm_cb->m_context_core_2A.clear();  
        pm_cb->m_context_core_3A.clear();  
        pm_cb->m_context_core_4A.clear();  

        pm_cb->m_context_shell_L2.reserve(m_unitball->m_unitball_dense.size());
        pm_cb->m_context_core_L2.reserve(m_unitball->m_unitball_dense.size());
        pm_cb->m_context_inner.reserve(m_unitball->m_unitball_dense.size());

        pm_cb->m_context_shell_1A.reserve(m_unitball->m_unitball_dense.size()); 
        pm_cb->m_context_shell_2A.reserve(m_unitball->m_unitball_dense.size()); 
        pm_cb->m_context_shell_3A.reserve(m_unitball->m_unitball_dense.size()); 
        pm_cb->m_context_shell_4A.reserve(m_unitball->m_unitball_dense.size()); 
        pm_cb->m_context_core_1A.reserve(m_unitball->m_unitball_dense.size());  
        pm_cb->m_context_core_2A.reserve(m_unitball->m_unitball_dense.size());  
        pm_cb->m_context_core_3A.reserve(m_unitball->m_unitball_dense.size());  
        pm_cb->m_context_core_4A.reserve(m_unitball->m_unitball_dense.size());  

        vector<double>  tmp_ray_start(3, 0); 
        vector<double>  tmp_unit(3, 0); 
        double          tmp_dist_sq;
        
        CONTEXT_RAY_INT_TYPE    tmp_context_shell_L2;
        CONTEXT_RAY_INT_TYPE    tmp_context_core_L2;
        CONTEXT_RAY_INT_TYPE    tmp_context_inner;

        CONTEXT_RAY_INT_TYPE    tmp_context_shell_1A;
        CONTEXT_RAY_INT_TYPE    tmp_context_shell_2A;
        CONTEXT_RAY_INT_TYPE    tmp_context_shell_3A;
        CONTEXT_RAY_INT_TYPE    tmp_context_shell_4A;
        CONTEXT_RAY_INT_TYPE    tmp_context_core_1A;
        CONTEXT_RAY_INT_TYPE    tmp_context_core_2A;
        CONTEXT_RAY_INT_TYPE    tmp_context_core_3A;
        CONTEXT_RAY_INT_TYPE    tmp_context_core_4A;

        size_t          i, j;
        size_t          ii, jj, kk; //index into m_mc_table
        int             ix, iy, iz;
        vector<double>  tmp_point(3, 0);

        tmp_ray_start[0] = pm_cb->m_ball_center[0];
        tmp_ray_start[1] = pm_cb->m_ball_center[1];
        tmp_ray_start[2] = pm_cb->m_ball_center[2];


        for (i = 0; i < m_unitball->m_unitball_dense.size(); i ++){
            tmp_context_shell_L2 = 0;
            tmp_context_core_L2 = 0;
            tmp_context_inner = 0;
            tmp_context_shell_1A = 0;
            tmp_context_shell_2A = 0;
            tmp_context_shell_3A = 0;
            tmp_context_shell_4A = 0;
            tmp_context_core_1A = 0;
            tmp_context_core_2A = 0;
            tmp_context_core_3A = 0;
            tmp_context_core_4A = 0;

            tmp_unit[0]  = m_unitball->m_unitball_dense[i][0];
            tmp_unit[1]  = m_unitball->m_unitball_dense[i][1];
            tmp_unit[2]  = m_unitball->m_unitball_dense[i][2];

            //tmp_unit[0]  = pm_cb->m_rotation_matrix[0][0] * m_unitball->m_unitball_dense[i][0];
            //tmp_unit[0] += pm_cb->m_rotation_matrix[0][1] * m_unitball->m_unitball_dense[i][1];
            //tmp_unit[0] += pm_cb->m_rotation_matrix[0][2] * m_unitball->m_unitball_dense[i][2];
            //tmp_unit[1]  = pm_cb->m_rotation_matrix[1][0] * m_unitball->m_unitball_dense[i][0];
            //tmp_unit[1] += pm_cb->m_rotation_matrix[1][1] * m_unitball->m_unitball_dense[i][1];
            //tmp_unit[1] += pm_cb->m_rotation_matrix[1][2] * m_unitball->m_unitball_dense[i][2];
            //tmp_unit[2]  = pm_cb->m_rotation_matrix[2][0] * m_unitball->m_unitball_dense[i][0];
            //tmp_unit[2] += pm_cb->m_rotation_matrix[2][1] * m_unitball->m_unitball_dense[i][1];
            //tmp_unit[2] += pm_cb->m_rotation_matrix[2][2] * m_unitball->m_unitball_dense[i][2];

            for (j = 0; j < CONTEXT_RAY_INT_SIZE; j ++){
                tmp_point = tmp_ray_start + (((j + 0.5) * g_sample_region_radius) / (CONTEXT_RAY_INT_SIZE * 1.0)) * tmp_unit;

                if (tmp_point[0] - m_min_x < 0){
                    continue;
                }
                ix = static_cast<int>((tmp_point[0] - m_min_x) * m_step_multiplier);
                if (ix > m_size_x - 1){
                    continue;
                }

                if (tmp_point[1] - m_min_y < 0){
                    continue;
                }
                iy = static_cast<int>((tmp_point[1] - m_min_y) * m_step_multiplier);
                if (iy > m_size_y - 1){
                    continue;
                }

                if (tmp_point[2] - m_min_z < 0){
                    continue;
                }
                iz = static_cast<int>((tmp_point[2] - m_min_z) * m_step_multiplier);
                if (iz > m_size_z - 1){
                    continue;
                }

                ii = mc_get_x(m_grid[ix][iy][iz]);
                jj = mc_get_y(m_grid[ix][iy][iz]);
                kk = mc_get_z(m_grid[ix][iy][iz]);

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_SES || 
                    (m_grid[ix][iy][iz] & RDL_MASK) == RDL_SHELL_L2 || 
                    (m_grid[ix][iy][iz] & RDL_MASK) == RDL_SHELL_L2_SURF){
                    
                    tmp_context_shell_L2 |= (0x1u << j);

                    ii = mc_get_x(m_grid[ix][iy][iz]);
                    jj = mc_get_y(m_grid[ix][iy][iz]);
                    kk = mc_get_z(m_grid[ix][iy][iz]);

                    tmp_dist_sq = (ii*ii + jj*jj + kk*kk) * (m_step_len*m_step_len);
                    if (tmp_dist_sq >= 0 && tmp_dist_sq < 1){
                        tmp_context_shell_1A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 1 && tmp_dist_sq < 4){
                        tmp_context_shell_2A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 4 && tmp_dist_sq < 9){
                        tmp_context_shell_3A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 9){
                        tmp_context_shell_4A |= (0x1u << j);
                    }
                }

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_CORE_L2){
                    tmp_context_core_L2 |= (0x1u << j);

                    ii = mc_get_x(m_grid[ix][iy][iz]);
                    jj = mc_get_y(m_grid[ix][iy][iz]);
                    kk = mc_get_z(m_grid[ix][iy][iz]);

                    tmp_dist_sq = (ii*ii + jj*jj + kk*kk) * (m_step_len*m_step_len);
                    if (tmp_dist_sq >= 0 && tmp_dist_sq < 1){
                        tmp_context_core_1A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 1 && tmp_dist_sq < 4){
                        tmp_context_core_2A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 4 && tmp_dist_sq < 9){
                        tmp_context_core_3A |= (0x1u << j);
                    }
                    else if (tmp_dist_sq >= 9){
                        tmp_context_core_4A |= (0x1u << j);
                    }
                }

                if ((m_grid[ix][iy][iz] & RDL_MASK) == RDL_INNER){
                    tmp_context_inner |= (0x1u << j);
                }
            }
            pm_cb->m_context_shell_L2.push_back(tmp_context_shell_L2);
            pm_cb->m_context_core_L2.push_back(tmp_context_core_L2);
            pm_cb->m_context_inner.push_back(tmp_context_inner);
            pm_cb->m_context_shell_1A.push_back(tmp_context_shell_1A); 
            pm_cb->m_context_shell_2A.push_back(tmp_context_shell_2A); 
            pm_cb->m_context_shell_3A.push_back(tmp_context_shell_3A); 
            pm_cb->m_context_shell_4A.push_back(tmp_context_shell_4A); 
            pm_cb->m_context_core_1A.push_back(tmp_context_core_1A); 
            pm_cb->m_context_core_2A.push_back(tmp_context_core_2A);
            pm_cb->m_context_core_3A.push_back(tmp_context_core_3A);
            pm_cb->m_context_core_4A.push_back(tmp_context_core_4A);
        }
}

void 
MyGridSES_SPD :: compute_CB_buried_verts(MyContextBall* pm_cb) {
        vector<unsigned int>    vec_local_verts_R;
        vector<double>          tmp_center(3, 0);
        tmp_center[0] = pm_cb->m_ball_center[0];
        tmp_center[1] = pm_cb->m_ball_center[1];
        tmp_center[2] = pm_cb->m_ball_center[2];
        m_surface->query_local(tmp_center, g_sample_region_radius, vec_local_verts_R);

        unsigned int    shift_dist;
        double          vert_dist;
        vector<double>  tmp_vector(3, 0);
        vector<double>  tmp_vector_rotated(3, 0);
        unsigned int    tmp_ray_id;

        CONTEXT_RAY_INT_TYPE    tmp_vert_mask;
        MyBuriedVert*   tmp_buried_vert;

        size_t i;

        pm_cb->m_buried_verts.clear();
        pm_cb->m_buried_verts.reserve(vec_local_verts_R.size());

        for (i = 0; i < vec_local_verts_R.size(); i ++){
            tmp_vector = m_surface->m_vert[vec_local_verts_R[i]]->point() - tmp_center;
            vert_dist = norm_3d(tmp_vector);
            if (vert_dist < 0.2) continue;
            if (vert_dist == g_sample_region_radius){
                shift_dist = CONTEXT_RAY_INT_SIZE;
            }
            else{
                shift_dist = static_cast<unsigned int>((vert_dist * 1.0 * CONTEXT_RAY_INT_SIZE)/g_sample_region_radius);
                if ((vert_dist * 1.0 * CONTEXT_RAY_INT_SIZE)/g_sample_region_radius - shift_dist > 0.5){
                    shift_dist ++;
                }
            }

            if (shift_dist == 0){
                tmp_vert_mask = 0;
            }
            else{
                tmp_vert_mask = (0x1u << (shift_dist - 1));
            }

            normalize_unit_3d(tmp_vector);

            //use M^T to rotate the vert into the normalized space
            tmp_vector_rotated[0]  = pm_cb->m_rotation_matrix[0][0] * tmp_vector[0];
            tmp_vector_rotated[0] += pm_cb->m_rotation_matrix[1][0] * tmp_vector[1];
            tmp_vector_rotated[0] += pm_cb->m_rotation_matrix[2][0] * tmp_vector[2];
            tmp_vector_rotated[1]  = pm_cb->m_rotation_matrix[0][1] * tmp_vector[0];
            tmp_vector_rotated[1] += pm_cb->m_rotation_matrix[1][1] * tmp_vector[1];
            tmp_vector_rotated[1] += pm_cb->m_rotation_matrix[2][1] * tmp_vector[2];
            tmp_vector_rotated[2]  = pm_cb->m_rotation_matrix[0][2] * tmp_vector[0];
            tmp_vector_rotated[2] += pm_cb->m_rotation_matrix[1][2] * tmp_vector[1];
            tmp_vector_rotated[2] += pm_cb->m_rotation_matrix[2][2] * tmp_vector[2];

            //now tmp_vector_rotated is in the normalized reference frame.

            tmp_ray_id = m_unitball->m_dense.query_nn(tmp_vector_rotated);

            tmp_buried_vert = new MyBuriedVert;
            tmp_buried_vert->m_dist = vert_dist;
            tmp_buried_vert->m_xyz = m_surface->m_vert[vec_local_verts_R[i]]->point();
            tmp_buried_vert->m_ray_id = tmp_ray_id;
            tmp_buried_vert->m_vert_mask = tmp_vert_mask;
            tmp_buried_vert->m_vert_id = vec_local_verts_R[i];
            tmp_buried_vert->m_area = m_surface->m_vert[vec_local_verts_R[i]]->m_area;
            pm_cb->m_buried_verts.push_back(tmp_buried_vert);
        }
}

void 
MyGridSES_SPD :: compute_CB_buried_verts_dockinfo(MyContextBall* pm_cb) {
        vector<unsigned int>    vec_local_verts_R;
        vector<double>          tmp_center(3, 0);
        tmp_center[0] = pm_cb->m_ball_center[0];
        tmp_center[1] = pm_cb->m_ball_center[1];
        tmp_center[2] = pm_cb->m_ball_center[2];
        m_surface->query_local(tmp_center, g_sample_region_radius, vec_local_verts_R);

        unsigned int    shift_dist;
        double          vert_dist;
        vector<double>  tmp_vector(3, 0);
        vector<double>  tmp_vector_rotated(3, 0);
        unsigned int    tmp_ray_id;

        CONTEXT_RAY_INT_TYPE    tmp_vert_mask;
        MyBuriedVert*   tmp_buried_vert;

        size_t i;

        pm_cb->m_buried_verts.clear();
        pm_cb->m_buried_verts.reserve(vec_local_verts_R.size());

        for (i = 0; i < vec_local_verts_R.size(); i ++){
            tmp_vector = m_surface->m_vert[vec_local_verts_R[i]]->point() - tmp_center;
            vert_dist = norm_3d(tmp_vector);
            if (vert_dist < 0.2) continue;
            if (vert_dist == g_sample_region_radius){
                shift_dist = CONTEXT_RAY_INT_SIZE;
            }
            else{
                shift_dist = static_cast<unsigned int>((vert_dist * 1.0 * CONTEXT_RAY_INT_SIZE)/g_sample_region_radius);
                if ((vert_dist * 1.0 * CONTEXT_RAY_INT_SIZE)/g_sample_region_radius - shift_dist > 0.5){
                    shift_dist ++;
                }
            }

            if (shift_dist == 0){
                tmp_vert_mask = 0;
            }
            else{
                tmp_vert_mask = (0x1u << (shift_dist - 1));
            }

            normalize_unit_3d(tmp_vector);

            tmp_vector_rotated[0]  = tmp_vector[0];
            tmp_vector_rotated[1]  = tmp_vector[1];
            tmp_vector_rotated[2]  = tmp_vector[2];

            //tmp_vector_rotated[0]  = pm_cb->m_rotation_matrix[0][0] * tmp_vector[0];
            //tmp_vector_rotated[0] += pm_cb->m_rotation_matrix[1][0] * tmp_vector[1];
            //tmp_vector_rotated[0] += pm_cb->m_rotation_matrix[2][0] * tmp_vector[2];
            //tmp_vector_rotated[1]  = pm_cb->m_rotation_matrix[0][1] * tmp_vector[0];
            //tmp_vector_rotated[1] += pm_cb->m_rotation_matrix[1][1] * tmp_vector[1];
            //tmp_vector_rotated[1] += pm_cb->m_rotation_matrix[2][1] * tmp_vector[2];
            //tmp_vector_rotated[2]  = pm_cb->m_rotation_matrix[0][2] * tmp_vector[0];
            //tmp_vector_rotated[2] += pm_cb->m_rotation_matrix[1][2] * tmp_vector[1];
            //tmp_vector_rotated[2] += pm_cb->m_rotation_matrix[2][2] * tmp_vector[2];

            tmp_ray_id = m_unitball->m_dense.query_nn(tmp_vector_rotated);

            tmp_buried_vert = new MyBuriedVert;
            tmp_buried_vert->m_dist = vert_dist;
            tmp_buried_vert->m_xyz = m_surface->m_vert[vec_local_verts_R[i]]->point();
            tmp_buried_vert->m_ray_id = tmp_ray_id;
            tmp_buried_vert->m_vert_mask = tmp_vert_mask;
            tmp_buried_vert->m_vert_id = vec_local_verts_R[i];
            tmp_buried_vert->m_area = m_surface->m_vert[vec_local_verts_R[i]]->m_area;
            pm_cb->m_buried_verts.push_back(tmp_buried_vert);
        }
}

void 
MyGridSES_SPD :: init(double pm_min_x, double pm_min_y, double pm_min_z, double pm_max_x, double pm_max_y, double pm_max_z) {
        m_critical_belt.clear();
        m_critical_concave.clear();
        m_critical_convex.clear();
        m_critical_points.clear();

        m_critical_belt_area.clear();
        m_critical_concave_area.clear();
        m_critical_convex_area.clear();
        m_critical_points_area.clear();

        mc_init();

        double time_begin = 0, time_end = 0;
        struct timeb tp;

        ftime(&tp);
        time_begin = tp.time + (tp.millitm * 1.0) / 1000;
        //allocate_grid, init each cell to RDL_OUTER
        allocate_grid(pm_min_x, pm_min_y, pm_min_z, pm_max_x, pm_max_y, pm_max_z);
        ftime(&tp);
        time_end = tp.time + (tp.millitm * 1.0) / 1000;
        cout << "\tINFO: time to allocate_grid " << time_end - time_begin << " seconds." << endl; 

        ftime(&tp);
        time_begin = tp.time + (tp.millitm * 1.0) / 1000;
        //set_surface_cells, set RDL_SES
        set_surface_cells(m_surface, RDL_SES);
        ftime(&tp);
        time_end = tp.time + (tp.millitm * 1.0) / 1000;
        cout << "\tINFO: time to set_surface_cells(m_surface, RDL_SES) " << time_end - time_begin << " seconds." << endl; 

        ftime(&tp);
        time_begin = tp.time + (tp.millitm * 1.0) / 1000;
        //set_inner_cells to RDL_INNER
        set_inner_cells(RDL_INNER, RDL_OUTER);
        ftime(&tp);
        time_end = tp.time + (tp.millitm * 1.0) / 1000;
        cout << "\tINFO: time to set_inner_cells(RDL_INNER, RDL_OUTER) " << time_end - time_begin << " seconds." << endl; 

        ftime(&tp);
        time_begin = tp.time + (tp.millitm * 1.0) / 1000;
        //marching cubes
        mc_looper();
        ftime(&tp);
        time_end = tp.time + (tp.millitm * 1.0) / 1000;
        cout << "\tINFO: time to mc_looper(); " << time_end - time_begin << " seconds." << endl; 


        size_t i, j, k;
        size_t num_inner = 0;
        size_t num_outer = 0;
        size_t num_ses = 0;
        size_t num_shell_L2_surf = 0;
        size_t num_shell_L2 = 0;
        size_t num_core_L2_surf = 0;
        size_t num_core_L2 = 0;


        for (i = 0; i < m_grid.size(); i ++){
            for (j = 0; j < m_grid[i].size(); j ++){
                for (k = 0; k < m_grid[i][j].size(); k ++){
                    if ((m_grid[i][j][k] & RDL_MASK) == RDL_OUTER) num_outer ++;
                    if ((m_grid[i][j][k] & RDL_MASK) == RDL_INNER) num_inner ++;
                    if ((m_grid[i][j][k] & RDL_MASK) == RDL_SES) num_ses ++;
                    if ((m_grid[i][j][k] & RDL_MASK) == RDL_SHELL_L2) num_shell_L2 ++;
                    if ((m_grid[i][j][k] & RDL_MASK) == RDL_SHELL_L2_SURF) num_shell_L2_surf ++;
                    if ((m_grid[i][j][k] & RDL_MASK) == RDL_CORE_L2) num_core_L2 ++;
                    if ((m_grid[i][j][k] & RDL_MASK) == RDL_CORE_L2_SURF) num_core_L2_surf ++;
                }
            }
        }

        cout << "num_inner = " << num_inner << endl;
        cout << "num_outer = " << num_outer << endl;
        cout << "num_ses = " << num_ses << endl;
        cout << "num_shell_L2 = " << num_shell_L2 << endl;
        cout << "num_shell_L2_surf = " << num_shell_L2_surf << endl;
        cout << "num_core_L2 = " << num_core_L2 << endl;
        cout << "num_core_L2_surf = " << num_core_L2_surf << endl;

}   

void 
MyGridSES_SPD :: allocate_grid(double pm_min_x, double pm_min_y, double pm_min_z, double pm_max_x, double pm_max_y, double pm_max_z) {
        m_min_x = pm_min_x;
        m_min_y = pm_min_y;
        m_min_z = pm_min_z;
        m_max_x = pm_max_x;
        m_max_y = pm_max_y;
        m_max_z = pm_max_z;

        int x_len, y_len, z_len;

        x_len = static_cast<int>((m_max_x - m_min_x) * m_step_multiplier + 1);
        y_len = static_cast<int>((m_max_y - m_min_y) * m_step_multiplier + 1);
        z_len = static_cast<int>((m_max_z - m_min_z) * m_step_multiplier + 1);

        m_size_x = x_len;
        m_size_y = y_len;
        m_size_z = z_len;

        vector<unsigned int>            z_block(z_len, RDL_OUTER);
        vector<vector<unsigned int> >   yz_block(y_len, z_block);
        m_grid.insert(m_grid.end(), x_len, yz_block);

        cout << "size_grid_x = " << m_max_x - m_min_x << endl;
        cout << "size_grid_y = " << m_max_y - m_min_y << endl;
        cout << "size_grid_z = " << m_max_z - m_min_z << endl;
        //cout << "size_grid = " << m_size_x * m_size_y * m_size_z << endl;

        g_sizeof_grid = m_size_x * m_size_y * m_size_z * sizeof(unsigned int);
}

void 
MyGridSES_SPD :: set_surface_cells(MySES* pm_surface, unsigned int pm_mark) {
        set<vector<double> >    tmp_trilet;
        set<vector<double> >::iterator it_dense_point;
        int ix, iy, iz;
        size_t i;
        vector<double>          tmp_point(3, 0);
        vector<vector<double> > tmp_triangle(3, tmp_point);

        for (i = 0; i < pm_surface->m_triangle.size(); i ++){
            //cout << "i = " << i << " in set_surface_cells(MySES* pm_surface, unsigned int pm_mark)" << endl;
            tmp_point = pm_surface->m_vert[pm_surface->m_triangle[i][0]]->point();
            tmp_triangle[0] = tmp_point;
            tmp_point = pm_surface->m_vert[pm_surface->m_triangle[i][1]]->point();
            tmp_triangle[1] = tmp_point;
            tmp_point = pm_surface->m_vert[pm_surface->m_triangle[i][2]]->point();
            tmp_triangle[2] = tmp_point;

            tmp_trilet.clear();
            compute_trianglet(tmp_triangle, m_step_len*0.95, tmp_trilet);
            for (it_dense_point = tmp_trilet.begin(); it_dense_point != tmp_trilet.end(); it_dense_point ++){
                ix = static_cast<int>(((*it_dense_point)[0] - m_min_x) * m_step_multiplier);
                iy = static_cast<int>(((*it_dense_point)[1] - m_min_y) * m_step_multiplier);
                iz = static_cast<int>(((*it_dense_point)[2] - m_min_z) * m_step_multiplier);
                m_grid[ix][iy][iz] = pm_mark;
            }
        }
}

void 
MyGridSES_SPD :: mc_cells(unsigned int pm_mark_init, unsigned int pm_mark, size_t pm_i, size_t pm_j, size_t pm_k) {
        size_t c_i, c_j, c_k;    //current cell's i, j, k, the indices into the mc_table
        size_t ii, jj, kk;       //indices into m_mc_table
        size_t nb_i, nb_j, nb_k; //neighbor cell's i, j, k, the indices into the mc_table
        size_t i, j, k;          //indices into m_grid

        vector<size_t>  tmp_cell(3, 0);
        vector<size_t>  cur_cell(3, 0);

        list<vector<size_t> > tmp_task_list;
        tmp_cell[0] = pm_i;
        tmp_cell[1] = pm_j;
        tmp_cell[2] = pm_k;
        tmp_task_list.push_back(tmp_cell);

        while(tmp_task_list.size() > 0){
            cur_cell = tmp_task_list.front();
            tmp_task_list.pop_front();

            if (cur_cell[0] > 0){
                c_i = mc_get_x(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                if (c_i >= m_mc_size - 1) continue;
                c_j = mc_get_y(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                c_k = mc_get_z(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);

                i = cur_cell[0]-1;
                j = cur_cell[1];
                k = cur_cell[2];

                ii = c_i+1;
                jj = c_j;
                kk = c_k;

                if ((m_grid[i][j][k] & RDL_MASK) == pm_mark){
                    nb_i = mc_get_x(m_grid[i][j][k]);
                    nb_j = mc_get_y(m_grid[i][j][k]);
                    nb_k = mc_get_z(m_grid[i][j][k]);

                    if (m_mc_table[ii][jj][kk] < m_mc_table[nb_i][nb_j][nb_k]){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
                else if ((m_grid[i][j][k] & RDL_MASK) == pm_mark_init){
                    if (m_mc_table[ii][jj][kk] < m_mc_localR){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
            }//if (cur_cell[0] > 0)

            if (cur_cell[0] < m_grid.size() - 1){
                c_i = mc_get_x(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                if (c_i >= m_mc_size - 1) continue;
                c_j = mc_get_y(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                c_k = mc_get_z(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);

                i = cur_cell[0]+1;
                j = cur_cell[1];
                k = cur_cell[2];

                ii = c_i+1;
                jj = c_j;
                kk = c_k;

                if ((m_grid[i][j][k] & RDL_MASK) == pm_mark){
                    nb_i = mc_get_x(m_grid[i][j][k]);
                    nb_j = mc_get_y(m_grid[i][j][k]);
                    nb_k = mc_get_z(m_grid[i][j][k]);

                    if (m_mc_table[ii][jj][kk] < m_mc_table[nb_i][nb_j][nb_k]){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
                else if ((m_grid[i][j][k] & RDL_MASK) == pm_mark_init){
                    if (m_mc_table[ii][jj][kk] < m_mc_localR){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
            }//if (cur_cell[0] < m_grid.size() - 1)


            if (cur_cell[1] > 0){
                c_j = mc_get_y(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                if (c_j >= m_mc_size - 1) continue;
                c_i = mc_get_x(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                c_k = mc_get_z(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);

                i = cur_cell[0];
                j = cur_cell[1]-1;
                k = cur_cell[2];

                ii = c_i;
                jj = c_j+1;
                kk = c_k;

                if ((m_grid[i][j][k] & RDL_MASK) == pm_mark){
                    nb_i = mc_get_x(m_grid[i][j][k]);
                    nb_j = mc_get_y(m_grid[i][j][k]);
                    nb_k = mc_get_z(m_grid[i][j][k]);

                    if (m_mc_table[ii][jj][kk] < m_mc_table[nb_i][nb_j][nb_k]){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
                else if ((m_grid[i][j][k] & RDL_MASK) == pm_mark_init){
                    if (m_mc_table[ii][jj][kk] < m_mc_localR){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
            }//if (cur_cell[1] > 0)


            if (cur_cell[1] < m_grid.size() - 1){
                c_j = mc_get_y(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                if (c_j >= m_mc_size - 1) continue;
                c_i = mc_get_x(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                c_k = mc_get_z(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);

                i = cur_cell[0];
                j = cur_cell[1]+1;
                k = cur_cell[2];

                ii = c_i;
                jj = c_j+1;
                kk = c_k;

                if ((m_grid[i][j][k] & RDL_MASK) == pm_mark){
                    nb_i = mc_get_x(m_grid[i][j][k]);
                    nb_j = mc_get_y(m_grid[i][j][k]);
                    nb_k = mc_get_z(m_grid[i][j][k]);

                    if (m_mc_table[ii][jj][kk] < m_mc_table[nb_i][nb_j][nb_k]){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
                else if ((m_grid[i][j][k] & RDL_MASK) == pm_mark_init){
                    if (m_mc_table[ii][jj][kk] < m_mc_localR){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
            }//if (cur_cell[1] < m_grid.size() - 1)

            if (cur_cell[2] > 0){
                c_k = mc_get_z(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                if (c_k >= m_mc_size - 1) continue;
                c_j = mc_get_y(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                c_i = mc_get_x(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);


                i = cur_cell[0];
                j = cur_cell[1];
                k = cur_cell[2]-1;

                ii = c_i;
                jj = c_j;
                kk = c_k+1;

                if ((m_grid[i][j][k] & RDL_MASK) == pm_mark){
                    nb_i = mc_get_x(m_grid[i][j][k]);
                    nb_j = mc_get_y(m_grid[i][j][k]);
                    nb_k = mc_get_z(m_grid[i][j][k]);

                    if (m_mc_table[ii][jj][kk] < m_mc_table[nb_i][nb_j][nb_k]){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
                else if ((m_grid[i][j][k] & RDL_MASK) == pm_mark_init){
                    if (m_mc_table[ii][jj][kk] < m_mc_localR){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
            }//if (cur_cell[2] > 0)


            if (cur_cell[2] < m_grid.size() - 1){
                c_k = mc_get_z(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                if (c_k >= m_mc_size - 1) continue;
                c_j = mc_get_y(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);
                c_i = mc_get_x(m_grid[cur_cell[0]][cur_cell[1]][cur_cell[2]]);

                i = cur_cell[0];
                j = cur_cell[1];
                k = cur_cell[2]+1;

                ii = c_i;
                jj = c_j;
                kk = c_k+1;

                if ((m_grid[i][j][k] & RDL_MASK) == pm_mark){
                    nb_i = mc_get_x(m_grid[i][j][k]);
                    nb_j = mc_get_y(m_grid[i][j][k]);
                    nb_k = mc_get_z(m_grid[i][j][k]);

                    if (m_mc_table[ii][jj][kk] < m_mc_table[nb_i][nb_j][nb_k]){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
                else if ((m_grid[i][j][k] & RDL_MASK) == pm_mark_init){
                    if (m_mc_table[ii][jj][kk] < m_mc_localR){
                        m_grid[i][j][k] = (pm_mark & RDL_MASK) | (ii << 16) | (jj << 8) | kk;
                        tmp_cell[0] = i;
                        tmp_cell[1] = j;
                        tmp_cell[2] = k;
                        tmp_task_list.push_front(tmp_cell);
                    }
                }
            }//if (cur_cell[2] < m_grid.size() - 1)
        }//while(tmp_task_list.size() > 0)
}

void 
MyGridSES_SPD :: mc_looper() {
        size_t i, j, k;
        for (i = 0; i < m_grid.size(); i ++){
            for (j = 0; j < m_grid[i].size(); j ++){
                for (k = 0; k < m_grid[i][j].size(); k ++){
                    if ((m_grid[i][j][k] & RDL_MASK) != RDL_SES) continue;

                    if (i > 0){
                        if((m_grid[i-1][j][k] & RDL_MASK) == RDL_OUTER || (m_grid[i-1][j][k] & RDL_MASK) == RDL_SHELL_L2){
                            m_grid[i-1][j][k] = RDL_SHELL_L2 | 0x00010000;
                            mc_cells(RDL_OUTER, RDL_SHELL_L2, i-1, j, k);
                        }
                        else if ((m_grid[i-1][j][k] & RDL_MASK) == RDL_INNER || (m_grid[i-1][j][k] & RDL_MASK) == RDL_CORE_L2){
                            m_grid[i-1][j][k] = RDL_CORE_L2 | 0x00010000;
                            mc_cells(RDL_INNER, RDL_CORE_L2, i-1, j, k);
                        }
                    }
                    if (i < m_grid.size() - 1){
                        if((m_grid[i+1][j][k] & RDL_MASK) == RDL_OUTER || (m_grid[i+1][j][k] & RDL_MASK) == RDL_SHELL_L2){
                            m_grid[i+1][j][k] = RDL_SHELL_L2 | 0x00010000;
                            mc_cells(RDL_OUTER, RDL_SHELL_L2, i+1, j, k);
                        }
                        else if ((m_grid[i+1][j][k] & RDL_MASK) == RDL_INNER || (m_grid[i+1][j][k] & RDL_MASK) == RDL_CORE_L2){
                            m_grid[i+1][j][k] = RDL_CORE_L2 | 0x00010000;
                            mc_cells(RDL_INNER, RDL_CORE_L2, i+1, j, k);
                        }
                    }

                    if (j > 0){
                        if((m_grid[i][j-1][k] & RDL_MASK) == RDL_OUTER || (m_grid[i][j-1][k] & RDL_MASK) == RDL_SHELL_L2){
                            m_grid[i][j-1][k] = RDL_SHELL_L2 | 0x00000100;
                            mc_cells(RDL_OUTER, RDL_SHELL_L2, i, j-1, k);
                        }
                        else if ((m_grid[i][j-1][k] & RDL_MASK) == RDL_INNER || (m_grid[i][j-1][k] & RDL_MASK) == RDL_CORE_L2){
                            m_grid[i][j-1][k] = RDL_CORE_L2 | 0x00000100;
                            mc_cells(RDL_INNER, RDL_CORE_L2, i, j-1, k);
                        }
                    }
                    if (j < m_grid[i].size() - 1){
                        if((m_grid[i][j+1][k] & RDL_MASK) == RDL_OUTER || (m_grid[i][j+1][k] & RDL_MASK) == RDL_SHELL_L2){
                            m_grid[i][j+1][k] = RDL_SHELL_L2 | 0x00000100;
                            mc_cells(RDL_OUTER, RDL_SHELL_L2, i, j+1, k);
                        }
                        else if ((m_grid[i][j+1][k] & RDL_MASK) == RDL_INNER || (m_grid[i][j+1][k] & RDL_MASK) == RDL_CORE_L2){
                            m_grid[i][j+1][k] = RDL_CORE_L2 | 0x00000100;
                            mc_cells(RDL_INNER, RDL_CORE_L2, i, j+1, k);
                        }
                    }

                    if (k > 0){
                        if((m_grid[i][j][k-1] & RDL_MASK) == RDL_OUTER || (m_grid[i][j][k-1] & RDL_MASK) == RDL_SHELL_L2){
                            m_grid[i][j][k-1] = RDL_SHELL_L2 | 0x00000001;
                            mc_cells(RDL_OUTER, RDL_SHELL_L2, i, j, k-1);
                        }
                        else if ((m_grid[i][j][k-1] & RDL_MASK) == RDL_INNER || (m_grid[i][j][k-1] & RDL_MASK) == RDL_CORE_L2){
                            m_grid[i][j][k-1] = RDL_CORE_L2 | 0x00000001;
                            mc_cells(RDL_INNER, RDL_CORE_L2, i, j, k-1);
                        }
                    }
                    if (k < m_grid[i][j].size() - 1){
                        if((m_grid[i][j][k+1] & RDL_MASK) == RDL_OUTER || (m_grid[i][j][k+1] & RDL_MASK) == RDL_SHELL_L2){
                            m_grid[i][j][k+1] = RDL_SHELL_L2 | 0x00000001;
                            mc_cells(RDL_OUTER, RDL_SHELL_L2, i, j, k+1);
                        }
                        else if ((m_grid[i][j][k+1] & RDL_MASK) == RDL_INNER || (m_grid[i][j][k+1] & RDL_MASK) == RDL_CORE_L2){
                            m_grid[i][j][k+1] = RDL_CORE_L2 | 0x00000001;
                            mc_cells(RDL_INNER, RDL_CORE_L2, i, j, k+1);
                        }
                    }
                }
            }
        }
    }

//go from any atom center as seed, mark pm_overriden with pm_mark
void 
MyGridSES_SPD :: set_inner_cells(unsigned int pm_mark, unsigned int pm_overriden) {
        vector<double>  tmp_seed(3, 0);
        int ix, iy, iz;
        size_t i;
		int num_components = 0;
        for (i = 0; i < m_surface->m_atom.size(); i ++){
            tmp_seed = m_surface->m_atom[i]->point();
            ix = static_cast<int>((tmp_seed[0] - m_min_x) * m_step_multiplier);
            iy = static_cast<int>((tmp_seed[1] - m_min_y) * m_step_multiplier);
            iz = static_cast<int>((tmp_seed[2] - m_min_z) * m_step_multiplier);
            if ((m_grid[ix][iy][iz] & RDL_MASK) == pm_overriden){
				num_components ++;
                mark_cells(ix, iy, iz, pm_mark, pm_overriden);
            }
        }

		cout << "#COMPONENTS: " << num_components << endl;
}

//check pm_check's neighbors, if pm_overriden then mark with pm_mark
void 
MyGridSES_SPD :: set_shell_cells(unsigned int pm_check, unsigned int pm_mark, unsigned int pm_overriden) {
        size_t i, j, k;
        for (i = 0; i < m_grid.size(); i ++){
            for (j = 0; j < m_grid[i].size(); j ++){
                for (k = 0; k < m_grid[i][j].size(); k ++){
                    if (m_grid[i][j][k] == pm_check){
                        if (i > 0){
                            if(m_grid[i-1][j][k] == pm_overriden){
                                mark_cells(i-1, j, k, pm_mark, pm_overriden);
                            }
                        }
                        if (i < m_grid.size() - 1){
                            if(m_grid[i+1][j][k] == pm_overriden){
                                mark_cells(i+1, j, k, pm_mark, pm_overriden);
                            }
                        }

                        if (j > 0){
                            if(m_grid[i][j-1][k] == pm_overriden){
                                mark_cells(i, j-1, k, pm_mark, pm_overriden);
                            }
                        }
                        if (j < m_grid[i].size() - 1){
                            if(m_grid[i][j+1][k] == pm_overriden){
                                mark_cells(i, j+1, k, pm_mark, pm_overriden);
                            }
                        }

                        if (k > 0){
                            if(m_grid[i][j][k-1] == pm_overriden){
                                mark_cells(i, j, k-1, pm_mark, pm_overriden);
                            }
                        }
                        if (k < m_grid[i][j].size() - 1){
                            if(m_grid[i][j][k+1] == pm_overriden){
                                mark_cells(i, j, k+1, pm_mark, pm_overriden);
                            }
                        }
                    }
                }
            }
        }
}

void 
MyGridSES_SPD :: mark_cells(int pm_ix, int pm_iy, int pm_iz, unsigned int pm_mark, unsigned int pm_overriden) {
        int ix, iy, iz;

        vector<vector<int> >    tmp_queue;
        vector<int>             tmp_cell(3, 0);
        tmp_cell[0] = pm_ix;
        tmp_cell[1] = pm_iy;
        tmp_cell[2] = pm_iz;
        tmp_queue.push_back(tmp_cell);

        while(tmp_queue.size() > 0){
            tmp_cell = tmp_queue.back();
            tmp_queue.pop_back();
            ix = tmp_cell[0];
            iy = tmp_cell[1];
            iz = tmp_cell[2];
            m_grid[ix][iy][iz] = pm_mark;

            if (ix > 0){
                if((m_grid[ix-1][iy][iz] & RDL_MASK) == pm_overriden){
                    tmp_cell[0] = ix - 1;
                    tmp_cell[1] = iy;
                    tmp_cell[2] = iz;
                    tmp_queue.push_back(tmp_cell);
                }
            }
            if (ix < m_size_x - 1){
                if((m_grid[ix+1][iy][iz] & RDL_MASK) == pm_overriden){
                    tmp_cell[0] = ix + 1;
                    tmp_cell[1] = iy;
                    tmp_cell[2] = iz;
                    tmp_queue.push_back(tmp_cell);
                }
            }

            if (iy > 0){
                if((m_grid[ix][iy-1][iz] & RDL_MASK) == pm_overriden){
                    tmp_cell[0] = ix;
                    tmp_cell[1] = iy - 1;
                    tmp_cell[2] = iz;
                    tmp_queue.push_back(tmp_cell);
                }
            }
            if (iy < m_size_y - 1){
                if((m_grid[ix][iy+1][iz] & RDL_MASK) == pm_overriden){
                    tmp_cell[0] = ix;
                    tmp_cell[1] = iy + 1;
                    tmp_cell[2] = iz;
                    tmp_queue.push_back(tmp_cell);
                }
            }

            if (iz > 0){
                if((m_grid[ix][iy][iz-1] & RDL_MASK) == pm_overriden){
                    tmp_cell[0] = ix;
                    tmp_cell[1] = iy;
                    tmp_cell[2] = iz - 1;
                    tmp_queue.push_back(tmp_cell);
                }
            }
            if (iz < m_size_z - 1){
                if((m_grid[ix][iy][iz+1] & RDL_MASK) == pm_overriden){
                    tmp_cell[0] = ix;
                    tmp_cell[1] = iy;
                    tmp_cell[2] = iz + 1;
                    tmp_queue.push_back(tmp_cell);
                }
            }
        }//end of while
}

void 
MyGridSES_SPD :: compute_trianglet(vector<vector<double> >& pm_triangle, double pm_threshold, set<vector<double> >& rtn_points) {
        rtn_points.clear();

        vector<vector<vector<double> > > tmp_queue;
        vector<vector<double> >          tmp_triangle;
        vector<double>  tmp_v1(3, 0);
        vector<double>  tmp_v2(3, 0);
        vector<double>  tmp_v3(3, 0);
        vector<double>  tmp_v_new(3, 0);
        double          tmp_edge_len;

        vector<pair<double, int> > tmp_sort;

        tmp_triangle = pm_triangle;
        tmp_queue.push_back(pm_triangle);

        while(tmp_queue.size() > 0){
            tmp_triangle = tmp_queue.back();
            tmp_queue.pop_back();
            tmp_sort.clear();
            tmp_v1 = tmp_triangle[0];
            tmp_v2 = tmp_triangle[1];
            tmp_v3 = tmp_triangle[2];

            tmp_edge_len = rdl_vector_ssd(tmp_v1, tmp_v2);
            tmp_sort.push_back(pair<double, int>(tmp_edge_len, 1));
            tmp_edge_len = rdl_vector_ssd(tmp_v2, tmp_v3);
            tmp_sort.push_back(pair<double, int>(tmp_edge_len, 2));
            tmp_edge_len = rdl_vector_ssd(tmp_v3, tmp_v1);
            tmp_sort.push_back(pair<double, int>(tmp_edge_len, 3));

            sort(tmp_sort.begin(), tmp_sort.end(), less_than_univ<double, int>());

            if (tmp_sort[2].first < pm_threshold){
                rtn_points.insert(tmp_v1);
                rtn_points.insert(tmp_v2);
                rtn_points.insert(tmp_v3);
            }
            else{
                if (tmp_sort[2].second == 1){
                    tmp_v_new[0] = (tmp_v1[0] + tmp_v2[0])/2;
                    tmp_v_new[1] = (tmp_v1[1] + tmp_v2[1])/2;
                    tmp_v_new[2] = (tmp_v1[2] + tmp_v2[2])/2;
                    tmp_triangle.clear();
                    tmp_triangle.push_back(tmp_v3);
                    tmp_triangle.push_back(tmp_v1);
                    tmp_triangle.push_back(tmp_v_new);
                    tmp_queue.push_back(tmp_triangle);
                    tmp_triangle.clear();
                    tmp_triangle.push_back(tmp_v3);
                    tmp_triangle.push_back(tmp_v2);
                    tmp_triangle.push_back(tmp_v_new);
                    tmp_queue.push_back(tmp_triangle);
                }
                else if (tmp_sort[2].second == 2){
                    tmp_v_new[0] = (tmp_v3[0] + tmp_v2[0])/2;
                    tmp_v_new[1] = (tmp_v3[1] + tmp_v2[1])/2;
                    tmp_v_new[2] = (tmp_v3[2] + tmp_v2[2])/2;
                    tmp_triangle.clear();
                    tmp_triangle.push_back(tmp_v1);
                    tmp_triangle.push_back(tmp_v2);
                    tmp_triangle.push_back(tmp_v_new);
                    tmp_queue.push_back(tmp_triangle);
                    tmp_triangle.clear();
                    tmp_triangle.push_back(tmp_v1);
                    tmp_triangle.push_back(tmp_v3);
                    tmp_triangle.push_back(tmp_v_new);
                    tmp_queue.push_back(tmp_triangle);
                }
                else{
                    tmp_v_new[0] = (tmp_v3[0] + tmp_v1[0])/2;
                    tmp_v_new[1] = (tmp_v3[1] + tmp_v1[1])/2;
                    tmp_v_new[2] = (tmp_v3[2] + tmp_v1[2])/2;
                    tmp_triangle.clear();
                    tmp_triangle.push_back(tmp_v2);
                    tmp_triangle.push_back(tmp_v1);
                    tmp_triangle.push_back(tmp_v_new);
                    tmp_queue.push_back(tmp_triangle);
                    tmp_triangle.clear();
                    tmp_triangle.push_back(tmp_v2);
                    tmp_triangle.push_back(tmp_v3);
                    tmp_triangle.push_back(tmp_v_new);
                    tmp_queue.push_back(tmp_triangle);
                }
            }
        }
        return;
}

void    
MyGridSES_Dense_SPD :: sload_CB(string pm_PDB_code) {
        string in_filename = g_dir_database_DB + pm_PDB_code + g_cb_suffix_DB + ".txt";
        std::ifstream ifs(in_filename.c_str(), std::ios::binary);
        if (!ifs.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << endl;
            exit(1);
        }

        boost::archive::text_iarchive la(ifs);
        sload_context_balls(la);

        size_t i;
        for (i = 0; i < m_context_balls.size(); i ++){
            m_context_balls[i]->m_unitball_template = &(m_unitball->m_unitball_dense);
        }
}

void 
MyGridSES_Dense_SPD :: ssave_CB(string pm_PDB_code) {
  string out_filename = g_dir_database_DB + pm_PDB_code + g_cb_suffix_DB + ".txt";
  std::ofstream ofs(out_filename.c_str(), std::ios::binary);
  if (!ofs.is_open()){
    cerr << "ERROR: cannot open file " << out_filename << endl;
    exit(1);
  }

  boost::archive::text_oarchive la(ofs);
  ssave_context_balls(la);
}

void 
MyGridSES_Dense_SPD:: compute_context_balls_criticalPointBased() {
  create_context_balls(&(m_unitball->m_unitball_dense));
  MyGridSES_SPD::compute_context_balls_criticalPointBased();
}

void 
MyGridSES_Dense_SPD :: compute_context_balls_criticalPointBased_dockinfo() {
  create_context_balls(&(m_unitball->m_unitball_dense));
  MyGridSES_SPD::compute_context_balls_criticalPointBased_dockinfo();
}

void 
MyGridSES_Dense_SPD :: init(double pm_min_x, double pm_min_y, double pm_min_z, double pm_max_x, double pm_max_y, double pm_max_z) {
  MyGridSES_SPD::init(pm_min_x, pm_min_y, pm_min_z, pm_max_x, pm_max_y, pm_max_z);
}


void
MyGridSES_Coarse_SPD :: sload_CB(string pm_PDB_code) {
        string in_filename = g_dir_database_Query + pm_PDB_code + g_cb_suffix_Query + ".txt";
        std::ifstream ifs(in_filename.c_str(), std::ios::binary);
        if (!ifs.is_open()){
            cerr << "ERROR: cannot open file " << in_filename << endl;
            exit(1);
        }

        boost::archive::text_iarchive la(ifs);
        sload_context_balls(la);

        size_t i;
        for (i = 0; i < m_context_balls.size(); i ++){
            m_context_balls[i]->m_unitball_template = &(m_unitball->m_unitball_coarse);
        }
}

void 
MyGridSES_Coarse_SPD :: ssave_CB(string pm_PDB_code) {
  string out_filename = g_dir_database_Query + pm_PDB_code + g_cb_suffix_Query + ".txt";
  std::ofstream ofs(out_filename.c_str(), std::ios::binary);
  if (!ofs.is_open()){
    cerr << "ERROR: cannot open file " << out_filename << endl;
    exit(1);
  }

  boost::archive::text_oarchive la(ofs);
  ssave_context_balls(la);
}

void 
MyGridSES_Coarse_SPD :: compute_context_balls_criticalPointBased() {
  create_context_balls(&(m_unitball->m_unitball_coarse));
  MyGridSES_SPD::compute_context_balls_criticalPointBased();
}

void 
MyGridSES_Coarse_SPD :: compute_context_balls_criticalPointBased_dockinfo() {
  create_context_balls(&(m_unitball->m_unitball_dense));
  MyGridSES_SPD::compute_context_balls_criticalPointBased_dockinfo();
}

void 
MyGridSES_Coarse_SPD :: init(double pm_min_x, double pm_min_y, double pm_min_z, double pm_max_x, double pm_max_y, double pm_max_z) {
  MyGridSES_SPD::init(pm_min_x, pm_min_y, pm_min_z, pm_max_x, pm_max_y, pm_max_z);
}
