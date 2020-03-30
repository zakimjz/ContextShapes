#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>

using namespace std;

std::string rdl_trim(std::string& s, const std::string& drop = " ")
{
	std::string r=s.erase(s.find_last_not_of(drop)+1);
	return r.erase(0,r.find_first_not_of(drop));
}

void rotate_PDB(string pm_fn_R_PDB, 
				string pm_fn_L_PDB, 
				string pm_fn_out_PDB, 
				vector<double>& pm_T0, 
				vector<vector<double> >& pm_R, 
				vector<double>& pm_T1, 
				int pm_sn)
{
    ifstream inR(pm_fn_R_PDB.c_str());
	ifstream inL(pm_fn_L_PDB.c_str());

	char	buffer[10];
	sprintf(buffer, "%d", pm_sn);
	string out_filename = pm_fn_out_PDB + ".out." + buffer + ".pdb";
	
    ofstream out(out_filename.c_str());

    if (!inR.is_open()){
        cerr << "ERROR: cannot open file " << pm_fn_R_PDB << endl;
        exit(1);
    }

    if (!inL.is_open()){
        cerr << "ERROR: cannot open file " << pm_fn_L_PDB << endl;
        exit(1);
    }

	if (!out.is_open()){
        cerr << "ERROR: cannot open file " << out_filename << endl;
        exit(1);
    }


	vector<double>	tmp_old_point(3, 0);
	vector<double>	tmp_new_point(3, 0);
	string str_line;
	string str_line_new;
	
    while(!inR.eof()){
		getline(inR, str_line);
		rdl_trim(str_line);
		if (str_line.length() > 4) {
			out << str_line << endl;
		}
	}

    while(!inL.eof()){
		getline(inL, str_line);
		rdl_trim(str_line);
		if (str_line.length() < 54) {
			out << str_line << endl;
			continue;
		}
		if (str_line.substr(0, 4) == "ATOM" || str_line.substr(0, 4) == "atom" ||
            str_line.substr(0, 6) == "HETATM" || str_line.substr(0, 6) == "hetatm"){
			str_line_new = str_line.substr(0, 30);
			tmp_old_point[0] = atof(str_line.substr(30, 8).c_str());
			tmp_old_point[1] = atof(str_line.substr(38, 8).c_str());
			tmp_old_point[2] = atof(str_line.substr(46, 8).c_str());
            tmp_new_point[0] = pm_R[0][0] * (tmp_old_point[0] + pm_T0[0]) + pm_R[0][1] * (tmp_old_point[1] + pm_T0[1]) + pm_R[0][2] * (tmp_old_point[2] + pm_T0[2]) + pm_T1[0];
            tmp_new_point[1] = pm_R[1][0] * (tmp_old_point[0] + pm_T0[0]) + pm_R[1][1] * (tmp_old_point[1] + pm_T0[1]) + pm_R[1][2] * (tmp_old_point[2] + pm_T0[2]) + pm_T1[1];
            tmp_new_point[2] = pm_R[2][0] * (tmp_old_point[0] + pm_T0[0]) + pm_R[2][1] * (tmp_old_point[1] + pm_T0[1]) + pm_R[2][2] * (tmp_old_point[2] + pm_T0[2]) + pm_T1[2];
			sprintf(buffer, "%8.3f", tmp_new_point[0]);
			str_line_new = str_line_new + buffer;
			sprintf(buffer, "%8.3f", tmp_new_point[1]);
			str_line_new = str_line_new + buffer;
			sprintf(buffer, "%8.3f", tmp_new_point[2]);
			str_line_new = str_line_new + buffer;
			str_line_new = str_line_new + str_line.substr(54);
			out << str_line_new << endl;
		}
		//else {
		//	out << str_line << endl;
		//}
    }    
    inR.close(); 
	inL.close(); 
	out.close();
}

int main(int argc, char *argv[])
{
    if (argc < 6){
        cerr << "Usage: " << argv[0] << " <fn_R_PDB> <fn_L_PDB> <param_filename> <#rotations> <fn_out_prefix>" << endl;
        exit(1);
    }	

    vector<double>  T0(3, 0);
    vector<double>  T1(3, 0);
    vector<vector<double> > R(3, vector<double>(3, 0));

    ifstream in(argv[3]);
    int num_rotation = atol(argv[4]);

    if (!in.is_open()){
        cerr << "ERROR: cannot open file " << argv[3] << endl;
        exit(1);
    }
    
    for (int i = 1; i <= num_rotation; i ++){
	  	if (in.eof()) break;
	    in >> T0[0] >> T0[1] >> T0[2];
	    in >> R[0][0] >> R[0][1] >> R[0][2]
	       >> R[1][0] >> R[1][1] >> R[1][2]
	       >> R[2][0] >> R[2][1] >> R[2][2];
	    in >> T1[0] >> T1[1] >> T1[2];
	
	    rotate_PDB(argv[1], argv[2], argv[5], T0, R, T1, i);
    }

    return 0;
}
