#ifndef _MY_MATCHPAIR_H_
#define _MY_MATCHPAIR_H_

#include "contextBall.h"

struct MyMatchPair{
    MyContextBall*  m_ball_db;
    MyContextBall*  m_ball_query;
	unsigned short  m_pose_id;
	double          m_score;
	double          m_BOV;
	vector<double>          m_BOV_v;
    vector<double>          m_buried_area_v;  
                                            //[ 0] buried_area_inner_inDB
                                            //[ 1] buried_area_core_4A_inDB
                                            //[ 2] buried_area_core_3A_inDB
                                            //[ 3] buried_area_core_2A_inDB
                                            //[ 4] buried_area_core_1A_inDB
                                            //[ 5] buried_area_shell_1A_inDB
                                            //[ 6] buried_area_shell_2A_inDB
                                            //[ 7] buried_area_shell_3A_inDB
                                            //[ 8] buried_area_shell_4A_inDB
                                            //[ 9] buried_area_inner_inQuery
                                            //[10] buried_area_core_4A_inQuery
                                            //[11] buried_area_core_3A_inQuery
                                            //[12] buried_area_core_2A_inQuery
                                            //[13] buried_area_core_1A_inQuery
                                            //[14] buried_area_shell_1A_inQuery
                                            //[15] buried_area_shell_2A_inQuery
                                            //[16] buried_area_shell_3A_inQuery
                                            //[17] buried_area_shell_4A_inQuery

	double          m_dist;
	double          m_rmsd_ca;
    double          m_rmsd_local;
    double          m_rmsd_all;

	vector<double>                  m_translate_query_to_origin;    //[3]
    vector<vector<double> >         m_rotate_query;                 //[3x3]
    vector<double>                  m_translate_query_to_db;        //[3]
};

#endif
