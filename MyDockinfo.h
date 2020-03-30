#ifndef _MYDOCKINFO_H_
#define _MYDOCKINFO_H_

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
#include <unistd.h>

#include "contextBall.h"
#include "util.h"
#include "gridSES.h"

using namespace std;

/* weight counter */
extern double weight_in_1_16bits [0x1u << 16] ;
extern double weight_in_2_16bits [0x1u << 16] ;

double iterated_weight_16bits (CONTEXT_RAY_INT_TYPE n, unsigned int b, unsigned int e);

void compute_weight_in_16bits();

inline double precomputed32_weight (CONTEXT_RAY_INT_TYPE n);

template <class T>
T norm_3d(vector<T>&);

void test_weight_table_32();

struct MyDockinfo
{
  MyGridSES_Dense_SPD*       m_receptor; 
  MyGridSES_Coarse_SPD*       m_ligand; 
  MyUnitBall*          m_unitball;

  double  compute_rmsd_local(MyContextBall* pm_CB_ligand, 
                            vector<double>& V_translate_ligand_to_origin, 
                            vector<vector<double> >& M_rotate_ligand,
                            vector<double>& V_translate_ligand_to_receptor);

  void compute_BSA_Grid(MyGridSES_SPD* pm_GridSES_Cells, MyContextBall* pm_CB_Cells, 
		     MyGridSES_SPD* pm_GridSES_Grid, MyContextBall* pm_CB_Grid, 
			 vector<double>& pm_BSA_inGrid);

  void compute_BSA_CB(MyContextBall* pm_CB_receptor, 
                      MyContextBall* pm_CB_ligand, 
                      vector<double>& pm_BSA_CB); 

  void compute_BOV_CB(MyContextBall* pm_CB_receptor,
                      MyContextBall* pm_CB_ligand, 
		      vector<double>& pm_BOV_CB);

  void compute_BOV_Grid(MyGridSES_SPD* pm_GridSES_Cells, MyContextBall* pm_CB_Cells, 
                        MyGridSES_SPD* pm_GridSES_Grid, MyContextBall* pm_CB_Grid, 
		        vector<double>& pm_BOV_inGrid);


  void match_dockinfo();

	
  ~MyDockinfo();

  MyDockinfo(MyGridSES_Dense_SPD* A, MyGridSES_Coarse_SPD* B, MyUnitBall* unitball);

};

#endif
