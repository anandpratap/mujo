	#ifndef __INCL_INCLUDE
	#define __INCL_INCLUDE
	
	#include <stdio.h>
	#include <iostream>
	#include <fstream>
	#include <iomanip>
	#include <cmath>
	#include <cassert>
	#include <cstdlib>

	#define INDEX(I, J, NI, NJ) ((I)*(NJ) + (J))
	#define NQ 4
	#define INDEXQ(I, J, K, NI, NJ) ((I)*(NJ)*(NQ) + (J)*(NQ) + K)

	#define PI 3.14159265358979323846
	
	
	#define GAMMA 1.4
	#define RHOINF 1.0
	#define MACHINF 0.5
	#define PINF 0.7142857142857143

	#include "adolc/adolc.h"
	#include "adolc/adolc_sparse.h"
	#include "petscksp.h"
	
	typedef double real_t;
	typedef adouble areal_t;

	#include "io.h"
	#include "metrics.h"
	#include "utils.h"
	#include "bc.h"
	#include "flux.h"
	#include "recon.h"
	#include "residual.h"
	#include "petsc.h"
	
	struct Grid{
		int nic, njc;
		int ni, nj;
		real_t *x, *y, *xc, *yc, *vol, *ds_xi, *ds_eta, *nxi_x, *nxi_y, *neta_x, *neta_y, *dt;
		real_t *q;
	public:
		Grid(int nii, int njj){
			ni = nii;
			nj = njj;
			nic = ni - 1;
			njc = nj - 1;

			x = new real_t[ni*nj];
			y = new real_t[ni*nj];
			
			xc = new real_t[nic*njc];
			yc = new real_t[nic*njc];

			vol = new real_t[nic*njc];
			dt = new real_t[nic*njc];
			ds_xi = new real_t[nic*njc];
			ds_eta = new real_t[nic*njc];

			nxi_x = new real_t[ni*njc];
			nxi_y = new real_t[ni*njc];

			neta_x = new real_t[nic*nj];
			neta_y = new real_t[nic*nj];
			q = new real_t[nic*njc*NQ];

		}
		~Grid(){
			delete[] x;
			delete[] y;
			delete[] xc;
			delete[] yc;
			delete[] vol;
			delete[] ds_xi;
			delete[] ds_eta;
			
			delete[] nxi_x;
			delete[] nxi_y;
			delete[] neta_x;
			delete[] neta_y;

			delete[] dt;
			delete[] q;
		}
	};


	
	#endif
