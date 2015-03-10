	#ifndef __INCL_BC
	#define __INCL_BC
	#include "include.h"
	template<class T>
	void bc(int ni, int nj, struct Grid *grid, T *q_ad, T rho[ni+1][nj+1], T u[ni+1][nj+1], T v[ni+1][nj+1], T p[ni+1][nj+1]){
		int nic = grid->nic;
		int njc = grid->njc;
		int nip = ni + 1;
		int njp = nj + 1;

		int nvar = njc*nic*NQ;

		// set interior points from the solution
		for(int i=0; i<nic; i++){
			for(int j=0; j<njc; j++){
				primvars(&q_ad[INDEXQ(i, j, 0, nic, njc)], &rho[i+1][j+1], &u[i+1][j+1], &v[i+1][j+1], &p[i+1][j+1]);
			
			}
		}
		
		// left and right
		for(int j=0; j<njp; j++){
			rho[0][j] = RHOINF;
			u[0][j] = RHOINF*MACHINF;
			v[0][j] = 0.0;
			p[0][j] = PINF;
			rho[ni][j] = RHOINF;
			u[ni][j] = RHOINF*MACHINF;
			v[ni][j] = 0.0;
			p[ni][j] = PINF;	
		}

		// top
		for(int i=0; i<nip; i++){
			rho[i][nj] = RHOINF;
			u[i][nj] = RHOINF*MACHINF;
			v[i][nj] = 0.0;
			p[i][nj] = PINF;
		}

		T un;
		T ds;
		real_t local_neta_x, local_neta_y;
		T uh, vh;
		int j1 = 41;
		int nb = 200;
		for(int i=j1; i<=j1+nb; i++){

			local_neta_x = grid->neta_x[INDEX(i-1, 0, nic, nj)];
			local_neta_y = grid->neta_y[INDEX(i-1, 0, nic, nj)];

			ds = local_neta_x*local_neta_x + local_neta_y*local_neta_y;
			
			p[i][0] = 2.0*p[i][1] - 1.0*p[i][2];
			rho[i][0] = 2.0*rho[i][1] - 1.0*rho[i][2];

			//p[i][0] = p[i][1];
			//rho[i][0] = rho[i][1];
			uh = 2.0*u[i][1] - 1.0*u[i][2];
			vh = 2.0*v[i][1] - 1.0*v[i][2];
			un = uh*local_neta_x + vh*local_neta_y;
			u[i][0] = uh - un*local_neta_x/ds;
			v[i][0] = vh - un*local_neta_y/ds;
		}


		for(int i=1; i<j1; i++){
			rho[i][0] = rho[nic-i+1][1];
			u[i][0] = u[nic-i+1][1];
			v[i][0] = v[nic-i+1][1];
			p[i][0] = p[nic-i+1][1];


			rho[nic-i+1][0] = rho[i][1];
			u[nic-i+1][0] = u[i][1];
			v[nic-i+1][0] = v[i][1];
			p[nic-i+1][0] = p[i][1];
		}

	}

	#endif