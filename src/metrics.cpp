	#include "include.h"
	void metrics(struct Grid *grid){
		int ni = grid->ni;
		int nj = grid->nj;
		int nic = grid->nic;
		int njc = grid->njc;

		real_t reta_x[nic][nj], reta_y[nic][nj];
		real_t rxi_x[ni][njc], rxi_y[ni][njc];

		for(int i=0; i<nic; i++){
			for(int j=0; j<nj; j++){
				reta_x[i][j] = grid->x[INDEX(i+1,j,ni,nj)] - grid->x[INDEX(i,j,ni,nj)];
				reta_y[i][j] = grid->y[INDEX(i+1,j,ni,nj)] - grid->y[INDEX(i,j,ni,nj)];

				grid->neta_x[INDEX(i, j, nic, nj)] = -reta_y[i][j];
				grid->neta_y[INDEX(i, j, nic, nj)] = reta_x[i][j];
		
			}
		}

		for(int i=0; i<ni; i++){
			for(int j=0; j<njc; j++){
				rxi_x[i][j] = grid->x[INDEX(i,j+1,ni,nj)] - grid->x[INDEX(i,j,ni,nj)];
				rxi_y[i][j] = grid->y[INDEX(i,j+1,ni,nj)] - grid->y[INDEX(i,j,ni,nj)];

				grid->nxi_x[INDEX(i, j, ni, njc)] = rxi_y[i][j];
				grid->nxi_y[INDEX(i, j, ni, njc)] = -rxi_x[i][j];

			}
		}

		for(int i=0; i<nic; i++){
			for(int j=0; j<njc; j++){
				grid->vol[INDEX(i, j, nic, njc)] = 0.5*(reta_x[i][j]*rxi_y[i][j] - rxi_x[i][j]*reta_y[i][j] \
					+ reta_x[i][j+1]*rxi_y[i+1][j] - rxi_x[i+1][j]*reta_y[i][j+1]);

				grid->xc[INDEX(i, j, nic, njc)] = 0.25*(grid->x[INDEX(i, j, ni, nj)] + grid->x[INDEX(i+1, j, ni, nj)]\
					+ grid->x[INDEX(i, j+1, ni, nj)] + grid->x[INDEX(i+1, j+1, ni, nj)]);
				grid->yc[INDEX(i, j, nic, njc)] = 0.25*(grid->y[INDEX(i, j, ni, nj)] + grid->y[INDEX(i+1, j, ni, nj)]\
					+ grid->y[INDEX(i, j+1, ni, nj)] + grid->y[INDEX(i+1, j+1, ni, nj)]);

				grid->ds_eta[INDEX(i, j, nic, njc)] = sqrt(std::pow(grid->neta_x[INDEX(i, j, nic, nj)],2) + std::pow(grid->neta_y[INDEX(i, j, nic, nj)],2));
				grid->ds_xi[INDEX(i, j, nic, njc)] = sqrt(std::pow(grid->nxi_x[INDEX(i, j, ni, njc)],2) + std::pow(grid->nxi_y[INDEX(i, j, ni, njc)],2));
			}
		}

	}
