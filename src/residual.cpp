	#include "include.h"
	void residual(struct Grid *grid, real_t *res){
		int ni = grid->ni;
		int nj = grid->nj;
		int nic = grid->nic;
		int njc = grid->njc;
    int nip = ni + 1;
    int njp = nj + 1;

    int nvar = njc*nic*NQ;
		real_t *q = grid->q;

    int order = 1;
		trace_on(1);

		areal_t *res_ad = new areal_t[nvar];
		areal_t *q_ad = new areal_t[nvar];

		res_ad[0:nvar] = 0.0;
		q_ad[0:nvar] = 0.0;

    for(int i=0; i<nvar; i++){
      q_ad[i] <<= q[i];
    }

    areal_t rho[nip][njp], u[nip][njp], v[nip][njp], p[nip][njp];

    bc(ni, nj, grid, q_ad, rho, u, v, p);
  
    real_t eig, dsmin;
    for(int i = 0; i < nic; i++) {
      for(int j = 0; j < njc; j++) {
        grid->dt[INDEX(i,j,nic,njc)] = min<real_t>(grid->ds_xi[INDEX(i, j, nic, njc)], grid->ds_xi[INDEX(i, j, nic, njc)])/(sqrt(GAMMA*p[i+1][j+1]/rho[i+1][j+1]) + \
          sqrt(u[i+1][j+1]*u[i+1][j+1] + v[i+1][j+1]*v[i+1][j+1])).value();
      }
    }
  

    areal_t rlft_xi[ni][njc], ulft_xi[ni][njc], vlft_xi[ni][njc], plft_xi[ni][njc];
    areal_t rrht_xi[ni][njc], urht_xi[ni][njc], vrht_xi[ni][njc], prht_xi[ni][njc];
    areal_t flux_xi[ni][njc][4];
 
    grad_xi(ni, nj, rho, u, v, p, rlft_xi, ulft_xi, vlft_xi, plft_xi, rrht_xi, urht_xi, vrht_xi, prht_xi, order);

    for(int i = 0; i < ni; i++) {
      for(int j = 0; j < njc; j++) {
        roeflux(grid->nxi_x[INDEX(i, j, ni, njc)], grid->nxi_y[INDEX(i, j, ni, njc)], rlft_xi[i][j], ulft_xi[i][j], vlft_xi[i][j], plft_xi[i][j], \
          rrht_xi[i][j], urht_xi[i][j], vrht_xi[i][j], prht_xi[i][j], flux_xi[i][j]);
      }
    }
  
    areal_t rlft_eta[nic][nj], ulft_eta[nic][nj], vlft_eta[nic][nj], plft_eta[nic][nj];
    areal_t rrht_eta[nic][nj], urht_eta[nic][nj], vrht_eta[nic][nj], prht_eta[nic][nj];
    areal_t flux_eta[nic][nj][4];

    grad_eta(ni, nj, rho, u, v, p, rlft_eta, ulft_eta, vlft_eta, plft_eta, rrht_eta, urht_eta, vrht_eta, prht_eta, order);
    for(int i = 0; i < nic; i++) {
      for(int j = 0; j < nj; j++) {
        roeflux(grid->neta_x[INDEX(i, j, nic, nj)], grid->neta_y[INDEX(i, j, nic, nj)], rlft_eta[i][j], ulft_eta[i][j], vlft_eta[i][j], plft_eta[i][j], \
          rrht_eta[i][j], urht_eta[i][j], vrht_eta[i][j], prht_eta[i][j], flux_eta[i][j]);
      }
    }

    for(int i=0; i<nic; i++){
      for(int j=0; j<njc; j++){
        for(int n=0; n<NQ; n++){
          res_ad[INDEXQ(i, j, n, nic, njc)] = -(flux_xi[i+1][j][n] - flux_xi[i][j][n]) - (flux_eta[i][j+1][n] - flux_eta[i][j][n]);
        }
      }
    }

    for(int i=0; i<nvar; i++){
      res_ad[i] >>= res[i];
    }

   
    delete[] res_ad;
    delete[] q_ad;
  
    trace_off();
  }