	#include "include.h"
	

	void initialize(struct Grid *grid){
		int nic = grid->nic;
		int njc = grid->njc;
		real_t *q = grid->q;

		for(int i=0; i<nic; i++){
			for(int j=0; j<njc; j++){
				q[INDEXQ(i, j, 0, nic, njc)] = RHOINF;
				q[INDEXQ(i, j, 1, nic, njc)] = RHOINF*MACHINF;
				q[INDEXQ(i, j, 2, nic, njc)] = 0.0;
				q[INDEXQ(i, j, 3, nic, njc)] = PINF/(GAMMA-1) + 0.5*RHOINF*(MACHINF*MACHINF);
			}
		}

	}

	void calc_l2norm(int ni, int nj, real_t *res, real_t *l2norm){
		*l2norm = 0.0;
		for(int i=0; i<ni-1; i++){
			for(int j=0; j<nj-1; j++){
				*l2norm += res[INDEXQ(i, j, 0, ni-1, nj-1)]*res[INDEXQ(i, j, 0, ni-1, nj-1)];
			//	cout<<res[INDEXQ(i, j, 0, ni-1, nj-1)]<<"\n";
			}
		}
		*l2norm = sqrt(*l2norm)/(ni-1)/(nj-1);
	}

	