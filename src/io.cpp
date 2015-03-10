	#include "include.h"
	
	void read_grid(struct Grid *grid){
		int ni = grid->ni;
		int nj = grid->nj;
		std::ifstream infile("grid.unf2");
		infile >> std::fixed >> std::setprecision(20);
		for(int j=0; j<nj; j++){
			for(int i=0; i<ni; i++){  
				infile>>grid->x[INDEX(i, j, ni, nj)]>>grid->y[INDEX(i, j, ni, nj)];
			}
		}

	}	

	void write_solution(struct Grid *grid){
		int nic = grid->nic;
		int njc = grid->njc;
		real_t *q = grid->q;
		real_t *xc = grid->xc;
		real_t *yc = grid->yc;

		real_t rho, u, v, p;


		std::ofstream outfile;
		outfile.open("data.dat");
		char buffer [500];
		outfile<<"title = \"Solution\""<<"\n";
		outfile<<"variables = \"x\" \"y\" \"rho\" \"u\" \"v\" \"p\""<<"\n"; 
		outfile<<"zone i="<<nic<<", j="<<njc<<", f=point\n";
		for(int j=0; j<njc; j++){
			for(int i=0; i<nic; i++){
				primvars(&q[INDEXQ(i, j, 0, nic, njc)], &rho, &u, &v, &p);
				sprintf(buffer, "%8.14E %8.14E %8.14E %8.14E %8.14E %8.14E\n", xc[INDEX(i, j, nic, njc)], yc[INDEX(i, j, nic, njc)], rho, u, v, p);
				outfile<<buffer;
			}
		}

		outfile.close();
		printf("Data dumped..\n");
	}

	void log_residual(int iter, real_t l2norm){
		std::ofstream outfile;
		
		if(iter == 0){
			outfile.open("residual.log", std::ofstream::out);
		}
		else{
			outfile.open("residual.log", std::ofstream::app);
		}
		outfile<<std::setprecision(20);
		outfile<<iter<<" "<<l2norm<<"\n";
		outfile.close();
	}
