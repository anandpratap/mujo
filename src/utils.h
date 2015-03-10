	#ifndef __INCL_UTILS
	#define __INCL_UTILS

	#include "include.h"


	template<class T>
	void primvars(T *q, T *rho, T *u, T *v, T *p){
		*rho = q[0];
		*u = q[1]/q[0];
		*v = q[2]/q[0];
		*p = (q[3] - 0.5*(q[1]*q[1] + q[2]*q[2])/q[0])*(GAMMA-1);
	}

	template <class T>
	T min(T a, T b){
		if(a > b){
			return b;
		}	
		else{
			return a;
		}
	}
	
	void calc_l2norm(int ni, int nj, real_t *res, real_t *l2norm);

	void initialize(struct Grid *grid);
	
	#endif