	#ifndef _INCL_RECON
	#define _INCL_RECON
	#include "include.h"
	
	template <class Tad>
	void first_order_xi(int ni, int nj, Tad q[ni+1][nj+1], Tad ql[ni][nj-1], Tad qr[ni][nj-1]){
		int njm = nj-1;
		ql[0:ni][0:njm] = q[0:ni][1:njm];
		qr[0:ni][0:njm] = q[1:ni][1:njm];
	}
	
	template <class Tad>
	void first_order_eta(int ni, int nj, Tad q[ni+1][nj+1], Tad ql[ni-1][nj], Tad qr[ni-1][nj]){
		int nim = ni-1;
		ql[0:nim][0:nj] = q[1:nim][0:nj];
		qr[0:nim][0:nj] = q[1:nim][1:nj];
	}
	
	
	template <class Tad>
	void second_order_xi(int ni, int nj, Tad q[ni+1][nj+1], Tad ql[ni][nj-1], Tad qr[ni][nj-1]){
		int njm = nj-1;
		ql[0:ni][0:njm] = q[0:ni][1:njm];
		qr[0:ni][0:njm] = q[1:ni][1:njm];
		
		int nim = ni-1;
		real_t thm = 2.0/3.0;
		real_t thp = 4.0/3.0;
		real_t eps = pow(10.0/nim, 3);
		Tad f2[ni][njm], a1[nim][njm], a2[nim][njm], f3qt[nim][njm];
		f2[:][:] = q[1:ni][1:njm] - q[0:ni][1:njm];
		a1[:][:] = 3.0*f2[1:nim][:]*f2[0:nim][:];
		a2[:][:] = 2*(f2[1:nim][:] - f2[0:nim][:])*(f2[1:nim][:] - f2[0:nim][:]) + a1[:][:];
		f3qt[:][:] = 0.25*(a1[:][:] + eps)/(a2[:][:] + eps);
		ql[1:nim][:] = ql[1:nim][:] + f3qt[:][:]*(thm*f2[0:nim][:] + thp*f2[1:nim][:]);
		qr[0:nim][:] = qr[0:nim][:] - f3qt[:][:]*(thp*f2[0:nim][:] + thm*f2[1:nim][:]);
	}
	
	template <class Tad>
	void second_order_eta(int ni, int nj, Tad q[ni+1][nj+1], Tad ql[ni-1][nj], Tad qr[ni-1][nj]){
		int nim = ni-1;
		ql[0:nim][0:nj] = q[1:nim][0:nj];
		qr[0:nim][0:nj] = q[1:nim][1:nj];
		
		int njm = nj-1;
		real_t thm = 2.0/3.0;
		real_t thp = 4.0/3.0;
		real_t eps = pow(10.0/njm, 3);
		Tad f2[nim][nj], a1[nim][njm], a2[nim][njm], f3qt[nim][njm];
		f2[:][:] = q[1:nim][1:nj] - q[1:nim][0:nj];
		a1[:][:] = 3.0*f2[:][1:njm]*f2[:][0:njm];
		a2[:][:] = 2*(f2[:][1:njm] - f2[:][0:njm])*(f2[:][1:njm] - f2[:][0:njm]) + a1[:][:];
		f3qt[:][:] = 0.25*(a1[:][:] + eps)/(a2[:][:] + eps);
		ql[:][1:njm] = ql[:][1:njm] + f3qt[:][:]*(thm*f2[:][0:njm] + thp*f2[:][1:njm]);
		qr[:][0:njm] = qr[:][0:njm] - f3qt[:][:]*(thp*f2[:][0:njm] + thm*f2[:][1:njm]);
		
		
	}
	
	template <class Tad>
	void grad_xi(int ni, int nj, Tad rho[ni+1][nj+1], Tad u[ni+1][nj+1], Tad v[ni+1][nj+1], Tad p[ni+1][nj+1],\
		Tad rlft[ni][nj-1], Tad ulft[ni][nj-1], Tad vlft[ni][nj-1], Tad plft[ni][nj-1], \
		Tad rrht[ni][nj-1], Tad urht[ni][nj-1], Tad vrht[ni][nj-1], Tad prht[ni][nj-1], int order = 2){
		
		if(order == 1){
			first_order_xi(ni, nj, rho, rlft, rrht);
			first_order_xi(ni, nj, u, ulft, urht);
			first_order_xi(ni, nj, v, vlft, vrht);
			first_order_xi(ni, nj, p, plft, prht);
		}
		else if(order == 2){
			second_order_xi(ni, nj, rho, rlft, rrht);
			second_order_xi(ni, nj, u, ulft, urht);
			second_order_xi(ni, nj, v, vlft, vrht);
			second_order_xi(ni, nj, p, plft, prht);
		}
		else{
			printf("order not implemented");
			exit(-1);
		}
	}
	
	template <class Tad>
	void grad_eta(int ni, int nj, Tad rho[ni+1][nj+1], Tad u[ni+1][nj+1], Tad v[ni+1][nj+1], Tad p[ni+1][nj+1],	\
		Tad rlft[ni-1][nj], Tad ulft[ni-1][nj], Tad vlft[ni-1][nj], Tad plft[ni-1][nj], \
		Tad rrht[ni-1][nj], Tad urht[ni-1][nj], Tad vrht[ni-1][nj], Tad prht[ni-1][nj], int order = 2){
		
		if(order == 1){
			first_order_eta(ni, nj, rho, rlft, rrht);
			first_order_eta(ni, nj, u, ulft, urht);
			first_order_eta(ni, nj, v, vlft, vrht);
			first_order_eta(ni, nj, p, plft, prht);
		}
		else if(order == 2){
			second_order_eta(ni, nj, rho, rlft, rrht);
			second_order_eta(ni, nj, u, ulft, urht);
			second_order_eta(ni, nj, v, vlft, vrht);
			second_order_eta(ni, nj, p, plft, prht);
		}
		else{
			printf("order not implemented");
			exit(-1);
		}
	}
	
	
	
	
	#endif

