	#ifndef _INCL_FLUX
	#define _INCL_FLUX
	#include "include.h"

	template <class T, class Tad>
	void roeflux(T nx, T ny,		\
		Tad rlft, Tad ulft, Tad vlft, Tad plft, \
		Tad rrht, Tad urht, Tad vrht, Tad prht, \
		Tad f[4]){
		T gm1 = GAMMA - 1.0;
		Tad ogm1 = 1/gm1;
		
		Tad rlfti = 1/rlft;
		Tad rulft = rlft*ulft;
		Tad rvlft = rlft*vlft;
		Tad uvl = 0.5f*(ulft*ulft + vlft*vlft);
		Tad elft = plft*ogm1 + rlft*uvl;
		Tad hlft = (elft + plft)*rlfti;
		
		
		Tad rrhti = 1/rrht;
		Tad rurht = rrht*urht;
		Tad rvrht = rrht*vrht;
		Tad uvr = 0.5f*(urht*urht + vrht*vrht);
		Tad erht = prht*ogm1 + rrht*uvr;
		Tad hrht = (erht + prht)*rrhti;
		
		Tad rat = sqrt(rrht*rlfti);
		Tad rati = 1/(rat+1);
		Tad rav = rat*rlft;
		Tad uav = (rat*urht + ulft)*rati;
		Tad vav = (rat*vrht + vlft)*rati;
		Tad hav = (rat*hrht + hlft)*rati;
		Tad uv = 0.5f*(uav*uav + vav*vav);
		Tad cav = sqrt(gm1*(hav - uv));
		
		Tad aq1 = rrht - rlft;
		Tad aq2 = urht - ulft;
		Tad aq3 = vrht - vlft;
		Tad aq4 = prht - plft;
		
		Tad dr = sqrt(nx*nx + ny*ny);
		Tad r1 = nx/dr;
		Tad r2 = ny/dr;
		

		Tad uu = r1*uav + r2*vav;
		Tad c2 = cav*cav;
		Tad c2i = 1/c2;
		Tad auu = fabs(uu);
		Tad aupc = fabs(uu + cav);
		Tad aumc = fabs(uu - cav);

		Tad uulft = r1*ulft + r2*vlft;
		Tad uurht = r1*urht + r2*vrht;
		Tad rcav = rav*cav;
		Tad aquu = uurht - uulft;
		Tad c2ih = 0.5f*c2i;
		Tad ruuav = auu*rav;
		
		Tad b1, b2, b3, b4, b5, b6, b7;
		b1 = auu*(aq1 - c2i*aq4);
		b2 = c2ih*aupc*(aq4 + rcav*aquu);
		b3 = c2ih*aumc*(aq4 - rcav*aquu);
		b4 = b1 + b2 + b3;
		b5 = cav*(b2 - b3);
		b6 = ruuav*(aq2 - r1*aquu);
		b7 = ruuav*(aq3 - r2*aquu);
		
		aq1 = b4;
		aq2 = uav*b4 + r1*b5 + b6;
		aq3 = vav*b4 + r2*b5 + b7;
		aq4 = hav*b4 + uu*b5 + uav*b6 + vav*b7 - c2*b1*ogm1;
		
		Tad aj = 0.5f*dr;
		Tad plar = plft + prht;
		Tad eplft = elft + plft;
		Tad eprht = erht + prht;
		f[0] = aj*(rlft*uulft + rrht*uurht - aq1);
		f[1] = aj*(rulft*uulft + rurht*uurht + r1*plar - aq2);
		f[2] = aj*(rvlft*uulft + rvrht*uurht + r2*plar - aq3);
		f[3] = aj*(eplft*uulft + eprht*uurht - aq4);
		
	}
	
	
#endif
