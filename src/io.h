	#ifndef __INCL_IO
	#define __INCL_IO

	void read_grid(struct Grid *grid);
	void write_solution(struct Grid *grid);
	void log_residual(int iter, real_t l2norm);
	#endif