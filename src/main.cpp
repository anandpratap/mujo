  #include "include.h"


  int main(int argc, char **argv){
    int ni = 281;
    int nj = 51;
    struct Grid grid(ni, nj);
    read_grid(&grid);
    metrics(&grid);
    initialize(&grid);
    write_solution(&grid);
    double *res = new real_t[(ni-1)*(nj-1)*NQ];
    double *dq = NULL;
    int nvar = (ni-1)*(nj-1)*NQ;
    int iter = 0;

    PetscSolver psolver(argc, argv, nvar);

    for(int t=0; t<100; t++){
     residual(&grid, res);
     unsigned int *rind  = NULL;    
     unsigned int *cind  = NULL;        
     double *values = NULL;       
     int nnz;
     int options[4];
     PetscErrorCode ierr;
  
     options[0] = 0;      
     options[1] = 0;      
     options[2] = 0;      
     options[3] = 1;
     real_t *q = grid.q;
     sparse_jac(1, nvar, nvar, 0, q, &nnz, &rind, &cind, &values, options);
     ierr = petsc_set_matrix(&psolver, &grid, nnz, rind, cind, values);CHKERRQ(ierr);
     ierr = petsc_set_vector(&psolver, &grid, res);CHKERRQ(ierr);
     ierr = petsc_solve(&psolver, &grid, &dq);CHKERRQ(ierr);
     for(int i=0; i<nvar; i++){
      grid.q[i] = grid.q[i] + dq[i];
    }
    free(rind); rind=NULL;
    free(cind); cind=NULL;
    free(values); values=NULL;

    double l2norm;
    calc_l2norm(ni, nj, res, &l2norm);
    cout<<"l2norm: "<<l2norm<<" iter:"<<iter<<"\n";
    write_solution(&grid);
    log_residual(iter, l2norm);
    iter += 1;
}
delete[] res;
return 0;
}
