  #include "include.h"
 PetscErrorCode petsc_set_matrix(struct PetscSolver *psolver, struct Grid *grid, int nnz, unsigned int *rind, unsigned int *cind, real_t *values){
  int ni = grid->ni;
  int nj = grid->nj;
  int i_idx, j_idx;
  int ii, jj;
  real_t cfl = 30.0;
  PetscScalar val;
  PetscErrorCode ierr;
  for(int i=0; i<nnz; i++){
    i_idx = rind[i];
    j_idx = cind[i];
    val = -values[i];

    if(i_idx == j_idx){
     div_t divresult;
     divresult = div (i_idx, (nj-1)*NQ);
     ii = divresult.quot;
     divresult = div (divresult.rem, NQ);
     jj = divresult.quot;
     val += grid->vol[INDEX(ii, jj, ni-1, nj-1)]/(grid->dt[INDEX(ii, jj, ni-1, nj-1)] * cfl);
   }
   ierr   = MatSetValue(psolver->A,i_idx,j_idx,val,INSERT_VALUES);CHKERRQ(ierr);
 }

 ierr = MatAssemblyBegin(psolver->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
 ierr = MatAssemblyEnd(psolver->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
 cout<<"Matrix set....\n";
 return ierr;
}
 
  PetscErrorCode petsc_set_vector(struct PetscSolver *psolver, struct Grid *grid, real_t *res){

    int ni = grid->ni;
    int nj = grid->nj;
    int nvar = (ni-1)*(nj-1)*NQ;

    PetscScalar val;
    PetscErrorCode ierr;
    
    for (int i=0; i<nvar; i++){
      val=res[i];
      ierr = VecSetValue(psolver->b, i, val, INSERT_VALUES);CHKERRQ(ierr);
    }

    cout<<"Vector set....\n";
    return ierr;
  }

  PetscErrorCode petsc_solve(struct PetscSolver *psolver, struct Grid *grid, real_t **dq){
   cout<<"Now solving....\n";
    PetscErrorCode ierr;
    ierr = KSPSolve(psolver->ksp, psolver->b, psolver->x);CHKERRQ(ierr);
    ierr = VecGetArray(psolver->x, dq);CHKERRQ(ierr);
    cout<<"Solver iteration complete....\n";
    return ierr;
  }
  
