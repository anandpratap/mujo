	#ifndef __INCL__PETSC
	#define __INCL__PETSC
	#include "include.h"
	
	#define CONFIG_PETSC_TOL 1e-12
	#define CONFIG_PETSC_MAXITER 1000
	
	struct PetscSolver{
		KSP ksp;
		Vec x, b;
		Mat A;
		PC  pc;
		/*PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,600,600,&viewer);
		PetscObjectSetName((PetscObject)viewer,"Line graph Plot");
		PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
		VecView(x,viewer);
*/
	public:
		PetscSolver(int argc, char **argv, int nvar){
			PetscErrorCode ierr;
			ierr = init_petsc(argc, argv, nvar);
			
		}		

		~PetscSolver(){
			VecDestroy(&b);
			VecDestroy(&x);
			MatDestroy(&A);
			PetscFinalize();
		}
		PetscErrorCode init_petsc(int argc, char **argv, int nvar){
			static char help[] = "";
			PetscReal norm, tol=1.e-6; 
			PetscErrorCode ierr;
			PetscInt var, n=nvar;
			PetscMPIInt size;
			PetscScalar one = 1.0;
			PetscBool nonzeroguess = PETSC_FALSE;
			PetscInitialize(&argc,&argv,(char*)0,help);
			ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
			cout<<"size: "<<size<<"\n";
			if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
			ierr = PetscOptionsGetInt(NULL,"-n",&n,NULL);CHKERRQ(ierr);
			ierr = PetscOptionsGetBool(NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(ierr);


			ierr = VecCreate(PETSC_COMM_WORLD, &x);CHKERRQ(ierr);
			ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
			ierr = VecSetSizes(x, PETSC_DECIDE, nvar);CHKERRQ(ierr);
			ierr = VecSetFromOptions(x);CHKERRQ(ierr);
			ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
			
			MatCreateSeqAIJ(PETSC_COMM_WORLD, nvar, nvar, 36, NULL, &A);
      

			ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
			ierr = KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN);CHKERRQ(ierr);
			ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
			ierr = PCSetType(pc, PCJACOBI);CHKERRQ(ierr);
			ierr = KSPSetTolerances(ksp, CONFIG_PETSC_TOL, PETSC_DEFAULT, PETSC_DEFAULT, CONFIG_PETSC_MAXITER);CHKERRQ(ierr);
			ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
			return ierr;
		}
	};


	PetscErrorCode petsc_set_matrix(struct PetscSolver *psolver, struct Grid *grid, int nnz, unsigned int *rind, unsigned int *cind, real_t *values);
	PetscErrorCode petsc_set_vector(struct PetscSolver *psolver, struct Grid *grid, real_t *res);
	PetscErrorCode petsc_solve(struct PetscSolver *psolver, struct Grid *grid, real_t **dq);

	#endif