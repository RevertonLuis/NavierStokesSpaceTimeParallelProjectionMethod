MODULE Variaveis_Solvers_V

   USE ClassMatrizA

   IMPLICIT NONE

   !----------------------GERAL------------------------------------------------------
   INTEGER :: Nx, Ny, Tn
   DOUBLE PRECISION :: tol
   DOUBLE PRECISION :: norma_L2
   DOUBLE PRECISION :: Lx, Ly
   DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: normas_L2
   DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: u, f
   DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: bp
   INTEGER :: verbose
   !----------------------GERAL------------------------------------------------------
   
   
   !----------------------MULTIGRID GAUSS-SEIDEL-------------------------------------
   INTEGER :: ite_down, ite_up
   INTEGER :: max_ciclos_V
   INTEGER :: Levels
   INTEGER, DIMENSION(:, :), ALLOCATABLE :: lv_if
   INTEGER, DIMENSION(:),   ALLOCATABLE :: lv_nx, lv_ny
   !----------------------MULTIGRID GAUSS-SEIDEL-------------------------------------
 
   
   TYPE ( MatrizA ), DIMENSION(:), ALLOCATABLE :: A
      
END MODULE Variaveis_Solvers_V
