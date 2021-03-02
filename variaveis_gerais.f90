MODULE VariaveisGerais

    USE ClassVelocidadeU
    USE ClassVelocidadeV
    USE CLassPressao
    
    IMPLICIT NONE
   
    INTEGER :: Nx, Ny, Nt, Experimento, Tn
    DOUBLE PRECISION :: Lx, Ly, Tf, Ti, hx, hy, ht, Re, Qui, CFL, Xf, Xi, Yf, Yi
   
    INTEGER :: Metodo_P, Metodo_U, Metodo_V, Extrapolacao_U, Extrapolacao_V
    DOUBLE PRECISION :: Tolerancia_P, Tolerancia_U, Tolerancia_V
   
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xu, xv, yu, yv
   
    ! Variaveis relacionadas ao paralelismo
    INTEGER :: Threads
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: threads_if

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: normas_L2_u, normas_L2_v, normas_L2_p 


    ! Descricao dos modelos
    CHARACTER(100), DIMENSION(6) :: Metodos
    CHARACTER(100), DIMENSION(6) :: Metodos_Siglas
   
    TYPE( VelocidadeU )      :: u
    TYPE( VelocidadeV )      :: v
    TYPE( Pressao )          :: p

    INTEGER :: iteracoes_externas, ite_sol

    ! Variaveis relacionadas
    ! escrita/leitura e arquivos de
    ! entrada/saida do programa
    INTEGER :: unit_arq_geral
    CHARACTER(LEN=100) :: arq_geral
   
END MODULE VariaveisGerais
