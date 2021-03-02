MODULE FontesSubrotinas

   CONTAINS

   
                  
      SUBROUTINE ComplementosU_SW(t, t0, pu, Nx, dx, dy, dt, Re, u, u0, v, v0, Dij, Aij, Aij0, Cij, pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pu, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dij, Aij, Aij0, Cij
         
         
         INTEGER :: pv, i, j
         DOUBLE PRECISION :: u_extrapolado, u_extrapolado0 
         
         j = (pu-1)/Nx + 1
         i = pu - (j-1)*Nx
         pv = pu + j - 1
         pp = pv
            
         
         u_extrapolado  = extrapola_u_S(t,  pu, i, Nx, u)
         u_extrapolado0 = extrapola_u_S(t0, pu, i, Nx, u0)
         
         ! Difusao
         Dij = ( u(pu+1) - 2.0d0*u(pu) + cc_u_w(t,j) )/(dx**2) + ( u(pu+Nx) - 2.0d0*u(pu) + u_extrapolado )/(dy**2) 
         
         ! Adveccao: Aij = u*du/dx + v*du/dy
         Aij  = u (pu)*(u (pu+1) - cc_u_w(t ,j) )/(2.0d0*dx) + ((v (pv) + v (pv+1) + cc_v_s(t , i) + cc_v_s(t, i+1))/4.0d0)*(u (pu+Nx) - u_extrapolado )/(2.0d0*dy)
         Aij0 = u0(pu)*(u0(pu+1) - cc_u_w(t0,j) )/(2.0d0*dx) + ((v0(pv) + v0(pv+1) + cc_v_s(t0, i) + cc_v_s(t0,i+1))/4.0d0)*(u0(pu+Nx) - u_extrapolado0)/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de U"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao U entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_u)?
         Cij = (dt/(2.0d0*Re))*( cc_u_w(t+1,j)/(dx**2) + (16.0d0/5.0)*(cc_u_s(t+1,i)/(dy**2)) )
         
      END SUBROUTINE ComplementosU_SW



      SUBROUTINE ComplementosU_S(t, t0, pu, Nx, dx, dy, dt, Re, u, u0, v, v0, Dij, Aij, Aij0, Cij, pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pu, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dij, Aij, Aij0, Cij 
         
         INTEGER :: pv, i, j
         DOUBLE PRECISION :: u_extrapolado, u_extrapolado0
                  
         j = (pu-1)/Nx + 1
         i = pu - (j-1)*Nx
         pv = pu + j - 1
         pp = pv
            
         u_extrapolado  = extrapola_u_S(t,  pu, i, Nx, u)
         u_extrapolado0 = extrapola_u_S(t0, pu, i, Nx, u0)
         
         ! Difusao
         Dij = ( u(pu+1) - 2.0d0*u(pu) + u(pu-1) )/(dx**2) + ( u(pu+Nx) - 2.0d0*u(pu) + u_extrapolado )/(dy**2) 
         
         ! Adveccao: Aij = u*du/dx + v*du/dy
         Aij  = u (pu)*(u (pu+1) - u (pu-1) )/(2.0d0*dx) + ((v (pv) + v (pv+1) + cc_v_s(t,  i) + cc_v_s(t,  i+1))/4.0d0)*(u (pu+Nx) - u_extrapolado )/(2.0d0*dy)
         Aij0 = u0(pu)*(u0(pu+1) - u0(pu-1) )/(2.0d0*dx) + ((v0(pv) + v0(pv+1) + cc_v_s(t0, i) + cc_v_s(t0, i+1))/4.0d0)*(u0(pu+Nx) - u_extrapolado0)/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de U"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao U entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_u)?
         Cij = (dt/(2.0d0*Re))*( (16.0d0/5.0)*(cc_u_s(t+1,i)/(dy**2)) )
         
         
      END SUBROUTINE ComplementosU_S
      
      
      SUBROUTINE ComplementosU_SE(t, t0, pu, Nx, dx, dy, dt, Re, u, u0, v, v0, Dij, Aij, Aij0, Cij, pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pu, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dij, Aij, Aij0, Cij 
         
         INTEGER :: pv, i, j
         DOUBLE PRECISION :: u_extrapolado, u_extrapolado0
         
         j = (pu-1)/Nx + 1
         i = pu - (j-1)*Nx
         pv = pu + j - 1
         pp = pv
            
         ! Usando interpolacao cubica
         u_extrapolado  = extrapola_u_S(t,  pu, i, Nx, u)
         u_extrapolado0 = extrapola_u_S(t0, pu, i, Nx, u0)
         
         ! Difusao
         Dij = ( cc_u_e(t,j) - 2.0d0*u(pu) + u(pu-1) )/(dx**2) + ( u(pu+Nx) - 2.0d0*u(pu) + u_extrapolado )/(dy**2) 
         
         ! Adveccao: Aij = u*du/dx + v*du/dy
         Aij  = u (pu)*(cc_u_e(t, j) - u (pu-1) )/(2.0d0*dx) + ((v (pv) + v (pv+1) + cc_v_s(t, i) + cc_v_s(t, i+1))/4.0d0)*(u (pu+Nx) - u_extrapolado )/(2.0d0*dy)
         Aij0 = u0(pu)*(cc_u_e(t0,j) - u0(pu-1) )/(2.0d0*dx) + ((v0(pv) + v0(pv+1) + cc_v_s(t0,i) + cc_v_s(t0,i+1))/4.0d0)*(u0(pu+Nx) - u_extrapolado0)/(2.0d0*dy)
         
         
         ! Termo oriundo do "lado esquerdo da equacao de U"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao U entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_u)?
         Cij = (dt/(2.0d0*Re))*( cc_u_e(t+1,j)/(dx**2) + (16.0d0/5.0)*(cc_u_s(t+1,i)/(dy**2)) )
         
      END SUBROUTINE ComplementosU_SE
      
      
      SUBROUTINE ComplementosU_W(t, t0, pu, Nx, dx, dy, dt, Re, u, u0, v, v0, Dij, Aij, Aij0, Cij, pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pu, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dij, Aij, Aij0, Cij 
         
         INTEGER :: pv, i, j
         
         j = (pu-1)/Nx + 1
         i = pu - (j-1)*Nx
         pv = pu + j - 1
         pp = pv
            
         ! Difusao
         Dij = ( u(pu+1) - 2.0d0*u(pu) + cc_u_w(t,j) )/(dx**2) + ( u(pu+Nx) - 2.0d0*u(pu) + u(pu-Nx) )/(dy**2) 
         
         ! Adveccao: Aij = u*du/dx + v*du/dy
         Aij  = u (pu)*(u (pu+1) - cc_u_w(t,  j) )/(2.0d0*dx) + ((v (pv) + v (pv+1) + v (pv-(Nx+1)) + v (pv+1-(Nx+1)))/4.0d0)*(u (pu+Nx) - u (pu-Nx))/(2.0d0*dy)
         Aij0 = u0(pu)*(u0(pu+1) - cc_u_w(t0, j) )/(2.0d0*dx) + ((v0(pv) + v0(pv+1) + v0(pv-(Nx+1)) + v0(pv+1-(Nx+1)))/4.0d0)*(u0(pu+Nx) - u0(pu-Nx))/(2.0d0*dy)
         
         
         ! Termo oriundo do "lado esquerdo da equacao de U"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao U entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_u)?
         Cij = (dt/(2.0d0*Re))*( cc_u_w(t+1,j)/(dx**2) )
         
      END SUBROUTINE ComplementosU_W
      
      
      SUBROUTINE ComplementosU_C(t, t0, pu, Nx, dx, dy, dt, Re, u, u0, v, v0, Dij, Aij, Aij0, Cij, pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pu, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dij, Aij, Aij0, Cij 
         
         INTEGER :: pv, i, j
         
         j = (pu-1)/Nx + 1
         i = pu - (j-1)*Nx
         pv = pu + j - 1
         pp = pv
            
         ! Difusao
         Dij = ( u(pu+1) - 2.0d0*u(pu) + u(pu-1) )/(dx**2) + ( u(pu+Nx) - 2.0d0*u(pu) + u(pu-Nx) )/(dy**2) 
         
         ! Adveccao: Aij = u*du/dx + v*du/dy
         Aij  = u (pu)*(u (pu+1) - u (pu-1) )/(2.0d0*dx) + ((v (pv) + v (pv+1) + v (pv-(Nx+1)) + v (pv+1-(Nx+1)))/4.0d0)*(u (pu+Nx) - u (pu-Nx))/(2.0d0*dy)
         Aij0 = u0(pu)*(u0(pu+1) - u0(pu-1) )/(2.0d0*dx) + ((v0(pv) + v0(pv+1) + v0(pv-(Nx+1)) + v0(pv+1-(Nx+1)))/4.0d0)*(u0(pu+Nx) - u0(pu-Nx))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de U"
         Cij = 0.0d0
                  
      END SUBROUTINE ComplementosU_C

      
      SUBROUTINE ComplementosU_E(t, t0, pu, Nx, dx, dy, dt, Re, u, u0, v, v0, Dij, Aij, Aij0, Cij, pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pu, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dij, Aij, Aij0, Cij 
         
         INTEGER :: pv, i, j
         
         j = (pu-1)/Nx + 1
         i = pu - (j-1)*Nx
         pv = pu + j - 1
         pp = pv
            
         ! Difusao
         Dij = ( cc_u_e(t, j) - 2.0d0*u(pu) + u(pu-1) )/(dx**2) + ( u(pu+Nx) - 2.0d0*u(pu) + u(pu-Nx) )/(dy**2) 
         
         ! Adveccao: Aij = u*du/dx + v*du/dy
         Aij  = u (pu)*(cc_u_e(t,  j) - u (pu-1) )/(2.0d0*dx) + ((v (pv) + v (pv+1) + v (pv-(Nx+1)) + v (pv+1-(Nx+1)))/4.0d0)*(u (pu+Nx) - u (pu-Nx))/(2.0d0*dy)
         Aij0 = u0(pu)*(cc_u_e(t0, j) - u0(pu-1) )/(2.0d0*dx) + ((v0(pv) + v0(pv+1) + v0(pv-(Nx+1)) + v0(pv+1-(Nx+1)))/4.0d0)*(u0(pu+Nx) - u0(pu-Nx))/(2.0d0*dy)
         
                  
         ! Termo oriundo do "lado esquerdo da equacao de U"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao U entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_u)?
         Cij = (dt/(2.0d0*Re))*( cc_u_e(t+1,j)/(dx**2) )
         
      END SUBROUTINE ComplementosU_E
      
      
      SUBROUTINE ComplementosU_NW(t, t0, pu, Nx, dx, dy, dt, Re, u, u0, v, v0, Dij, Aij, Aij0, Cij, pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pu, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dij, Aij, Aij0, Cij 
         
         INTEGER :: pv, i, j
         DOUBLE PRECISION :: u_extrapolado, u_extrapolado0
         
         j = (pu-1)/Nx + 1
         i = pu - (j-1)*Nx
         pv = pu + j - 1
         pp = pv
            
         ! Usando interpolacao cubica
         u_extrapolado  = extrapola_u_N(t,  pu, i, Nx, u)
         u_extrapolado0 = extrapola_u_N(t0, pu, i, Nx, u0)
         
         ! Difusao
         Dij = ( u(pu+1) - 2.0d0*u(pu) + cc_u_w(t, j) )/(dx**2) + ( u_extrapolado - 2.0d0*u(pu) + u(pu-Nx) )/(dy**2) 
         
         ! Adveccao: Aij = u*du/dx + v*du/dy
         Aij  = u (pu)*(u (pu+1) - cc_u_w(t,  j) )/(2.0d0*dx) + ((cc_v_n(t, i) + cc_v_n(t ,i+1) + v (pv-(Nx+1)) + v (pv+1-(Nx+1)))/4.0d0)*(u_extrapolado  - u (pu-Nx))/(2.0d0*dy)
         Aij0 = u0(pu)*(u0(pu+1) - cc_u_w(t0, j) )/(2.0d0*dx) + ((cc_v_n(t0,i) + cc_v_n(t0,i+1) + v0(pv-(Nx+1)) + v0(pv+1-(Nx+1)))/4.0d0)*(u_extrapolado0 - u0(pu-Nx))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de U"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao U entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_u)?
         Cij = (dt/(2.0d0*Re))*( cc_u_w(t+1,j)/(dx**2) + (16.0d0/5.0)*(cc_u_n(t+1,i)/(dy**2)) )
         
      END SUBROUTINE ComplementosU_NW
      
      
      SUBROUTINE ComplementosU_N(t, t0, pu, Nx, dx, dy, dt, Re, u, u0, v, v0, Dij, Aij, Aij0, Cij, pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pu, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dij, Aij, Aij0, Cij 
         
         INTEGER :: pv, i, j
         DOUBLE PRECISION :: u_extrapolado, u_extrapolado0
         
         j = (pu-1)/Nx + 1
         i = pu - (j-1)*Nx
         pv = pu + j - 1
         pp = pv
         
         ! Usando interpolacao cubica
         u_extrapolado  = extrapola_u_N(t,  pu, i, Nx, u)
         u_extrapolado0 = extrapola_u_N(t0, pu, i, Nx, u0)
         
         ! Difusao
         Dij = ( u(pu+1) - 2.0d0*u(pu) + u(pu-1) )/(dx**2) + ( u_extrapolado - 2.0d0*u(pu) + u(pu-Nx) )/(dy**2) 
         
         ! Adveccao: Aij = u*du/dx + v*du/dy
         Aij  = u (pu)*(u (pu+1) - u (pu-1) )/(2.0d0*dx) + ((cc_v_n(t, i) + cc_v_n(t, i+1) + v (pv-(Nx+1)) + v (pv+1-(Nx+1)))/4.0d0)*(u_extrapolado  - u (pu-Nx))/(2.0d0*dy)
         Aij0 = u0(pu)*(u0(pu+1) - u0(pu-1) )/(2.0d0*dx) + ((cc_v_n(t0,i) + cc_v_n(t0,i+1) + v0(pv-(Nx+1)) + v0(pv+1-(Nx+1)))/4.0d0)*(u_extrapolado0 - u0(pu-Nx))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de U"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao U entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_u)?
         Cij = (dt/(2.0d0*Re))*( (16.0d0/5.0)*(cc_u_n(t+1,i)/(dy**2)) )
                  
      END SUBROUTINE ComplementosU_N
      
      
      SUBROUTINE ComplementosU_NE(t, t0, pu, Nx, dx, dy, dt, Re, u, u0, v, v0, Dij, Aij, Aij0, Cij, pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pu, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dij, Aij, Aij0, Cij 
         
         INTEGER :: pv, i, j
         DOUBLE PRECISION :: u_extrapolado, u_extrapolado0
         
         j = (pu-1)/Nx + 1
         i = pu - (j-1)*Nx
         pv = pu + j - 1
         pp = pv
            
         ! Usando interpolacao cubica
         u_extrapolado  = extrapola_u_N(t,  pu, i, Nx, u)
         u_extrapolado0 = extrapola_u_N(t0, pu, i, Nx, u0)
         
         ! Difusao
         Dij = ( cc_u_e(t,j) - 2.0d0*u(pu) + u(pu-1) )/(dx**2) + ( u_extrapolado - 2.0d0*u(pu) + u(pu-Nx) )/(dy**2) 
         
         ! Adveccao: Aij = u*du/dx + v*du/dy
         Aij  = u (pu)*(cc_u_e(t, j) - u (pu-1) )/(2.0d0*dx) + ((cc_v_n(t, i) + cc_v_n(t, i+1) + v (pv-(Nx+1)) + v (pv+1-(Nx+1)))/4.0d0)*(u_extrapolado  - u (pu-Nx))/(2.0d0*dy)
         Aij0 = u0(pu)*(cc_u_e(t0,j) - u0(pu-1) )/(2.0d0*dx) + ((cc_v_n(t0,i) + cc_v_n(t0,i+1) + v0(pv-(Nx+1)) + v0(pv+1-(Nx+1)))/4.0d0)*(u_extrapolado0 - u0(pu-Nx))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de U"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao U entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_u)?
         Cij = (dt/(2.0d0*Re))*( cc_u_e(t+1,j)/(dx**2) + (16.0d0/5.0)*(cc_u_n(t+1,i)/(dy**2)) )
         
      END SUBROUTINE ComplementosU_NE
      
      
      
      
!-------------------------------------------------------------------------------------------------------------------------------      
      
      
      
      
      SUBROUTINE ComplementosV_SW(t, t0, pv, Nx, dx, dy, dt, Re, u, u0, v, v0, Dvij, Avij, Avij0, Cij,  pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pv, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dvij, Avij, Avij0, Cij 
         
         
         INTEGER :: pu, i, j
         DOUBLE PRECISION :: v_extrapolado, v_extrapolado0
         
         j = (pv-1)/Nx + 1
         i = pv - (j-1)*Nx
         pu = pv - j + 1
         pp = pv
            
         ! Usando interpolacao cubica
         v_extrapolado  = extrapola_v_W(t,  pv, j, Nx, v)
         v_extrapolado0 = extrapola_v_W(t0, pv, j, Nx, v0)
         
         ! Difusao
         Dvij = ( v(pv+1) - 2.0d0*v(pv) + v_extrapolado )/(dx**2) + ( v(pv+Nx) - 2.0d0*v(pv) + cc_v_s(t,i) )/(dy**2) 
         
         ! Adveccao: Avij = u*dv/dx + v*dv/dy
         Avij  = ( (u (pu) + cc_u_w(t, j) + u (pu+(Nx-1)) + cc_u_w(t, j+1))/4.0d0)*(v (pv+1) - v_extrapolado )/(2.0d0*dx) + v (pv)*(v (pv+Nx) - cc_v_s(t ,i))/(2.0d0*dy)
         Avij0 = ( (u0(pu) + cc_u_w(t0,j) + u0(pu+(Nx-1)) + cc_u_w(t0,j+1))/4.0d0)*(v0(pv+1) - v_extrapolado0)/(2.0d0*dx) + v0(pv)*(v0(pv+Nx) - cc_v_s(t0,i))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de V"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao V entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_v)?
         Cij = (dt/(2.0d0*Re))*( (16.0d0/5.0)*cc_v_w(t+1,j)/(dx**2) + cc_v_s(t+1, i)/(dy**2) )   
         
         
      END SUBROUTINE ComplementosV_SW
      
      
      SUBROUTINE ComplementosV_S(t, t0, pv, Nx, dx, dy, dt, Re, u, u0, v, v0, Dvij, Avij, Avij0, Cij,  pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pv, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dvij, Avij, Avij0, Cij 
         
         
         INTEGER :: pu, i, j
         
         j = (pv-1)/Nx + 1
         i = pv - (j-1)*Nx
         pu = pv - j + 1
         pp = pv
            
         ! Difusao
         Dvij = ( v(pv+1) - 2.0d0*v(pv) + v(pv-1) )/(dx**2) + ( v(pv+Nx) - 2.0d0*v(pv) + cc_v_s(t,i) )/(dy**2) 
         
         ! Adveccao: Avij = u*dv/dx + v*dv/dy
         Avij  = ( (u (pu) + u (pu-1) + u (pu+(Nx-1)) + u (pu+(Nx-1)-1))/4.0d0)*(v (pv+1) - v (pv-1) )/(2.0d0*dx) + v (pv)*(v (pv+Nx) - cc_v_s(t ,i))/(2.0d0*dy)
         Avij0 = ( (u0(pu) + u0(pu-1) + u0(pu+(Nx-1)) + u0(pu+(Nx-1)-1))/4.0d0)*(v0(pv+1) - v0(pv-1) )/(2.0d0*dx) + v0(pv)*(v0(pv+Nx) - cc_v_s(t0,i))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de V"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao V entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_v)?
         Cij = (dt/(2.0d0*Re))*( cc_v_s(t+1, i)/(dy**2) )   
                  
         
      END SUBROUTINE ComplementosV_S
      
      
      SUBROUTINE ComplementosV_SE(t, t0, pv, Nx, dx, dy, dt, Re, u, u0, v, v0, Dvij, Avij, Avij0, Cij,  pp )
      
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pv, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dvij, Avij, Avij0, Cij 
         
         
         INTEGER :: pu, i, j
         DOUBLE PRECISION :: v_extrapolado, v_extrapolado0 
         
         j = (pv-1)/Nx + 1
         i = pv - (j-1)*Nx
         pu = pv - j + 1
         pp = pv
            
         v_extrapolado  = extrapola_v_E(t,  pv, j, Nx, v)
         v_extrapolado0 = extrapola_v_E(t0, pv, j, Nx, v0)
         
         ! Difusao
         Dvij = ( v_extrapolado - 2.0d0*v(pv) + v(pv-1) )/(dx**2) + ( v(pv+Nx) - 2.0d0*v(pv) + cc_v_s(t,i) )/(dy**2) 
         
         ! Adveccao: Avij = u*dv/dx + v*dv/dy
         Avij  = ( (cc_u_e(t, j) + u (pu-1) + cc_u_e(t, j+1) + u (pu+(Nx-1)-1))/4.0d0)*(v_extrapolado  - v (pv-1))/(2.0d0*dx) + v (pv)*(v (pv+Nx) - cc_v_s(t ,i))/(2.0d0*dy)
         Avij0 = ( (cc_u_e(t0,j) + u0(pu-1) + cc_u_e(t0,j+1) + u0(pu+(Nx-1)-1))/4.0d0)*(v_extrapolado0 - v0(pv-1))/(2.0d0*dx) + v0(pv)*(v0(pv+Nx) - cc_v_s(t0,i))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de V"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao V entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_v)?
         Cij = (dt/(2.0d0*Re))*( (16.0d0/5.0)*cc_v_e(t+1,j)/(dx**2) + cc_v_s(t+1, i)/(dy**2) )            
         
      END SUBROUTINE ComplementosV_SE
      
      
      SUBROUTINE ComplementosV_W(t, t0, pv, Nx, dx, dy, dt, Re, u, u0, v, v0, Dvij, Avij, Avij0, Cij,  pp )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pv, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dvij, Avij, Avij0, Cij 
         
         INTEGER :: pu, i, j
         DOUBLE PRECISION :: v_extrapolado, v_extrapolado0 
         
         j = (pv-1)/Nx + 1
         i = pv - (j-1)*Nx
         pu = pv - j + 1
         pp = pv
            
         
         ! Usando interpolacao cubica
         v_extrapolado  = extrapola_v_W(t,  pv, j, Nx, v)
         v_extrapolado0 = extrapola_v_W(t0, pv, j, Nx, v0)
         
         ! Difusao
         Dvij = ( v(pv+1) - 2.0d0*v(pv) + v_extrapolado )/(dx**2) + ( v(pv+Nx) - 2.0d0*v(pv) + v(pv-Nx) )/(dy**2) 
         
         ! Adveccao: Avij = u*dv/dx + v*dv/dy
         Avij  = ( (u (pu) + cc_u_w(t, j) + u (pu+(Nx-1)) + cc_u_w(t, j+1))/4.0d0)*(v (pv+1) - v_extrapolado )/(2.0d0*dx) + v (pv)*(v (pv+Nx) - v (pv-Nx))/(2.0d0*dy)
         Avij0 = ( (u0(pu) + cc_u_w(t0,j) + u0(pu+(Nx-1)) + cc_u_w(t0,j+1))/4.0d0)*(v0(pv+1) - v_extrapolado0)/(2.0d0*dx) + v0(pv)*(v0(pv+Nx) - v0(pv-Nx))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de V"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao V entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_v)?
         Cij = (dt/(2.0d0*Re))*( (16.0d0/5.0)*cc_v_w(t+1,j)/(dx**2) )   
                  
         
      END SUBROUTINE ComplementosV_W
      
      
      SUBROUTINE ComplementosV_C(t, t0, pv, Nx, dx, dy, dt, Re, u, u0, v, v0, Dvij, Avij, Avij0, Cij,  pp )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pv, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dvij, Avij, Avij0, Cij 
         
         INTEGER :: pu, i, j
         
         j = (pv-1)/Nx + 1
         i = pv - (j-1)*Nx
         pu = pv - j + 1
         pp = pv
            
         ! Difusao
         Dvij = ( v(pv+1) - 2.0d0*v(pv) + v(pv-1) )/(dx**2) + ( v(pv+Nx) - 2.0d0*v(pv) + v(pv-Nx) )/(dy**2) 
         
         ! Adveccao: Avij = u*dv/dx + v*dv/dy
         Avij  = ( (u (pu) + u (pu-1) + u (pu+(Nx-1)) + u (pu+(Nx-1)-1))/4.0d0)*(v (pv+1) - v (pv-1) )/(2.0d0*dx) + v (pv)*(v (pv+Nx) - v (pv-Nx))/(2.0d0*dy)
         Avij0 = ( (u0(pu) + u0(pu-1) + u0(pu+(Nx-1)) + u0(pu+(Nx-1)-1))/4.0d0)*(v0(pv+1) - v0(pv-1) )/(2.0d0*dx) + v0(pv)*(v0(pv+Nx) - v0(pv-Nx))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de V"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao V entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_v)?
         Cij = 0.0d0
         
      END SUBROUTINE ComplementosV_C
      
      
      SUBROUTINE ComplementosV_E(t, t0, pv, Nx, dx, dy, dt, Re, u, u0, v, v0, Dvij, Avij, Avij0, Cij,  pp )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pv, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dvij, Avij, Avij0, Cij 
         
         INTEGER :: pu, i, j
         DOUBLE PRECISION :: v_extrapolado, v_extrapolado0 
         
         j = (pv-1)/Nx + 1
         i = pv - (j-1)*Nx
         pu = pv - j + 1
         pp = pv
            
         
         ! Usando interpolacao cubica
         v_extrapolado  = extrapola_v_E(t,  pv, j, Nx, v)
         v_extrapolado0 = extrapola_v_E(t0, pv, j, Nx, v0)
         
         ! Difusao
         Dvij = ( v_extrapolado - 2.0d0*v(pv) + v(pv-1) )/(dx**2) + ( v(pv+Nx) - 2.0d0*v(pv) + v(pv-Nx) )/(dy**2) 
         
         ! Adveccao: Avij = u*dv/dx + v*dv/dy
         Avij  = ( (cc_u_e(t , j) + u (pu-1) + cc_u_e(t, j+1) + u (pu+(Nx-1)-1))/4.0d0)*(v_extrapolado  - v (pv-1) )/(2.0d0*dx) + v (pv)*(v (pv+Nx) - v (pv-Nx))/(2.0d0*dy)
         Avij0 = ( (cc_u_e(t0, j) + u0(pu-1) + cc_u_e(t0,j+1) + u0(pu+(Nx-1)-1))/4.0d0)*(v_extrapolado0 - v0(pv-1) )/(2.0d0*dx) + v0(pv)*(v0(pv+Nx) - v0(pv-Nx))/(2.0d0*dy)
                  
         ! Termo oriundo do "lado esquerdo da equacao de V"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao V entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_v)?
         Cij = (dt/(2.0d0*Re))*( (16.0d0/5.0)*cc_v_e(t+1,j)/(dx**2) )   
         
         
      END SUBROUTINE ComplementosV_E
      
      SUBROUTINE ComplementosV_NW(t, t0, pv, Nx, dx, dy, dt, Re, u, u0, v, v0, Dvij, Avij, Avij0, Cij,  pp )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pv, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dvij, Avij, Avij0, Cij 
         
         INTEGER :: pu, i, j
         DOUBLE PRECISION :: v_extrapolado, v_extrapolado0
         
         j = (pv-1)/Nx + 1
         i = pv - (j-1)*Nx
         pu = pv - j + 1
         pp = pv
            
         ! Usando interpolacao cubica
         v_extrapolado  = extrapola_v_W(t,  pv, j, Nx, v)
         v_extrapolado0 = extrapola_v_W(t0, pv, j, Nx, v0)
         
         ! Difusao
         Dvij = ( v(pv+1) - 2.0d0*v(pv) + v_extrapolado )/(dx**2) + ( cc_v_n(t,i) - 2.0d0*v(pv) + v(pv-Nx) )/(dy**2) 
         
         ! Adveccao: Avij = u*dv/dx + v*dv/dy
         Avij  = ( (u (pu) + cc_u_w(t, j) + u (pu+(Nx-1)) + cc_u_w(t, j+1))/4.0d0)*(v (pv+1) - v_extrapolado )/(2.0d0*dx) + v (pv)*(cc_v_n(t ,i) - v (pv-Nx))/(2.0d0*dy)
         Avij0 = ( (u0(pu) + cc_u_w(t0,j) + u0(pu+(Nx-1)) + cc_u_w(t0,j+1))/4.0d0)*(v0(pv+1) - v_extrapolado0)/(2.0d0*dx) + v0(pv)*(cc_v_n(t0,i) - v0(pv-Nx))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de V"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao V entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_v)?
         Cij = (dt/(2.0d0*Re))*( (16.0d0/5.0)*cc_v_w(t+1,j)/(dx**2) + cc_v_n(t+1, i)/(dy**2) )   
                  
         
      END SUBROUTINE ComplementosV_NW
      
      SUBROUTINE ComplementosV_N(t, t0, pv, Nx, dx, dy, dt, Re, u, u0, v, v0, Dvij, Avij, Avij0, Cij,  pp )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pv, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dvij, Avij, Avij0, Cij 
         
         
         INTEGER :: pu, i, j
         
         j = (pv-1)/Nx + 1
         i = pv - (j-1)*Nx
         pu = pv - j + 1
         pp = pv
            
         
         ! Difusao
         Dvij = ( v(pv+1) - 2.0d0*v(pv) + v(pv-1) )/(dx**2) + ( cc_v_n(t,i) - 2.0d0*v(pv) + v(pv-Nx) )/(dy**2) 
         
         ! Adveccao: Avij = u*dv/dx + v*dv/dy
         Avij  = ( (u (pu) + u (pu-1) + u (pu+(Nx-1)) + u (pu+(Nx-1)-1))/4.0d0)*(v (pv+1) - v (pv-1) )/(2.0d0*dx) + v (pv)*(cc_v_n(t ,i) - v (pv-Nx))/(2.0d0*dy)
         Avij0 = ( (u0(pu) + u0(pu-1) + u0(pu+(Nx-1)) + u0(pu+(Nx-1)-1))/4.0d0)*(v0(pv+1) - v0(pv-1) )/(2.0d0*dx) + v0(pv)*(cc_v_n(t0,i) - v0(pv-Nx))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de V"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao V entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_v)?
         Cij = (dt/(2.0d0*Re))*( cc_v_n(t+1, i)/(dy**2) )   
                  
         
      END SUBROUTINE ComplementosV_N
      
      SUBROUTINE ComplementosV_NE(t, t0, pv, Nx, dx, dy, dt, Re, u, u0, v, v0, Dvij, Avij, Avij0, Cij,  pp )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, t0, pv, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: Dvij, Avij, Avij0, Cij 
         
         INTEGER :: pu, i, j
         DOUBLE PRECISION :: v_extrapolado, v_extrapolado0
         
         j = (pv-1)/Nx + 1
         i = pv - (j-1)*Nx
         pu = pv - j + 1
         pp = pv
            
         ! Usando interpolacao cubica
         v_extrapolado  = extrapola_v_E(t,  pv, j, Nx, v)
         v_extrapolado0 = extrapola_v_E(t0, pv, j, Nx, v0)
         
         ! Difusao
         Dvij = ( v_extrapolado - 2.0d0*v(pv) + v(pv-1) )/(dx**2) + ( cc_v_n(t,i) - 2.0d0*v(pv) + v(pv-Nx) )/(dy**2) 
         
         ! Adveccao: Avij = u*dv/dx + v*dv/dy
         Avij  = ( (cc_u_e(t, j) + u (pu-1) + cc_u_e(t, j+1) + u (pu+(Nx-1)-1))/4.0d0)*(v_extrapolado  - v (pv-1))/(2.0d0*dx) + v (pv)*(cc_v_n(t ,i) - v (pv-Nx))/(2.0d0*dy)
         Avij0 = ( (cc_u_e(t0,j) + u0(pu-1) + cc_u_e(t0,j+1) + u0(pu+(Nx-1)-1))/4.0d0)*(v_extrapolado0 - v0(pv-1))/(2.0d0*dx) + v0(pv)*(cc_v_n(t0,i) - v0(pv-Nx))/(2.0d0*dy)
         
         ! Termo oriundo do "lado esquerdo da equacao de V"
         ! OBS: Como o complemento Cij veio do lado esquerdo da equacao V entao ele representa a solucao no passo de tempo + 1
         ! OBS: como o complemento Cij depende da extrapolacao utilizada mover para as funcoes de extrapolacao (extrapola_v)?
         Cij = (dt/(2.0d0*Re))*( (16.0d0/5.0)*cc_v_e(t+1,j)/(dx**2) + cc_v_n(t+1, i)/(dy**2) )                     
         
      END SUBROUTINE ComplementosV_NE
      
      
!-------------------------------------------------------------------------------------------------------------------------------      

     DOUBLE PRECISION FUNCTION StencilP_SW(t, p, Nx, dx, dy, dt, ut, vt )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, p, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ut, vt
         
         INTEGER :: pv, pu, i, j
         
         j = (p-1)/Nx + 1
         i = p - (j-1)*Nx
         pu = p - j + 1
         pv = p
             
         StencilP_SW = -(1.0d0/dt)*((ut(pu) - cc_u_w(t+1, j))/dx + (vt(pv) - cc_v_s(t+1, i))/dy)  
                  
      END FUNCTION StencilP_SW
      

      DOUBLE PRECISION FUNCTION StencilP_S(t, p, Nx, dx, dy, dt, ut, vt )

         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, p, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ut, vt
         
         INTEGER :: pv, pu, i, j
         
         j = (p-1)/Nx + 1
         i = p - (j-1)*Nx
         pu = p - j + 1
         pv = p

         StencilP_S = -(1.0d0/dt)*((ut(pu) - ut(pu-1))/dx + (vt(pv) - cc_v_s(t+1, i))/dy)
         
      END FUNCTION StencilP_S


      DOUBLE PRECISION FUNCTION StencilP_SE(t, p, Nx, dx, dy, dt, ut, vt )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, p, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ut, vt
         
         INTEGER :: pv, pu, i, j
         
         j = (p-1)/Nx + 1
         i = p - (j-1)*Nx
         pu = p - j + 1
         pv = p
             
         StencilP_SE = -(1.0d0/dt)*((cc_u_e(t+1, j) - ut(pu-1))/dx + (vt(pv) - cc_v_s(t+1, i))/dy)  
                  
      END FUNCTION StencilP_SE


      DOUBLE PRECISION FUNCTION StencilP_W(t, p, Nx, dx, dy, dt, ut, vt )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, p, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ut, vt
         
         INTEGER :: pv, pu, i, j
         
         j = (p-1)/Nx + 1
         i = p - (j-1)*Nx
         pu = p - j + 1
         pv = p
             
         StencilP_W = -(1.0d0/dt)*((ut(pu) - cc_u_w(t+1, j))/dx + (vt(pv) - vt(pv-Nx))/dy)  
                  
      END FUNCTION StencilP_W
      

      DOUBLE PRECISION FUNCTION StencilP_C(t, p, Nx, dx, dy, dt, ut, vt )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, p, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ut, vt
         
         INTEGER :: pv, pu, i, j
         
         j = (p-1)/Nx + 1
         i = p - (j-1)*Nx
         pu = p - j + 1
         pv = p
             
         StencilP_C = -(1.0d0/dt)*((ut(pu) - ut(pu-1))/dx + (vt(pv) - vt(pv-Nx))/dy)  
                  
      END FUNCTION StencilP_C 

      
      DOUBLE PRECISION FUNCTION StencilP_E(t, p, Nx, dx, dy, dt, ut, vt )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, p, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ut, vt
         
         INTEGER :: pv, pu, i, j
         
         j = (p-1)/Nx + 1
         i = p - (j-1)*Nx
         pu = p - j + 1
         pv = p
             
         StencilP_E = -(1.0d0/dt)*((cc_u_e(t+1,j) - ut(pu-1))/dx + (vt(pv) - vt(pv-Nx))/dy)  
                  
      END FUNCTION StencilP_E
 
      
      DOUBLE PRECISION FUNCTION StencilP_NW(t, p, Nx, dx, dy, dt, ut, vt )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, p, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ut, vt
         
         INTEGER :: pv, pu, i, j
         
         j = (p-1)/Nx + 1
         i = p - (j-1)*Nx
         pu = p - j + 1
         pv = p
             
         StencilP_NW = -(1.0d0/dt)*((ut(pu) - cc_u_w(t+1, j))/dx + (cc_v_n(t+1, i) - vt(pv-Nx))/dy)  
                  
      END FUNCTION StencilP_NW



      DOUBLE PRECISION FUNCTION StencilP_N(t, p, Nx, dx, dy, dt, ut, vt )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, p, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ut, vt
         
         INTEGER :: pv, pu, i, j
         
         j = (p-1)/Nx + 1
         i = p - (j-1)*Nx
         pu = p - j + 1
         pv = p
             
         StencilP_N = -(1.0d0/dt)*((ut(pu) - ut(pu-1))/dx + (cc_v_n(t+1, i) - vt(pv-Nx))/dy)  
                  
      END FUNCTION StencilP_N
    

      DOUBLE PRECISION FUNCTION StencilP_NE(t, p, Nx, dx, dy, dt, ut, vt )
         USE FuncoesAlias
      
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, p, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ut, vt
         
         INTEGER :: pv, pu, i, j
         
         j = (p-1)/Nx + 1
         i = p - (j-1)*Nx
         pu = p - j + 1
         pv = p
             
         StencilP_NE = -(1.0d0/dt)*((cc_u_e(t+1, j) - ut(pu-1))/dx + (cc_v_n(t+1, i) - vt(pv-Nx))/dy)  
                  
      END FUNCTION StencilP_NE

END MODULE FontesSubrotinas
