MODULE ExtrapolacoesDeUeV

   CONTAINS

   
      DOUBLE PRECISION FUNCTION Extrapola_u_S_Linear(t, pu, i, Nx, u )
      
         USE FuncoesAlias
         
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, pu, i, Nx
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
         
         Extrapola_u_S_Linear = 2.0d0*cc_u_s(t, i) - u(pu)
         
      END FUNCTION Extrapola_u_S_Linear

               
      DOUBLE PRECISION FUNCTION Extrapola_u_S_Cubica(t, pu, i, Nx, u )
      
         USE FuncoesAlias
         
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, pu, i, Nx
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
         
         Extrapola_u_S_Cubica = (16.0d0/5.0d0)*cc_u_s(t, i) - 3.0d0*u(pu) + u(pu+Nx) - (1.0d0/5.0d0)*u(pu+2*Nx) 
         
      END FUNCTION Extrapola_u_S_Cubica


      DOUBLE PRECISION FUNCTION Extrapola_u_N_Linear(t, pu, i, Nx, u )
      
         USE FuncoesAlias
         
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, pu, i, Nx
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
         
         Extrapola_u_N_Linear = 2.0d0*cc_u_n(t, i) - u(pu)
         
      END FUNCTION Extrapola_u_N_Linear

               
      DOUBLE PRECISION FUNCTION Extrapola_u_N_Cubica(t, pu, i, Nx, u )
      
         USE FuncoesAlias
         
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, pu, i, Nx
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
         
         Extrapola_u_N_Cubica = (16.0d0/5.0d0)*cc_u_N(t, i) - 3.0d0*u(pu) + u(pu-Nx) - (1.0d0/5.0d0)*u(pu-2*Nx) 
         
      END FUNCTION Extrapola_u_N_Cubica
      
      
      
      DOUBLE PRECISION FUNCTION Extrapola_v_W_Linear(t, pv, j, Nx, v )
      
         USE FuncoesAlias
         
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, pv, j, Nx
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: v
         
         Extrapola_v_W_Linear = 2.0d0*cc_v_w(t, j) - v(pv)
         
      END FUNCTION Extrapola_v_W_Linear

               
      DOUBLE PRECISION FUNCTION Extrapola_v_W_Cubica(t, pv, j, Nx, v )
      
         USE FuncoesAlias
         
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, pv, j, Nx
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: v
         
         Extrapola_v_W_Cubica = (16.0d0/5.0d0)*cc_v_w(t, j) - 3.0d0*v(pv) + v(pv+1) - (1.0d0/5.0d0)*v(pv+2) 
         
      END FUNCTION Extrapola_v_W_Cubica
      
      
      DOUBLE PRECISION FUNCTION Extrapola_v_E_Linear(t, pv, j, Nx, v )
      
         USE FuncoesAlias
         
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, pv, j, Nx
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: v
         
         Extrapola_v_E_Linear = 2.0d0*cc_v_e(t, j) - v(pv)
         
      END FUNCTION Extrapola_v_E_Linear

               
      DOUBLE PRECISION FUNCTION Extrapola_v_E_Cubica(t, pv, j, Nx, v )
      
         USE FuncoesAlias
         
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: t, pv, j, Nx
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: v
         
         Extrapola_v_E_Cubica = (16.0d0/5.0d0)*cc_v_e(t, j) - 3.0d0*v(pv) + v(pv-1) - (1.0d0/5.0d0)*v(pv-2) 
         
      END FUNCTION Extrapola_v_E_Cubica
      
      
END MODULE ExtrapolacoesDeUeV
