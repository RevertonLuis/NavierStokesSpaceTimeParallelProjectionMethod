MODULE ExperimentosNumericos

   ! Modulo que implementa os problemas com Solucao (ou benchmarks)
   
   IMPLICIT  NONE

   CONTAINS
   !-------------------- FUNCOES DO EXPERIMENTO 5 -----------------------
   DOUBLE PRECISION FUNCTION CC_U_W_experimento5(t, j) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, j
    
      valor = -DCOS(Xi)*DSIN( yu(j) )*DEXP(-2.0d0*(ht*(t-1) + Ti))
         
   END FUNCTION CC_U_W_experimento5
      
   
   DOUBLE PRECISION FUNCTION CC_U_E_experimento5(t, j) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, j
    
      valor = -DCOS(Xf)*DSIN( yu(j) )*DEXP(-2.0d0*(ht*(t-1) + Ti))
         
   END FUNCTION CC_U_E_experimento5
      

   DOUBLE PRECISION FUNCTION CC_U_S_experimento5(t, i) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, i
    
      valor = -DCOS( xu(i) )*DSIN( Yi )*DEXP(-2.0d0*(ht*(t-1)+Ti))
         
   END FUNCTION CC_U_S_experimento5
      

   DOUBLE PRECISION FUNCTION CC_U_N_experimento5(t, i) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, i
    
      valor = -DCOS( xu(i) )*DSIN( Yf )*DEXP(-2.0d0*(ht*(t-1)+Ti))
         
   END FUNCTION CC_U_N_experimento5

   DOUBLE PRECISION FUNCTION CC_V_W_experimento5(t, j) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, j
    
      valor = DSIN(Xi)*DCOS( yv(j) )*DEXP(-2.0d0*(ht*(t-1)+Ti))
         
   END FUNCTION CC_V_W_experimento5


   DOUBLE PRECISION FUNCTION CC_V_E_experimento5(t, j) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, j
    
      valor = DSIN(Xf)*DCOS( yv(j) )*DEXP(-2.0d0*(ht*(t-1)+Ti))
         
   END FUNCTION CC_V_E_experimento5

   
   DOUBLE PRECISION FUNCTION CC_V_S_experimento5(t, i) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, i
    
      valor = DSIN( xv(i) )*DCOS(Yi)*DEXP(-2.0d0*(ht*(t-1)+Ti))
         
   END FUNCTION CC_V_S_experimento5


   DOUBLE PRECISION FUNCTION CC_V_N_experimento5(t, i) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, i
    
      valor = DSIN( xv(i) )*DCOS(Yf)*DEXP(-2.0d0*(ht*(t-1)+Ti))
         
   END FUNCTION CC_V_N_experimento5
    
   DOUBLE PRECISION FUNCTION PressaoAnaliticaExperimento5(t, volume) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, volume
      INTEGER :: i, j
      
      j = FLOOR ( DBLE(volume-1)/Nx ) + 1
      i = volume - (j - 1)*Nx
      valor = -(1.0d0/(4.0d0*Re))*( DCOS(2.0d0*xv(i)) + DCOS(2.0d0*yu(j)) )*DEXP(-4.0d0 * (ht*(t-1.5d0) + Ti) )   
      
   END FUNCTION PressaoAnaliticaExperimento5
   
   DOUBLE PRECISION FUNCTION PressaoInicialExperimento5(t, volume) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, volume
      INTEGER :: i, j
      
      j = FLOOR ( DBLE(volume-1)/Nx ) + 1
      i = volume - (j - 1)*Nx
      valor = -(1.0d0/(4.0d0*Re))*( DCOS(2.0d0*xv(i)) + DCOS(2.0d0*yu(j)) )
      
   END FUNCTION PressaoInicialExperimento5

   DOUBLE PRECISION FUNCTION VelocidadeU_InicialExperimento5(t, volume) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, volume
      INTEGER :: i, j
      
      j = FLOOR ( DBLE(volume-1)/(Nx-1) ) + 1
      i = volume - (j - 1)*(Nx-1)
      valor = -DCOS(xu(i))*DSIN(yu(j))
      
   END FUNCTION VelocidadeU_InicialExperimento5
   
   DOUBLE PRECISION FUNCTION VelocidadeV_InicialExperimento5(t, volume) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, volume
      INTEGER :: i, j
      
      j = FLOOR ( DBLE(volume-1)/Nx ) + 1
      i = volume - (j - 1)*Nx
      valor = DSIN(xv(i))*DCOS(yv(j)) 
      
   END FUNCTION VelocidadeV_InicialExperimento5

   DOUBLE PRECISION FUNCTION VelocidadeU_AnaliticaExperimento5(t, volume) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, volume
      INTEGER :: i, j

      j = FLOOR ( DBLE(volume-1)/(Nx-1) ) + 1
      i = volume - (j - 1)*(Nx-1)
      valor = -DCOS(xu(i))*DSIN(yu(j))*DEXP(-2.0d0*(ht*(t-1)+Ti))

   END FUNCTION VelocidadeU_AnaliticaExperimento5

   DOUBLE PRECISION FUNCTION VelocidadeV_AnaliticaExperimento5(t, volume) RESULT(valor)
      USE VariaveisGerais
      INTEGER, INTENT(IN) :: t, volume
      INTEGER :: i, j

      j = FLOOR ( DBLE(volume-1)/(Nx) ) + 1
      i = volume - (j - 1)*(Nx)
      valor = DSIN(xv(i))*DCOS(yv(j))*DEXP(-2.0d0*(ht*(t-1)+Ti))

   END FUNCTION VelocidadeV_AnaliticaExperimento5
 
   !-------------------- FUNCOES DO EXPERIMENTO 5 -----------------------
   
END MODULE ExperimentosNumericos
