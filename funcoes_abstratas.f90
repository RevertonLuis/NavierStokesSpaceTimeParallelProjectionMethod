MODULE FuncoesAbstratas

   IMPLICIT NONE

   !CONTAINS

   ABSTRACT INTERFACE

      DOUBLE PRECISION FUNCTION FuncaoTempoVolume (t, volume)
         INTEGER, INTENT(IN) :: t, volume
      END FUNCTION FuncaoTempoVolume


      DOUBLE PRECISION FUNCTION FuncaoTempoVolume2 (t, volume, indice, Nx, vetor)
         INTEGER, INTENT(IN) :: t, volume, Nx, indice
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vetor
      END FUNCTION FuncaoTempoVolume2


      DOUBLE PRECISION FUNCTION FuncaoTempoVolume3 (t, volume, Nx, dx, dy, dt, vetor1, vetor2 )
         INTEGER, INTENT (IN) :: volume, Nx, t
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vetor1
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vetor2
      END FUNCTION FuncaoTempoVolume3

      DOUBLE PRECISION FUNCTION FuncaoResiduoVolume (volume, Nx, dx, dy, vetor, vetor_bp )
         INTEGER, INTENT (IN) :: volume, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy 
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vetor_bp
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vetor
      END FUNCTION FuncaoResiduoVolume
      
      SUBROUTINE SubrotinaTempoVolume (t, t0, volume, Nx, dx, dy, dt, Re, u, u0, v, v0, D, A, A0, C, pp )
         INTEGER, INTENT (IN) :: t, t0, volume, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, u0, v, v0
         INTEGER,          INTENT(INOUT) :: pp
         DOUBLE PRECISION, INTENT(INOUT) :: D, A, A0, C
      END SUBROUTINE SubrotinaTempoVolume
   
                  
      SUBROUTINE SubrotinaTempoVolume2 (volume, Nx, dx, dy, dt, Re, vetor, vetor_bp )
         INTEGER, INTENT (IN) :: volume, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy, dt, Re
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: vetor_bp
         DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: vetor
      END SUBROUTINE SubrotinaTempoVolume2


      SUBROUTINE SubrotinaTempoVolume3 (volume, Nx, dx, dy, vetor, vetor_bp )
         INTEGER, INTENT (IN) :: volume, Nx
         DOUBLE PRECISION, INTENT(IN) :: dx, dy 
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: vetor_bp
         DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: vetor
      END SUBROUTINE SubrotinaTempoVolume3 

      SUBROUTINE SubrotinaProlongacaoVolume (volume_fina, volume_grossa, Nx_fina, Nx_grossa, vetor_fina, vetor_grossa )
         INTEGER, INTENT (IN) :: volume_fina, volume_grossa, Nx_fina, Nx_grossa
         DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: vetor_fina
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: vetor_grossa         
      END SUBROUTINE SubrotinaProlongacaoVolume
      
      
      
   END INTERFACE

END MODULE FuncoesAbstratas
