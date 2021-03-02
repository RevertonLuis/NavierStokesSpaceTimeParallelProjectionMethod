MODULE ClassMatrizA

   ! Modulo que implementa a classe coeficientes
   
   IMPLICIT  NONE
   
   PRIVATE

   ! Os stencils da matriz A
   TYPE, PRIVATE :: stencil
      INTEGER,          DIMENSION(:), ALLOCATABLE :: p ! As posicoes do stencil     -1, +1, +Nx, -Nx,  0, +2Nx,  -2, etc...
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: c ! Os coeficientes do stencil aw, ae,  an,  as, ap,  ann, aww, etc...
   END TYPE stencil
      
   TYPE, PUBLIC :: MatrizA
   
      ! Atritubos Externos serao atualizados por outra rotina
      ! Atributos Internos serao atualizados na classe
   
      !---------Atributos Externos-----------------------
      INTEGER :: Nx, Ny
      DOUBLE PRECISION :: hx, hy, ht, Re
      !--------------------------------------------------
      
      
      !--------Atributos Internos------------------------
      ! Para o problema estudado sao 9 stencils na matriz para U, V e P
      TYPE( stencil ), DIMENSION(9) :: s       ! Os stencils da matriz A
      INTEGER, DIMENSION(:), ALLOCATABLE :: v  ! Cada volume V vai utilizar um stencils s
      !--------------------------------------------------

      CONTAINS
     
         ! Metodos da Classe
         procedure :: init                ! Metodo que inicia a classe
         
         !--------------------------------------------------------------------------
         ! Obs: este metodo cria um vetor do mesmo tamanho que a matriz
         ! se estiver utilizando muita memoria pode-se calcular o stencil
         ! utilizado pelo volume utilizando o metodo calcula_stencil. Desta
         ! forma cada vez que o volume e considerado seu stencil e calculado
         ! em tempo real, ou seja, sera trocado memoria por processamento.
         procedure :: atribui_stencils    ! Metodo que ira atribuir a cada volume 
                                          ! os stencils, utiliza mais memoria mas calcula os 
                                          ! stencils uma unica vez
                  
         procedure :: calcula_stencil     ! Metodo que calcula o stencil utilizado pelo volume 
                                          ! em tempo real, nao utiliza memoria extra mas calcula 
                                          ! os stencils toda vez que o volume for considerado
         !-----------------------------------------------------------------------------
                  
         procedure :: stencils_u  ! Metodo que calcula os stencils para U
         procedure :: stencils_v  ! Metodo que calcula os stencils para V
         procedure :: stencils_p  ! Metodo que calcula os stencils para P
         
   END TYPE MatrizA

   CONTAINS

      SUBROUTINE init(self, Nx, Ny, hx, hy, ht, Re)

         !------------ Metodo que inicia as variaveis gerais da malha ---------------------

         CLASS( MatrizA ), INTENT(INOUT) :: self
         
         INTEGER         , INTENT(IN)      :: Nx, Ny
         DOUBLE PRECISION, INTENT(IN)      :: hx, hy, ht, Re

         ! Atualiando os atributos da classe
         self % Nx = Nx
         self % Ny = Ny
         self % hx = hx
         self % hy = hy
         self % Re = Re
         self % ht = ht
                  
      END SUBROUTINE init

      
      SUBROUTINE atribui_stencils( self )
         
         CLASS (MatrizA), INTENT(INOUT) :: self
         
         INTEGER :: p, Nx, Ny

         ! Apenas para ficar mais clean o codigo abaixo
         Nx = self % Nx
         Ny = self % Ny
         
         ! Reservando o vetor de stencils
         ALLOCATE( self % v( Nx * Ny ) )
                 
         ! Atribuindo o stencil central para todos os volumes,
         ! desta forma so sera necessario corrigir os contornos
         self % v(:) = 5
         
         ! Posicao do stencil, sao nove posicoes
         ! 1 = SW (South West)
         ! 2 = S  (South)
         ! 3 = SE (South East)
         ! 4 = W  (West)
         ! 5 = C  (Central)
         ! 6 = E  (East)
         ! 7 = NW (North West)
         ! 8 = N  (North)
         ! 9 = NE (North East)
                  
         ! Canto SOUTH WEST recebe stencil 1
         p = 1
         self % v(p) = 1
         
         ! Contorno SOUTH recebe stencil 2
         DO p = 2, Nx - 1
            self % v(p) = 2
         END DO
                           
         ! Canto SOUTH EAST  recebe stencil 3
         p = Nx
         self % v(p) = 3
         
         ! Contorno WEST recebe stencil 4
         DO p = 1 + Nx, Nx *( Ny - 2 ) + 1, Nx
            self % v(p) = 4
         END DO
         
         ! Volumes no CENTRO recebem stencil 5, mas ja foram considerados
         
         ! Contorno EAST recebe stencil 6
         DO p = 2 * Nx, Nx * (Ny - 1), Nx
            self % v(p) = 6
         END DO
         
         ! Canto NORTH WEST recebe stencil 7
         p = Nx *( Ny - 1 ) + 1
         self % v(p) = 7
         
         ! Contorno NORTH recebe stencil 8
         DO p = Nx *( Ny - 1 ) + 2, Nx * Ny - 1
            self % v(p) = 8
         END DO
                  
         ! Canto NORTH EAST recebe stencil 9
         p = Nx * Ny
         self % v(p) = 9
         
      END SUBROUTINE atribui_stencils

      INTEGER FUNCTION calcula_stencil( self, volume )

         CLASS (MatrizA), INTENT(INOUT) :: self
         
         INTEGER, INTENT(IN) :: volume
         INTEGER             :: i, j, Nx, Ny
         INTEGER             :: s1, s2, s3, s4, s5, s6, s7, s8, s9
         
         j = (volume-1)/(self % Nx) + 1
         i = volume - (j-1)*(self % Nx)
            
         s2 = 1/j
         s4 = 1/i
         s6 = i/(self % Nx)
         s8 = j/(self % Ny)
         s1 = (s2 + s4)/2
         s3 = (s2 + s6)/2
         s5 = (5 - (s2 + s4 + s6 + s8))/5
         s7 = (s4 + s8)/2
         s9 = (s6 + s8)/2
         calcula_stencil = ABS( 2*s2 + 4*s4 + 6*s6 + 8*s8 - (7*s1 + 11*s3 + 19*s7 + 23*s9) + 5*s5 )

      END FUNCTION calcula_stencil

! STENCILS U ------------------------------------------------------------------------------------------------------------------------------------------      
                  
      SUBROUTINE stencils_u( self, alfa2 )
         
         CLASS (MatrizA), INTENT(INOUT) :: self
         
         INTEGER :: Nx, Ny
         DOUBLE PRECISION :: Re, ht, hx, hy
         DOUBLE PRECISION, INTENT(IN) :: alfa2
         
         ! Apenas para ficar mais clean o codigo abaixo          
         hx = self % hx
         hy = self % hy
         ht = self % ht
         Re = self % Re
         Nx = self % Nx
         
         ! Sao nove stencils
         ! 1 = SW (South West)
         ! 2 = S  (South)
         ! 3 = SE (South East)
         ! 4 = W  (West)
         ! 5 = C  (Central)
         ! 6 = E  (East)
         ! 7 = NW (North West)
         ! 8 = N  (North)
         ! 9 = NE (North East)
         
         !------- 1 = SW (South West) -------
         IF ( .NOT. ALLOCATED( self % s(1) % c ) )    ALLOCATE( self % s(1) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(1) % p ) )    ALLOCATE( self % s(1) % p(4) )
         
         ! ap (sempre na primeira posicao dos coeficientes)
         self % s(1) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 5.0d0 * ht/(Re*hy**2)
         self % s(1) % p(1) = 0
         
         ! ae
         self % s(1) % c(2) = ht/(Re*hx**2)
         self % s(1) % p(2) = 1

         ! an
         self % s(1) % c(3) = 2.0d0 * ht/(Re*hy**2)
         self % s(1) % p(3) = Nx

         !ann 
         self % s(1) % c(4) = -(1.0d0/5.0d0)*(ht/(Re*hy**2))
         self % s(1) % p(4) = 2*Nx
         !------- 1 = SW (South West) -------
         
         !------- 2 = S (South) -------
         IF ( .NOT. ALLOCATED( self % s(2) % c ) )    ALLOCATE( self % s(2) % c(5) )
         IF ( .NOT. ALLOCATED( self % s(2) % p ) )    ALLOCATE( self % s(2) % p(5) )

         ! ap (sempre na primeira posicao dos coeficientes)
         self % s(2) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 5.0d0 * ht/(Re*hy**2)
         self % s(2) % p(1) = 0

         ! aw
         self % s(2) % c(2) = ht/(Re*hx**2)
         self % s(2) % p(2) = -1

         ! ae
         self % s(2) % c(3) = ht/(Re*hx**2)
         self % s(2) % p(3) = 1

         ! an
         self % s(2) % c(4) = 2.0d0 * ht/(Re*hy**2)
         self % s(2) % p(4) = Nx

         ! ann
         self % s(2) % c(5) = -(1.0d0/5.0d0)*(ht/(Re*hy**2))
         self % s(2) % p(5) = 2*Nx
         !------- 2 = S (South) ------- 
         
         !------- 3 = SE (South East) -------
         IF ( .NOT. ALLOCATED( self % s(3) % c ) )    ALLOCATE( self % s(3) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(3) % p ) )    ALLOCATE( self % s(3) % p(4) )
         
         ! ap (sempre na primeira posicao dos coeficientes)
         self % s(3) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 5.0d0 * ht/(Re*hy**2)        
         self % s(3) % p(1) = 0
   
         ! aw
         self % s(3) % c(2) = ht/(Re*hx**2)
         self % s(3) % p(2) = -1
         
         ! an
         self % s(3) % c(3) = 2.0d0 * ht/(Re*hy**2)
         self % s(3) % p(3) = Nx

         ! ann
         self % s(3) % c(4) = -(1.0d0/5.0d0)*(ht/(Re*hy**2))
         self % s(3) % p(4) = 2*Nx
         !------- 3 = SE (South East) -------
       
         !------- 4 = W (West) -------
         IF ( .NOT. ALLOCATED( self % s(4) % c ) )    ALLOCATE( self % s(4) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(4) % p ) )    ALLOCATE( self % s(4) % p(4) )
         
         ! ap (sempre na primeira posicao dos coeficientes) 
         self % s(4) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 2.0d0 * ht/(Re*hy**2) 
         self % s(4) % p(1) = 0
         
         ! ae
         self % s(4) % c(2) = ht/(Re*hx**2)
         self % s(4) % p(2) = 1

         ! as
         self % s(4) % c(3) = ht/(Re*hy**2)
         self % s(4) % p(3) = -Nx

         ! an
         self % s(4) % c(4) = ht/(Re*hy**2)
         self % s(4) % p(4) = Nx
         !------- 4 = W (West) -------  
         
         !------- 5 = C (Central) ------- 
         IF ( .NOT. ALLOCATED( self % s(5) % c ) )    ALLOCATE( self % s(5) % c(5) )
         IF ( .NOT. ALLOCATED( self % s(5) % p ) )    ALLOCATE( self % s(5) % p(5) )
          
         ! ap (sempre na primeira posicao dos coeficientes) 
         self % s(5) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 2.0d0 * ht/(Re*hy**2)  
         self % s(5) % p(1) = 0
        
         ! aw
         self % s(5) % c(2) = ht/(Re*hx**2)
         self % s(5) % p(2) = -1

         ! ae
         self % s(5) % c(3) = ht/(Re*hx**2)
         self % s(5) % p(3) = 1

         ! as
         self % s(5) % c(4) = ht/(Re*hy**2)
         self % s(5) % p(4) = -Nx
        
         ! an
         self % s(5) % c(5) = ht/(Re*hy**2)
         self % s(5) % p(5) = Nx
         !------- 5 = C (Central) ------- 
  
         !------- 6 = E (East) -------
         IF ( .NOT. ALLOCATED( self % s(6) % c ) )    ALLOCATE( self % s(6) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(6) % p ) )    ALLOCATE( self % s(6) % p(4) )
         
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(6) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 2.0d0 * ht/(Re*hy**2) 
         self % s(6) % p(1) = 0

         ! aw
         self % s(6) % c(2) = ht/(Re*hx**2)
         self % s(6) % p(2) = -1

         ! as
         self % s(6) % c(3) = ht/(Re*hy**2)
         self % s(6) % p(3) = -Nx

         ! an
         self % s(6) % c(4) = ht/(Re*hy**2)
         self % s(6) % p(4) = Nx
         !------- 6 = E (East) ------- 

         !------- 7 = NW (North West) -------
         IF ( .NOT. ALLOCATED( self % s(7) % c ) )    ALLOCATE( self % s(7) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(7) % p ) )    ALLOCATE( self % s(7) % p(4) )
         
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(7) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 5.0d0 * ht/(Re*hy**2)
         self % s(7) % p(1) = 0

         ! ae
         self % s(7) % c(2) = ht/(Re*hx**2)
         self % s(7) % p(2) = 1

         ! as
         self % s(7) % c(3) = 2.0d0 * ht/(Re*hy**2)
         self % s(7) % p(3) = -Nx

         ! ass
         self % s(7) % c(4) = -(1.0d0/5.0d0)*(ht/(Re*hy**2)) 
         self % s(7) % p(4) = -2*Nx
         !------- 7 = NW (North West) -------
 

         !------- 8 = N (North) -------
         IF ( .NOT. ALLOCATED( self % s(8) % c ) )    ALLOCATE( self % s(8) % c(5) )
         IF ( .NOT. ALLOCATED( self % s(8) % p ) )    ALLOCATE( self % s(8) % p(5) )
         
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(8) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 5.0d0 * ht/(Re*hy**2)
         self % s(8) % p(1) = 0
         
         ! aw
         self % s(8) % c(2) =  ht/(Re*hx**2)
         self % s(8) % p(2) = -1

         ! ae
         self % s(8) % c(3) = ht/(Re*hx**2)
         self % s(8) % p(3) = 1
 
         ! as
         self % s(8) % c(4) = 2.0d0 * ht/(Re*hy**2)
         self % s(8) % p(4) = -Nx

         ! ass
         self % s(8) % c(5) = -(1.0d0/5.0d0)*(ht/(Re*hy**2))
         self % s(8) % p(5) = -2*Nx
         !------- 8 = N (North) -------

         !------- 9 = NE (North East) ------- 
         IF ( .NOT. ALLOCATED( self % s(9) % c ) )    ALLOCATE( self % s(9) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(9) % p ) )    ALLOCATE( self % s(9) % p(4) )
         
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(9) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 5.0d0 * ht/(Re*hy**2)
         self % s(9) % p(1) = 0

         ! aw
         self % s(9) % c(2) =  ht/(Re*hx**2)
         self % s(9) % p(2) = -1
         
         ! as
         self % s(9) % c(3) = 2.0d0 * ht/(Re*hy**2)
         self % s(9) % p(3) = -Nx

         ! ass
         self % s(9) % c(4) = -(1.0d0/5.0d0)*(ht/(Re*hy**2))
         self % s(9) % p(4) = -2*Nx
         !------- 9 = NE (North East) ------- 

      END SUBROUTINE stencils_u

! STENCILS V ------------------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE stencils_v( self, alfa2 )
         
         CLASS (MatrizA), INTENT(INOUT) :: self
         
         INTEGER :: Nx, Ny
         DOUBLE PRECISION :: Re, ht, hx, hy
         DOUBLE PRECISION, INTENT(IN) :: alfa2
         
         ! Apenas para ficar mais clean o codigo abaixo          
         hx = self % hx
         hy = self % hy
         ht = self % ht
         Re = self % Re
         Nx = self % Nx
         
         ! Sao nove stencils
         ! 1 = SW (South West)
         ! 2 = S  (South)
         ! 3 = SE (South East)
         ! 4 = W  (West)
         ! 5 = C  (Central)
         ! 6 = E  (East)
         ! 7 = NW (North West)
         ! 8 = N  (North)
         ! 9 = NE (North East)
         
         !------- 1 = SW (South West) -------
         IF ( .NOT. ALLOCATED( self % s(1) % c ) )    ALLOCATE( self % s(1) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(1) % p ) )    ALLOCATE( self % s(1) % p(4) )
         
         ! ap (sempre na primeira posicao dos coeficientes)
         self % s(1) % c(1) = alfa2 + 5.0d0*(ht/(Re*hx**2)) + 2.0d0 * ht/(Re*hy**2)
         self % s(1) % p(1) = 0
         
         ! ae
         self % s(1) % c(2) = 2.0d0 * ht/(Re*hx**2)
         self % s(1) % p(2) = 1

         ! an
         self % s(1) % c(3) = ht/(Re*hy**2)
         self % s(1) % p(3) = Nx

         !aee 
         self % s(1) % c(4) = -(1.0d0/5.0d0)*(ht/(Re*hx**2))
         self % s(1) % p(4) = 2
         !------- 1 = SW (South West) -------
         
         !------- 2 = S (South) -------
         IF ( .NOT. ALLOCATED( self % s(2) % c ) )    ALLOCATE( self % s(2) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(2) % p ) )    ALLOCATE( self % s(2) % p(4) )
         
         ! ap (sempre na primeira posicao dos coeficientes)
         self % s(2) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 2.0d0 * ht/(Re*hy**2) 
         self % s(2) % p(1) = 0

         ! aw
         self % s(2) % c(2) = ht/(Re*hx**2)
         self % s(2) % p(2) = -1

         ! ae
         self % s(2) % c(3) = ht/(Re*hx**2)
         self % s(2) % p(3) = 1

         ! an
         self % s(2) % c(4) = ht/(Re*hy**2)
         self % s(2) % p(4) = Nx

         !------- 2 = S (South) ------- 
         
         !------- 3 = SE (South East) -------
         IF ( .NOT. ALLOCATED( self % s(3) % c ) )    ALLOCATE( self % s(3) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(3) % p ) )    ALLOCATE( self % s(3) % p(4) )
         
         ! ap (sempre na primeira posicao dos coeficientes)
         self % s(3) % c(1) = alfa2 + 5.0d0*(ht/(Re*hx**2)) + 2.0d0 * ht/(Re*hy**2)
         self % s(3) % p(1) = 0
   
         ! aw
         self % s(3) % c(2) = 2.0d0 * ht/(Re*hx**2)
         self % s(3) % p(2) = -1
         
         ! an
         self % s(3) % c(3) = ht/(Re*hy**2)
         self % s(3) % p(3) = Nx

         ! aww
         self % s(3) % c(4) = -(1.0d0/5.0d0)*(ht/(Re*hx**2))
         self % s(3) % p(4) = -2
         !------- 3 = SE (South East) -------
       
         !------- 4 = W (West) -------
         IF ( .NOT. ALLOCATED( self % s(4) % c ) )    ALLOCATE( self % s(4) % c(5) )
         IF ( .NOT. ALLOCATED( self % s(4) % p ) )    ALLOCATE( self % s(4) % p(5) )
         
         ! ap (sempre na primeira posicao dos coeficientes) 
         self % s(4) % c(1) = alfa2 + 5.0d0*(ht/(Re*hx**2)) + 2.0d0 * ht/(Re*hy**2)
         self % s(4) % p(1) = 0
         
         ! ae
         self % s(4) % c(2) = 2.0d0 * ht/(Re*hx**2)
         self % s(4) % p(2) = 1

         ! as
         self % s(4) % c(3) = ht/(Re*hy**2)
         self % s(4) % p(3) = -Nx

         ! an
         self % s(4) % c(4) = ht/(Re*hy**2)
         self % s(4) % p(4) = Nx
         
         ! aee
         self % s(4) % c(5) = -(1.0d0/5.0d0)*(ht/(Re*hx**2))
         self % s(4) % p(5) = 2
         !------- 4 = W (West) -------  
         
         !------- 5 = C (Central) ------- 
         IF ( .NOT. ALLOCATED( self % s(5) % c ) )    ALLOCATE( self % s(5) % c(5) )
         IF ( .NOT. ALLOCATED( self % s(5) % p ) )    ALLOCATE( self % s(5) % p(5) )
                  
         ! ap (sempre na primeira posicao dos coeficientes) 
         self % s(5) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 2.0d0 * ht/(Re*hy**2) 
         self % s(5) % p(1) = 0
        
         ! aw
         self % s(5) % c(2) = ht/(Re*hx**2)
         self % s(5) % p(2) = -1

         ! ae
         self % s(5) % c(3) = ht/(Re*hx**2)
         self % s(5) % p(3) = 1

         ! as
         self % s(5) % c(4) = ht/(Re*hy**2)
         self % s(5) % p(4) = -Nx
        
         ! an
         self % s(5) % c(5) = ht/(Re*hy**2)
         self % s(5) % p(5) = Nx
         !------- 5 = C (Central) ------- 
  
         !------- 6 = E (East) -------
         IF ( .NOT. ALLOCATED( self % s(6) % c ) )    ALLOCATE( self % s(6) % c(5) )
         IF ( .NOT. ALLOCATED( self % s(6) % p ) )    ALLOCATE( self % s(6) % p(5) )
                  
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(6) % c(1) = alfa2 + 5.0d0*(ht/(Re*hx**2)) + 2.0d0 * ht/(Re*hy**2)
         self % s(6) % p(1) = 0

         ! aw
         self % s(6) % c(2) = 2.0d0 * ht/(Re*hx**2)
         self % s(6) % p(2) = -1

         ! as
         self % s(6) % c(3) = ht/(Re*hy**2)
         self % s(6) % p(3) = -Nx

         ! an
         self % s(6) % c(4) = ht/(Re*hy**2)
         self % s(6) % p(4) = Nx
         
         ! aww
         self % s(6) % c(5) = -(1.0d0/5.0d0)*(ht/(Re*hx**2))
         self % s(6) % p(5) = -2
         !------- 6 = E (East) ------- 

         !------- 7 = NW (North West) -------
         IF ( .NOT. ALLOCATED( self % s(7) % c ) )    ALLOCATE( self % s(7) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(7) % p ) )    ALLOCATE( self % s(7) % p(4) )
               
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(7) % c(1) = alfa2 + 5.0d0*(ht/(Re*hx**2)) + 2.0d0 * ht/(Re*hy**2)
         self % s(7) % p(1) = 0

         ! ae
         self % s(7) % c(2) =  2.0d0 * ht/(Re*hx**2)
         self % s(7) % p(2) = 1

         ! as
         self % s(7) % c(3) =  ht/(Re*hy**2)
         self % s(7) % p(3) = -Nx

         ! aee
         self % s(7) % c(4) = -(1.0d0/5.0d0)*(ht/(Re*hx**2))
         self % s(7) % p(4) = 2
         !------- 7 = NW (North West) -------
 

         !------- 8 = N (North) -------
         IF ( .NOT. ALLOCATED( self % s(8) % c ) )    ALLOCATE( self % s(8) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(8) % p ) )    ALLOCATE( self % s(8) % p(4) )
                  
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(8) % c(1) = alfa2 + 2.0d0 * ht/(Re*hx**2) + 2.0d0 * ht/(Re*hy**2) 
         self % s(8) % p(1) = 0
         
         ! aw
         self % s(8) % c(2) = ht/(Re*hx**2)
         self % s(8) % p(2) = -1

         ! ae
         self % s(8) % c(3) =  ht/(Re*hx**2)
         self % s(8) % p(3) = 1
 
         ! as
         self % s(8) % c(4) =  ht/(Re*hy**2)
         self % s(8) % p(4) = -Nx

         !------- 8 = N (North) -------

         !------- 9 = NE (North East) ------- 
         IF ( .NOT. ALLOCATED( self % s(9) % c ) )    ALLOCATE( self % s(9) % c(4) )
         IF ( .NOT. ALLOCATED( self % s(9) % p ) )    ALLOCATE( self % s(9) % p(4) )
                  
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(9) % c(1) =  alfa2 + 5.0d0*(ht/(Re*hx**2)) + 2.0d0 * ht/(Re*hy**2)
         self % s(9) % p(1) = 0

         ! aw
         self % s(9) % c(2) = 2.0d0 * ht/(Re*hx**2)
         self % s(9) % p(2) = -1
         
         ! as
         self % s(9) % c(3) = ht/(Re*hy**2)
         self % s(9) % p(3) = -Nx

         ! aww
         self % s(9) % c(4) = -(1.0d0/5.0d0)*(ht/(Re*hx**2))
         self % s(9) % p(4) = -2
         !------- 9 = NE (North East) ------- 

      END SUBROUTINE stencils_v
      
! STENCILS P ------------------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE stencils_p( self )
         
         CLASS (MatrizA), INTENT(INOUT) :: self
         
         INTEGER :: Nx, Ny
         DOUBLE PRECISION :: hx, hy
         
         ! Apenas para ficar mais clean o codigo abaixo          
         hx = self % hx
         hy = self % hy
         Nx = self % Nx
         
         ! Sao nove stencils
         ! 1 = SW (South West)
         ! 2 = S  (South)
         ! 3 = SE (South East)
         ! 4 = W  (West)
         ! 5 = C  (Central)
         ! 6 = E  (East)
         ! 7 = NW (North West)
         ! 8 = N  (North)
         ! 9 = NE (North East)
         
         !------- 1 = SW (South West) -------
         ALLOCATE( self % s(1) % c(3) )
         ALLOCATE( self % s(1) % p(3) )
         
         ! ap (sempre na primeira posicao dos coeficientes)
         self % s(1) % c(1) = 1.0d0/(hx**2) + 1.0d0/(hy**2)
         self % s(1) % p(1) = 0
         
         ! ae
         self % s(1) % c(2) = 1.0d0/(hx**2)
         self % s(1) % p(2) = 1

         ! an
         self % s(1) % c(3) = 1.0d0/(hy**2)
         self % s(1) % p(3) = Nx

         !------- 1 = SW (South West) -------
         
         !------- 2 = S (South) -------
         ALLOCATE( self % s(2) % c(4) )
         ALLOCATE( self % s(2) % p(4) )

         ! ap (sempre na primeira posicao dos coeficientes)
         self % s(2) % c(1) = 2.0d0/(hx**2) + 1.0d0/(hy**2)
         self % s(2) % p(1) = 0

         ! aw
         self % s(2) % c(2) = 1.0d0/(hx**2)
         self % s(2) % p(2) = -1

         ! ae
         self % s(2) % c(3) =  1.0d0/(hx**2)
         self % s(2) % p(3) = 1

         ! an
         self % s(2) % c(4) = 1.0d0/(hy**2)
         self % s(2) % p(4) = Nx

         !------- 2 = S (South) ------- 
         
         !------- 3 = SE (South East) -------
         ALLOCATE( self % s(3) % c(3) ) 
         ALLOCATE( self % s(3) % p(3) ) 

         ! ap (sempre na primeira posicao dos coeficientes)
         self % s(3) % c(1) = 1.0d0/(hx**2) + 1.0d0/(hy**2)
         self % s(3) % p(1) = 0
   
         ! aw
         self % s(3) % c(2) = 1.0d0/(hx**2)
         self % s(3) % p(2) = -1
         
         ! an
         self % s(3) % c(3) =  1.0d0/(hy**2)
         self % s(3) % p(3) = Nx
         !------- 3 = SE (South East) -------
       
         !------- 4 = W (West) -------
         ALLOCATE( self % s(4) % c(4) ) 
         ALLOCATE( self % s(4) % p(4) ) 

         ! ap (sempre na primeira posicao dos coeficientes) 
         self % s(4) % c(1) = 1.0d0/(hx**2) + 2.0d0/(hy**2)
         self % s(4) % p(1) = 0
         
         ! ae
         self % s(4) % c(2) = 1.0d0/(hx**2)
         self % s(4) % p(2) = 1

         ! as
         self % s(4) % c(3) = 1.0d0/(hy**2) 
         self % s(4) % p(3) = -Nx

         ! an
         self % s(4) % c(4) = 1.0d0/(hy**2)
         self % s(4) % p(4) = Nx
         !------- 4 = W (West) -------  
         
         !------- 5 = C (Central) ------- 
         ALLOCATE( self % s(5) % c(5) )
         ALLOCATE( self % s(5) % p(5) )
          
         ! ap (sempre na primeira posicao dos coeficientes) 
         self % s(5) % c(1) = 2.0d0/(hx**2) + 2.0d0/(hy**2)
         self % s(5) % p(1) = 0
        
         ! aw
         self % s(5) % c(2) = 1.0d0/(hx**2)
         self % s(5) % p(2) = -1

         ! ae
         self % s(5) % c(3) = 1.0d0/(hx**2)
         self % s(5) % p(3) = 1

         ! as
         self % s(5) % c(4) = 1.0d0/(hy**2) 
         self % s(5) % p(4) = -Nx
        
         ! an
         self % s(5) % c(5) = 1.0d0/(hy**2) 
         self % s(5) % p(5) = Nx
         !------- 5 = C (Central) ------- 
  
         !------- 6 = E (East) -------
         ALLOCATE( self % s(6) % c(4) ) 
         ALLOCATE( self % s(6) % p(4) ) 
         
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(6) % c(1) = 1.0d0/(hx**2) + 2.0d0/(hy**2)
         self % s(6) % p(1) = 0

         ! aw
         self % s(6) % c(2) = 1.0d0/(hx**2)
         self % s(6) % p(2) = -1

         ! as
         self % s(6) % c(3) = 1.0d0/(hy**2) 
         self % s(6) % p(3) = -Nx

         ! an
         self % s(6) % c(4) = 1.0d0/(hy**2)
         self % s(6) % p(4) = Nx
         !------- 6 = E (East) ------- 

         !------- 7 = NW (North West) -------
         ALLOCATE( self % s(7) % c(3) )
         ALLOCATE( self % s(7) % p(3) )
      
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(7) % c(1) = 1.0d0/(hx**2) + 1.0d0/(hy**2)
         self % s(7) % p(1) = 0

         ! ae
         self % s(7) % c(2) = 1.0d0/(hx**2)
         self % s(7) % p(2) = 1

         ! as
         self % s(7) % c(3) =  1.0d0/(hy**2) 
         self % s(7) % p(3) = -Nx

         !------- 7 = NW (North West) -------
 

         !------- 8 = N (North) -------
         ALLOCATE( self % s(8) % c(4) )
         ALLOCATE( self % s(8) % p(4) )
         
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(8) % c(1) = 2.0d0/(hx**2) + 1.0d0/(hy**2)
         self % s(8) % p(1) = 0
         
         ! aw
         self % s(8) % c(2) = 1.0d0/(hx**2)
         self % s(8) % p(2) = -1

         ! ae
         self % s(8) % c(3) = 1.0d0/(hx**2)
         self % s(8) % p(3) = 1
 
         ! as
         self % s(8) % c(4) = 1.0d0/(hy**2)
         self % s(8) % p(4) = -Nx
         !------- 8 = N (North) -------

         !------- 9 = NE (North East) ------- 
         ALLOCATE( self % s(9) % c(3) )
         ALLOCATE( self % s(9) % p(3) )
         
         ! ap (sempre na primeira posicao dos coeficientes)  
         self % s(9) % c(1) = 1.0d0/(hx**2) + 1.0d0/(hy**2)
         self % s(9) % p(1) = 0

         ! aw
         self % s(9) % c(2) = 1.0d0/(hx**2)
         self % s(9) % p(2) = -1
         
         ! as
         self % s(9) % c(3) = 1.0d0/(hy**2)
         self % s(9) % p(3) = -Nx
         !------- 9 = NE (North East) ------- 

      END SUBROUTINE stencils_p

END MODULE ClassMatrizA
