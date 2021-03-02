MODULE ClassVelocidadeV

   USE ClassArrayFonteComplementosV
   USE ClassMatrizA
   
   ! Modulo que implementa a classe para a velocidade V
   
   IMPLICIT  NONE
  
   PRIVATE
  
   TYPE, PUBLIC :: VelocidadeV
   
      ! Atritubos Externos serao atualizados por outra rotina
      ! Atributos Internos serao atualizados na classe
   
      !---------Atributos Externos-----------------------
      INTEGER           :: Tn                                            ! Numero de passos de tempo resolvidos em paralelo
      INTEGER           :: Nx, Ny                                      ! Nx e Ny da malha
      DOUBLE PRECISION  :: Lx, Ly                                      ! Comprimento do dominio
      DOUBLE PRECISION  :: ht                                          ! Refinamento temporal (necessario por causa dos coeficientes)
      DOUBLE PRECiSION  :: Re                                          ! Numero de Reynolds (necessario por causa dos coeficientes)
      DOUBLE PRECISION  :: Tolerancia                                  ! Tolerancia no criterio de parada dos solvers
      !--------------------------------------------------
      
      
      !--------Atributos Internos------------------------
      DOUBLE PRECISION                                 :: hx, hy       ! hx e hy da malha
      
      ! ATENCAO: DEFINIR UMA VETOR DE VETORES E RESERVAR VARIAVEIS v0, v1, v2, v3, ..., vt dinamicamente (de acordo com o a ordem do metodo de projecao no futuro)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE      :: v   ! Vetores da velocidade V
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE      :: v_bp, v_ex, residuo
      DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE      :: v_ref
      DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE      :: normas_L2

      ! A matriz A de V
      TYPE ( MatrizA ) :: A
      
      TYPE (ArrayFonteComplementosV), DIMENSION(9) :: FonteComplementosV 
      !--------------------------------------------------

      CONTAINS
     
         ! Metodos da Classe
         procedure :: init                    ! Metodo que inicia a classe
         procedure :: aplica_condicao_inicial ! 
         procedure :: calcula_termo_fonte
         procedure :: single_grid
         procedure :: calcula_residuo
         procedure :: calcula_norma_L2 
         procedure :: escreve_matriz
         procedure :: escreve_comparacao
     
   END TYPE VelocidadeV

   CONTAINS

      SUBROUTINE init(self, Nx, Ny, Tn, Lx, Ly, ht, Re, Tolerancia)
         
         USE Subrotinas_MG_GS  
       
         !------------ Metodo que inicia as variaveis gerais da malha ---------------------

         CLASS( VelocidadeV ), INTENT(INOUT) :: self
         
         INTEGER         , INTENT(IN) :: Nx, Ny, Tn
         DOUBLE PRECISION, INTENT(IN) :: Lx, Ly, ht, Re, Tolerancia
         
         INTEGER :: stencil
         
         ! Atualiando os atributos da classe
         self % Tn = Tn
         self % Lx = Lx
         self % Ly = Ly
         self % Re = Re
         self % ht = ht
         self % Tolerancia = Tolerancia
         
         ! Calculando hx e hy (o mesmo valor para u, v e p)
         self % hx = self % Lx / Nx
         self % hy = self % Ly / Ny
         
         ! Atencao: Nx e Ny depende do tipo de variavel considerada (u, v ou p)
         ! entao aqui sera "corrigido" Nx e Ny para a velocidade U
         ! Apos esta correcao vai aparecer no codigo self % Nx e self % Ny (corrigidos) e nao
         ! Nx-1 e Ny
         self % Nx = Nx 
         self % Ny = Ny - 1
         
         ! Reservando os vetores v0, v1, vt
         ! Lembrando que a dimensao real e v0( Nx*(Ny-1) )
         ALLOCATE( self % v( self % Nx * self % Ny, self % Tn ), &
                   self % v_bp( self % Nx * self % Ny, self % Tn ), &
                   self % v_ex( self % Nx * self % Ny, self % Tn ), &
                   self % residuo( self % Nx * self % Ny, self % Tn ), &
                   self % normas_L2( self % Tn ), &
                   self % v_ref( self % Nx * self % Ny ) )
                   
         ! Iniciando os vetores
         self % v(:, :) = 0.0d0
         self % v_bp(:, :) = 0.0d0
         self % v_ex(:, :) = 0.0d0
         self % v_ref(:) = 0.0d0
         self % residuo(:, :) = 0.0d0
         self % normas_L2(:) = 0.0d0
                  
         !-------------------------- matriz A de V --------------------------------
         ! Inicializando a matriz A de U
         CALL self % A % init(self % Nx, self % Ny, &
                              self % hx, self % hy, self % ht, & 
                              self % Re)
                              
         ! Atribuindo um stencil a cada volume de A                     
         CALL self % A % atribui_stencils()
         
         ! Calculandos os coeficientes e posicoes dos stencils
         CALL self % A % stencils_v(1.0d0)
         !-------------------------- matriz A de V --------------------------------                
         
         ! Inicializando o array de fontes de U
         DO stencil = 1, 9
            CALL self % FonteComplementosV(stencil) % init(stencil)
         END DO
         
         CALL calcula_exata(self)
         
         CALL Atualiza_Variaveis_MG_GS_V(self % Lx, self % Ly, &
                                         self % Nx, self % Ny, &
                                         self % Tolerancia, &
                                         3, 3, &
                                         50, &
                                         self % ht, &
                                         self % Re, self % Tn )
         
      END SUBROUTINE init

      
      SUBROUTINE aplica_condicao_inicial( self )
         USE FuncoesAlias

         CLASS (VelocidadeV), INTENT(INOUT) :: self
         
         INTEGER :: volume
         
         DO volume = 1, self % Nx * self % Ny
            
            self % v(volume, :) = v_inicial(0, volume)
            
         END DO

      END SUBROUTINE aplica_condicao_inicial
      
      SUBROUTINE calcula_exata( self )
         USE FuncoesAlias

         CLASS (VelocidadeV), INTENT(INOUT) :: self
         
         INTEGER :: volume, t
         
         DO t = 1, self % Tn
         
             DO volume = 1, self % Nx * self % Ny
            
                 self % v_ex(volume, t) = v_analitica(t, volume)
            
             END DO
             
         END DO

      END SUBROUTINE calcula_exata
      
      SUBROUTINE calcula_termo_fonte( self, t, t0, u, pressao, correcao)
         
         CLASS (VelocidadeV) :: self
         INTEGER :: t, t0
         DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: u, pressao, correcao
         
         INTEGER :: volume, pp, stencil
         DOUBLE PRECISION :: Dv, Av, Av0, Complemento
 
         DO volume = 1, self % Nx * self % Ny
         
            ! Stencil do volume
            stencil = self % A % v(volume)
            !stencil = self % A % calcula_stencil(volume)
            
            ! Calculando a Adveccao e Difusao
            CALL self % FonteComplementosV(stencil) % CalculaComplementos(t, t0, &
                                                                          volume, &
                                                                          self % Nx, &
                                                                          self % hx, self % hy, self % ht, &
                                                                          self % Re, &
                                                                          u(:, t), u(:, t0), &
                                                                          self % v(:, t), self % v(:, t0), &
                                                                          Dv, Av, Av0, Complemento, pp )
            
            ! Aqui esta o calculo do termo fonte
            self % v_bp(volume, t) =  self % ht*( -1.5d0*Av + 0.5d0*Av0 - &
                                      ( pressao(pp + self % Nx, t) - pressao(pp, t) )/(self % hx) - &
                                      ( correcao(pp + self % Nx, t) - correcao(pp, t) )/(self % hx) ) &
                                      + 2.0d0 * Complemento !+ Sx(t, pu) )               
             
         END DO
      END SUBROUTINE calcula_termo_fonte

      
      SUBROUTINE single_grid( self, t, rb)
         
         USE GS
   
         CLASS (VelocidadeV) :: self
         
         INTEGER, INTENT(IN) :: t, rb
         
         ! Resolvendo com Gauss Seidel
         CALL GaussSeidel_UV( self % A, self % v, self % v_bp, t, rb)
         
      END SUBROUTINE single_grid
      

      SUBROUTINE calcula_residuo( self, t )
   
          USE Residuo
        
          CLASS (VelocidadeV) :: self
         
          INTEGER, INTENT(IN) :: t
         
          ! Calcula residuo
          CALL calcula_residuo_UV( self % A, self % v, self % v_bp, self % residuo, t )
         
      END SUBROUTINE calcula_residuo

       
      SUBROUTINE calcula_norma_L2( self, t )
   
          USE Residuo
        
          CLASS (VelocidadeV) :: self
         
          INTEGER, INTENT(IN) :: t
         
          ! Calcula residuo
          self % normas_L2(t) = SUM(self % residuo(:, t) * self % residuo(:, t))
         
      END SUBROUTINE calcula_norma_L2




      SUBROUTINE escreve_matriz(self, nome, vetor)
   
         CLASS (VelocidadeV) :: self
      
         CHARACTER(len=*), INTENT(IN) :: nome
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vetor
         
         INTEGER :: pv, i, j
         CHARACTER(len=(self % Nx*23)) :: linha
         CHARACTER(len=23) :: valor
      
         OPEN(UNIT=41, FILE=TRIM(ADJUSTL(nome)), STATUS='REPLACE', ACTION='WRITE')
      
         ! Para pontos internos 
         DO j = 1, self % Ny
            DO i = 1, self % Nx
               
               pv = i + (j-1)* self % Nx

               WRITE(valor, "(1PE22.15)") vetor(pv)
               linha((i-1)*23+1:(i-1)*23+23) = valor
             
            END DO
            
            WRITE(41,1) linha

         END DO         
         CLOSE(41)
                
         1 FORMAT( A )
      
      END SUBROUTINE escreve_matriz
      
      
      SUBROUTINE escreve_comparacao(self, tempo, nome, vetor)

         USE FuncoesAlias

         CLASS (VelocidadeV) :: self
      
         CHARACTER(len=*), INTENT(IN) :: nome
         INTEGER, INTENT(IN) :: tempo
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vetor
         
         INTEGER :: pv, i, j, p_Linf
         DOUBLE PRECISION :: Linf
       
         1 FORMAT( 1PE22.15, 3X, 1PE22.15, 3X, 1PE22.15, 3X, I10, 3X, I10 )
         2 FORMAT( 1PE22.15, 1X, "= ht", /,  &
                   1PE22.15, 1X, "= Erro Maximo", /,  &
                   1PE22.15, 1X, "= Analitica (onde o erro aconteceu)", /, &
                   1PE22.15, 1X, "=  Numerica (onde o erro aconteceu)" )

  
         OPEN(UNIT=41, FILE=TRIM(ADJUSTL(nome)), STATUS='REPLACE', ACTION='WRITE')
   
         Linf   = 0.0d0
         p_Linf  = 0
         
         ! Para pontos internos 
         DO j = 1, self % Ny
            DO i = 1, self % Nx
               
               pv = i + (j-1)* self % Nx

               !WRITE(41, 1) v_analitica(tempo, pv), vetor(pv), v_analitica(tempo, pv) - vetor(pv)
                       
               IF ( ABS(v_analitica(tempo, pv) - vetor(pv)) > Linf ) THEN
                  Linf = ABS(v_analitica(tempo, pv) - vetor(pv))
                  p_Linf  = pv
               END IF

            END DO
            
         END DO     

         WRITE(41, 2) self % ht, Linf, v_analitica(tempo, p_Linf), vetor(p_Linf)
         
         CLOSE(41)
        
      END SUBROUTINE escreve_comparacao

       

END MODULE ClassVelocidadeV
