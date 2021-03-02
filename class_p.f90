MODULE ClassPressao

   USE ClassArrayFonteP
   USE ClassMatrizA
      
   ! Modulo que implementa a classe para Pressao
   
   IMPLICIT  NONE
  
   PRIVATE
     
   TYPE, PUBLIC :: Pressao
   
      ! Atritubos Externos serao atualizados por outra rotina
      ! Atributos Internos serao atualizados na classe
   
      !---------Atributos Externos----------------------- 
      INTEGER           :: Tn, Threads                                   ! Numero de passos de tempo resolvidos em paralelo
      INTEGER           :: Nx, Ny                                        ! Nx e Ny da malha
      DOUBLE PRECISION  :: Lx, Ly                                        ! Comprimento do dominio
      DOUBLE PRECISION  :: ht                                            ! Refinamento temporal (necessario por causa dos coeficientes)
      DOUBLE PRECISION  :: Re                                            ! Numero de Reynolds (necessario por causa dos coeficientes)
      DOUBLE PRECISION  :: Qui
      DOUBLE PRECISION  :: Tolerancia                                    ! Tolerancia no criterio de parada dos solvers
      !--------------------------------------------------
      
      
      !--------Atributos Internos------------------------
      DOUBLE PRECISION                                 :: integral_de_p  ! A integral da pressao          
      DOUBLE PRECISION                                 :: hx, hy         ! hx e hy da malha
      DOUBLE PRECISION                                 :: C_Normalizacao ! Constante de Normalizacao
      
      ! ATENCAO: DEFINIR UMA VETOR DE VETORES E RESERVAR VARIAVEIS v0, v1, v2, v3, ..., vt dinamicamente (de acordo com o a ordem do metodo de projecao no futuro)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE      :: p, q, p_bp, p_ex, residuo       ! Vetores da pressao
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE         :: p_ref
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: threads_if
      DOUBLE PRECISION, DIMENSION(:),    ALLOCATABLE      :: normas_L2
       

      ! A matriz A de V
      TYPE ( MatrizA ) :: A
      
      TYPE (ArrayFonteP),  DIMENSION(9) :: FonteP
      !--------------------------------------------------

      CONTAINS
     
         ! Metodos da Classe
         procedure :: init                      ! Metodo que inicia a classe
         procedure :: calcula_integral_da_pressao  ! Metodo que calcula a integral da pressao para garantir a unicidade da solucao
         procedure :: corrige_pressao_com_a_integral
         procedure :: aplica_condicao_inicial
         procedure :: normaliza_pressao
         procedure :: calcula_termo_fonte
         procedure :: single_grid
         procedure :: multi_grid
         procedure :: atualiza_pressao
         procedure :: escreve_matriz
         procedure :: escreve_comparacao
         procedure :: calcula_subdominios
         procedure :: calcula_residuo
         procedure :: calcula_norma_L2

   END TYPE Pressao

   CONTAINS

      SUBROUTINE init(self, Nx, Ny, Tn, Lx, Ly, ht, Re, Qui, Tolerancia, Threads)
         
         USE Subrotinas_MG_GS
       
         !------------ Metodo que inicia as variaveis gerais da malha ---------------------

         CLASS( Pressao ), INTENT(INOUT) :: self
         
         INTEGER         , INTENT(IN) :: Nx, Ny, Tn, Threads
         DOUBLE PRECISION, INTENT(IN) :: Lx, Ly, ht, Re, Qui, Tolerancia
         
         INTEGER :: stencil

         ! Atualiando os atributos da classe
         self % Tn = Tn
         self % Lx  = Lx
         self % Ly  = Ly
         self % Re  = Re
         self % Qui = Qui
         self % ht  = ht
         self % Tolerancia = Tolerancia
         self % Threads = Threads
         
         ! Calculando hx e hy (o mesmo valor para u, v e p)
         self % hx = self % Lx / Nx
         self % hy = self % Ly / Ny
         
         self % Nx = Nx 
         self % Ny = Ny 
         
         ! Reservando os vetores p, q
         ALLOCATE( self % p( self % Nx * self % Ny, self % Tn ), &
                   self % q( self % Nx * self % Ny, self % Tn ), &
                   self % p_bp( self % Nx * self % Ny, self % Tn ), &
                   self % p_ex( self % Nx * self % Ny, self % Tn ), &
                   self % threads_if(self % Threads, 2 ), &
                   self % residuo( self % Nx * self % Ny, self % Tn ), &
                   self % normas_L2( self % Tn ), &
                   self % p_ref( self % Nx * self % Ny ) )

         ! Iniciando os vetores
         self % p(:, :) = 0.0d0          
         self % q(:, :)  = 0.0d0
         self % p_bp(:, :)  = 0.0d0
         self % p_ex(:, :) = 0.0d0
         self % p_ref(:) = 0.0d0
         self % residuo(:, :) = 0.0d0   
         self % normas_L2( self % Tn ) = 0.d0       
         
         !-------------------------- matriz A de P --------------------------------
         ! Inicializando a matriz A de U
         CALL self % A % init(self % Nx, self % Ny, &
                              self % hx, self % hy, self % ht, & 
                              self % Re)
                              
         ! Atribuindo um stencil a cada volume de A                     
         CALL self % A % atribui_stencils()
         
         ! Calculandos os coeficientes e posicoes dos stencils
         CALL self % A % stencils_p()
         !-------------------------- matriz A de P --------------------------------     
         
         ! Inicializando o array de fontes de P
          DO stencil = 1, 9
            CALL self % FonteP(stencil) % init(stencil)
         END DO

         CALL calcula_exata(self)
         CALL calcula_subdominios(self)
         
         CALL Atualiza_Variaveis_MG_GS_P(self % Lx, self % Ly, &
                                         self % Nx, self % Ny, &
                                         self % Tolerancia, &
                                         3, 3, &
                                         50, &
                                         self % ht, &
                                         self % Re, self % Tn )
         
      END SUBROUTINE init

      SUBROUTINE calcula_subdominios(self)
      
          CLASS (Pressao), INTENT(INOUT) :: self
          
          INTEGER :: Nsub, tc, C, thread, it, ft
          
          ! Estimando incorretamente o numero de passos de tempo por thread
          Nsub = FLOOR(DBLE(self % Nx * self % Ny)/DBLE(self % Threads))
        
          ! tc e a thread de referencia
          tc = self % Threads + 1 - (self % Nx* self % Ny - self % Threads * Nsub)
       
          ! Calculando o inicio e fim de cada thread, isto e, quais passos de tempo cada thread vai assumir
          DO thread = 1, self % Threads
              IF (thread >= tc) THEN
                  C = 1    
              ELSE 
                  C = 0
              END IF
           
              it = (thread -1) * Nsub + 1 + C * (thread - tc)
              ft = thread * Nsub + C * (thread - tc + 1)
           
              ! Salvando estas informacoes
              self % threads_if(thread, 1) = it
              self % threads_if(thread, 2) = ft
                            
          END DO
          
      END SUBROUTINE calcula_subdominios
      
      SUBROUTINE calcula_integral_da_pressao( self, t )
   
         CLASS (Pressao), INTENT(INOUT) :: self
         
         INTEGER, INTENT(IN) :: t
         
         ! Regra do Retangulo
         self % integral_de_p = self % hx* self % hy * SUM( self % p(:, t) )
      
         WRITE(*,*) "Integral da pressao: ", self % integral_de_p
   
         
      END SUBROUTINE calcula_integral_da_pressao

      
      SUBROUTINE aplica_condicao_inicial( self )
         USE FuncoesAlias

         CLASS (Pressao), INTENT(INOUT) :: self
         
         INTEGER :: volume
         
         DO volume = 1, self % Nx * self % Ny
            
            self % p(volume, :) = pressao_inicial(0, volume)
            
         END DO

      END SUBROUTINE aplica_condicao_inicial
      

      SUBROUTINE normaliza_pressao( self )
         CLASS (Pressao) :: self
         
         !self % C_Normalizacao = SUM(self % p1 * self % hx * self % hy) / SUM((self % p1 * 0.0d0+1.0d0)* self % hx * self %hy)
         
         self % C_Normalizacao = SUM(self % p(:,1)) / (self % Nx * self % Ny)
         
         self % p(:,:) = self % p(:,:) - self % C_Normalizacao
         
      END SUBROUTINE normaliza_pressao
      
      SUBROUTINE calcula_exata( self )
         USE FuncoesAlias

         CLASS (Pressao), INTENT(INOUT) :: self
         
         INTEGER :: volume, t
         DOUBLE PRECISION :: C
         
         DO t = 1, self % Tn
         
             DO volume = 1, self % Nx * self % Ny
            
                 self % p_ex(volume, t) = pressao_analitica(t, volume)
            
             END DO
             
            C = SUM(self % p_ex(:,t)) / (self % Nx * self % Ny)
         
            self % p_ex(:,t) = self % p_ex(:,t) - C
             
         END DO

      END SUBROUTINE calcula_exata

      SUBROUTINE calcula_termo_fonte( self, t, t0, ut, vt)
         
         CLASS(Pressao), INTENT(INOUT) :: self
         
         INTEGER, INTENT(IN) :: t, t0
         DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: ut, vt
         
         INTEGER :: volume, stencil
         
         DO volume = 1, self % Nx * self % Ny
         
            ! Stencil do volume
            stencil = self % A % v(volume)
            !stencil = self % A % calcula_stencil(volume)
         
            self % p_bp(volume, t) = self % FonteP( stencil ) % FonteStencil(t0, &
                                                                             volume, &
                                                                             self % Nx, &
                                                                             self % hx, self % hy, self % ht, &
                                                                             ut(:, t), &
                                                                             vt(:, t))
            self % p_bp(volume, t) = self % p_bp(volume, t)
            !write(*,*) self % p_bp(volume), volume
           
         END DO

      END SUBROUTINE calcula_termo_fonte


      SUBROUTINE single_grid( self, t, iteracoes )
      
         USE GS
   
         CLASS (Pressao) :: self
         
         INTEGER, INTENT(IN) :: t, iteracoes
         
         ! Resolvendo com Gauss Seidel
         CALL GaussSeidel( self % A, self % q(:, t), self % p_bp(:, t), iteracoes )
         
      END SUBROUTINE single_grid

      SUBROUTINE calcula_residuo( self, t )
   
          USE Residuo
        
          CLASS (Pressao) :: self
         
          INTEGER, INTENT(IN) :: t
         
          ! Calcula residuo
          CALL calcula_residuo_P( self % A, self % q, self % p_bp, self % residuo, t )
         
      END SUBROUTINE calcula_residuo 

      SUBROUTINE calcula_norma_L2( self, t )
   
          USE Residuo
        
          CLASS (Pressao) :: self
         
          INTEGER, INTENT(IN) :: t
         
          ! Calcula residuo
          self % normas_L2(t) = SUM(self % residuo(:, t) * self % residuo(:, t))
         
      END SUBROUTINE calcula_norma_L2  

      
      
      SUBROUTINE multi_grid( self, t, ciclos_V )
      
         USE Subrotinas_MG_GS
         
         CLASS (Pressao) :: self
         
         INTEGER, INTENT(IN) :: t, ciclos_V
         
         CALL MG_GS_P( self % q(:, t), self % p_bp(:, t), ciclos_V, t )
         
      END SUBROUTINE multi_grid


      
      
      
      SUBROUTINE corrige_pressao_com_a_integral(self, t)
         
         CLASS (Pressao) :: self
    
         INTEGER, INTENT(IN) :: t
         DOUBLE PRECISION :: integral_q, integral_q_antes, integral_q_depois

         !write(*,*)
         integral_q_antes = self % hx* self % hy * SUM( self % q(:, t) ) ! Calcula integral
         !write(*,*) "Integral de Q ANTES", integral_q_antes
      
         !---------------- Para descobrir o fator da integral (VERIFICAR SE ESTA PARTE NAO E EQUIVALENTE AO CP/C1 NORMALIZACAO)
         self % q(:, t) = self % q(:, t) - integral_q_antes               ! Atualiza vetor com integral
         integral_q_depois = self % hx* self % hy * SUM( self % q(:, t) ) ! Recalcula a integral
         self % q(:, t) = self % q(:, t) + integral_q_antes               ! Remove o valor da integral, ou seja, volta pra q original
         !--------------------------------------------------------------
      
         self % q(:, t) = self % q(:, t) - integral_q_antes/(1.0d0 - integral_q_depois/integral_q_antes)  ! Agora sim esta garantido integral nula
         
         !integral_q = self % hx* self % hy * SUM( self % q(:) ) ! Calcula integral
         !write(*,*) "Integral de Q DEPOIS", integral_q
         
      END SUBROUTINE corrige_pressao_com_a_integral
   

      SUBROUTINE atualiza_pressao( self, t, t0, thread )

         CLASS (Pressao) :: self
         
         INTEGER, INTENT(IN) :: t, t0, thread
         INTEGER :: volume
         
         !self % ps(:) = 2 * self % p1(:) - self % p0(:)
         !self % p0(:) = self % p1(:)
         
         DO volume = self % threads_if(thread, 1), self % threads_if(thread, 2)
         
            self % p(volume, t) = self % p(volume, t0) +  &
                                  self % q(volume, t) + &
                                  self % Qui*(self % ht/(1.0d0 * self % Re)) * self % p_bp(volume, t)
            
         END DO
         
         !self % ps(:) = 2 * self % p1(:) - self % p0(:)
         
      END SUBROUTINE atualiza_pressao


      SUBROUTINE escreve_matriz(self, nome, vetor)
   
         CLASS (Pressao) :: self
      
         CHARACTER(len=*), INTENT(IN) :: nome
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vetor
         
         INTEGER :: pp, i, j
         CHARACTER(len=(self % Nx*23)) :: linha
         CHARACTER(len=23) :: valor
      
         OPEN(UNIT=41, FILE=TRIM(ADJUSTL(nome)), STATUS='REPLACE', ACTION='WRITE')
      
         ! Para pontos internos 
         DO j = 1, self % Ny
            DO i = 1, self % Nx
               
               pp = i + (j-1)* self % Nx

               WRITE(valor, "(1PE22.15)") vetor(pp)
               linha((i-1)*23+1:(i-1)*23+23) = valor
             
            END DO
            
            WRITE(41,1) linha

         END DO         
         CLOSE(41)
                
         1 FORMAT( A )
      
      END SUBROUTINE escreve_matriz

      SUBROUTINE escreve_comparacao(self, tempo, nome, vetor)

         USE FuncoesAlias

         CLASS (Pressao) :: self
      
         CHARACTER(len=*), INTENT(IN) :: nome
         INTEGER, INTENT(IN) :: tempo
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vetor
         
         INTEGER :: pp, i, j, p_Linf
         DOUBLE PRECISION :: C1, CP, Linf
       
         1 FORMAT( 1PE22.15, 3X, 1PE22.15, 3X, 1PE22.15, 3X, I10, 3X, I10 )
         2 FORMAT( 1PE22.15, 1X, "= ht", /,  &
                   1PE22.15, 1X, "= Erro Maximo", /,  &
                   1PE22.15, 1X, "= Analitica (onde o erro aconteceu)", /, &
                   1PE22.15, 1X, "=  Numerica (onde o erro aconteceu)" )


        ! Obtendo as constantes de normalizacao
        CP = 0.0d0
        C1 = 0.0d0
        DO pp = 1, self % Nx * self % Ny

            CP = CP +  pressao_analitica(tempo, pp) * self % hx * self % hy
            C1 = C1 + self % hx * self % hy

        END DO

  
         OPEN(UNIT=41, FILE=TRIM(ADJUSTL(nome)), STATUS='REPLACE', ACTION='WRITE')
   
         Linf   = 0.0d0
         p_Linf  = 0
         
         ! Para pontos internos 
         DO j = 1, self % Ny
            DO i = 1, self % Nx
               
               pp = i + (j-1)* self % Nx

               !WRITE(41, 1) pressao_analitica(tempo, pp) - CP/C1, vetor(pp), pressao_analitica(tempo, pp) - CP/C1 - vetor(pp)
                       
               IF ( ABS(pressao_analitica(tempo, pp) - CP/C1 - vetor(pp)) > Linf ) THEN
                  Linf = ABS(pressao_analitica(tempo, pp) - CP/C1 - vetor(pp))
                  p_Linf  = pp
               END IF

            END DO
            
         END DO     

         WRITE(41, 2) self % ht, Linf, pressao_analitica(tempo, p_Linf) - CP/C1, vetor(p_Linf)
         
         CLOSE(41)
        
      END SUBROUTINE escreve_comparacao

       

END MODULE ClassPressao
