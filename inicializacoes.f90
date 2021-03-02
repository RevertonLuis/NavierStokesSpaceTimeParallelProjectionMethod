MODULE Inicializacoes
   
   IMPLICIT NONE
   
   CONTAINS
      
   SUBROUTINE inicializa_variaveis
      
      USE VariaveisGerais
      USE NavierStokesInOut
      
      INTEGER :: i, j 
      
      DOUBLE PRECISION :: pi, Tnt, ht1
      
      ! Lendo arquivo de configuracoes gerais
      CALL read_configuracoes_gerais
            
      pi=4.D0*DATAN(1.D0)
      Xf = pi/4.
      Xi = -pi/4.
      Yf = pi/4.
      Yi = -pi/4.
      
      Lx = Xf - Xi
      Ly = Yf - Yi
      
      hx = Lx/Nx
      hy = Ly/Ny
      
      !------------------ Tratamento da CFL 
      ht1 = CFL/(1.0d0/hx + 1.0d0/hy)
      
      !Tn => number of iterations
      Tnt = NINT( Tf/ht1 )
      ! Para este ht e garantido que vai convergir, mas e muito lento
      ht = (Tf - Ti)/Tnt
            
      !ht = (Tf - Ti)*ht1/Tf
      IF ( Experimento .eq. 5 ) THEN
         ht = ((Tf - Ti)/(2.0d0*Tf))*(Re/100.0d0)*(32.0d0/Nx)*26*hx
      ELSE
         ht = ((Tf - Ti)/(2.0d0*Tf))*(100.0d0/(2.0d0*Re))*(32.0d0/Nx)*16*hx          
      END IF
      
      ! O tamanho dos vetores xp e yp Nx e Ny respectivamente
      IF ( .NOT. ALLOCATED( xu ) )     ALLOCATE( xu ( Nx - 1 ) )
      IF ( .NOT. ALLOCATED( yu ) )     ALLOCATE( yu ( Ny ) )
      IF ( .NOT. ALLOCATED( xv ) )     ALLOCATE( xv ( Nx ) )
      IF ( .NOT. ALLOCATED( yv ) )     ALLOCATE( yv ( Ny - 1 ) )
      
      xv(1) = hx/2.0d0                      + Xi 
      DO i = 1, Nx-1
         
         xu(i)   = hx*i                     + Xi
         xv(i+1) = hx*i + hx/2.0d0          + Xi
         
      END DO
     
       
      yu(1) = hy/2.0d0                        + Yi
      DO j = 1, Ny-1
      
         yu(j+1) = hy*j + hy/2.0d0            + Yi
         yv(j)   = hy*j                       + Yi
         
      END DO      
 
      ! Criando o Alias das funcoes pressao_analitica, u_analitica, v_analitica, condicoes de contorno e condicoes iniciais
      CALL CriaExperimentoAlias(Experimento, Extrapolacao_U, Extrapolacao_V)

      CALL atualiza_variaveis_paralelismo()
      
      CALL u % init( Nx, Ny, Tn, Lx, Ly, ht, Re, Tolerancia_U) 
      CALL u % aplica_condicao_inicial()
       
      CALL v % init( Nx, Ny, Tn, Lx, Ly, ht, Re, Tolerancia_V) 
      CALL v % aplica_condicao_inicial()
      
      CALL p % init( Nx, Ny, Tn, Lx, Ly, ht, Re, Qui, Tolerancia_P, Threads) 
      CALL p % aplica_condicao_inicial()
      CALL p % normaliza_pressao() 
   
      CALL escreve_cabecalho
      
   END SUBROUTINE inicializa_variaveis

   SUBROUTINE CriaExperimentoAlias(Experimento, Extrapolacao_U, Extrapolacao_V)
      
      USE FuncoesAlias
      USE ExperimentosNumericos
      USE ExtrapolacoesDeUeV
   
      INTEGER, INTENT(IN) :: Experimento, Extrapolacao_U, Extrapolacao_V
      
    
      ! Aplicando as condicoes de contorno e calculando as solucoes analiticas para o problema estudado
      SELECT CASE (Experimento)
      
         CASE (5)
            
            cc_u_w => CC_U_W_Experimento5
            cc_u_e => CC_U_E_Experimento5
            cc_u_s => CC_U_S_Experimento5
            cc_u_n => CC_U_N_Experimento5
            
            cc_v_w => CC_V_W_Experimento5
            cc_v_e => CC_V_E_Experimento5
            cc_v_s => CC_V_S_Experimento5
            cc_v_n => CC_V_N_Experimento5
           
            pressao_analitica => PressaoAnaliticaExperimento5
            u_analitica => VelocidadeU_AnaliticaExperimento5
            v_analitica => VelocidadeV_AnaliticaExperimento5

            pressao_inicial   => PressaoInicialExperimento5
            u_inicial         => VelocidadeU_InicialExperimento5
            v_inicial         => VelocidadeV_InicialExperimento5
            
      END SELECT
      
      extrapola_u_S => Extrapola_u_S_Cubica
      extrapola_u_N => Extrapola_u_N_Cubica
      extrapola_v_W => Extrapola_v_W_Cubica
      extrapola_v_E => Extrapola_v_E_Cubica
   !
   !   
   !    SELECT CASE (Extrapolacao)
   !   
   !      CASE (1)
   !         
   !          extrapola_u_S => Extrapola_u_S_Linear
   !          extrapola_u_N => Extrapola_u_N_Linear
   !          extrapola_v_W => Extrapola_v_W_Linear
   !          extrapola_v_E => Extrapola_v_E_Linear
   !          
   !      CASE (2)
   !         
   !          extrapola_u_S => Extrapola_u_S_Cubica
   !          extrapola_u_N => Extrapola_u_N_Cubica
   !          extrapola_v_W => Extrapola_v_W_Cubica
   !          extrapola_v_E => Extrapola_v_E_Cubica
   !                       
   !   END SELECT
        
   END SUBROUTINE CriaExperimentoAlias
   
   
   SUBROUTINE atualiza_variaveis_paralelismo()
   
       USE OMP_LIB
       USE VariaveisGerais
   
       INTEGER :: Tne
       INTEGER :: threads_max, C, tc, it, ft, ptt, thread
       
       ! Atualizando variaveis relacionadas ao paralelismo 
       ! Obtendo o numero maximo de threads do sistema
       threads_max = OMP_get_max_threads()
      
       ! Comparando com o numero de threads que o usuario quer utilizar
       IF (Threads > threads_max) THEN
          WRITE(*,*) ""
          WRITE(*,*) "ATENCAO: numero de threads fornecidas excede o numero de threads disponivel no sistema"
          WRITE(*,*) "Utilizando o numero maximo de threads disponivel:", threads_max
          WRITE(*,*) "Consulte o arquivo " // TRIM(ADJUSTL(arq_geral)) // " para ver os parametros utilizados pelo programa" 
          WRITE(*,*) ""
          Threads = threads_max
       END IF
       
       ALLOCATE( threads_if( Threads, 2) )
       
       Tne = Tn - 1
       
       !------------------------ DISSERTACAO PAGINAS 77-79 ----------------------
       ! Estimando incorretamente o numero de passos de tempo por thread
       ptt = FLOOR(DBLE(Tne)/DBLE(Threads))
       
       ! tc e a thread de referencia
       tc = Threads + 1 - (Tne - Threads * ptt)
       
       ! Calculando o inicio e fim de cada thread, isto e, quais passos de tempo cada thread vai assumir
       DO thread = 1, Threads
           IF (thread >= tc) THEN
               C = 1    
           ELSE 
               C = 0
           END IF
           
           it = (thread -1) * ptt + 1 + C * (thread - tc)
           ft = thread * ptt + C * (thread - tc + 1)
           
           ! Salvando estas informacoes
           threads_if(thread, 1) = it
           threads_if(thread, 2) = ft
           
       END DO
       !------------------------ DISSERTACAO PAGINAS 77-79 ----------------------
          
   END SUBROUTINE atualiza_variaveis_paralelismo

END MODULE Inicializacoes
