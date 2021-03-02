PROGRAM NavierStokes
   
    ! Para usar o openmp
    USE OMP_LIB
    USE VariaveisGerais
    USE Inicializacoes
    USE MG_GS_U
    USE MG_GS_V
    
    IMPLICIT NONE
    
    INTEGER :: t, ite_ext, ite_s, ite_uvp, flag_invalido, ite_ud
    INTEGER :: thread, thread_unit
    CHARACTER(len=100) :: thread_str, ite_sol_str, ite_ext_str, Nx_str, threads_str
    
    ! Variaveis para benchmark
    ! Tempos inicial e final de processamento
    DOUBLE PRECISION :: t1, t2, tt1, tt2, ttt1, ttt2
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: normas_L2_mg_u, normas_L2_mg_v, normas_L2_mg_p
    
    INTEGER :: ciclo, level, Levels
    
    11 FORMAT( 1X, "Multigrid U => Ciclo:", 1X, I3, 3X, "Parada:", 1X, 1PE16.10, 3X, "Norma L2:", 1X,  1PE16.10, 3X, "(Norma L2)**(1/ciclo)", 1X, 1PE16.10)
    12 FORMAT( 1X, "Multigrid V => Ciclo:", 1X, I3, 3X, "Parada:", 1X, 1PE16.10, 3X, "Norma L2:", 1X,  1PE16.10, 3X, "(Norma L2)**(1/ciclo)", 1X, 1PE16.10)
    
    CALL inicializa_variaveis
    
    Levels = LOG( DBLE(Nx))/LOG( DBLE(2) ) - 1
    flag_invalido = 0
    ite_ud = 6
    
    WRITE(ite_ext_str, "(I12)") iteracoes_externas
    WRITE(ite_sol_str, "(I12)") ite_sol
    WRITE(Nx_str, "(I12)") Nx
    WRITE(threads_str, "(I12)") Threads
    
    ALLOCATE(normas_L2_u(iteracoes_externas), &
             normas_L2_v(iteracoes_externas), &
             normas_L2_p(iteracoes_externas))  
             
    ALLOCATE(normas_L2_mg_u(ite_sol), &
             normas_L2_mg_v(ite_sol), &
             normas_L2_mg_p(ite_sol))
                 
    ! Calculando a integral de p para garantir unicidade da solucao da pressao
    !CALL p % calcula_integral_da_pressao(1)
    
    ! Obtendo o tempo inicial
    t1 = OMP_GET_WTIME()
    
    ! Forcando o sistema a utilizar o numero de threads fornecido pelo usuario
    CALL OMP_SET_NUM_THREADS(Threads)
    
    !------------PARALELISMO COMECA AQUI---------------------------------------
    !$OMP PARALLEL PRIVATE(thread, ite_ext, t, thread_unit, thread_str, tt1, tt2, ttt1, ttt2, ite_s, ciclo, level) 
    
       ! Descobrindo o numero id do processador (0-threads-1)
       thread = OMP_get_thread_num() + 1 ! +1 para que a id do processador 0 seja 1
    
       !WRITE(thread_str, "(I12)") thread
       !thread_unit = thread * 100 
       !OPEN(UNIT=thread_unit, FILE="./resultados/log_thread_" // TRIM(ADJUSTL(thread_str)) // ".txt", STATUS='REPLACE', ACTION='WRITE' )
       
       DO ite_ext = 1, iteracoes_externas
 
           !WRITE(thread_unit,*) "##############################################################"
           !WRITE(thread_unit,*) "Iteracao Externa:", ite_ext 
           !WRITE(thread_unit,*) "##############################################################"
    
           ttt1 = OMP_GET_WTIME() 

           DO t = threads_if( thread, 1), threads_if( thread, 2)
           
               ! Tempo de processamento dos termos fontes
               tt1 = OMP_GET_WTIME() 
    
               ! Calculando os termos fontes
               IF ( t .eq. 1 ) THEN
                   CALL u % calcula_termo_fonte(1, 1, v % v, p % p, p % q)
                   CALL v % calcula_termo_fonte(1, 1, u % u, p % p, p % q)
               ELSE
                   CALL u % calcula_termo_fonte(t, t-1, v % v, p % p, p % q)
                   CALL v % calcula_termo_fonte(t, t-1, u % u, p % p, p % q)
               END IF

               ! Tempo de processamento dos termos fontes
               tt2 = OMP_GET_WTIME() 
     
               !WRITE(thread_unit,*) "t =", t, "Fontes UV:", tt2 - tt1 

           END DO
           !$OMP BARRIER

           !---------------------- MULTIGRID U ----------------------
           DO t = threads_if( thread, 1), threads_if( thread, 2)
               CALL atualiza_uo( t, u % u(:, t) )                   
           END DO
           !$OMP BARRIER  
           
           DO ciclo = 1, ite_sol
               DO t = threads_if( thread, 1), threads_if( thread, 2)
                   CALL atualiza_bp_u(t, u % u_bp(:, t))
               END DO
               !$OMP BARRIER  
               
               DO t = threads_if( thread, 1), threads_if( thread, 2)
                   CALL calcula_residuo_u(t, 1)
               END DO
               !$OMP BARRIER  
               
               DO t = threads_if( thread, 1), threads_if( thread, 2)
                   CALL calcula_norma_L2_u(ciclo, t+1)
               END DO
               !$OMP BARRIER  
               
               !!$OMP SINGLE
               !! Somando o valor da norma L2
               !!CALL soma_norma_L2_u(ciclo, normas_L2_mg_u)
               !!WRITE(*, 11) ciclo, normas_L2_mg_u(ciclo)/normas_L2_mg_u(1), normas_L2_mg_u(ciclo), (normas_L2_mg_u(ciclo)/normas_L2_mg_u(1))**(1.0d0/ciclo)
               !!$OMP END SINGLE 
                     
               DO level = 1, Levels - 1
               
                   ! So atualiza bp = residuo, ou seja, apenas para level > 1 (o primeiro e bp = fonte)
                   IF ( level > 1 ) THEN

                       DO t = threads_if( thread, 1), threads_if( thread, 2)
                           CALL atualiza_bp_u2(t, level)
                       END DO
                       !$OMP BARRIER     
                   END IF
                   
                   ! Gauss-Seidel RB aqui
                   DO ite_s = 1, ite_ud
                       DO t = threads_if( thread, 1), threads_if( thread, 2)
                           ! Resolvendo todos os reds do dominio
                           CALL GS_RB_U(t, level, 1)
                   
                       END DO
                       !$OMP BARRIER
                       
                       DO t = threads_if( thread, 1), threads_if( thread, 2)
                           ! Resolvendo todos os blacks do dominio
                           CALL GS_RB_U(t, level, -1)
                       END DO        
                       !$OMP BARRIER  
                   END DO
               
                   DO t = threads_if( thread, 1), threads_if( thread, 2)
                       CALL calcula_residuo_u(t, level)
                   END DO     
                   !$OMP BARRIER                                   
                   DO t = threads_if( thread, 1), threads_if( thread, 2)
                       CALL restricao_u ( t+1, level )
                   END DO                                      
                   !$OMP BARRIER  
               END DO
                              
               ! Falta conta no ultimo level
               
               DO level = Levels, 2, -1

                   DO t = threads_if( thread, 1), threads_if( thread, 2)
                       CALL prolongacao_u ( t+1, level )
                   END DO
                   !$OMP BARRIER  
                   
                   ! Limpando o vetor para evitar erros no proximo ciclo
                   !DO t = threads_if( thread, 1), threads_if( thread, 2)
                   !   CALL limpa_u( t, level )
                   !END DO                                
                   !$OMP BARRIER  
                   
                   ! Gauss-Seidel RB aqui
                   DO ite_s = 1, ite_ud
                       DO t = threads_if( thread, 1), threads_if( thread, 2)
                           ! Resolvendo todos os reds do dominio
                           CALL GS_RB_U(t, level-1, 1)
                   
                       END DO
                       !$OMP BARRIER
                       
                       DO t = threads_if( thread, 1), threads_if( thread, 2)
                           ! Resolvendo todos os blacks do dominio
                           CALL GS_RB_U(t, level-1, -1)
                       END DO        
                       !$OMP BARRIER  
                   END DO
                                                       
               END DO              
                        
           END DO
           
           DO t = threads_if( thread, 1), threads_if( thread, 2)
               CALL atualiza_u2( t+1, u % u(:, t+1) )                   
           END DO
           !$OMP BARRIER  
           !---------------------- MULTIGRID U ----------------------
             
           !!$OMP SINGLE         
           !WRITE(*,*) ""
           !!$OMP END SINGLE         
           
           !---------------------- MULTIGRID V ----------------------
           DO t = threads_if( thread, 1), threads_if( thread, 2)
               CALL atualiza_vo( t, v % v(:, t) )                   
           END DO
           !$OMP BARRIER  
           
           DO ciclo = 1, ite_sol
               DO t = threads_if( thread, 1), threads_if( thread, 2)
                   CALL atualiza_bp_v(t, v % v_bp(:, t))
               END DO
               !$OMP BARRIER  
               
               DO t = threads_if( thread, 1), threads_if( thread, 2)
                   CALL calcula_residuo_v(t, 1)
               END DO
               !$OMP BARRIER  
               
               DO t = threads_if( thread, 1), threads_if( thread, 2)
                   CALL calcula_norma_L2_v(ciclo, t+1)
               END DO
               !$OMP BARRIER  
               
               !!$OMP SINGLE
               !! Somando o valor da norma L2
               !!CALL soma_norma_L2_v(ciclo, normas_L2_mg_v)
               !!WRITE(*, 12) ciclo, normas_L2_mg_v(ciclo)/normas_L2_mg_v(1), normas_L2_mg_v(ciclo), (normas_L2_mg_v(ciclo)/normas_L2_mg_v(1))**(1.0d0/ciclo)
               !!$OMP END SINGLE 
               
               DO level = 1, Levels - 1
               
                   ! So atualiza bp = residuo, ou seja, apenas para level > 1 (o primeiro e bp = fonte)
                   IF ( level > 1 ) THEN

                       DO t = threads_if( thread, 1), threads_if( thread, 2)
                           CALL atualiza_bp_v2(t, level)
                       END DO
                       !$OMP BARRIER     
                   END IF
                   
                   ! Gauss-Seidel RB aqui
                   DO ite_s = 1, ite_ud
                       DO t = threads_if( thread, 1), threads_if( thread, 2)
                           ! Resolvendo todos os reds do dominio
                           CALL GS_RB_V(t, level, 1)
                   
                       END DO
                       !$OMP BARRIER
                       
                       DO t = threads_if( thread, 1), threads_if( thread, 2)
                           ! Resolvendo todos os blacks do dominio
                           CALL GS_RB_V(t, level, -1)
                       END DO        
                       !$OMP BARRIER  
                   END DO
               
                   DO t = threads_if( thread, 1), threads_if( thread, 2)
                       CALL calcula_residuo_v(t, level)
                   END DO     
                   !$OMP BARRIER                                   
                   DO t = threads_if( thread, 1), threads_if( thread, 2)
                       CALL restricao_v ( t+1, level )
                   END DO                                      
                   !$OMP BARRIER  
               END DO
                              
               ! Falta conta no ultimo level
               
               DO level = Levels, 2, -1

                   DO t = threads_if( thread, 1), threads_if( thread, 2)
                       CALL prolongacao_v ( t+1, level )
                   END DO
                   !$OMP BARRIER  
                   
                   ! Limpando o vetor para evitar erros no proximo ciclo
                   !DO t = threads_if( thread, 1), threads_if( thread, 2)
                   !   CALL limpa_u( t, level )
                   !END DO                                
                   !$OMP BARRIER  
                   
                   ! Gauss-Seidel RB aqui
                   DO ite_s = 1, ite_ud
                       DO t = threads_if( thread, 1), threads_if( thread, 2)
                           ! Resolvendo todos os reds do dominio
                           CALL GS_RB_V(t, level-1, 1)
                   
                       END DO
                       !$OMP BARRIER
                       
                       DO t = threads_if( thread, 1), threads_if( thread, 2)
                           ! Resolvendo todos os blacks do dominio
                           CALL GS_RB_V(t, level-1, -1)
                       END DO        
                       !$OMP BARRIER  
                   END DO
                                                       
               END DO              
                        
           END DO
           
           DO t = threads_if( thread, 1), threads_if( thread, 2)
               CALL atualiza_v2( t+1, v % v(:, t+1) )                   
           END DO
           !$OMP BARRIER  
           !---------------------- MULTIGRID V ----------------------
           
           !write(*,*)  u % u(1:7, 2)
           !write(*,*)  u % u(8:14, 2)
           !write(*,*)  u % u(15:21, 2)
           !write(*,*)  u % u(22:28, 2)
           !write(*,*)  u % u(29:35, 2)
           !write(*,*)  u % u(36:42, 2)
           !write(*,*)  u % u(43:49, 2)
           !write(*,*)  u % u(50:56, 2)
           
           !write(*,*)  v % v(1:8, 2)
           !write(*,*)  v % v(9:16, 2)
           !write(*,*)  v % v(17:24, 2)
           !write(*,*)  v % v(25:32, 2)
           !write(*,*)  v % v(33:40, 2)
           !write(*,*)  v % v(41:48, 2)
           !write(*,*)  v % v(49:56, 2)
           
           !call exit()
           
                     
           ! Calcula o Residuo
           !DO t = threads_if( thread, 1), threads_if( thread, 2)
           !    ! Calculando todos os residuos do dominio
           !    CALL u % calcula_residuo(t)
           !    CALL v % calcula_residuo(t) 
           !    CALL u % calcula_norma_L2(t)
           !    CALL v % calcula_norma_L2(t)                            
           !END DO
           !$OMP BARRIER 

           !$OMP SINGLE
           ! Somando o valor da norma L2
           !normas_L2_u(ite_ext) = DSQRT(SUM(u % normas_L2(:)))
           !normas_L2_v(ite_ext) = DSQRT(SUM(v % normas_L2(:)))
           !WRITE(*,*) "Parada U:", normas_L2_u(ite_ext)/normas_L2_u(1), "Norma L2 de U:", normas_L2_u(ite_ext), ite_ext, MAXVAL(ABS(u % u(:, :) - u % u_ex(:, :)))
           !WRITE(*,*) "Parada V:", normas_L2_v(ite_ext)/normas_L2_v(1), "Norma L2 de V:", normas_L2_v(ite_ext), ite_ext, MAXVAL(ABS(v % v(:, :) - v % v_ex(:, :)))
           normas_L2_u(ite_ext) = maxval(abs(u % u(:, Tn) - u % u_ref(:) ))
           normas_L2_v(ite_ext) = maxval(abs(v % v(:, Tn) - v % v_ref(:) ))
           u % u_ref(:) = u % u(:, Tn)
           v % v_ref(:) = v % v(:, Tn)
           
           !$OMP END SINGLE 
           
           
           tt2 = OMP_GET_WTIME() 
           !WRITE(thread_unit, *) "Solver UV:", tt2 - tt1 
           
           DO t = threads_if( thread, 1), threads_if( thread, 2)
               
               tt1 = OMP_GET_WTIME() 
               CALL p % calcula_termo_fonte(t+1, t, u % u, v % v)
               tt2 = OMP_GET_WTIME()
               !WRITE(thread_unit,*) "t = ", t, "Fontes P:", tt2 - tt1 
               
               tt1 = OMP_GET_WTIME()
               !CALL p % single_grid(t+1, ite_sol*50)
               CALL p % multi_grid(t+1, ite_sol)
               tt2 = OMP_GET_WTIME()
               !WRITE(thread_unit,*) "t = ", t, "Solver P:", tt2 - tt1  
               
               !CALL p % calcula_residuo(t)
               !CALL p % calcula_norma_L2(t)
               !write(*,*) DSQRT((p % normas_L2( t )))
               
               CALL p % corrige_pressao_com_a_integral(t+1)
                
           END DO
           !$OMP BARRIER
           
           !write(*,*)  p % p_bp(1:8, 2)
           !write(*,*)  p % p_bp(9:16, 2)
           !write(*,*)  p % p_bp(17:24, 2)
           !write(*,*)  p % p_bp(25:32, 2)
           !write(*,*)  p % p_bp(33:40, 2)
           !write(*,*)  p % p_bp(41:48, 2)
           !write(*,*)  p % p_bp(49:56, 2)
           !write(*,*)  p % p_bp(57:64, 2)
           !call exit()
           
           !write(*,*)  p % q(1:8, 2)
           !write(*,*)  p % q(9:16, 2)
           !write(*,*)  p % q(17:24, 2)
           !write(*,*)  p % q(25:32, 2)
           !write(*,*)  p % q(33:40, 2)
           !write(*,*)  p % q(41:48, 2)
           !write(*,*)  p % q(49:56, 2)
           !write(*,*)  p % q(57:64, 2)
           !call exit()
           
           tt1 = OMP_GET_WTIME()  
           DO t = 1, Tn-1
               CALL p % atualiza_pressao(t+1, t, thread)
           END DO
           tt2 = OMP_GET_WTIME() 
           !WRITE(thread_unit,*) "Atualizacao de P:", tt2 - tt1  
           
           !write(*,*)  p % p(1:8, 2)
           !write(*,*)  p % p(9:16, 2)
           !write(*,*)  p % p(17:24, 2)
           !write(*,*)  p % p(25:32, 2)
           !write(*,*)  p % p(33:40, 2)
           !write(*,*)  p % p(41:48, 2)
           !write(*,*)  p % p(49:56, 2)
           !write(*,*)  p % p(57:64, 2)
           !call exit()
           
           
           !$OMP BARRIER
           
           !$OMP SINGLE
           ! Somando o valor da norma L2
           !normas_L2_p(ite_ext) = DSQRT((p % normas_L2( threads_if( Threads, 2) )))
           !WRITE(*,*) "Parada P:", normas_L2_p(ite_ext)/normas_L2_p(1), "Norma L2 de P:", normas_L2_p(ite_ext), ite_ext, maxval(abs(p % p(:, Tn) - p % p_ref(:)))
           normas_L2_p(ite_ext) = maxval(abs(p % p(:, Tn) - p % p_ref(:)))
           p % p_ref(:) = p % p(:, Tn)
           ite_uvp = ite_ext
           
           !$OMP END SINGLE 
           
           !DO t = threads_if( thread, 1), threads_if( thread, 2)
           !    WRITE(thread_unit,*) "#---------------------------------------------------------"
           !    WRITE(thread_unit,*) "Pressao, t =", t, "Linf", MAXVAL(ABS(p % p(:, t+1) - p % p_ex(:, t+1)))
           !    WRITE(thread_unit,*) "Veloc u, t =", t, "Linf", MAXVAL(ABS(u % u(:, t+1) - u % u_ex(:, t+1)))
           !    WRITE(thread_unit,*) "Veloc v, t =", t, "Linf", MAXVAL(ABS(v % v(:, t+1) - v % v_ex(:, t+1)))
           !    WRITE(thread_unit,*) "#---------------------------------------------------------"                                  
           !END DO
          
           ttt2 = OMP_GET_WTIME() 
           !WRITE(thread_unit,*) "Processamento total da thread:", ttt2 - ttt1

           IF (normas_L2_u(ite_ext) < Tolerancia_U .and. &
               normas_L2_v(ite_ext) < Tolerancia_V) THEN !.and. normas_L2_p(ite_ext) < Tolerancia_P) THEN
               ite_uvp = ite_ext
               EXIT
           END IF
           
           ! Verificando os casos NaN
           IF (normas_L2_u(ite_ext)/normas_L2_u(1) .ne. normas_L2_u(ite_ext)/normas_L2_u(1) .or. &
               normas_L2_v(ite_ext)/normas_L2_v(1) .ne. normas_L2_v(ite_ext)/normas_L2_v(1) .or. &
               normas_L2_p(ite_ext)/normas_L2_p(1) .ne. normas_L2_p(ite_ext)/normas_L2_p(1)) THEN
               ite_uvp = ite_ext
               WRITE(*,*) "Caso NaN  Nx: " // TRIM(ADJUSTL(Nx_str)) // " ite_solver: " // TRIM(ADJUSTL(ite_sol_str))
               EXIT
           END IF
           
           
       END DO 
    
       !CLOSE(thread_unit)
       
    !$OMP END PARALLEL
    
    ! Obtendo o tempo final
    t2 = OMP_GET_WTIME()
    IF (flag_invalido .eq. 0) THEN
        OPEN(UNIT=1, FILE="./resultados/convergencia_UVP_ite_ext_" // TRIM(ADJUSTL(ite_ext_str)) // &
                                                       "_ite_sol_" // TRIM(ADJUSTL(ite_sol_str)) // &
                                                       "_Nx_Ny_"   // TRIM(ADJUSTL(Nx_str)) // &
                                                       "_Threads_" // TRIM(ADJUSTL(Threads_str)) // &
                                                       ".txt", STATUS='REPLACE', ACTION='WRITE' )
    ELSE
       OPEN(UNIT=1, FILE="./resultados/INVALIDO_convergencia_UVP_ite_ext_" // TRIM(ADJUSTL(ite_ext_str)) // &
                                                               "_ite_sol_" // TRIM(ADJUSTL(ite_sol_str)) // &
                                                               "_Nx_Ny_"   // TRIM(ADJUSTL(Nx_str)) // &
                                                               "_Threads_" // TRIM(ADJUSTL(Threads_str)) // &
                                                               ".txt", STATUS='REPLACE', ACTION='WRITE' )
    END IF
    
    WRITE(1, *) "ite_ext, NormaL2(ite_uvp)/NormaL2(1): U, V, P e NormaL2(ite_uvp): U, V, P"
    DO ite_ext = 1, ite_uvp
        WRITE(1, *) ite_ext, & 
                    normas_L2_u(ite_ext), &
                    normas_L2_v(ite_ext), &
                    normas_L2_p(ite_ext)
    END DO 
    
    WRITE(1,*) "Tempo de processamento TOTAL", t2-t1
    
    CLOSE(1)
    
     CALL u % escreve_comparacao(Tn, "./resultados/resultados_U_ite_ext_" // TRIM(ADJUSTL(ite_ext_str)) // &
                                    "_ite_sol_" // TRIM(ADJUSTL(ite_sol_str)) // &
                                    "_Nx_Ny_"   // TRIM(ADJUSTL(Nx_str)) // &
                                    "_Threads_" // TRIM(ADJUSTL(Threads_str)) // &
                                    ".txt", u%u(:, Tn))
    CALL v % escreve_comparacao(Tn, "./resultados/resultados_V_ite_ext_" // TRIM(ADJUSTL(ite_ext_str)) // &
                                    "_ite_sol_" // TRIM(ADJUSTL(ite_sol_str)) // &
                                    "_Nx_Ny_"   // TRIM(ADJUSTL(Nx_str)) // &
                                    "_Threads_" // TRIM(ADJUSTL(Threads_str)) // &
                                    ".txt", v%v(:, Tn))
    CALL p % escreve_comparacao(Tn, "./resultados/resultados_P_ite_ext_" // TRIM(ADJUSTL(ite_ext_str)) // &
                                    "_ite_sol_" // TRIM(ADJUSTL(ite_sol_str)) // &
                                    "_Nx_Ny_"   // TRIM(ADJUSTL(Nx_str)) // &
                                    "_Threads_" // TRIM(ADJUSTL(Threads_str)) // &
                                    ".txt", p%p(:, Tn))
    
    
END PROGRAM NavierStokes
