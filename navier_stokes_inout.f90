MODULE NavierStokesInOut

    IMPLICIT NONE
   
    CONTAINS
 
    SUBROUTINE read_configuracoes_gerais
      
        USE, intrinsic :: ISO_FORTRAN_ENV
        USE VariaveisGerais
     
        IMPLICIT NONE

        LOGICAL :: fim_do_arquivo = .FALSE.
        CHARACTER(LEN=100) :: linha
        INTEGER :: ReadStatus, linha_comprimento_real, fim, argumento
        
        OPEN(1, FILE='configuracoes_gerais.txt', STATUS='old', ACTION='read')
            
        DO WHILE (fim_do_arquivo .NEQV. .TRUE.)

            READ (1, '(A100)', IOSTAT=ReadStatus ) linha
         
            ! Calcula o comprimento real da linha
            linha_comprimento_real = len(TRIM(ADJUSTL(linha)))
         
            ! Encontra a referencia da variavel
            DO fim = 1, linha_comprimento_real
            
                ! Lendo a variavel Nx
                IF ( linha(1:fim) .eq. "Nx" ) THEN
                    READ(linha(fim+1:), *) Nx
                    EXIT
                END IF
            
                ! Lendo a variavel Ny
                IF ( linha(1:fim) .eq. "Ny" ) THEN
                    READ(linha(fim+1:), *) Ny
                    EXIT
                END IF
                
                ! Lendo a variavel Xf
                IF ( linha(1:fim) .eq. "Xf" ) THEN
                    READ(linha(fim+1:), *) Xf
                    EXIT
                END IF
            
                ! Lendo a variavel Xi
                IF ( linha(1:fim) .eq. "Xi" ) THEN
                    READ(linha(fim+1:), *) Xi
                    EXIT
                END IF
            
                ! Lendo a variavel Yf
                IF ( linha(1:fim) .eq. "Yf" ) THEN
                    READ(linha(fim+1:), *) Yf
                    EXIT
                END IF
            
                ! Lendo a variavel Yi
                IF ( linha(1:fim) .eq. "Yi" ) THEN
                    READ(linha(fim+1:), *) Yi
                    EXIT
                END IF
                
                ! Lendo a variavel Ti
                IF ( linha(1:fim) .eq. "Ti" ) THEN
                    READ(linha(fim+1:), *) Ti
                    EXIT
                END IF
                
                ! Lendo a variavel Tf
                IF ( linha(1:fim) .eq. "Tf" ) THEN
                    READ(linha(fim+1:), *) Tf
                    EXIT
                END IF
                
                ! Lendo a variavel CFL
                IF ( linha(1:fim) .eq. "CFL" ) THEN
                    READ(linha(fim+1:), *) CFL
                    EXIT
                END IF
                
                ! Lendo a variavel Re
                IF ( linha(1:fim) .eq. "Re" ) THEN
                    READ(linha(fim+1:), *) Re
                    EXIT
                END IF
                
                ! Lendo a variavel Experimento
                IF ( linha(1:fim) .eq. "Experimento" ) THEN
                    READ(linha(fim+1:), *) Experimento
                    EXIT
                END IF
                
                ! Lendo a variavel Qui
                IF ( linha(1:fim) .eq. "Qui" ) THEN
                    READ(linha(fim+1:), *) Qui
                    EXIT
                END IF
                
                ! Lendo a variavel Tn
                IF ( linha(1:fim) .eq. "Tn" ) THEN
                    READ(linha(fim+1:), *) Tn
                    EXIT
                END IF
                
                ! Lendo a variavel extrapolacao_u
                IF ( linha(1:fim) .eq. "Extrapolacao_U" ) THEN
                    READ(linha(fim+1:), *) Extrapolacao_U
                    EXIT
                END IF
                
                ! Lendo a variavel extrapolacao_v
                IF ( linha(1:fim) .eq. "Extrapolacao_V" ) THEN
                    READ(linha(fim+1:), *) Extrapolacao_V
                    EXIT
                END IF
                
                ! Lendo Metodos
                IF ( linha(1:fim) .eq. "Metodo_P" ) THEN
                    READ(linha(fim+1:), *) Metodo_P
                    EXIT
                END IF 

                IF ( linha(1:fim) .eq. "Metodo_U" ) THEN
                    READ(linha(fim+1:), *) Metodo_U
                    EXIT
                END IF 

                IF ( linha(1:fim) .eq. "Metodo_V" ) THEN
                    READ(linha(fim+1:), *) Metodo_V
                    EXIT
                END IF
                
                ! Lendo Tolerancias
                IF ( linha(1:fim) .eq. "Tolerancia_P" ) THEN
                    READ(linha(fim+1:), *) Tolerancia_P
                    EXIT
                END IF 

                IF ( linha(1:fim) .eq. "Tolerancia_U" ) THEN
                    READ(linha(fim+1:), *) Tolerancia_U
                    EXIT
                END IF 
                
                IF ( linha(1:fim) .eq. "Tolerancia_V" ) THEN
                    READ(linha(fim+1:), *) Tolerancia_V
                    EXIT
                END IF 
            
                IF ( linha(1:fim) .eq. "Threads" ) THEN
                    READ(linha(fim+1:), *) Threads
                    EXIT
                END IF 
                

                IF ( linha(1:fim) .eq. "Iteracoes_Externas" ) THEN
                    READ(linha(fim+1:), *) iteracoes_externas
                    EXIT
                END IF 

                IF ( linha(1:fim) .eq. "Iteracoes_Solver" ) THEN
                    READ(linha(fim+1:), *) ite_sol
                    EXIT
                END IF 

                !----------
                 IF ( linha(1:fim) .eq. "unit_arq_geral" ) THEN
                    READ(linha(fim+1:), *) unit_arq_geral
                    EXIT
                END IF 
                
                IF ( linha(1:fim) .eq. "arq_geral" ) THEN
                    READ(linha(fim+1:), '(a)') arq_geral
                    EXIT
                END IF 
                !---------
                    
            END DO
            
            IF ( ReadStatus /= 0 ) THEN
                IF ( ReadStatus == IOSTAT_END ) then
                    fim_do_arquivo = .TRUE. 
                ELSE
                    WRITE ( *, '( / "Error on read: ", I0 )' )  ReadStatus
                    fim_do_arquivo = .TRUE.
                END IF
            END IF
            
        END DO 
      
        CLOSE(1)
        
    END SUBROUTINE read_configuracoes_gerais
    
    SUBROUTINE escreve_cabecalho
    
        USE VariaveisGerais
        
        CHARACTER(len=20) :: parte1
        CHARACTER(len=100) :: parte2
        
        Metodos(1) = "MG_GS       => Multigrid com solver Gauss-Seidel"
        Metodos(2) = "CG          => Gradiente Conjugado"
        Metodos(3) = "PCG_MG_GS   => Gradiente Conjugado Precondicionado com MG_GS"
        Metodos(4) = "ILU         => Metodo ILU"
        Metodos(5) = "MG_ILU      => Multigrid com solver ILU"
        Metodos(6) = "PCG_MG_ILU  => Gradiente Conjugado Precondicionado com MG_ILU"
       
        Metodos_Siglas(1) = "MG_GS"
        Metodos_Siglas(2) = "CG"
        Metodos_Siglas(3) = "PCG_MG_GS"
        Metodos_Siglas(4) = "ILU"
        Metodos_Siglas(5) = "MG_ILU"
        Metodos_Siglas(6) = "PCG_MG_ILU"
       
        1 FORMAT( 1X, A13, 1X, I10, 3X, A100 )
        2 FORMAT( 1X, A13, 1X, 1PE10.1, 3X, A100 )
        
        OPEN(UNIT=unit_arq_geral, FILE=TRIM(ADJUSTL(arq_geral)), STATUS='REPLACE', ACTION='WRITE' )
        
        parte1 = "Threads"
        parte2 = "Numero de Threads utilizadas"
        WRITE(unit_arq_geral, 1) parte1, Threads, parte2
       
        parte1 = "Tn"
        parte2 = "Numero de passos de tempo (ATENCAO: apenas para testes)"
        WRITE(unit_arq_geral, 1) parte1, Tn, parte2

        parte1 = "iteracoes_externas"
        parte2 = "Iteracoes Externas do algoritmo"
        WRITE(unit_arq_geral, 1) parte1, iteracoes_externas, parte2

        parte1 = "ite_sol"
        parte2 = "Iteracoes do Solver"
        WRITE(unit_arq_geral, 1) parte1, ite_sol, parte2
        
        parte1 = "Nx"
        parte2 = "Numero de volumes na direcao x (SEM FICTICIOS)"
        WRITE(unit_arq_geral, 1) parte1, Nx, parte2
        
        parte1 = "Ny"
        parte2 = "Numero de volumes na direcao y (SEM FICTICIOS)"
        WRITE(unit_arq_geral, 1) parte1, Ny, parte2
        
        parte1 = "Metodo para P"
        WRITE(unit_arq_geral, 1) parte1, Metodo_P, Metodos(Metodo_P)
        
        parte1 = "Metodo para U"
        WRITE(unit_arq_geral, 1) parte1, Metodo_U, Metodos(Metodo_U)
        
        parte1 = "Metodo para V"
        WRITE(unit_arq_geral, 1) parte1, Metodo_V, Metodos(Metodo_V)
        
        parte1 = "Tolerancia P"
        parte2 = "Tolerancia utilizada no criterio de parada da solucao de P"
        WRITE(unit_arq_geral, 2) parte1, Tolerancia_P, parte2
        
        parte1 = "Tolerancia U"
        parte2 = "Tolerancia utilizada no criterio de parada da solucao de U"
        WRITE(unit_arq_geral, 2) parte1, Tolerancia_U, parte2
        
        parte1 = "Tolerancia V"
        parte2 = "Tolerancia utilizada no criterio de parada da solucao de V"
        WRITE(unit_arq_geral, 2) parte1, Tolerancia_V, parte2
        
        parte1 = "Xf"
        parte2 = "Xf ponto x do fim do dominio"
        WRITE(unit_arq_geral, 2) parte1,  Xf, parte2
        
        parte1 = "Xi"
        parte2 = "Xi ponto x do fim do dominio"
        WRITE(unit_arq_geral, 2) parte1,  Xi, parte2
        
        parte1 = "Yf"
        parte2 = "Yf ponto y do fim do dominio"
        WRITE(unit_arq_geral, 2) parte1,  Yf, parte2
        
        parte1 = "Yi"
        parte2 = "Yi ponto Y do fim do dominio"
        WRITE(unit_arq_geral, 2) parte1,  Yi, parte2
        
        parte1 = "Tf"
        parte2 = "Tempo em que se quer conhecer a solucao de u, v e p"
        WRITE(unit_arq_geral, 2) parte1,  Tf, parte2
        
        parte1 = "Ti"
        parte2 = "Tempo inicial da simulacao"
        WRITE(unit_arq_geral, 2) parte1,  Ti, parte2
        
        parte1 = "CFL"
        parte2 = "Condicao CFL (garante convergencia com menos passos de tempo)"
        WRITE(unit_arq_geral, 2) parte1, CFL, parte2
        
        parte1 = "Re"
        parte2 = "Numero de Reynolds"
        WRITE(unit_arq_geral, 2) parte1, Re, parte2
               
        parte1 = "Qui"
        parte2 = "Qui=1=>Rotacional, Qui=0=>Padrao"
        WRITE(unit_arq_geral, 2) parte1, Qui, parte2

        parte1 = "Experimento"
        parte2 = "Experimento sendo realizado"
        WRITE(unit_arq_geral, 1) parte1, Experimento, parte2
        
        WRITE(unit_arq_geral, *) "Experimento=5=> Pearson, Experimento=4=> Cavidade"
                 
        CLOSE(unit_arq_geral)
        
    END SUBROUTINE escreve_cabecalho

END MODULE NavierStokesInOut
