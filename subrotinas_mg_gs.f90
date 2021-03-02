MODULE Subrotinas_MG_GS
   
   CONTAINS
   
   SUBROUTINE Atualiza_Variaveis_MG_GS_P( Lx_o, Ly_o, &
                                          Nx_o, Ny_o, &
                                          tol_o, &
                                          ite_down_o, ite_up_o, &
                                          max_ciclos_V_o, ht, Re, Tn_o )

      USE Variaveis_Solvers_P
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: Nx_o, Ny_o, ite_down_o, ite_up_o, max_ciclos_V_o, Tn_o
      DOUBLE PRECISION, INTENT(IN) :: Lx_o, Ly_o, tol_o, ht, Re
      
      INTEGER :: level
      DOUBLE PRECISION :: hx, hy
      
      Lx = Lx_o
      Ly = Ly_o
      Nx = Nx_o
      Ny = Ny_o
      tol = tol_o
      ite_down = ite_down_o
      ite_up   = ite_up_o
      max_ciclos_V = max_ciclos_V_o
      Tn = Tn_o
      
      ! Calculando o numero de niveis
      Levels = LOG( DBLE(Nx))/LOG( DBLE(2) )

      ! lv_if => duas dimensoes 1: inicio da malha, 2: fim da malha
      ALLOCATE( lv_if(Levels, 2), lv_nx(Levels), lv_ny(Levels), normas_L2(max_ciclos_V, Tn) )

      ! Inicio da primeira malha comeca em 1
      lv_if(1,1) = 1
      
      ! Fim da primeira malha termina em Nx*Ny
      lv_if(1,2) = Nx*Ny

      ! Para cada nivel de malha
      DO level = 1, Levels
         
         ! Calcula o Nx e Ny de cada nivel de malha
         lv_nx(level) = Nx/(2**(level-1))
         lv_ny(level) = Ny/(2**(level-1))

         ! O inicio e fim da primeira malha ja foram calculados
         IF ( level > 1 ) THEN

            ! Calculando o indice que define o fim do vetor
            lv_if(level,2) = lv_if(level-1,2) + lv_nx(level)*lv_ny(level)

            ! Calculando o indice que define o inicio do vetor
            lv_if(level,1) = lv_if(level-1,2) + 1

         END IF
         
         !write(*,*) "Para P", level, lv_if(level, 1), lv_if(level,2)       

      END DO

      ! Alocando as variaveis
      ! lv_if( Levels, 2 ) => comprimento total do array com todas as malhas
      ALLOCATE( u( lv_if( Levels, 2 ), Tn ) , &
                f( lv_if( Levels, 2 ), Tn ) , &
                bp( lv_if( Levels, 2 ), Tn ) &
              )
      
     ! Para evitar lixo numerico na memoria         
     u(:, :) = 0.0d0
     f(:, :) = 0.0d0
     bp(:, :) = 0.0d0           
      
     !-------------------------- matriz A de P --------------------------------
     ALLOCATE( A(Levels) )
      
     ! Inicializando a matriz A de U
     DO level = 1, Levels
              
         hx  = Lx/DBLE( lv_nx(level) )
         hy  = Ly/DBLE( lv_ny(level) )
    
         CALL A(level) % init( lv_nx(level), lv_ny(level), &
                               hx, hy, ht, &
                               Re)
                              
         ! Atribuindo um stencil a cada volume de A                     
         CALL A(level) % atribui_stencils()
         
         ! Calculandos os coeficientes e posicoes dos stencils
         CALL A(level) % stencils_p()
         !-------------------------- matriz A de P --------------------------------     
     END DO
     
                      
   END SUBROUTINE Atualiza_Variaveis_MG_GS_P
   
   SUBROUTINE RESTRICAO_P (fina, grossa, Nxf, Nxg, Nyg)
      IMPLICIT NONE
        
      INTEGER, INTENT(IN)                           :: Nxf, Nxg, Nyg
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: fina
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: grossa
      
      ! Declaracao de variaveis interna a esta rotina
      INTEGER :: pg, pf
      
      DO pg = 1, Nxg * Nyg
         
         ! Deduzido de if = 2ig - 1, jg = 2jg - 1 e Nf = 2Nxg
         pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * Nxg - 1
         
         ! Fazendo a restricao dos nos internos de acordo com o stencil escolhido 
         grossa(pg) = &
        ( fina(pf) + fina(pf + 1 ) + fina(pf + Nxf) + fina(pf + Nxf + 1) )/4.0d0
      
      END DO
      
   END SUBROUTINE RESTRICAO_P
   
   SUBROUTINE Prolongacao_P (fina, grossa, Nxf, Nxg, Nyg )
      IMPLICIT NONE
        
      INTEGER, INTENT(IN)                            :: Nxf, Nxg, Nyg
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)     :: grossa
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)  :: fina

      ! Declaracao de variaveis interna a esta rotina
      INTEGER :: p, pf, pg
      
      ! Volumes no CENTRO recebem stencil 5, mas ja foram considerados
      DO p = 1, (Nxg-2)*(Nyg-2)
         
         pg = p + 2*FLOOR( DBLE(p-1)/DBLE(Nxg-2) ) + Nxg + 1
         pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * Nxg - 1    
         
         fina(pf)           = fina(pf)           + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg - 1 ) )/16.0d0
         fina(pf+1)         = fina(pf+1)         + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg + 1 ) )/16.0d0
         fina(pf + Nxf)     = fina(pf + Nxf)     + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg - 1 ) )/16.0d0
         fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg + 1 ) )/16.0d0
                  
      END DO
      
      ! Canto SOUTH WEST
      pg = 1
      pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * Nxg - 1
      ! Usando Neumann dp/dt = 0 entao p-1 = p e p-nx = p
      fina(pf)           = fina(pf)           + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg       ) + grossa(pg           ) )/16.0d0
      fina(pf+1)         = fina(pf+1)         + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg       ) + grossa(pg       + 1 ) )/16.0d0
      fina(pf + Nxf)     = fina(pf + Nxf)     + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg     ) )/16.0d0
      fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg + 1 ) )/16.0d0
      
      ! Contorno SOUTH recebe stencil 2
      DO pg = 2, Nxg - 1         
         pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * Nxg - 1         
         ! Usando Neumann dp/dt = 0 entao p - Nx = p
         fina(pf)           = fina(pf)           + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg       ) + grossa(pg       - 1 ) )/16.0d0
         fina(pf+1)         = fina(pf+1)         + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg       ) + grossa(pg       + 1 ) )/16.0d0
         fina(pf + Nxf)     = fina(pf + Nxf)     + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg - 1 ) )/16.0d0
         fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg + 1 ) )/16.0d0
         
      END DO
                           
      ! Canto SOUTH EAST  recebe stencil 3
      pg = Nxg
      pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * Nxg - 1
      ! Usando Neumann dp/dt = 0 entao p+1 = p e p-nx = p
      fina(pf)           = fina(pf)           + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg       ) + grossa(pg       - 1 ) )/16.0d0
      fina(pf+1)         = fina(pf+1)         + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg       ) + grossa(pg           ) )/16.0d0
      fina(pf + Nxf)     = fina(pf + Nxf)     + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg - 1 ) )/16.0d0
      fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg + 1 ) )/16.0d0
      
      ! Contorno WEST recebe stencil 4
      DO pg = 1 + Nxg, Nxg *( Nyg - 2 ) + 1, Nxg
          pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * Nxg - 1
          ! Usando Neumann dp/dt = 0 entao p-1 = p
          fina(pf)           = fina(pf)           + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg     ) )/16.0d0
          fina(pf+1)         = fina(pf+1)         + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg + 1 ) )/16.0d0
          fina(pf + Nxf)     = fina(pf + Nxf)     + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg     ) )/16.0d0
          fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg + 1 ) )/16.0d0
                            
      END DO
         
      ! Contorno EAST recebe stencil 6
      DO pg = 2 * Nxg, Nxg * (Nyg - 1), Nxg
            pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * Nxg - 1
            ! Usando Neumann dp/dt = 0 entao p+1 = p
            fina(pf)           = fina(pf)           + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg - 1 ) )/16.0d0
            fina(pf+1)         = fina(pf+1)         + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg     ) )/16.0d0
            fina(pf + Nxf)     = fina(pf + Nxf)     + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg - 1 ) )/16.0d0
            fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg + Nxg ) + grossa(pg + Nxg     ) )/16.0d0
            
      END DO
         
     ! Canto NORTH WEST recebe stencil 7
     pg = Nxg *( Nyg - 1 ) + 1
     pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * Nxg - 1
     ! Usando Neumann dp/dt = 0 entao p-1 = p, p+nx = p
     fina(pf)           = fina(pf)           + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg     ) )/16.0d0
     fina(pf+1)         = fina(pf+1)         + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg + 1 ) )/16.0d0
     fina(pf + Nxf)     = fina(pf + Nxf)     + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg       ) + grossa(pg           ) )/16.0d0
     fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg       ) + grossa(pg       + 1 ) )/16.0d0
     
     ! Contorno NORTH recebe stencil 8
     DO pg = Nxg *( Nyg - 1 ) + 2, Nxg * Nyg - 1
         pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * Nxg - 1
         ! Usando Neumann dp/dt = 0 entao p+nx = p
         fina(pf)           = fina(pf)           + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg - 1 ) )/16.0d0
         fina(pf+1)         = fina(pf+1)         + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg + 1 ) )/16.0d0
         fina(pf + Nxf)     = fina(pf + Nxf)     + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg       ) + grossa(pg       - 1 ) )/16.0d0
         fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + 3.0d0 * grossa(pg       ) + grossa(pg       + 1 ) )/16.0d0
         
     END DO
                  
     ! Canto NORTH EAST recebe stencil 9
     pg = Nxg * Nyg
     pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * Nxg - 1
     ! Usando Neumann dp/dt = 0 entao p+1 =p e p+nx = p
     fina(pf)           = fina(pf)           + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg - 1 ) )/16.0d0
     fina(pf+1)         = fina(pf+1)         + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg - Nxg ) + grossa(pg - Nxg     ) )/16.0d0
     fina(pf + Nxf)     = fina(pf + Nxf)     + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + 3.0d0 * grossa(pg       ) + grossa(pg       - 1 ) )/16.0d0
     fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 9.0d0 * grossa(pg) + 3.0d0 * grossa(pg    ) + 3.0d0 * grossa(pg       ) + grossa(pg           ) )/16.0d0
      
   END SUBROUTINE Prolongacao_P
   
   
      
   SUBROUTINE calcula_residuo_p(A, u, b, r)
      
        USE ClassMatrizA
   
       IMPLICIT NONE
       
       TYPE ( MatrizA ), INTENT(IN) :: A
       DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: r
       DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u, b
       
       INTEGER :: p, stencil, ite, se ! se = elemento do stencil
       DOUBLE PRECISION :: ax, ux, sm, correcao ! Soma dos componentes do stencil
       
       
       DO p = 1, A % Nx * A % Ny
               
           ! Obtendo o stencil do volume
           stencil = A % v(p)
               
           ! Para cada coeficiente do stencil da matriz
           ! Lembrando que o primeiro e sempre ap e Gauss-Seidel
           ! e definido como 
           ! u(p) = (aw * u(p-1) + ae * u(p+1) + as * u(p-Nx) ... + b(p)) / ap
           ! O codigo abaixo soma aw * u(p-1) + ae * u(p+1) + as * u(p-Nx) ...
           sm = 0.0d0
           DO se = 2, SIZE(A % s(stencil) % p) !(comeca em 2 porque nao considera o ap)
              
               ! ax =     aw,     ae,       as,      an,     aww,        ann, etc...
               ! ux = u(p-1), u(p+1),  u(p-Nx), u(p+Nx),  u(p-2),  u(p+2*Nx), etc...
               ! Logo
               ! ax * ux = aw * u(p-1), ae * u(p+1), as * u(p-Nx), etc...
                   
               ax = A % s(stencil) % c(se)
               ux = u(p + A % s(stencil) % p(se))
                   
               sm = sm +  ax * ux
                  
           END DO 
               
           ! Obtendo o ap
           ax = A % s(stencil) % c(1)
                                 
           ! Somando b(p) e subtraindo ap * up para compor o residuo
           r(p) = sm + b(p) - ax * u(p)
               
       END DO
       
   END SUBROUTINE calcula_residuo_p
   
  
   SUBROUTINE gs_p( A, u, b, iteracoes )
      
       USE ClassMatrizA
   
       IMPLICIT NONE
       
       TYPE ( MatrizA ), INTENT(IN) :: A
       DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u
       DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: b
       INTEGER, INTENT(IN) :: iteracoes
       
       INTEGER :: p, stencil, ite, se ! se = elemento do stencil
       DOUBLE PRECISION :: ax, ux, sm ! Soma dos componentes do stencil
       
       DO ite = 1, iteracoes
       
           DO p = 1, A % Nx * A % Ny
               
               ! Obtendo o stencil do volume
               stencil = A % v(p)
               
               ! Para cada coeficiente do stencil da matriz
               ! Lembrando que o primeiro e sempre ap e Gauss-Seidel
               ! e definido como 
               ! u(p) = (aw * u(p-1) + ae * u(p+1) + as * u(p-Nx) ... + b(p)) / ap
               ! O codigo abaixo soma aw * u(p-1) + ae * u(p+1) + as * u(p-Nx) ...
               sm = 0.0d0
               DO se = 2, SIZE(A % s(stencil) % p) !(comeca em 2 porque nao considera o ap)
               
                   ! ax =     aw,     ae,       as,      an,     aww,        ann, etc...
                   ! ux = u(p-1), u(p+1),  u(p-Nx), u(p+Nx),  u(p-2),  u(p+2*Nx), etc...
                   ! Logo
                   ! ax * ux = aw * u(p-1), ae * u(p+1), as * u(p-Nx), etc...
                   
                   ax = A % s(stencil) % c(se)
                   ux = u(p + A % s(stencil) % p(se))
                   
                   sm = sm +  ax * ux
                   
                   !if (A % s(stencil) % p(se) .eq. -1) then
                   !write(*,*) p, ax, A % s(stencil) % p(se)
                   !end if
               !    write(*,*) p, ax
                   
               END DO 
               
               ! Obtendo o ap
               ax = A % s(stencil) % c(1)
               
               ! Somando b(p) e divindo por ap para completar o Gauss-Seidel
               u(p) = (sm + b(p))/ax
               
              ! write(*,*) p, ax, b(p), u(p)
               
           END DO
           !write(*,*)
           !call exit()
          
       END DO
       
   END SUBROUTINE gs_p

   
   SUBROUTINE MG_GS_P( u_o, f_o, ciclos_V, t)
      
      USE Variaveis_Solvers_P
      
      INTEGER, INTENT(IN) :: ciclos_V, t
      
      ! u_o = estimativa inicial (pressao) 
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: u_o
      
      ! termo fonte ( div (ut, vt) )
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: f_o 
      
      INTEGER :: level, ciclo, iteracao
      
      11 FORMAT( 1X, "Multigrid P => Ciclo:", 1X, I3, 3X, "Parada:", 1X, 1PE16.10, 3X, "Norma L2:", 1X,  1PE16.10, 3X, "(Norma L2)**(1/ciclo)", 1X, 1PE16.10)
        
      ! Atualizando a primeira malha do vetor u, ou seja, 
      !  u (lv_if(1,1):lv_if(1,2)) = u_o(:) => estimativa inicial e a pressao 
      !  calculada em um passo de tempo anterior
      u(lv_if(1,1):lv_if(1,2), t) = u_o(:)
      
      DO ciclo = 1, 1 !max_ciclos_V
         
         bp(lv_if(1,1):lv_if(1,2), t) = f_o(:) 
         
         CALL calcula_residuo_p ( A(1), &
                                  u(lv_if(1,1):lv_if(1,2), t), &
                                  bp(lv_if(1,1):lv_if(1,2), t), &
                                  f(lv_if(1,1):lv_if(1,2), t) )
         norma_L2 = DSQRT(SUM(f(lv_if(1,1):lv_if(1,2), t) * f(lv_if(1,1):lv_if(1,2), t)))
         normas_L2(ciclo, t) = norma_L2
         !WRITE(*, 11) ciclo, normas_L2(ciclo)/normas_L2(1), normas_L2(ciclo), (normas_L2(ciclo)/normas_L2(1))**(1.0d0/ciclo)
        
         IF ( normas_L2(ciclo, t)/normas_L2(1, t) < tol ) THEN

            EXIT

         END IF

         DO level = 1, Levels - 1

            ! So atualiza bp = residuo, ou seja, apenas para level > 1 (o primeiro e bp = fonte)
            IF ( level > 1 ) THEN

               bp(lv_if(level,1):lv_if(level,2), t) = f(lv_if(level,1):lv_if(level,2), t)

            END IF
                        
            CALL gs_p ( A(level), u(lv_if(level,1):lv_if(level,2), t), bp(lv_if(level,1):lv_if(level,2), t), ite_down )
            
            CALL calcula_residuo_p ( A(level), &
                                     u(lv_if(level,1):lv_if(level,2), t), &
                                     bp(lv_if(level,1):lv_if(level,2), t), &
                                     f(lv_if(level,1):lv_if(level,2), t) )
                                  
            CALL Restricao_P ( f(lv_if(level,1):lv_if(level,2), t),               &
                               f(lv_if(level+1, 1):lv_if(level+1, 2), t),         &
                               lv_nx(level),                                   &
                               lv_nx(level+1), lv_ny(level+1) )
         
         END DO
      
         level = Levels
         bp(lv_if(level,1):lv_if(level,2), t) = f(lv_if(level,1):lv_if(level,2), t)
         CALL gs_p ( A(level), u(lv_if(level,1):lv_if(level,2), t), bp(lv_if(level,1):lv_if(level,2), t), ite_down )                  
         
         DO level = Levels, 2, -1

             CALL Prolongacao_P ( u(lv_if(level-1,1):lv_if(level-1,2), t), &
                                  u(lv_if(level, 1):lv_if(level, 2), t),   &
                                  lv_nx(level-1),                       &
                                  lv_nx(level), lv_ny(level) )
            
             ! Limpando o vetor para evitar erros no proximo ciclo
             u(lv_if(level, 1):lv_if(level, 2), t) = 0.0d0                                       
             CALL gs_p ( A(level-1), u(lv_if(level-1,1):lv_if(level-1,2), t), &
                                    bp(lv_if(level-1,1):lv_if(level-1,2), t), ite_up ) 
                                    
         END DO             
      END DO
      
      u_o(:) = u(lv_if(1,1):lv_if(1,2), t)
      
   END SUBROUTINE MG_GS_P 
  
   
   SUBROUTINE Atualiza_Variaveis_MG_GS_U( Lx_o, Ly_o, &
                                          Nx_o, Ny_o, &
                                          tol_o, &
                                          ite_down_o, ite_up_o, &
                                          max_ciclos_V_o, ht, Re, Tn_o )

      USE Variaveis_Solvers_U
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: Nx_o, Ny_o, ite_down_o, ite_up_o, max_ciclos_V_o, Tn_o
      DOUBLE PRECISION, INTENT(IN) :: Lx_o, Ly_o, tol_o, ht, Re
      
      INTEGER :: level
      DOUBLE PRECISION :: hx, hy
      
      Lx = Lx_o
      Ly = Ly_o
      Nx = Nx_o
      Ny = Ny_o
      tol = tol_o
      ite_down = ite_down_o
      ite_up   = ite_up_o
      max_ciclos_V = max_ciclos_V_o
      Tn = Tn_o
      
      ! Calculando o numero de niveis
      Levels = LOG( DBLE(Nx+1))/LOG( DBLE(2) ) - 1
      
      ! lv_if => duas dimensoes 1: inicio da malha, 2: fim da malha
      ALLOCATE( lv_if(Levels, 2), lv_nx(Levels), lv_ny(Levels), normas_L2(max_ciclos_V, Tn) )

      ! Inicio da primeira malha comeca em 1
      lv_if(1,1) = 1
      
      ! Fim da primeira malha termina em Nx*Ny
      lv_if(1,2) = Nx*Ny

      ! Para cada nivel de malha
      DO level = 1, Levels
         
         ! Calcula o Nx e Ny de cada nivel de malha
         lv_nx(level) = (Nx+1)/(2**(level-1)) - 1
         lv_ny(level) = Ny/(2**(level-1))

         ! O inicio e fim da primeira malha ja foram calculados
         IF ( level > 1 ) THEN

            ! Calculando o indice que define o fim do vetor
            lv_if(level,2) = lv_if(level-1,2) + lv_nx(level)*lv_ny(level)

            ! Calculando o indice que define o inicio do vetor
            lv_if(level,1) = lv_if(level-1,2) + 1

         END IF
         
         !write(*,*) "Para P", level, lv_if(level, 1), lv_if(level,2)       

      END DO
      
      ! Alocando as variaveis
      ! lv_if( Levels, 2 ) => comprimento total do array com todas as malhas
      ALLOCATE( u( lv_if( Levels, 2 ), Tn ) , &
                f( lv_if( Levels, 2 ), Tn ) , &
                bp( lv_if( Levels, 2 ), Tn ) &
              )
      
     ! Para evitar lixo numerico na memoria         
     u(:, :) = 0.0d0
     f(:, :) = 0.0d0
     bp(:, :) = 0.0d0           
      
     !-------------------------- matriz A de P --------------------------------
     ALLOCATE( A(Levels) )
      
     ! Inicializando a matriz A de U
     DO level = 1, Levels
              
         hx  = Lx/DBLE( lv_nx(level) + 1 )
         hy  = Ly/DBLE( lv_ny(level) )
         
         CALL A(level) % init( lv_nx(level), lv_ny(level), &
                               hx, hy, ht, &
                               Re)
                              
         ! Atribuindo um stencil a cada volume de A                     
         CALL A(level) % atribui_stencils()
         
         ! Calculandos os coeficientes e posicoes dos stencils
         CALL A(level) % stencils_u(1.0d0)
         !-------------------------- matriz A de P --------------------------------     
     END DO
                      
   END SUBROUTINE Atualiza_Variaveis_MG_GS_U
   
   
   SUBROUTINE Atualiza_Variaveis_MG_GS_V( Lx_o, Ly_o, &
                                          Nx_o, Ny_o, &
                                          tol_o, &
                                          ite_down_o, ite_up_o, &
                                          max_ciclos_V_o, ht, Re, Tn_o )

      USE Variaveis_Solvers_V
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: Nx_o, Ny_o, ite_down_o, ite_up_o, max_ciclos_V_o, Tn_o
      DOUBLE PRECISION, INTENT(IN) :: Lx_o, Ly_o, tol_o, ht, Re
      
      INTEGER :: level
      DOUBLE PRECISION :: hx, hy
      
      Lx = Lx_o
      Ly = Ly_o
      Nx = Nx_o
      Ny = Ny_o
      tol = tol_o
      ite_down = ite_down_o
      ite_up   = ite_up_o
      max_ciclos_V = max_ciclos_V_o
      Tn = Tn_o
      
      ! Calculando o numero de niveis
      Levels = LOG( DBLE(Nx))/LOG( DBLE(2) ) - 1
      
      ! lv_if => duas dimensoes 1: inicio da malha, 2: fim da malha
      ALLOCATE( lv_if(Levels, 2), lv_nx(Levels), lv_ny(Levels), normas_L2(max_ciclos_V, Tn) )

      ! Inicio da primeira malha comeca em 1
      lv_if(1,1) = 1
      
      ! Fim da primeira malha termina em Nx*Ny
      lv_if(1,2) = Nx*Ny

      ! Para cada nivel de malha
      DO level = 1, Levels
         
         ! Calcula o Nx e Ny de cada nivel de malha
         lv_nx(level) = Nx/(2**(level-1))
         lv_ny(level) = (Ny+1)/(2**(level-1)) - 1

         ! O inicio e fim da primeira malha ja foram calculados
         IF ( level > 1 ) THEN

            ! Calculando o indice que define o fim do vetor
            lv_if(level,2) = lv_if(level-1,2) + lv_nx(level)*lv_ny(level)

            ! Calculando o indice que define o inicio do vetor
            lv_if(level,1) = lv_if(level-1,2) + 1

         END IF
         
         !write(*,*) "Para P", level, lv_if(level, 1), lv_if(level,2)       

      END DO
      
      ! Alocando as variaveis
      ! lv_if( Levels, 2 ) => comprimento total do array com todas as malhas
      ALLOCATE( u( lv_if( Levels, 2 ), Tn ) , &
                f( lv_if( Levels, 2 ), Tn ) , &
                bp( lv_if( Levels, 2 ), Tn ) &
              )
      
     ! Para evitar lixo numerico na memoria         
     u(:, :) = 0.0d0
     f(:, :) = 0.0d0
     bp(:, :) = 0.0d0           
      
     !-------------------------- matriz A de P --------------------------------
     ALLOCATE( A(Levels) )
      
     ! Inicializando a matriz A de U
     DO level = 1, Levels
              
         hx  = Lx/DBLE( lv_nx(level) )
         hy  = Ly/DBLE( lv_ny(level) + 1 )
         
         CALL A(level) % init( lv_nx(level), lv_ny(level), &
                               hx, hy, ht, &
                               Re)
                              
         ! Atribuindo um stencil a cada volume de A                     
         CALL A(level) % atribui_stencils()
         
         ! Calculandos os coeficientes e posicoes dos stencils
         CALL A(level) % stencils_v(1.0d0)
         !-------------------------- matriz A de P --------------------------------     
     END DO
                      
   END SUBROUTINE Atualiza_Variaveis_MG_GS_V
   
   
END MODULE Subrotinas_MG_GS
