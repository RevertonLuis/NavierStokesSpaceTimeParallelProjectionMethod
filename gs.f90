MODULE GS

   CONTAINS
   
   SUBROUTINE GaussSeidel_UV2( A, u, b, t, rb, fator)
      
       USE ClassMatrizA
   
       IMPLICIT NONE
       
       TYPE ( MatrizA ), INTENT(IN) :: A
       DOUBLE PRECISION, DIMENSION(:, :), INTENT(INOUT) :: u
       DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: b
       INTEGER, INTENT(IN) :: t, rb
       DOUBLE PRECISION, INTENT(IN) :: fator
       
       INTEGER :: p, pp, stencil, ite, se ! se = elemento do stencil
       DOUBLE PRECISION :: ax, ux, sm, correcao ! Soma dos componentes do stencil
       INTEGER :: linha, inicio, new_rb
       
       new_rb = (1 + rb)/2
       inicio = ABS(MOD(t,2) - new_rb )
       
              DO pp = 1 + inicio, A % Nx * A % Ny, 2
               
                   ! Calculando a linha de p
                   linha = (pp-1)/(A % Nx) + 1
                   
                   ! Fazendo a correcao do caso Nx par (black points sao pares ao inves de impares)
                   p = pp + MOD(linha + 1, 2)*MOD(A % Nx + 1, 2)*(1 - 2*inicio)
   
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
                       ux = u(p + A % s(stencil) % p(se), t+1)
                   
                       sm = sm +  ax * ux
                   
                   END DO 
               
                  ! Obtendo o ap
                  ax = A % s(stencil) % c(1)
                  sm  = sm + u(p, t) * fator
                                 
                  ! Somando b(p) e divindo por ap para completar o Gauss-Seidel
                  u(p, t+1) = (sm + b(p, t))/ax
               
              END DO
       
       
   END SUBROUTINE GaussSeidel_UV2

   
   
   SUBROUTINE GaussSeidel_UV( A, u, b, t, rb)
      
       USE ClassMatrizA
   
       IMPLICIT NONE
       
       TYPE ( MatrizA ), INTENT(IN) :: A
       DOUBLE PRECISION, DIMENSION(:, :), INTENT(INOUT) :: u
       DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: b
       INTEGER, INTENT(IN) :: t, rb
       
       INTEGER :: p, pp, stencil, ite, se ! se = elemento do stencil
       DOUBLE PRECISION :: ax, ux, sm, correcao ! Soma dos componentes do stencil
       INTEGER :: linha, inicio, new_rb
       
       new_rb = (1 + rb)/2
       inicio = ABS(MOD(t,2) - new_rb )
       
              DO pp = 1 + inicio, A % Nx * A % Ny, 2
               
                   ! Calculando a linha de p
                   linha = (pp-1)/(A % Nx) + 1
                   
                   ! Fazendo a correcao do caso Nx par (black points sao pares ao inves de impares)
                   p = pp + MOD(linha + 1, 2)*MOD(A % Nx + 1, 2)*(1 - 2*inicio)
   
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
                       ux = u(p + A % s(stencil) % p(se), t+1)
                   
                       sm = sm +  ax * ux
                   
                   END DO 
               
                  ! Obtendo o ap
                  ax = A % s(stencil) % c(1)
                  sm  = sm + u(p, t)
                                 
                  ! Somando b(p) e divindo por ap para completar o Gauss-Seidel
                  u(p, t+1) = (sm + b(p, t))/ax
               
              END DO
       
       
   END SUBROUTINE GaussSeidel_UV

   
   
   SUBROUTINE GaussSeidel( A, u, b, iteracoes )
      
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
                   
               END DO 
               
               ! Obtendo o ap
               ax = A % s(stencil) % c(1)
               
               ! Somando b(p) e divindo por ap para completar o Gauss-Seidel
               u(p) = (sm + b(p))/ax
               
           END DO
           !call exit()
          
       END DO
       
   END SUBROUTINE GaussSeidel

END MODULE GS
