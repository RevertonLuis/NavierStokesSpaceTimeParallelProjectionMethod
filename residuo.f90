MODULE Residuo

   CONTAINS
   
   SUBROUTINE calcula_residuo_UV( A, u, b, r, t)
      
       USE ClassMatrizA
   
       IMPLICIT NONE
       
       TYPE ( MatrizA ), INTENT(IN) :: A
       DOUBLE PRECISION, DIMENSION(:, :), INTENT(INOUT) :: r
       DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: u, b
       INTEGER, INTENT(IN) :: t
       
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
               ux = u(p + A % s(stencil) % p(se), t+1)
                   
               sm = sm +  ax * ux
                  
           END DO 
               
           ! Obtendo o ap
           ax = A % s(stencil) % c(1)
           sm  = sm + u(p, t)
                                 
           ! Somando b(p) e subtraindo ap * up para compor o residuo
           r(p, t+1) = sm + b(p, t) - ax * u(p, t+1)
               
       END DO
       
       
   END SUBROUTINE calcula_residuo_UV

   
   
   SUBROUTINE calcula_residuo_P(A, u, b, r, t)
      
        USE ClassMatrizA
   
       IMPLICIT NONE
       
       TYPE ( MatrizA ), INTENT(IN) :: A
       DOUBLE PRECISION, DIMENSION(:, :), INTENT(INOUT) :: r
       DOUBLE PRECISION, DIMENSION(:, :), INTENT(IN) :: u, b
       INTEGER, INTENT(IN) :: t
       
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
               ux = u(p + A % s(stencil) % p(se), t+1)
                   
               sm = sm +  ax * ux
                  
           END DO 
               
           ! Obtendo o ap
           ax = A % s(stencil) % c(1)
                                 
           ! Somando b(p) e subtraindo ap * up para compor o residuo
           r(p, t+1) = sm + b(p, t) - ax * u(p, t+1)
               
       END DO
       
   END SUBROUTINE calcula_residuo_P

END MODULE Residuo
