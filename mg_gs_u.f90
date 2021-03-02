MODULE MG_GS_U


    CONTAINS
    
    SUBROUTINE atualiza_u2(t, uo)
    
        USE Variaveis_Solvers_U
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t
        DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: uo
    
         uo(:) = u(lv_if(1,1):lv_if(1,2), t)
                
    END SUBROUTINE atualiza_u2
    
    
    SUBROUTINE atualiza_uo(t, uo)
    
        USE Variaveis_Solvers_U
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: uo
    
        u(lv_if(1,1):lv_if(1,2), t) = uo(:)
        
    END SUBROUTINE atualiza_uo
        
    
    SUBROUTINE atualiza_bp_u(t, f_o)
    
        USE Variaveis_Solvers_U
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: f_o
    
        bp(lv_if(1,1):lv_if(1,2), t) = f_o(:) 
        
    END SUBROUTINE atualiza_bp_u
    
    SUBROUTINE limpa_u(t, level)
    
        USE Variaveis_Solvers_U
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t, level
        
        u(lv_if(level, 1):lv_if(level, 2), t) = 0.0d0       
        
     END SUBROUTINE limpa_u
        
        
    SUBROUTINE atualiza_bp_u2(t, level)
    
        USE Variaveis_Solvers_U
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t, level
        
        bp(lv_if(level,1):lv_if(level,2), t) = f(lv_if(level,1):lv_if(level,2), t) 
        
    END SUBROUTINE atualiza_bp_u2
    
    SUBROUTINE calcula_residuo_u(t, level)
        
        USE Variaveis_Solvers_U
        USE Residuo
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: t, level
        
        CALL calcula_residuo_UV( A(level), &
                                 u(lv_if(level,1):lv_if(level,2), :), &
                                bp(lv_if(level,1):lv_if(level,2), :), &
                                 f(lv_if(level,1):lv_if(level,2), :), &
                                 t )
    
    END SUBROUTINE calcula_residuo_u
    
    SUBROUTINE calcula_norma_L2_u(ciclo, t)
    
        USE Variaveis_Solvers_U
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t, ciclo
    
        normas_L2(ciclo, t) = SUM(f(lv_if(1,1):lv_if(1,2), t) * f(lv_if(1,1):lv_if(1,2), t))
        
    END SUBROUTINE calcula_norma_L2_u


    SUBROUTINE soma_norma_L2_u(ciclo, normas_L2_u)
    
        USE Variaveis_Solvers_U
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ciclo
        DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: normas_L2_u
    
        normas_L2_u(ciclo) = DSQRT( SUM(normas_L2(ciclo, : )) )
        
    END SUBROUTINE soma_norma_L2_u


    SUBROUTINE GS_RB_U(t, level, rb)
    
        USE Variaveis_Solvers_U
        USE GS
        IMPLICIT NONE
            
        INTEGER, INTENT(IN) :: t, level, rb
        
        DOUBLE PRECISIOn :: fator
        
        fator = 1.0d0
        if (level > 1) then
            fator = 0.0d0
        end if
        fator = 1.0d0
             
        CALL GaussSeidel_UV2( A(level), &
                              u(lv_if(level,1):lv_if(level,2), :), &
                              bp(lv_if(level,1):lv_if(level,2), :), &
                              t, &
                              rb, fator )
    
    END SUBROUTINE GS_RB_U
    
    SUBROUTINE restringe_u (fina, grossa, Nxf, Nxg, Nyg )
      IMPLICIT NONE
        
      INTEGER, INTENT(IN)                           :: Nxf, Nxg, Nyg
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: fina
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: grossa
      
      INTEGER :: pg, pf
      
      DO pg = 1, Nxg * Nyg
      
          ! Deduzido de if = 2ig, jg = 2jg - 1 e Nf = 2Nxg + 1
          pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * (Nxg + 1)
         
         ! Fazendo a restricao dos nos internos de acordo com o stencil escolhido 
         !grossa(pg) = &
         !( 2.0d0*fina(pf)               + 2.0d0*fina(pf + Nxf) + &
         !       fina(pf - 1)           +       fina(pf + 1) + &
         !       fina(pf + Nxf - 1) +       fina(pf + Nxf + 1) )/8.0d0
      
         grossa(pg) = ( fina(pf) + fina(pf + Nxf) )/2.0d0 
                          
      END DO
      
    END SUBROUTINE restringe_u
    
    SUBROUTINE restricao_u(t, level)
    
        USE Variaveis_Solvers_U
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: t, level
        
        CALL restringe_u( f(lv_if(level,1):lv_if(level,2), t),               &
                          f(lv_if(level+1, 1):lv_if(level+1, 2), t),         &
                          lv_nx(level),                                      &
                          lv_nx(level+1), lv_ny(level+1) )
        
    
    END SUBROUTINE restricao_u
    
    
    SUBROUTINE prolonga_u (fina, grossa, Nxf, Nxg, Nyg )
      IMPLICIT NONE
        
      INTEGER, INTENT(IN)                            :: Nxf, Nxg, Nyg
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)     :: grossa
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)  :: fina

      ! Declaracao de variaveis interna a esta rotina
      INTEGER :: p, pf, pg
      
      ! Volumes no CENTRO recebem stencil 5, mas ja foram considerados
      DO p = 1, (Nxg-2)*(Nyg-2)
         
         pg = p + 2*FLOOR( DBLE(p-1)/DBLE(Nxg-2) ) + Nxg + 1
         ! Deduzido de if = 2ig, jg = 2jg - 1 e Nf = 2Nxg + 1
         pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * (Nxg + 1)
         
         ! 6 velocidades u na malha fina para serem atualizados:
         ! pf, pf - 1, pf + 1, pf + Nxf, pf + Nxf - 1, pf + Nxf + 1
         fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg) +         grossa(pg - Nxg) )/4.0d0
         fina(pf + Nxf)     = fina(pf + Nxf )    + ( 3.0d0 * grossa(pg) +         grossa(pg + Nxg) )/4.0d0
         fina(pf + 1)       = fina(pf + 1)       + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + grossa(pg - Nxg + 1) + grossa(pg - Nxg) )/8.0d0
         fina(pf - 1)       = fina(pf - 1)       + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + grossa(pg - Nxg - 1) + grossa(pg - Nxg) )/8.0d0
         fina(pf + Nxf - 1) = fina(pf + Nxf - 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + grossa(pg + Nxg - 1) + grossa(pg + Nxg) )/8.0d0
         fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + grossa(pg + Nxg + 1) + grossa(pg + Nxg) )/8.0d0
         
      END DO
      
      ! Canto SOUTH WEST
      pg = 1
      pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * (Nxg + 1)
      ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
      ! grossa(pg - 1)   = 0 (esta no contorno)
      ! grossa(pg + Nx - 1) = 0 (esta no contorno)
      ! grossa(pg - Nxg) = - 3.0d0 * grossa(pg) + grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg) 
      ! grossa(pg - Nxg - 1) = 0 (extrapolacao de valor no contorno que e 0)
      ! grossa(pg - Nxg + 1) = - 3.0d0 * grossa(pg + 1) + grossa(pg + 1 + Nxg) - (1.0d0/5.0d0) * grossa(pg + 1 + 2*Nxg) 
      fina(pf)           = fina(pf)           + (                              grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2 * Nxg) )/4.0d0
      fina(pf + Nxf)     = fina(pf + Nxf )    + ( 3.0d0 * grossa(pg) +         grossa(pg + Nxg) )/4.0d0
      fina(pf + 1)       = fina(pf + 1)       + ( grossa(pg + 1 + Nxg) - (1.0d0/5.0d0) * grossa(pg + 1 + 2*Nxg) + grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg)  )/8.0d0
      fina(pf - 1)       = fina(pf - 1)       + ( grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg)  )/8.0d0
      fina(pf + Nxf - 1) = fina(pf + Nxf - 1) + ( 3.0d0 * grossa(pg) + grossa(pg + Nxg) )/8.0d0
      fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + grossa(pg + Nxg + 1) + grossa(pg + Nxg) )/8.0d0
      
      ! Contorno SOUTH recebe stencil 2
      DO pg = 2, Nxg - 1         
          pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * (Nxg + 1)
          ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
          ! grossa(pg - Nxg) = - 3.0d0 * grossa(pg) + grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg) 
          ! grossa(pg - Nxg + 1) = - 3.0d0 * grossa(pg + 1) + grossa(pg + 1 + Nxg) - (1.0d0/5.0d0) * grossa(pg + 1 + 2*Nxg) 
          ! grossa(pg - Nxg - 1) = - 3.0d0 * grossa(pg - 1) + grossa(pg - 1 + Nxg) - (1.0d0/5.0d0) * grossa(pg - 1 + 2*Nxg) 
          fina(pf)           = fina(pf)           + ( grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg)  )/4.0d0
          fina(pf + Nxf)     = fina(pf + Nxf )    + ( 3.0d0 * grossa(pg) +         grossa(pg + Nxg) )/4.0d0
          fina(pf + 1)       = fina(pf + 1)       + ( grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg) + grossa(pg + 1 + Nxg) - (1.0d0/5.0d0) * grossa(pg + 1 + 2*Nxg)   )/8.0d0
          fina(pf - 1)       = fina(pf - 1)       + ( grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg) + grossa(pg - 1 + Nxg) - (1.0d0/5.0d0) * grossa(pg - 1 + 2*Nxg)  )/8.0d0
          fina(pf + Nxf - 1) = fina(pf + Nxf - 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + grossa(pg + Nxg - 1) + grossa(pg + Nxg) )/8.0d0
          fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + grossa(pg + Nxg + 1) + grossa(pg + Nxg) )/8.0d0          
      END DO  
      
      ! Canto SOUTH EAST  recebe stencil 3
      pg = Nxg
      pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * (Nxg + 1)
      ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
      ! grossa(pg + 1)   = 0 (esta no contorno)
      ! grossa(pg + Nx + 1) = 0 (esta no contorno)
      ! grossa(pg - Nxg) = - 3.0d0 * grossa(pg) + grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg) 
      ! grossa(pg - Nxg - 1) = - 3.0d0 * grossa(pg - 1) + grossa(pg - 1 + Nxg) - (1.0d0/5.0d0) * grossa(pg - 1 + 2*Nxg) 
      ! grossa(pg - Nxg + 1) = 0 (extrapolacao de valor no contorno que e 0)
      fina(pf)           = fina(pf)           + ( grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg) )/4.0d0
      fina(pf + Nxf)     = fina(pf + Nxf )    + ( 3.0d0 * grossa(pg) + grossa(pg + Nxg) )/4.0d0
      fina(pf + 1)       = fina(pf + 1)       + ( grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg) )/8.0d0
      fina(pf - 1)       = fina(pf - 1)       + ( grossa(pg - 1 + Nxg) - (1.0d0/5.0d0) * grossa(pg - 1 + 2*Nxg) + grossa(pg + Nxg) - (1.0d0/5.0d0) * grossa(pg + 2*Nxg) )/8.0d0
      fina(pf + Nxf - 1) = fina(pf + Nxf - 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + grossa(pg + Nxg - 1) + grossa(pg + Nxg) )/8.0d0
      fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + grossa(pg + Nxg) )/8.0d0
      
      ! Contorno WEST recebe stencil 4
      DO pg = 1 + Nxg, Nxg *( Nyg - 2 ) + 1, Nxg
          pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * (Nxg + 1)
          
          ! grossa(pg - 1)   = 0 (esta no contorno)
          ! grossa(pg + Nx - 1) = 0 (esta no contorno)
          ! grossa(pg - Nxg - 1) = 0 (esta no contorno)
          fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg) +         grossa(pg - Nxg) )/4.0d0
          fina(pf + Nxf)     = fina(pf + Nxf )    + ( 3.0d0 * grossa(pg) +         grossa(pg + Nxg) )/4.0d0
          fina(pf + 1)       = fina(pf + 1)       + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + grossa(pg - Nxg + 1) + grossa(pg - Nxg) )/8.0d0
          fina(pf - 1)       = fina(pf - 1)       + ( 3.0d0 * grossa(pg)                                                 + grossa(pg - Nxg) )/8.0d0
          fina(pf + Nxf - 1) = fina(pf + Nxf - 1) + ( 3.0d0 * grossa(pg)                                                 + grossa(pg + Nxg) )/8.0d0
          fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + grossa(pg + Nxg + 1) + grossa(pg + Nxg) )/8.0d0
          
      END DO
      
      ! Contorno EAST recebe stencil 6
      DO pg = 2 * Nxg, Nxg * (Nyg - 1), Nxg
          pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * (Nxg + 1)
          
          ! grossa(pg + 1)   = 0 (esta no contorno)
          ! grossa(pg + Nx + 1) = 0 (esta no contorno)
          ! grossa(pg - Nxg + 1) = 0 (esta no contorno)
      
          fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg) +         grossa(pg - Nxg) )/4.0d0
          fina(pf + Nxf)     = fina(pf + Nxf )    + ( 3.0d0 * grossa(pg) +         grossa(pg + Nxg) )/4.0d0
          fina(pf + 1)       = fina(pf + 1)       + ( 3.0d0 * grossa(pg)                                                 + grossa(pg - Nxg) )/8.0d0
          fina(pf - 1)       = fina(pf - 1)       + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + grossa(pg - Nxg - 1) + grossa(pg - Nxg) )/8.0d0
          fina(pf + Nxf - 1) = fina(pf + Nxf - 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + grossa(pg + Nxg - 1) + grossa(pg + Nxg) )/8.0d0
          fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg)                                                 + grossa(pg + Nxg) )/8.0d0
      
      END DO
      
      ! Canto NORTH WEST recebe stencil 7
      pg = Nxg *( Nyg - 1 ) + 1
      pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * (Nxg + 1)
      ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
      ! grossa(pg - 1)   = 0 (esta no contorno)
      ! grossa(pg - Nx - 1) = 0 (esta no contorno)
      ! grossa(pg + Nxg) = - 3.0d0 * grossa(pg) + grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg) 
      ! grossa(pg + Nxg - 1) = 0 (extrapolacao de valor no contorno que e 0)
      ! grossa(pg + Nxg + 1) = - 3.0d0 * grossa(pg + 1) + grossa(pg + 1 - Nxg) - (1.0d0/5.0d0) * grossa(pg + 1 - 2*Nxg)
      
      fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg) +         grossa(pg - Nxg) )/4.0d0
      fina(pf + Nxf)     = fina(pf + Nxf )    + ( grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg)  )/4.0d0
      fina(pf + 1)       = fina(pf + 1)       + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + grossa(pg - Nxg + 1) + grossa(pg - Nxg) )/8.0d0
      fina(pf - 1)       = fina(pf - 1)       + ( 3.0d0 * grossa(pg)                                                 + grossa(pg - Nxg) )/8.0d0
      fina(pf + Nxf - 1) = fina(pf + Nxf - 1) + ( grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg)  )/8.0d0
      fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( grossa(pg + 1 - Nxg) - (1.0d0/5.0d0) * grossa(pg + 1 - 2*Nxg) + grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg)  )/8.0d0 
      
      ! Contorno NORTH recebe stencil 8
      DO pg = Nxg *( Nyg - 1 ) + 2, Nxg * Nyg - 1
          pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * (Nxg + 1)
          
          ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
          ! grossa(pg + Nxg) = - 3.0d0 * grossa(pg) + grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg) 
          ! grossa(pg + Nxg + 1) = - 3.0d0 * grossa(pg + 1) + grossa(pg + 1 - Nxg) - (1.0d0/5.0d0) * grossa(pg + 1 - 2*Nxg) 
          ! grossa(pg + Nxg - 1) = - 3.0d0 * grossa(pg - 1) + grossa(pg - 1 - Nxg) - (1.0d0/5.0d0) * grossa(pg - 1 - 2*Nxg) 
          
          fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg) +         grossa(pg - Nxg) )/4.0d0
          fina(pf + Nxf)     = fina(pf + Nxf )    + ( grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg)  )/4.0d0
          fina(pf + 1)       = fina(pf + 1)       + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + 1) + grossa(pg - Nxg + 1) + grossa(pg - Nxg) )/8.0d0
          fina(pf - 1)       = fina(pf - 1)       + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + grossa(pg - Nxg - 1) + grossa(pg - Nxg) )/8.0d0
          fina(pf + Nxf - 1) = fina(pf + Nxf - 1) + ( grossa(pg - 1 - Nxg) - (1.0d0/5.0d0) * grossa(pg - 1 - 2*Nxg) + grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg) )/8.0d0
          fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( grossa(pg + 1 - Nxg) - (1.0d0/5.0d0) * grossa(pg + 1 - 2*Nxg) + grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg) )/8.0d0
      
      END DO
      
      ! Canto NORTH EAST recebe stencil 9
      pg = Nxg * Nyg
      pf = 2 * pg + 2 * FLOOR( DBLE(pg-1)/Nxg ) * (Nxg + 1)
      
      ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
      ! grossa(pg + 1)   = 0 (esta no contorno)
      ! grossa(pg - Nx + 1) = 0 (esta no contorno)
      ! grossa(pg + Nxg) = - 3.0d0 * grossa(pg) + grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg) 
      ! grossa(pg + Nxg + 1) = 0 (extrapolacao de valor no contorno que e 0)
      ! grossa(pg + Nxg - 1) = - 3.0d0 * grossa(pg - 1) + grossa(pg - 1 - Nxg) - (1.0d0/5.0d0) * grossa(pg - 1 - 2*Nxg)
      fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg) +         grossa(pg - Nxg) )/4.0d0
      fina(pf + Nxf)     = fina(pf + Nxf )    + ( grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg)  )/4.0d0
      fina(pf + 1)       = fina(pf + 1)       + ( 3.0d0 * grossa(pg) + grossa(pg - Nxg) )/8.0d0
      fina(pf - 1)       = fina(pf - 1)       + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - 1) + grossa(pg - Nxg - 1) + grossa(pg - Nxg) )/8.0d0
      fina(pf + Nxf - 1) = fina(pf + Nxf - 1) + ( grossa(pg - 1 - Nxg) - (1.0d0/5.0d0) * grossa(pg - 1 - 2*Nxg) + grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg)  )/8.0d0
      fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( grossa(pg - Nxg) - (1.0d0/5.0d0) * grossa(pg - 2*Nxg) )/8.0d0
      
   END SUBROUTINE prolonga_U
    
    
    
    
     SUBROUTINE prolongacao_u(t, level)
    
        USE Variaveis_Solvers_U
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: t, level
        
        CALL prolonga_u( u(lv_if(level-1,1):lv_if(level-1,2), t), &
                         u(lv_if(level, 1):lv_if(level, 2), t),   &
                         lv_nx(level-1),                       &
                         lv_nx(level), lv_ny(level) )
        
    
    END SUBROUTINE prolongacao_u
    
END MODULE MG_GS_U