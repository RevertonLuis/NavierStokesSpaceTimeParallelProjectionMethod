MODULE MG_GS_V


    CONTAINS
    
    SUBROUTINE atualiza_v2(t, uo)
    
        USE Variaveis_Solvers_V
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t
        DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: uo
    
         uo(:) = u(lv_if(1,1):lv_if(1,2), t)
                
    END SUBROUTINE atualiza_v2
    
    
    SUBROUTINE atualiza_vo(t, uo)
    
        USE Variaveis_Solvers_V
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: uo
    
        u(lv_if(1,1):lv_if(1,2), t) = uo(:)
        
    END SUBROUTINE atualiza_vo
        
    
    SUBROUTINE atualiza_bp_v(t, f_o)
    
        USE Variaveis_Solvers_V
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: f_o
    
        bp(lv_if(1,1):lv_if(1,2), t) = f_o(:) 
        
    END SUBROUTINE atualiza_bp_v
    
    SUBROUTINE limpa_v(t, level)
    
        USE Variaveis_Solvers_V
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t, level
        
        u(lv_if(level, 1):lv_if(level, 2), t) = 0.0d0       
        
     END SUBROUTINE limpa_v
        
        
    SUBROUTINE atualiza_bp_v2(t, level)
    
        USE Variaveis_Solvers_V
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t, level
        
        bp(lv_if(level,1):lv_if(level,2), t) = f(lv_if(level,1):lv_if(level,2), t) 
        
    END SUBROUTINE atualiza_bp_v2
    
    SUBROUTINE calcula_residuo_v(t, level)
        
        USE Variaveis_Solvers_V
        USE Residuo
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: t, level
        
        CALL calcula_residuo_UV( A(level), &
                                 u(lv_if(level,1):lv_if(level,2), :), &
                                bp(lv_if(level,1):lv_if(level,2), :), &
                                 f(lv_if(level,1):lv_if(level,2), :), &
                                 t )
    
    END SUBROUTINE calcula_residuo_v
    
    SUBROUTINE calcula_norma_L2_v(ciclo, t)
    
        USE Variaveis_Solvers_V
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: t, ciclo
    
        normas_L2(ciclo, t) = SUM(f(lv_if(1,1):lv_if(1,2), t) * f(lv_if(1,1):lv_if(1,2), t))
        
    END SUBROUTINE calcula_norma_L2_v


    SUBROUTINE soma_norma_L2_v(ciclo, normas_L2_v)
    
        USE Variaveis_Solvers_V
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ciclo
        DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: normas_L2_v
    
        normas_L2_v(ciclo) = DSQRT( SUM(normas_L2(ciclo, : )) )
        
    END SUBROUTINE soma_norma_L2_v


    SUBROUTINE GS_RB_V(t, level, rb)
    
        USE Variaveis_Solvers_V
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
    
    END SUBROUTINE GS_RB_V
    
    SUBROUTINE restringe_v (fina, grossa, Nxf, Nxg, Nyg )
      IMPLICIT NONE
        
      INTEGER, INTENT(IN)                           :: Nxf, Nxg, Nyg
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: fina
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: grossa
      
      INTEGER :: pg, pf
      
      DO pg = 1, Nxg * Nyg
      
          ! Deduzido de if = 2ig-1, jg = 2jg e Nf = 2Nxg
          pf = 2 * pg + 2 * (FLOOR( DBLE(pg-1)/Nxg ) + 1) * Nxg - 1
         
         ! Fazendo a restricao dos nos internos de acordo com o stencil escolhido 
         !grossa(pg) = &
         !( 2.0d0*fina(pf)               + 2.0d0*fina(pf + Nxf) + &
         !       fina(pf - 1)           +       fina(pf + 1) + &
         !       fina(pf + Nxf - 1) +       fina(pf + Nxf + 1) )/8.0d0
      
         grossa(pg) = ( fina(pf) + fina(pf + 1) )/2.0d0
                          
      END DO
      
    END SUBROUTINE restringe_v
    
    SUBROUTINE restricao_v(t, level)
    
        USE Variaveis_Solvers_V
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: t, level
        
        CALL restringe_v( f(lv_if(level,1):lv_if(level,2), t),               &
                          f(lv_if(level+1, 1):lv_if(level+1, 2), t),         &
                          lv_nx(level),                                      &
                          lv_nx(level+1), lv_ny(level+1) )
        
    
    END SUBROUTINE restricao_v
    
    
    SUBROUTINE prolonga_v (fina, grossa, Nxf, Nxg, Nyg )
      IMPLICIT NONE
        
      INTEGER, INTENT(IN)                            :: Nxf, Nxg, Nyg
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN)     :: grossa
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)  :: fina

      ! Declaracao de variaveis interna a esta rotina
      INTEGER :: p, pf, pg
      
      ! Volumes no CENTRO recebem stencil 5, mas ja foram considerados
      DO p = 1, (Nxg-2)*(Nyg-2)
         
         pg = p + 2*FLOOR( DBLE(p-1)/DBLE(Nxg-2) ) + Nxg + 1
         ! Deduzido de if = 2ig-1, jg = 2jg e Nf = 2Nxg
         pf = 2 * pg + 2 * (FLOOR( DBLE(pg-1)/Nxg ) + 1) * Nxg - 1
         
         ! 6 velocidades u na malha fina para serem atualizados:
         ! pf, pf - 1, pf + 1, pf + Nxf, pf + Nxf - 1, pf + Nxf + 1
         fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg)                            + grossa(pg - 1) )/4.0d0
         fina(pf + 1)       = fina(pf + 1 )      + ( 3.0d0 * grossa(pg)                            + grossa(pg + 1) )/4.0d0
         fina(pf + Nxf)     = fina(pf + Nxf)     + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + Nxg) + grossa(pg - 1) + grossa(pg + Nxg - 1) )/8.0d0
         fina(pf - Nxf)     = fina(pf - Nxf)     + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - Nxg) + grossa(pg - 1) + grossa(pg - Nxg - 1) )/8.0d0                
         fina(pf - Nxf + 1) = fina(pf - Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - Nxg) + grossa(pg + 1) + grossa(pg - Nxg + 1) )/8.0d0
         fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + Nxg) + grossa(pg + 1) + grossa(pg + Nxg + 1) )/8.0d0
         
      END DO
      
      ! Canto SOUTH WEST
      pg = 1
      pf = 2 * pg + 2 * (FLOOR( DBLE(pg-1)/Nxg ) + 1) * Nxg - 1
      ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
      ! grossa(pg - 1) =  -3.0d0 * grossa(pg) + grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2)  
      ! grossa(pg - Nxg) = 0 (no contorno)
      ! grossa(pg - Nxg - 1) = 0 (extrapolacao do contorno)
      ! grossa(pg - Nxg + 1) = 0 (esta no contorno)
      ! grossa(pg + Nxg - 1) = -3.0d0 * grossa(pg + Nxg) + grossa(pg + Nxg + 1) - (1.0d0/5.0d0) * grossa(pg + Nxg + 2)  
      
      fina(pf)           = fina(pf)           + ( grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2)   )/4.0d0
      fina(pf + 1)       = fina(pf + 1 )      + ( 3.0d0 * grossa(pg) + grossa(pg + 1) )/4.0d0
      fina(pf + Nxf)     = fina(pf + Nxf)     + ( grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2) + grossa(pg + Nxg + 1) - (1.0d0/5.0d0) * grossa(pg + Nxg + 2) )/8.0d0
      fina(pf - Nxf)     = fina(pf - Nxf)     + ( grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2) )/8.0d0                
      fina(pf - Nxf + 1) = fina(pf - Nxf + 1) + ( 3.0d0 * grossa(pg) + grossa(pg + 1) )/8.0d0
      fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + Nxg) + grossa(pg + 1) + grossa(pg + Nxg + 1) )/8.0d0
      
      ! Contorno SOUTH recebe stencil 2
      DO pg = 2, Nxg - 1         
          pf = 2 * pg + 2 * (FLOOR( DBLE(pg-1)/Nxg ) + 1) * Nxg - 1
          ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
          ! grossa(pg - Nxg) = 0 (no contorno)
          ! grossa(pg - Nxg - 1) = 0 (esta no contorno)
          ! grossa(pg - Nxg + 1) = 0 (esta no contorno)
          fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg)                            + grossa(pg - 1) )/4.0d0
          fina(pf + 1)       = fina(pf + 1 )      + ( 3.0d0 * grossa(pg)                            + grossa(pg + 1) )/4.0d0
          fina(pf + Nxf)     = fina(pf + Nxf)     + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + Nxg) + grossa(pg - 1) + grossa(pg + Nxg - 1) )/8.0d0
          fina(pf - Nxf)     = fina(pf - Nxf)     + ( 3.0d0 * grossa(pg) + grossa(pg - 1) )/8.0d0                
          fina(pf - Nxf + 1) = fina(pf - Nxf + 1) + ( 3.0d0 * grossa(pg) + grossa(pg + 1) )/8.0d0
          fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + Nxg) + grossa(pg + 1) + grossa(pg + Nxg + 1) )/8.0d0
       
      END DO  
      
      ! Canto SOUTH EAST  recebe stencil 3
      pg = Nxg
      pf = 2 * pg + 2 * (FLOOR( DBLE(pg-1)/Nxg ) + 1) * Nxg - 1
      ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
      ! grossa(pg + 1) =  -3.0d0 * grossa(pg) + grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2)  
      ! grossa(pg - Nxg) = 0 (no contorno)
      ! grossa(pg - Nxg - 1) = 0 (esta no contorno)
      ! grossa(pg - Nxg + 1) = 0 (extrapolacao do contorno)
      ! grossa(pg + Nxg + 1) = -3.0d0 * grossa(pg + Nxg) + grossa(pg + Nxg - 1) - (1.0d0/5.0d0) * grossa(pg + Nxg - 2)  
      fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg)                            + grossa(pg - 1) )/4.0d0
      fina(pf + 1)       = fina(pf + 1 )      + ( grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2)  )/4.0d0
      fina(pf + Nxf)     = fina(pf + Nxf)     + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + Nxg) + grossa(pg - 1) + grossa(pg + Nxg - 1) )/8.0d0
      fina(pf - Nxf)     = fina(pf - Nxf)     + ( 3.0d0 * grossa(pg) + grossa(pg - 1) )/8.0d0                
      fina(pf - Nxf + 1) = fina(pf - Nxf + 1) + ( grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2) )/8.0d0
      fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + grossa(pg + 1) + grossa(pg + Nxg - 1) - (1.0d0/5.0d0) * grossa(pg + Nxg - 2)   )/8.0d0
                 
      ! Contorno WEST recebe stencil 4
      DO pg = 1 + Nxg, Nxg *( Nyg - 2 ) + 1, Nxg
          pf = 2 * pg + 2 * (FLOOR( DBLE(pg-1)/Nxg ) + 1) * Nxg - 1
          ! grossa(pg - 1) =  -3.0d0 * grossa(pg) + grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2)  
          ! grossa(pg - Nxg - 1) = -3.0d0 * grossa(pg - Nxg) + grossa(pg - Nxg + 1) - (1.0d0/5.0d0) * grossa(pg - Nxg + 2)  
          ! grossa(pg + Nxg - 1) = -3.0d0 * grossa(pg + Nxg) + grossa(pg + Nxg + 1) - (1.0d0/5.0d0) * grossa(pg + Nxg + 2)  
          fina(pf)           = fina(pf)           + ( grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2)   )/4.0d0
          fina(pf + 1)       = fina(pf + 1 )      + ( 3.0d0 * grossa(pg)                            + grossa(pg + 1) )/4.0d0
          fina(pf + Nxf)     = fina(pf + Nxf)     + ( grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2) + grossa(pg + Nxg + 1) - (1.0d0/5.0d0) * grossa(pg + Nxg + 2)   )/8.0d0
          fina(pf - Nxf)     = fina(pf - Nxf)     + ( grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2) + grossa(pg - Nxg + 1) - (1.0d0/5.0d0) * grossa(pg - Nxg + 2)  )/8.0d0                
          fina(pf - Nxf + 1) = fina(pf - Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - Nxg) + grossa(pg + 1) + grossa(pg - Nxg + 1) )/8.0d0
          fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + Nxg) + grossa(pg + 1) + grossa(pg + Nxg + 1) )/8.0d0
          
      END DO
      
      ! Contorno EAST recebe stencil 6
      DO pg = 2 * Nxg, Nxg * (Nyg - 1), Nxg
          pf = 2 * pg + 2 * (FLOOR( DBLE(pg-1)/Nxg ) + 1) * Nxg - 1
          ! grossa(pg + 1) =  -3.0d0 * grossa(pg) + grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2)  
          ! grossa(pg - Nxg + 1) = -3.0d0 * grossa(pg - Nxg) + grossa(pg - Nxg - 1) - (1.0d0/5.0d0) * grossa(pg - Nxg - 2)  
          ! grossa(pg + Nxg + 1) = -3.0d0 * grossa(pg + Nxg) + grossa(pg + Nxg - 1) - (1.0d0/5.0d0) * grossa(pg + Nxg - 2)  
          fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg)                            + grossa(pg - 1) )/4.0d0
          fina(pf + 1)       = fina(pf + 1 )      + ( grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2)   )/4.0d0
          fina(pf + Nxf)     = fina(pf + Nxf)     + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg + Nxg) + grossa(pg - 1) + grossa(pg + Nxg - 1) )/8.0d0
          fina(pf - Nxf)     = fina(pf - Nxf)     + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - Nxg) + grossa(pg - 1) + grossa(pg - Nxg - 1) )/8.0d0                
          fina(pf - Nxf + 1) = fina(pf - Nxf + 1) + ( grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2) + grossa(pg - Nxg - 1) - (1.0d0/5.0d0) * grossa(pg - Nxg - 2) )/8.0d0
          fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2) + grossa(pg + Nxg - 1) - (1.0d0/5.0d0) * grossa(pg + Nxg - 2) )/8.0d0
          
      END DO
      
      ! Canto NORTH WEST recebe stencil 7
      pg = Nxg *( Nyg - 1 ) + 1
      pf = 2 * pg + 2 * (FLOOR( DBLE(pg-1)/Nxg ) + 1) * Nxg - 1
      ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
      ! grossa(pg - 1) =  -3.0d0 * grossa(pg) + grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2)  
      ! grossa(pg + Nxg) = 0 (no contorno)
      ! grossa(pg + Nxg - 1) = 0 (extrapolacao do contorno)
      ! grossa(pg + Nxg + 1) = 0 (esta no contorno)
      ! grossa(pg - Nxg - 1) = -3.0d0 * grossa(pg + Nxg) + grossa(pg + Nxg + 1) - (1.0d0/5.0d0) * grossa(pg + Nxg + 2)  
      fina(pf)           = fina(pf)           + ( grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2)   )/4.0d0
      fina(pf + 1)       = fina(pf + 1 )      + ( 3.0d0 * grossa(pg)                            + grossa(pg + 1) )/4.0d0
      fina(pf + Nxf)     = fina(pf + Nxf)     + ( grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2)   )/8.0d0
      fina(pf - Nxf)     = fina(pf - Nxf)     + ( grossa(pg + 1) - (1.0d0/5.0d0) * grossa(pg + 2) + grossa(pg + Nxg + 1) - (1.0d0/5.0d0) * grossa(pg + Nxg + 2) )/8.0d0                
      fina(pf - Nxf + 1) = fina(pf - Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - Nxg) + grossa(pg + 1) + grossa(pg - Nxg + 1) )/8.0d0
      fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + grossa(pg + 1) )/8.0d0
            
            
      ! Contorno NORTH recebe stencil 8
      DO pg = Nxg *( Nyg - 1 ) + 2, Nxg * Nyg - 1
          pf = 2 * pg + 2 * (FLOOR( DBLE(pg-1)/Nxg ) + 1) * Nxg - 1
          
          ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
          ! grossa(pg + Nxg) = 0 (no contorno)
          ! grossa(pg + Nxg - 1) = 0 (esta no contorno)
          ! grossa(pg + Nxg + 1) = 0 (esta no contorno)
          fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg)                            + grossa(pg - 1) )/4.0d0
          fina(pf + 1)       = fina(pf + 1 )      + ( 3.0d0 * grossa(pg)                            + grossa(pg + 1) )/4.0d0
          fina(pf + Nxf)     = fina(pf + Nxf)     + ( 3.0d0 * grossa(pg) + grossa(pg - 1) )/8.0d0
          fina(pf - Nxf)     = fina(pf - Nxf)     + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - Nxg) + grossa(pg - 1) + grossa(pg - Nxg - 1) )/8.0d0                
          fina(pf - Nxf + 1) = fina(pf - Nxf + 1) + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - Nxg) + grossa(pg + 1) + grossa(pg - Nxg + 1) )/8.0d0
          fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( 3.0d0 * grossa(pg) + grossa(pg + 1) )/8.0d0
      
      END DO
      
      ! Canto NORTH EAST recebe stencil 9
      pg = Nxg * Nyg
      pf = 2 * pg + 2 * (FLOOR( DBLE(pg-1)/Nxg ) + 1) * Nxg - 1
      
      ! Usando extrapolacao cubica e o fato de que o residuo e 0 nos contornos
      ! grossa(pg + 1) =  -3.0d0 * grossa(pg) + grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2)  
      ! grossa(pg + Nxg) = 0 (no contorno)
      ! grossa(pg + Nxg - 1) = 0 (extrapolacao do contorno)
      ! grossa(pg + Nxg + 1) = 0 (esta no contorno)
      ! grossa(pg - Nxg + 1) = -3.0d0 * grossa(pg - Nxg) + grossa(pg - Nxg - 1) - (1.0d0/5.0d0) * grossa(pg - Nxg - 2)
      fina(pf)           = fina(pf)           + ( 3.0d0 * grossa(pg)                            + grossa(pg - 1) )/4.0d0
      fina(pf + 1)       = fina(pf + 1 )      + ( grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2) )/4.0d0
      fina(pf + Nxf)     = fina(pf + Nxf)     + ( 3.0d0 * grossa(pg) + grossa(pg - 1) )/8.0d0
      fina(pf - Nxf)     = fina(pf - Nxf)     + ( 3.0d0 * grossa(pg) + 3.0d0 * grossa(pg - Nxg) + grossa(pg - 1) + grossa(pg - Nxg - 1) )/8.0d0                
      fina(pf - Nxf + 1) = fina(pf - Nxf + 1) + ( grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2) + grossa(pg - Nxg - 1) - (1.0d0/5.0d0) * grossa(pg - Nxg - 2) )/8.0d0
      fina(pf + Nxf + 1) = fina(pf + Nxf + 1) + ( grossa(pg - 1) - (1.0d0/5.0d0) * grossa(pg - 2) )/8.0d0  
      
   END SUBROUTINE prolonga_V
    
    
   SUBROUTINE prolongacao_v(t, level)
    
        USE Variaveis_Solvers_V
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: t, level
        
        CALL prolonga_v( u(lv_if(level-1,1):lv_if(level-1,2), t), &
                         u(lv_if(level, 1):lv_if(level, 2), t),   &
                         lv_nx(level-1),                       &
                         lv_nx(level), lv_ny(level) )
        
    
    END SUBROUTINE prolongacao_v
    
END MODULE MG_GS_V