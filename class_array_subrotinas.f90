MODULE ClassArrayFonteComplementosU
   
   USE FuncoesAbstratas

   IMPLICIT  NONE
   
   PRIVATE

   TYPE, PUBLIC :: ArrayFonteComplementosU

      PROCEDURE (SubrotinaTempoVolume), POINTER, NOPASS  :: CalculaComplementos

      CONTAINS
      
      procedure :: init
      
   END TYPE ArrayFonteComplementosU
   
   CONTAINS
   
      SUBROUTINE init(self, stencil)

         USE FontesSubrotinas

         CLASS ( ArrayFonteComplementosU ) :: self
         INTEGER, INTENT(IN) :: stencil
         
         SELECT CASE ( stencil )
         
            CASE ( 1 )
            
               self % CalculaComplementos => ComplementosU_SW
               
            CASE ( 2 )
               
               self % CalculaComplementos => ComplementosU_S

            CASE ( 3 )
            
               self % CalculaComplementos => ComplementosU_SE
                              
            CASE ( 4 )
            
               self % CalculaComplementos => ComplementosU_W
               
            CASE ( 5 )
            
               self % CalculaComplementos => ComplementosU_C
               
            CASE ( 6 )
            
               self % CalculaComplementos => ComplementosU_E
               
            CASE ( 7 )
            
               self % CalculaComplementos => ComplementosU_NW   
               
            CASE ( 8 )
            
               self % CalculaComplementos => ComplementosU_N   
               
            CASE ( 9 )
            
               self % CalculaComplementos => ComplementosU_NE   
                              
         END SELECT
      
      END SUBROUTINE init
     

END MODULE ClassArrayFonteComplementosU



MODULE ClassArrayFonteComplementosV
   
   USE FuncoesAbstratas

   IMPLICIT  NONE
   
   PRIVATE

   TYPE, PUBLIC :: ArrayFonteComplementosV

      PROCEDURE (SubrotinaTempoVolume), POINTER, NOPASS  :: CalculaComplementos

      CONTAINS
      
      procedure :: init
      
   END TYPE ArrayFonteComplementosV
   
   CONTAINS
   
      SUBROUTINE init(self, stencil)

         USE FontesSubrotinas

         CLASS ( ArrayFonteComplementosV ) :: self
         INTEGER, INTENT(IN) :: stencil
         
         SELECT CASE ( stencil )
         
            CASE ( 1 )
            
               self % CalculaComplementos => ComplementosV_SW
               
            CASE ( 2 )
               
               self % CalculaComplementos => ComplementosV_S

            CASE ( 3 )
            
               self % CalculaComplementos => ComplementosV_SE
                              
            CASE ( 4 )
            
               self % CalculaComplementos => ComplementosV_W
               
            CASE ( 5 )
            
               self % CalculaComplementos => ComplementosV_C
               
            CASE ( 6 )
            
               self % CalculaComplementos => ComplementosV_E
               
            CASE ( 7 )
            
               self % CalculaComplementos => ComplementosV_NW   
               
            CASE ( 8 )
            
               self % CalculaComplementos => ComplementosV_N   
               
            CASE ( 9 )
            
               self % CalculaComplementos => ComplementosV_NE   
                              
         END SELECT
      
      END SUBROUTINE init
     

END MODULE ClassArrayFonteComplementosV



MODULE ClassArrayFonteP
   
   USE FuncoesAbstratas

   IMPLICIT  NONE
   
   PRIVATE

   TYPE, PUBLIC :: ArrayFonteP

      PROCEDURE (FuncaoTempoVolume3), POINTER, NOPASS  :: FonteStencil

      CONTAINS
      
      procedure :: init
      
   END TYPE ArrayFonteP
   
   CONTAINS
   
      SUBROUTINE init(self, stencil)

         USE FontesSubrotinas

         CLASS ( ArrayFonteP ) :: self
         INTEGER, INTENT(IN) :: stencil
         
         SELECT CASE ( stencil )
         
            CASE ( 1 )
            
               self % FonteStencil => StencilP_SW
               
            CASE ( 2 )
               
               self % FonteStencil => StencilP_S

            CASE ( 3 )
            
               self % FonteStencil => StencilP_SE
                              
            CASE ( 4 )
            
               self % FonteStencil => StencilP_W
               
            CASE ( 5 )
            
               self % FonteStencil => StencilP_C
               
            CASE ( 6 )
            
               self % FonteStencil => StencilP_E
               
            CASE ( 7 )
            
               self % FonteStencil => StencilP_NW   
               
            CASE ( 8 )
            
               self % FonteStencil => StencilP_N   
               
            CASE ( 9 )
            
               self % FonteStencil => StencilP_NE   
                              
         END SELECT
      
      END SUBROUTINE init
     

END MODULE ClassArrayFonteP