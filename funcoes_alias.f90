MODULE FuncoesAlias

   USE FuncoesAbstratas
   
   PROCEDURE (FuncaoTempoVolume),  POINTER  :: cc_u_w, cc_u_e, cc_u_s, cc_u_n, cc_v_w, cc_v_e, cc_v_s, cc_v_n, &
                                               pressao_analitica, pressao_inicial, &
                                               u_analitica, u_inicial, &
                                               v_analitica, v_inicial
                                              
   PROCEDURE (FuncaoTempoVolume2), POINTER  :: extrapola_u_S, extrapola_u_N, extrapola_v_W, extrapola_v_E

END MODULE FuncoesAlias
