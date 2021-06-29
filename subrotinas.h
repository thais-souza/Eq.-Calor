//==================================================================================
// DECLARAÇÃO DAS SUBROTINAS
//==================================================================================
#include"subrotinas/ad.c"           // CONTÉM SUBROTINAS PARA ALOCAÇÃO DINÂMICA
#include"subrotinas/leitura.c"      // SUBROTINA PARA A LEITURA EXTERNA
#include"subrotinas/inicial.c"      // FUNÇÃO PARA O CÁLCULO DA SOLUÇÃO INICIAL
#include"subrotinas/ed2d.c"         // RESOLUÇÃO DO PROBLEMA PRIMAL E DUAL
#include"subrotinas/mrk.c"          // RESOLUÇÃO DAS EDPs VIA MÉTODO DE RUNGE-KUTTA
#include"subrotinas/ordem_mrk.c"    // VERIFICAÇÃO COMPUTACIONAL DA ORDEM DOS MRKs
#include"subrotinas/residuo.c"      // CALCULA OS RESÍDUOS DOS PROBLEMAS PRIMAL/DUAL
#include"subrotinas/fronteira.c"    // PARA TRATAMENTO DAS CONDIÇÕES DE FRONTEIRA
#include"subrotinas/grafico.c"      // GERAÇÃO DE GRÁFICOS VIA GNUPLOT
#include"subrotinas/analitica.c"    // SOLUÇÃO ANALÍTICA DE TESTE
#include"subrotinas/fxyt.c"         // FUNÇÃO TERMO FONTE
#include"subrotinas/alfa.c"         // FUNÇÃO COEFICIENTE DE DIFUSÃO
#include"subrotinas/ra.c"           // MÉTODO DA RAZÃO ÁUREA
#include"subrotinas/er.c"           // CALCULA O ERRO RELATIVO ENTRE DOIS VETORES
#include"subrotinas/mgc.c"          // MÉTODO DOS GRADIENTES CONJUGADOS
#include"subrotinas/integral.c"     // CALCULA AS INTEGRAIS TRIPLAS
#include"subrotinas/gradiente.c"    // CALCULA OS VETORES GRADIENTES
