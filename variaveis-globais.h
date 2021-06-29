//==================================================================================
// DECLARAÇÃO DE VARIÁVEIS GLOBAIS
//==================================================================================
int NX ;                     // NÚMERO DE DIVISÕES NO INTERVALO [0,LX] NO EIXO X
int NY ;                     // NÚMERO DE DIVISÕES NO INTERVALO [0,LY] NO EIXO Y
int NT ;                     // NÚMERO DE DIVISÕES NO INTERVALO [0,LT] NO EIXO Z
int ORDEM ;                  // ORDEM DO MÉTODO DE RUNGE-KUTTA
int GA ;                     // BANDEIRA DE DECISÃO PARA GERAÇÃO DE ANIMAÇÕES
int IGA ;                    // INTERVALO PARA GERAÇÃO DAS ANIMAÇÕES
int GRAD ;                   // ESCOLHA NA FORMA DE CALCULAR OS GRADIENTES
int MG ;                     // ESCOLHA NA FORMA DO MÉTODO DOS GRADIENTES CONJUGADOS
int F1 ;                     // CONDIÇÃO DE FRONTEIRA 1
int F2 ;                     // CONDIÇÃO DE FRONTEIRA 2
int F3 ;                     // CONDIÇÃO DE FRONTEIRA 3
int F4 ;                     // CONDIÇÃO DE FRONTEIRA 4
double LX ;                  // DIMENSÃO ESPACIAL NA DIREÇÃO DO EIXO X
double LY ;                  // DIMENSÃO ESPACIAL NA DIREÇÃO DO EIXO Y
double LT ;                  // DIMENSÃO TEMPORAL
double DX ;                  // ESPAÇAMENTO NO INTERVALO [0,LX] NO EIXO X
double DY ;                  // ESPAÇAMENTO NO INTERVALO [0,LY] NO EIXO Y
double DT ;                  // ESPAÇAMENTO NO INTERVALO [0,LT] NO EIXO Z
double W ;                   // CONSTANTE MULTIPLICATIVA PARA O CÁLCULO DO RESÍDUO
//==================================================================================
// FUNÇÕES
//==================================================================================
void LEITURA () ;            // LEITURA DOS DADOS EXTERNOS
void INICIALP () ;           // FUNÇÃO PARA GERAR A SOLUÇÃO INICIAL PRIMAL
void INICIALD () ;           // FUNÇÃO PARA GERAR A SOLUÇÃO INICIAL DUAL
void ED2DP () ;              // RESOLUÇÃO DA EQUAÇÃO DE DIFUSÃO DO PROBLEMA PRIMAL
void ED2DD () ;              // RESOLUÇÃO DA EQUAÇÃO DE DIFUSÃO DO PROBLEMA DUAL
void MRK20P () ;             // RESOLUÇÃO DA EQUAÇÃO DE DIFUSÃO VIA MRK DE ORDEM 2
void MRK30P () ;             // RESOLUÇÃO DA EQUAÇÃO DE DIFUSÃO VIA MRK DE ORDEM 3
void MRK40P () ;             // RESOLUÇÃO DA EQUAÇÃO DE DIFUSÃO VIA MRK DE ORDEM 4
void MRK20D () ;             // RESOLUÇÃO DO PROBLEMA DUAL VIA MRK DE ORDEM 2
void MRK30D () ;             // RESOLUÇÃO DO PROBLEMA DUAL VIA MRK DE ORDEM 3
void MRK40D () ;             // RESOLUÇÃO DO PROBLEMA DUAL VIA MRK DE ORDEM 4
void ORDEM_MRK () ;          // VERIFICAÇÃO DA ORDEM DOS MRKs
void RP () ;                 // RESÍDUO ASSOCIADO AO PROBLEMA PRIMAL
void RD () ;                 // RESÍDUO ASSOCIADO AO PROBLEMA DUAL/ADJUNTO
void FRONTEIRAP () ;         // TRATAMENTO DA SOLUÇÃO NA FRONTEIRA PRIMAL
void FRONTEIRAD () ;         // TRATAMENTO DA SOLUÇÃO NA FRONTEIRA DUAL
void GRAFICO () ;            // GERAÇÃO DE GRÁFICOS VIA GNUPLOT
void MGC () ;                // MÉTODO DOS GRADIENTES CONJUGADOS
void GRADDF () ;             // CALCULA O VETOR GRADIENTE VIA DIFERENÇA FINITA
void GRADMA () ;             // CALCULA O VETOR GRADIENTE VIA MÉTODO ADJUNTO
double RA () ;               // MÉTODO DA RAZÃO ÁUREA
double PHIXYT () ;           // FUNÇÃO SOLUÇÃO ANALÍTICA PARA TESTE
double PHIXYT0 () ;          // FUNÇÃO SOLUÇÃO INICIAL
double PHIXYT1 () ;          // FUNÇÃO SOLUÇÃO FRONTEIRA 1
double PHIXYT2 () ;          // FUNÇÃO SOLUÇÃO FRONTEIRA 2
double PHIXYT3 () ;          // FUNÇÃO SOLUÇÃO FRONTEIRA 3
double PHIXYT4 () ;          // FUNÇÃO SOLUÇÃO FRONTEIRA 4
double PSIXYT () ;           // FUNÇÃO SOLUÇÃO ANALÍTICA PARA TESTE
double PSIXYT0 () ;          // FUNÇÃO SOLUÇÃO INICIAL
double PSIXYT1 () ;          // FUNÇÃO SOLUÇÃO FRONTEIRA 1
double PSIXYT2 () ;          // FUNÇÃO SOLUÇÃO FRONTEIRA 2
double PSIXYT3 () ;          // FUNÇÃO SOLUÇÃO FRONTEIRA 3
double PSIXYT4 () ;          // FUNÇÃO SOLUÇÃO FRONTEIRA 4
double FXYT () ;             // FUNÇÃO TERMO FONTE
double ALFA () ;             // FUNÇÃO COEFICIENTE DE DIFUSÃO
double IP () ;               // CALCULA A INTEGRAL TRIPLA RELACIONADA AO FUNCIONAL
double ID () ;               // CALCULA A INTEGRAL TRIPLA RELACIONADA AO GRADIENTE
double ER () ;               // CALCULA O ERRO RELATIVO ENTRE DOIS VETORES
