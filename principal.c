///======================================================================================
// UNIVERSIDADE FEDERAL DE UBERLÂNDIA
// FACULDADE DE MATEMÁTICA
// PROJETO: IDENTIFICAÇÃO DE PARÂMETROS ENVOLVENDO A EQUAÇÃO DA DIFUSÃO BIDIMENSIONAL
// COM TERMO FONTE VIA MÉTODOS DE DIFERENÇAS FINITAS E ADJUNTOS
// ORIENTANDA: THAÍS BARBOSA CAETANO SOUZA
// ORIENTADOR: ALESSANDRO ALVES SANTANA
//======================================================================================
// PRÉ-PROCESSADORES
//======================================================================================
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<omp.h>
//======================================================================================
// CONSTANTES GLOBAIS
//======================================================================================
#define pi 3.141592653589793
//======================================================================================
// VARIÁVEIS GLOBAIS
//======================================================================================
#include"variaveis-globais.h"
//======================================================================================
// SUBROTINAS GLOBAIS
//======================================================================================
#include"subrotinas.h"
//======================================================================================
// PROGRAMA PRINCIPAL
//======================================================================================
int main ( int argc , char *argv [] )
{
  //====================================================================================
  // DECLARAÇÃO E INICIALIZAÇÃO DE VARIÁVEIS
  //====================================================================================
  int npa ;                    // NÚMERO DE PARÂMETROS RELACIONADOS À ALFA
  int npf ;                    // NÚMERO DE PARÂMETROS RELACIONADOS AO TERMO FONTE
  int i , j ;                  // VARIÁVEIS AUXILIARES
  double x , y ;               // VARIÁVEIS AUXILIARES PARA OS DESLOCAMENTOS ESPACIAIS
  double re ;                  // RAZÃO DE ESTABILIDADE
  double alfa ;                // ALFA PARA A ESTABILIDADE DOS MRKs
  double alfa_aux ;            // ALFA AUXILIAR PARA O CÁLCULO DA ESTABILIDADE DOS MRKs
  double beta ;                // BETA PARA A ESTABILIDADE DOS MRKs
  double *pa = NULL ;          // PARÂMETROS A SEREM IDENTIFICADOS EM ALFA
  double *pf = NULL ;          // PARÂMETROS A SEREM IDENTIFICADOS NO TERMO FONTE
  double *pra = NULL ;         // PARÂMETROS DE REFERÊNCIA RELACIONADOS À ALFA
  double *prf = NULL ;         // PARÂMETROS DE REFERÊNCIA RELACIONADOS AO TERMO FONTE
  double ***phir = NULL ;      // SOLUÇÃO DE REFERÊNCIA
  double ***phi = NULL ;       // SOLUÇÃO DO PROBLEMA PRIMAL
  double ***psi = NULL ;       // SOLUÇÃO DO PROBLEMA DUAL
  //====================================================================================
  // LEITURA DA ENTRADA
  //====================================================================================
  LEITURA ( argv , &npa , &pra , &pa , &npf , &prf , &pf ) ;
  //====================================================================================
  // DETERMINAÇÃO DOS INTERVALOS ESPACIAIS
  //====================================================================================
  DX = LX / NX ;
  DY = DX ;
  DT = LT / NT;
  NY = ( int ) ( LY / DY ) ;
  //====================================================================================
  // DETERMINAÇÃO DO ALFA PARA O CÁLCULO DO DX
  //====================================================================================
  /*for ( i = 0 ; i < npa ; i ++ )
    {
      pa [ i ] = 0.4 ;
    }
  alfa = ALFA ( npa , pa , 0.0 , 0.0 ) ;
  for ( j = 0 ; j <= NY ; j ++ )
    {
      y = j * DY ;
      for ( i = 0 ; i <= NX ; i ++ )
	{
	  x = i * DX ;
	  alfa_aux = ALFA ( npa , pa , x , y ) ;
	  if ( alfa_aux > alfa ) alfa = alfa_aux ;
	}
    }
  printf ( "ALFA: %lf\n" , alfa ) ;
  getchar ();
  //====================================================================================
  // DETERMINAÇÃO DO BETA PARA O CÁLCULO DO DX
  //====================================================================================
  beta = 0.0 ;
  if ( ORDEM == 2 ) beta = 2.0 ;
  else if ( ORDEM == 3 ) beta = 2.51 ;
  else if ( ORDEM == 4 ) beta = 2.78 ;
  //====================================================================================
  // DETERMINAÇÃO DO INTERVALO TEMPORAL
  //====================================================================================
  re = 0.24 ;
  printf ( "DX = %lf\nDY = %lf\n" , DX , DY ) ;
  DT = re * beta * pow ( DX , 2 ) * pow ( DY , 2 ) ;
  DT = DT / ( alfa * ( pow ( DX , 2 ) + pow ( DY , 2 ) ) ) ;
  NT = LT / DT ;*/
  W = 1.0 / pow ( DX , 2 ) ;
  //====================================================================================
  // VERIFICAÇÃO DA PARIDADE DO NT - NECESSÁRIA PARA A APLICAÇÃO DE 1/3 DE SIMPSON
  //====================================================================================
  /*if ( NT % 2 != 0 )
    {
      NT += 1 ;
      DT = LT / NT ;
    }
    DT = LT / NT ;*/
  printf ( "DT = %lf\n" , DT ) ;
  printf ( "NT = % d\n" , NT ) ;
  //getchar();
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA PARA AS MATRIZES
  //====================================================================================
  AMR3D ( &phir , NX + 1 , NY + 1 , NT + 1 ) ;
  AMR3D ( &phi , NX + 1 , NY + 1 , NT + 1 ) ;
  AMR3D ( &psi , NX + 1 , NY + 1 , NT + 1 ) ;
  //====================================================================================
  // RESOLUÇÃO ANALÍTICA DA EQUAÇÃO DIFERENCIAL DE REFERÊNCIA
  //====================================================================================
  ED2DP ( npa , pra , npf , prf , phir ) ;
  /*
  int na, nb;
  double da, db, a ,b;
  FILE *func;
  system ( "mkdir graficos" ) ;
  func = fopen ( "graficos/funcional.dat" , "w" );
  na = 40 ;
  nb = 40 ;
  da = ( 0.03 ) / na ;
  db = ( 0.03 ) / nb ;
  a = pra [ 0 ] - 0.015 ;
  b = pra [ 1 ] - 0.015 ;
  for ( j = 0 ; j <= nb ; j ++ )
    {
      pa [ 1 ] = b + db * j ;
      for ( i = 0 ; i <= na ; i ++ )
	{
	  printf ("Rodando %d\n", i);
	  pa [ 0 ] = a + da * i ;
	  fprintf ( func , "% lf\t% lf\t% e\n", pa [ 0 ] , pa [ 1 ] ,
		    IP ( npa , pa , npf , pf , phi , phir ) ) ;
	}
      fprintf ( func , "\n" ) ;
    }
    fclose ( func ) ;
  */
  //====================================================================================
  // APLICAÇÃO DO MÉTODO DE GRADIENTES CONJUGADOS
  //====================================================================================
  //time_t tempoi, tempof;
  //time (&tempoi);
  //MGC ( npa , pa , pra , npf , pf , prf , phi , phir , psi ) ;
  //time (&tempof);
  //printf ("\nTEMPO: %f\n\n", difftime (tempof, tempoi));
  //====================================================================================
  // GERAÇÃO DE GRÁFICOS
  //====================================================================================
  //if ( GA == 1 ) GRAFICO ( phir ) ;
  //====================================================================================
  // VERIFICAÇÃO DA ORDEM DOS MRKs
  //====================================================================================
  ORDEM_MRK ( phir ) ;
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA DOS ESPAÇOS ALOCADOS
  //====================================================================================
  LVR ( &pra ) ;
  LVR ( &prf ) ;
  LVR ( &pa ) ;
  LVR ( &pf ) ;
  LMR3D ( &phir , NX + 1 , NY + 1 ) ;
  LMR3D ( &phi , NX + 1 , NY + 1 ) ;
  LMR3D ( &psi , NX + 1 , NY + 1 ) ;
  //====================================================================================
  // FIM DO PROGRAMA
  //====================================================================================
  return ( 0 ) ;
}
