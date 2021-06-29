//======================================================================================
// SUBROTINAS PARA O CÁLCULO DOS VETORES GRADIENTES
//======================================================================================
// SUBROTINA PARA O CÁLCULO DOS VETORES GRADIENTES VIA MÉTODO DE DIFERENÇAS FINITAS
//======================================================================================
void GRADDF ( int npa , double *pa , int npf , double *pf, double ***phi ,
	      double ***phir , double *ga , double *gf )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;             // VARIÁVEL AUXILIAR
  double p1 , p2 ;    // VARIÁVEL AUXILIAR PARA O CÁLCULO DO VETOR GRADIENTE
  double p3 ;         // VARIÁVEL AUXILIAR PARA O CÁLCULO DO VETOR GRADIENTE
  double eps ;        // VALOR DE PERTURBAÇÃO NOS PARÂMETROS
  double norma ;      // VALOR DA NORMA DO VETOR CALCULADO
  FILE *grad ;        // PONTEIRO PARA O ARMAZENAMENTO DO GRADIENTE CALCULADO
  //====================================================================================
  // INICIALIZAÇÃO DE VARIÁVEIS
  //====================================================================================
  norma = 0.0 ;
  eps = 0.001 ;
  p1 = IP ( npa , pa , npf , pf , phi , phir ) ;
  //====================================================================================
  // CÁLCULO DO GRADIENTE DE ALFA
  //====================================================================================
  for ( i = 0 ; i < npa ; i ++ )
    {
      pa [ i ] += eps ; 
      p2 = IP ( npa , pa , npf , pf , phi , phir ) ;
      pa [ i ] += eps ; 
      p3 = IP ( npa , pa , npf , pf , phi , phir ) ;
      pa [ i ] -= 2.0 * eps ;
      ga [ i ] = ( - 3.0 * p1 + 4.0 * p2 - p3 ) / ( 2 * eps ) ;
    }
  //====================================================================================
  // CÁLCULO DO GRADIENTE DO TERMO FONTE
  //====================================================================================
  for ( i = 0 ; i < npf ; i ++ )
    {
      pf [ i ] += eps ; 
      p2 = IP ( npa , pa , npf , pf , phi , phir ) ;
      pf [ i ] += eps ; 
      p3 = IP ( npa , pa , npf , pf , phi , phir ) ;
      pf [ i ] -= 2.0 * eps ;
      gf [ i ] = ( - 3.0 * p1 + 4.0 * p2 - p3 ) / ( 2 * eps ) ;
    }
  //====================================================================================
  // CÁLCULO DA NORMA DO GRADIENTE DE ALFA E NORMALIZAÇÃO DO GRADIENTE DE ALFA
  //====================================================================================
  norma = 0.0 ;
  for ( i = 0 ; i < npa ; i ++ )
    {
      norma += pow ( ga [ i ] , 2 ) ;
    }
  norma = sqrt ( norma ) ;
  for ( i = 0 ; i < npa ; i ++ )
    {
      ga [ i ] /= norma ;
    }
  //====================================================================================
  // CÁLCULO DA NORMA DO GRADIENTE DO TERMO FONTE E NORMALIZAÇÃO DO GRADIENTE
  //====================================================================================
  norma = 0.0 ;
  for ( i = 0 ; i < npf ; i ++ )
    {
      norma += pow ( gf [ i ] , 2 ) ;
    }
  norma = sqrt ( norma ) ;
  for ( i = 0 ; i < npf ; i ++ )
    {
      gf [ i ] /= norma ;
    }
  //====================================================================================
  // ABERTURA DO ARQUIVO PARA O ARMAZENAMENTO DOS GRADIENTES
  //====================================================================================
  grad = fopen ( "graficos/direcao_busca.dat", "a" );
  //====================================================================================
  // VERIFICAÇÃO DA EXISTÊNCIA DA PASTA GRÁFICOS
  //====================================================================================
  if ( grad == NULL )
    {
      system ( "mkdir graficos" ) ;
      grad = fopen ( "graficos/direcao_busca.dat", "a" );
    }
  //====================================================================================
  // ARMAZENAMENTO DOS PARÂMETROS INICIAIS
  //====================================================================================
  for ( i = 0 ; i < npa ; i ++ )
    {
      fprintf ( grad, "% lf\t" , pa [ i ] ) ;
    }
  for ( i = 0 ; i < npf ; i ++ )
    {
      fprintf ( grad, "% lf\t" , pf [ i ] ) ;
    }
  //====================================================================================
  // ARMAZENAMENTO DOS GRADIENTES
  //====================================================================================
  for ( i = 0 ; i < npa ; i ++ )
    {
      fprintf ( grad, "% lf\t" , ( - ga [ i ] ) ) ;
    }
  for ( i = 0 ; i < npf ; i ++ )
    {
      fprintf ( grad, "% lf\t" , ( - gf [ i ] ) ) ;
    }
  fprintf ( grad , "\n" ) ;
  fclose ( grad ) ;
}

//======================================================================================
// SUBROTINA PARA O CÁLCULO DOS VETORES GRADIENTES VIA MÉTODO ADJUNTO
//======================================================================================
void GRADMA ( int npa , double *pa , int npf , double *pf , double ***phi ,
	      double ***phir , double ***psi , double *ga , double *gf )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;             // VARIÁVEL AUXILIAR
  double norma ;      // VALOR DA NORMA DO VETOR CALCULADO
  double *ca ;        // VETOR PARA O ISOLAMENTO DAS VARIÁVEIS EM ALFA
  double *cf ;        // VETOR PARA O ISOLAMENTO DAS VARIÁVEIS NO TERMO FONTE 
  FILE *grad ;        // PONTEIRO PARA O ARMAZENAMENTO DO GRADIENTE CALCULADO
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA
  //====================================================================================
  AVR ( &ca , npa ) ;
  AVR ( &cf , npf ) ;
  //====================================================================================
  // INICIALIZAÇÃO DE VARIÁVEIS
  //====================================================================================
  norma = 0.0 ;
  for ( i = 0 ; i < npa ; i ++ )
    {
      ca [ i ] = 0.0 ;
    }
  for ( i = 0 ; i < npf ; i ++ )
    {
      cf [ i ] = 0.0 ;
    }
  //====================================================================================
  // RESOLUÇÃO DA EQUAÇÃO DIFERENCIAL EM ANÁLISE
  //====================================================================================
  ED2DP ( npa , pa , npf , pf , phi ) ;
  ED2DD ( npa , pa , phi , phir , psi ) ;
  //====================================================================================
  // CÁLCULO DO GRADIENTE DE ALFA
  //====================================================================================
  for ( i = 0 ; i < npa ; i ++ )
    {
      ca [ i ] = 1.0 ;
      ga [ i ] = - ID ( npa , ca , npf , cf , phi , psi , 0 ) ;
      ca [ i ] = 0.0 ;
    }
  //====================================================================================
  // NORMALIZAÇÃO DO GRADIENTE DE ALFA
  //====================================================================================
  for ( i = 0 ; i < npa ; i ++ )
    {
      norma += pow ( ga [ i ] , 2 ) ;
    }
  norma = sqrt ( norma ) ;
  for ( i = 0 ;  i < npa ; i ++ )
    {
      ga [ i ] /= norma ;
    }
  //====================================================================================
  // CÁLCULO DO GRADIENTE DO TERMO FONTE
  //====================================================================================
  for ( i = 0 ; i < npf ; i ++ )
    {
      cf [ i ] = 1.0 ;
      gf [ i ] = - ID ( npa , ca , npf , cf , phi , psi , 1 ) ;
      cf [ i ] = 0.0 ;
    }
  //====================================================================================
  // NORMALIZAÇÃO DO GRADIENTE DE ALFA
  //====================================================================================
  norma = 0.0 ;
  for ( i = 0 ; i < npf ; i ++ )
    {
      norma += pow ( gf [ i ] , 2 ) ;
    }
  norma = sqrt ( norma ) ;
  for ( i = 0 ;  i < npa ; i ++ )
    {
      gf [ i ] /= norma ;
    }
  //====================================================================================
  // ABERTURA DO ARQUIVO PARA O ARMAZENAMENTO DOS GRADIENTES
  //====================================================================================
  grad = fopen ( "graficos/direcao_busca.dat", "a" );
  //====================================================================================
  // VERIFICAÇÃO DA EXISTÊNCIA DA PASTA GRÁFICOS
  //====================================================================================
  if ( grad == NULL )
    {
      system ( "mkdir graficos" ) ;
      grad = fopen ( "graficos/direcao_busca.dat", "a" );
    }
  //====================================================================================
  // ARMAZENAMENTO DOS PARÂMETROS INICIAIS
  //====================================================================================
  for ( i = 0 ; i < npa ; i ++ )
    {
      fprintf ( grad, "% lf\t" , pa [ i ] ) ;
    }
  for ( i = 0 ; i < npa ; i ++ )
    {
      fprintf ( grad, "% lf\t" , pf [ i ] ) ;
    }
  //====================================================================================
  // ARMAZENAMENTO DOS GRADIENTES
  //====================================================================================
  for ( i = 0 ; i < npa ; i ++ )
    {
      fprintf ( grad, "% lf\t" , ( - ga [ i ] ) ) ;
    }
  for ( i = 0 ; i < npf ; i ++ )
    {
      fprintf ( grad, "% lf\t" , ( - gf [ i ] ) ) ;
    }
  fprintf ( grad , "\n" ) ;
  fclose ( grad ) ;
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA
  //====================================================================================
  LVR ( &ca ) ;
  LVR ( &cf ) ;
}
