//======================================================================================
// SUBROTINA PARA A APLICAÇÃO DO MÉTODO DOS GRADIENTES CONJUGADOS
//======================================================================================
void MGC ( int npa , double *pa , double *pra , int npf , double *pf , double *prf ,
	   double ***phi , double ***phir , double ***psi )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;                     // VARIÁVEL AUXILIAR
  int intermax ;              // VARIÁVEL PARA A CONTAGEM DO NÚMERO DE ITERAÇÕES
  double b1a , b2a ;          // VARIÁVEIS INTERMEDIÁRIAS PARA O CÁLCULO DO BETA
  double b1f , b2f ;          // VARIÁVEIS INTERMEDIÁRIAS PARA O CÁLCULO DO BETA
  double betaa , betaf ;      // AUXILIARES PARA O CÁLCULO DA PRÓXIMA DIREÇÃO DE BUSCA
  double erro ;               // VARIÁVEL DE ERRO RELATIVO PARA A PARADA DO MGS
  double e1 , e2 ;            // VARIÁVEIS AUXILIARES PARA O CÁLCULO DO ERRO
  double tol ;                // TOLERÂNCIA PARA O ERRO RELATIVO
  double *pa_old ;            // PONTEIRO PARA O VETOR DE PARÂMETROS ANTERIOR (ALFA)
  double *pf_old ;            // PONTEIRO PARA O VETOR DE PARÂMETROS ANTERIOR (FONTE)
  double *ga_old , *ga_new ;  // PONTEIROS PARA VETORES GRADIENTE DO FUNCIONAL INTEGRAL
  double *gf_old , *gf_new ;  // PONTEIROS PARA VETORES GRADIENTE DO FUNCIONAL INTEGRAL
  double *gra_old , *gra_new ;
  double *grf_old , *grf_new ;
  double *sa_old , *sa_new ;  // PONTEIROS PARA A DETERMINAÇÃO DA DIREÇÃO DE BUSCA
  double *sf_old , *sf_new ;
  double *sra , *srf ;        // PONTEIROS PARA O ARMAZENAMENTO DA DIREÇÃO DE BUSCA
  FILE *pc ;                  // PONTEIRO PARA O ARMAZENAMENTO DOS PARÂMETROS CALCULADOS
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA PARA OS VETORES
  //====================================================================================
  AVR ( &pa_old , npa ) ; AVR ( &pf_old , npf ) ;
  AVR ( &ga_old , npa ) ; AVR ( &gf_old , npf ) ;
  AVR ( &ga_new , npa ) ; AVR ( &gf_new , npf ) ;
  AVR ( &gra_old , npa ) ; AVR ( &gra_new , npa ) ;
  AVR ( &grf_old , npf ) ; AVR ( &grf_new , npf ) ;
  AVR ( &sa_old , npa ) ; AVR ( &sf_old , npf ) ;
  AVR ( &sa_new , npa ) ; AVR ( &sf_new , npf ) ;
  AVR ( &sra , npa ) ; AVR ( &srf , npf ) ;
  //================================================================================
  // ABERTURA DO ARQUIVO PARA ARMAZENAMENTO DOS PARÂMETROS INICIAIS
  //================================================================================
  pc = fopen ( "graficos/parametros.dat", "w" );
  //====================================================================================
  // VERIFICAÇÃO DA EXISTÊNCIA DA PASTA GRÁFICOS
  //====================================================================================
  if ( pc == NULL )
    {
      system ( "mkdir graficos" ) ;
      pc = fopen ( "graficos/parametros.dat", "w" );
    }
  //====================================================================================
  // INICIALIZAÇÃO DE VARIÁVEIS
  //====================================================================================
  intermax = 0 ;
  erro = 1.0 ;
  tol = 1e-9 ;
  fprintf ( pc , "% d\t" , intermax ) ;
  for ( i = 0 ; i < npa ; i ++ )
    {
      pa [ i ] = pa_old [ i ] = 0.01 ;
      fprintf ( pc , "% lf\t" , pa [ i ] ) ;
    }
  for ( i = 0 ; i < npf ; i ++ )
    {
      pf [ i ] = pf_old [ i ] =  0.01 ;
      fprintf ( pc , "% lf\t" , pf [ i ] ) ;
    }
  fprintf ( pc , "% e\t\n" , erro ) ;
  if ( GRAD == 0 ) GRADDF ( npa , pa , npf , pf , phi , phir , ga_old , gf_old ) ;
  else GRADMA ( npa , pa , npf , pf , phi , phir , psi , ga_old , gf_old ) ;
  for ( i = 0 ; i < npa ; i ++ )
    {
      sa_old [ i ] = - ga_old [ i ] ;
    }
  for ( i = 0 ; i < npf ; i ++ )
    {
      sf_old [ i ] = - gf_old [ i ] ;
    }
  time_t tempo1, tempo2;
  //====================================================================================
  // APLICAÇÃO DO MÉTODO DE GRADIENTES CONJUGADOS
  //====================================================================================
  while ( erro > tol && intermax < 101 )
    {
      time (&tempo1);
      //================================================================================
      // INICIALIZAÇÃO DE VARIÁVEIS
      //================================================================================
      b1a = b2a = betaa = 0.0 ;
      b1f = b2f = betaf = 0.0 ;
      e1 = e2 = erro = 0.0 ;
      //================================================================================
      // PRIMEIRO PASSO DO MÉTODO
      //================================================================================
      RA ( npa , sa_old , npf , sf_old , phi , phir , pa , pf ) ;
      //================================================================================
      // SEGUNDO PASSO DO MÉTODO
      //================================================================================
      if ( GRAD == 0 ) GRADDF ( npa , pa , npf , pf , phi , phir , ga_new , gf_new ) ;
      else GRADMA ( npa , pa , npf , pf , phi , phir , psi , ga_new , gf_new ) ;
      //================================================================================
      // TERCEIRO PASSO DO MÉTODO
      //================================================================================
      //================================================================================
      // CÁLCULO DO BETA
      //================================================================================
      if ( MG == 0 )
	{
	  //============================================================================
	  // APLICAÇÃO DA FÓRMULA DE FLETCHER E REEVES ( ALFA )
	  //============================================================================
	  for ( i = 0 ; i < npa ; i ++ )
	    {
	      b1a += ga_new [ i ] * ga_new [ i ] ;
	      b2a += ga_old [ i ] * ga_old [ i ] ;
	    }
	  //============================================================================
	  // APLICAÇÃO DA FÓRMULA DE FLETCHER E REEVES ( TERMO FONTE )
	  //============================================================================
	  for ( i = 0 ; i < npf ; i ++ )
	    {
	      b1f += gf_new [ i ] * gf_new [ i ] ;
	      b2f += gf_old [ i ] * gf_old [ i ] ;
	    }
	}
      else if ( MG == 1 )
	{
	  //============================================================================
	  // APLICAÇÃO DA FÓRMULA DE POLAK E RIBIERE ( ALFA )
	  //============================================================================
	  for ( i = 0 ; i < npa ; i ++ )
	    {
	      b1a += ( ga_new [ i ] - ga_old [ i ] ) * ga_new [ i ] ;
	      b2a += ga_old [ i ] * ga_old [ i ] ;
	    }
	  //============================================================================
	  // APLICAÇÃO DA FÓRMULA DE POLAK E RIBIERE ( TERMO FONTE )
	  //============================================================================
	  for ( i = 0 ; i < npf ; i ++ )
	    {
	      b1f += ( gf_new [ i ] - gf_old [ i ] ) * gf_new [ i ] ;
	      b2f += gf_old [ i ] * gf_old [ i ] ;
	    }
	}
      betaa = b1a / b2a ; betaf = b1f / b2f ;
      //================================================================================
      // QUARTO PASSO DO MÉTODO
      //================================================================================
      for ( i = 0 ; i < npa ; i ++ )
	{
	  sa_new [ i ] = - ga_new [ i ] + betaa * sa_old [ i ] ;
	}
      for ( i = 0 ; i < npf ; i ++ )
	{
	  sf_new [ i ] = - gf_new [ i ] + betaf * sf_old [ i ] ;
	}
      //================================================================================
      // DETERMINAÇÃO DO ERRO RELATIVO
      //================================================================================
      //for ( i = 0 ; i < npa ; i ++ )
	//{
        //e1 += pow ( ( pa [ i ] - pa_old [ i ] ) , 2 ) ;
        //e2 += pow ( pa [ i ] , 2 ) ;
	//}
      //for ( i = 0 ; i < npf ; i ++ )
	//{
        //e1 += pow ( ( pf [ i ] - pf_old [ i ] ) , 2 ) ;
        //e2 += pow ( pf [ i ] , 2 ) ;
	//}
      //e1 = sqrt ( e1 ) ;
      //e2 = sqrt ( e2 ) ;
      //erro = e1 / e2 ;
      erro = IP ( npa , pa , npf , pf , phi , phir ) ;
      //================================================================================
      // ATUALIZAÇÃO DO NÚMERO DE ITERAÇÕES
      //================================================================================
      intermax ++ ;
      //================================================================================
      // ATUALIZAÇÃO DO VETOR DE PARÂMETROS
      //================================================================================
      time (&tempo2);
      fprintf ( pc , "% d\t" , intermax ) ;
      printf ( "Iteração: %d\n", intermax ) ;
      for ( i = 0 ; i < npa ; i ++ )
	{
	  pa_old [ i ] = pa [ i ] ;
	  ga_old [ i ] = ga_new [ i ] ;
	  sa_old [ i ] = sa_new [ i ] ;
	  fprintf ( pc , "% lf\t" , pa [ i ] ) ;
	  printf ( "PA: % lf\t" , pa [ i ] ) ;
	}
      for ( i = 0 ; i < npf ; i ++ )
	{
	  pf_old [ i ] = pf [ i ] ;
	  gf_old [ i ] = gf_new [ i ] ;
	  sf_old [ i ] = sf_new [ i ] ;
	  fprintf ( pc , "% lf\t" , pf [ i ] ) ;
	  printf ( "PF: % lf\t" , pf [ i ] ) ;
	}
      fprintf ( pc , "% e\t" , erro ) ;
      fprintf ( pc , "% f\t\n" , difftime (tempo2, tempo1) ) ;
      printf ("\nTEMPO: %f\n\n", difftime (tempo2, tempo1));
      printf ( "\nERRO: % e\n\n" , erro ) ;
    }
  fclose ( pc ) ;
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA DOS ESPAÇOS ALOCADOS
  //====================================================================================
  LVR ( &pa_old ) ; LVR ( &pf_old ) ;
  LVR ( &ga_old ) ; LVR ( &gf_old ) ;
  LVR ( &ga_new ) ; LVR ( &gf_new ) ;
  LVR ( &gra_new ) ; LVR ( &gra_old ) ;
  LVR ( &grf_new ) ; LVR ( &grf_old ) ;
  LVR ( &sa_old ) ; LVR ( &sf_old ) ;
  LVR ( &sa_new ) ; LVR ( &sf_new ) ;
  LVR ( &sra ) ; LVR ( &srf ) ;
}
