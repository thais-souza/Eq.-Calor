//======================================================================================
// SUBROTINA PARA A APLICAÇÃO DO MÉTODO DE RAZÃO ÁUREA
//======================================================================================
double RA ( int npa , double *sa , int npf , double *sf , double ***phi ,
	    double ***phir , double *pa , double *pf )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;                  // VARIÁVEL AUXILIAR
  double t ;               // CONSTANTE RAZÃO ÁUREA
  double x1 , x2 ;         // VALORES INTERNOS DO INTERVALO DE BUSCA
  double f1 , f2 ;         // VALOR DO FUNCIONAL SUBMETIDO AOS VALORES INTERNOS
  double a , b ;           // VALORES LIMITES DO INTERVALO DE BUSCA
  double tol ;             // TOLERÂNCIA DE ERRO NA APLICAÇÃO DO MÉTODO
  double alfa ;            // PARÂMETRO DE DIREÇÃO DE BUSCA
  double *px1a , *px2a ;   // PONTEIRO DE ARMAZENAMENTO DOS PARÂMETOS INTERNOS ( ALFA )
  double *px1f , *px2f ;   // PONTEIRO DE ARMAZENAMENTO DOS PARÂMETOS INTERNOS ( FONTE )
  double normaa ;          // CÁLCULO DA NORMA INFINITA RELATIVA AO ALFA
  double normaf ;          // CÁLCULO DA NORMA INFINITA RELATIVA AO TERMO FONTE
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA PARA OS VETORES
  //====================================================================================
  AVR ( &px1a , npa ) ;
  AVR ( &px2a , npa ) ;
  AVR ( &px1f , npf ) ;
  AVR ( &px2f , npf ) ;
  for ( i = 0 ; i < npa ; i ++ )
    {
      printf ( "%e\t" , sa [ i ] ) ;
    }
  for ( i = 0 ; i < npf ; i ++ )
    {
      printf ( "%e\t" , sf [ i ] ) ;
    }
  printf ( "\n" ) ;
  //====================================================================================
  // INICIALIZAÇÃO DE VARIÁVEIS
  //====================================================================================
  tol = 1e-6 ;
  normaa = pa [ 0 ] ;
  for ( i = 0 ; i < npa ; i ++ )
    {
      if ( pa [ i ] > normaa ) normaa = pa [ i ] ;
    }
  a = 0.0 ;
  b = 0.15 ;
  t = ( sqrt ( 5.0 ) - 1.0 ) / 2.0 ;
  x1 = a + ( 1.0 - t ) * ( b - a ) ;
  x2 = a + t * ( b - a ) ;
  for ( i = 0 ; i < npa ; i ++ )
    {
      px1a [ i ] = pa [ i ] + sa [ i ] * x1 ;
      px2a [ i ] = pa [ i ] + sa [ i ] * x2 ;
      printf ( "PAi: %lf\tPAf: %lf\n", pa [ i ] + sa [ i ] * a , pa [ i ] + sa [ i ] * b ) ;
    }
  for ( i = 0 ; i < npf ; i ++ )
    {
      px1f [ i ] = pf [ i ] ;
      px2f [ i ] = pf [ i ] ;
    }
  f1 = IP ( npa , px1a , npf , px1f , phi , phir ) ;
  f2 = IP ( npa , px2a , npf , px2f , phi , phir ) ;
  //====================================================================================
  // APLICAÇÃO DO MÉTODO DA RAZÃO ÁUREA PARA ALFA
  //====================================================================================
  while ( ( b - a ) > tol )
    {
      if ( f1 > f2 )
	{ 
	  a = x1 ;
	  x1 = x2 ;
	  f1 = f2 ;
	  x2 = a + t * ( b - a ) ;
	  for ( i = 0 ; i < npa ; i ++ )
	    {
	      px2a [ i ] = pa [ i ] + sa [ i ] * x2 ;
	    }
	  f2 = IP ( npa , px2a , npf , px2f , phi , phir ) ;
	}
      else
	{
	  b = x2 ;
	  x2 = x1 ;
	  f2 = f1 ;
	  x1 =  a + ( 1.0 - t ) * ( b - a ) ;
	  for ( i = 0 ; i < npa ; i ++ )
	    {
	      px1a [ i ] = pa [ i ] + sa [ i ] * x1 ;
	    }
	  f1 = IP ( npa , px1a , npf , px1f , phi , phir ) ;
	}
    }
  //====================================================================================
  // DETERMINAÇÃO DO ALFA FINAL
  //====================================================================================
  alfa = ( x1 + x2 ) / 2.0 ;
  //====================================================================================
  // DETERMINAÇÃO DO VETOR DE PARÂMETROS - MINIMIZA O FUNCIONAL NA DIREÇÃO DE BUSCA
  //====================================================================================
  for ( i = 0 ; i < npa ; i ++ )
    {
      pa [ i ] = pa [ i ] + sa [ i ] * alfa ;
    }
  //====================================================================================
  // INICIALIZAÇÃO DE VARIÁVEIS
  //====================================================================================
  normaf = pf [ 0 ] ;
  for ( i = 0 ; i < npf ; i ++ )
    {
      if ( pf [ i ] > normaf ) normaf = pf [ i ] ;
    }
  a = 0.0 ;
  b = 0.15 ;
  t = ( sqrt ( 5.0 ) - 1.0 ) / 2.0 ;
  x1 = a + ( 1.0 - t ) * ( b - a ) ;
  x2 = a + t * ( b - a ) ;
  for ( i = 0 ; i < npf ; i ++ )
    {
      px1f [ i ] = pf [ i ] + sf [ i ] * x1 ;
      px2f [ i ] = pf [ i ] + sf [ i ] * x2 ;
      printf ( "PFi: %lf\tPFf: %lf\n", pf [ i ] + sf [ i ] * a , pf [ i ] + sf [ i ] * b ) ;
    }
  for ( i = 0 ; i < npa ; i ++ )
    {
      px1a [ i ] = pa [ i ] ;
      px2a [ i ] = pa [ i ] ;
    }
  f1 = IP ( npa , px1a , npf , px1f , phi , phir ) ;
  f2 = IP ( npa , px2a , npf , px2f , phi , phir ) ;
  //====================================================================================
  // APLICAÇÃO DO MÉTODO DA RAZÃO ÁUREA PARA O TERMO FONTE
  //====================================================================================
  while ( ( b - a ) > tol )
    {
      if ( f1 > f2 )
	{ 
	  a = x1 ;
	  x1 = x2 ;
	  f1 = f2 ;
	  x2 = a + t * ( b - a ) ;
	  for ( i = 0 ; i < npf ; i ++ )
	    {
	      px2f [ i ] = pf [ i ] + sf [ i ] * x2 ;
	    }
	  f2 = IP ( npa , px2a , npf , px2f , phi , phir ) ;
	}
      else
	{
	  b = x2 ;
	  x2 = x1 ;
	  f2 = f1 ;
	  x1 =  a + ( 1.0 - t ) * ( b - a ) ;
	  for ( i = 0 ; i < npf ; i ++ )
	    {
	      px1f [ i ] = pf [ i ] + sf [ i ] * x1 ;
	    }
	  f1 = IP ( npa , px1a , npf , px1f , phi , phir ) ;
	}
    }
  //====================================================================================
  // DETERMINAÇÃO DO ALFA FINAL
  //====================================================================================
  alfa = ( x1 + x2 ) / 2.0 ;
  //====================================================================================
  // DETERMINAÇÃO DO VETOR DE PARÂMETROS - MINIMIZA O FUNCIONAL NA DIREÇÃO DE BUSCA
  //====================================================================================
  for ( i = 0 ; i < npf ; i ++ )
    {
      pf [ i ] = pf [ i ] + sf [ i ] * alfa ;
    }
  //====================================================================================
  // LIBERAÇÃO DAS MEMÓRIAS ALOCADAS
  //====================================================================================
  LVR ( &px1a ) ;
  LVR ( &px2a ) ;
  LVR ( &px1f ) ;
  LVR ( &px2f ) ;
  //====================================================================================
  // RETORNO DO VALOR DE ALFA ENCONTRADO
  //====================================================================================
  return ( alfa ) ;
}
