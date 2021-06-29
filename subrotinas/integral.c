//======================================================================================
// SUBROTINA PARA O CÁLCULO DA INTEGRAL TRIPLA VIA REGRA 1/3 DE SIMPSON
//======================================================================================
// SUBROTINA PARA O CÁLCULO DO FUNCIONAL
//======================================================================================
double IP ( int npa , double *pa , int npf , double *pf , double ***phi ,
	    double ***phir )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int nx , ny , nt ;         // NÚMERO DE DIVISÕES NOS INTERVALOS DE X, Y E T
  int i , j , k ;            // VARIÁVEIS AUXILIARES PARA OS DESLOCAMENTOS EM X, Y E T
  double som1 , som2 ;       // VARIÁVEIS AUXILIARES PARA O CÁLCULO DOS SOMATÓRIOS
  double integral ;          // RESULTADO DA INTEGRAL TRIPLA
  double ***integrando ;     // PONTEIRO PARA O ARMAZENAMENTO DO INTEGRANDO
  double *fyt ;              // PONTEIRO PARA O ARMAZENAMENTO DA INTEGRAL EM X
  double *ft ;               // PONTEIRO PARA O ARMAZENAMENTO DA INTEGRAL EM X E Y
  //====================================================================================
  // INICIALIZAÇÃO DE VARIÁVEIS
  //====================================================================================
  nx = NX / 2 ;
  ny = NY / 2 ;
  nt = NT / 2 ;
  som1 = 0.0 ;
  som2 = 0.0 ;
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA
  //====================================================================================
  AMR3D ( &integrando , NX + 1 , NY + 1 , NT + 1 ) ;
  AVR ( &fyt , NY + 1 ) ;
  AVR ( &ft , NT + 1 ) ;
  //====================================================================================
  // RESOLUÇÃO DA EQUAÇÃO DIFERENCIAL EM ANÁLISE
  //====================================================================================
  ED2DP ( npa , pa , npf , pf , phi ) ;
  //====================================================================================
  // CÁLCULO DO INTEGRANDO
  //====================================================================================
  for ( k = 0 ; k <= NT ; k ++ )
    {
      for ( i = 0 ; i <= NX ; i ++ )
	{
	  for ( j = 0 ; j <= NY ; j ++ )
	    {
	      integrando [ i ] [ j ] [ k ] = pow ( ( phir [ i ] [ j ] [ k ] -
				   phi [ i ] [ j ] [ k ] ) , 2 ) ;
	    }
	}
    }
  //====================================================================================
  // CÁLCULO DAS INTEGRAIS EM RELAÇÃO A X E Y
  //====================================================================================
  for ( k = 0 ; k <= NT ; k ++ )
    {
      for ( j = 0 ; j <= NY ; j ++ )
	{
	  som1 = 0.0 ;
	  som2 = 0.0 ;
	  for ( i = 0 ; i <= ( nx - 1 ) ; i ++ )
	    {
	      som1 += integrando [ 2 * i + 1 ] [ j ] [ k ] ;
	    }
	  for ( i = 1 ; i <= ( nx - 1 ) ; i ++ )
	    {
	      som2 += integrando [ 2 * i ] [ j ] [ k ] ;
	    }
	  fyt [ j ] = ( DX / 3.0 ) * ( integrando [ 0 ] [ j ] [ k ] +
				       integrando [ NX ] [ j ] [ k ] +
				       4 * som1 + 2 * som2 ) ;
	}
      som1 = 0.0 ;
      som2 = 0.0 ;
      for ( j = 0 ; j <= ( ny - 1 ) ; j ++ ) 
	{
	  som1 += fyt [ 2 * j + 1 ] ;
	}
      for ( j = 1 ; j <= ( ny - 1 ) ; j ++ )
	{
	  som2 += fyt [ 2 * j ] ;
	}
      ft [ k ] = ( DY / 3.0 ) * ( fyt [ 0 ] + fyt [ NY ] + 4 * som1 + 2 * som2 ) ;
    }
  //====================================================================================
  // CÁLCULO DA INTEGRAL EM RELAÇÃO A T
  //====================================================================================
  som1 = 0.0 ;
  som2 = 0.0 ;
  for ( k = 0 ; k <= ( nt - 1 ) ; k ++ )
    {
      som1 += ft [ 2 * k + 1 ] ;
    }
  for ( k = 1 ; k <= ( nt - 1 ) ; k ++ )
    {
      som2 += ft [ 2 * k ] ;
    }
  integral = ( DT / 3.0 ) * ( ft [ 0 ] + ft [ NT ] + 4 * som1 + 2 * som2 ) ;
  integral /= 2.0 ;
  printf ( "%e\n", ( integral ) ) ;
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA
  //====================================================================================
  LMR3D ( &integrando , NX + 1 , NY + 1 ) ;
  LVR ( &fyt ) ;
  LVR ( &ft ) ;
  //====================================================================================
  // RETORNO DO RESULTADO
  //====================================================================================
  return (  ( integral ) ) ;
}

//======================================================================================
// SUBROTINA PARA O CÁLCULO DA INTEGRAL RELATIVA AO PROBLEMA DUAL
//======================================================================================
double ID ( int npa , double *ca , int npf , double *cf , double ***phi ,
	    double ***psi , int id )
{
  //====================================================================================
  // SE ID = 0, INTEGRAL RELATIVA AO GRADIENTE EM ALFA. SE ID = 1, INTEGRAL RELATIVA
  // AO GRADIENTE NO TERMO FONTE
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int nx , ny , nt ;         // NÚMERO DE DIVISÕES NOS INTERVALOS DE X, Y E T
  int i , j , k ;            // VARIÁVEIS AUXILIARES PARA OS DESLOCAMENTOS EM X, Y E T
  double x , y , t ;         // VARIÁVEIS REAIS PARA OS DESLOCAMENTOS EM X, Y E T
  double som1 , som2 ;       // VARIÁVEIS AUXILIARES PARA O CÁLCULO DOS SOMATÓRIOS
  double integral ;          // RESULTADO DA INTEGRAL TRIPLA
  double int1 , int2 ;       // VARIÁVEIS AUXILIARES PARA O CÁLCULO DO INTEGRANDO
  double ***integrando ;     // PONTEIRO PARA O ARMAZENAMENTO DO INTEGRANDO
  double *fyt ;              // PONTEIRO PARA O ARMAZENAMENTO DA INTEGRAL EM X
  double *ft ;               // PONTEIRO PARA O ARMAZENAMENTO DA INTEGRAL EM X E Y
  //====================================================================================
  // INICIALIZAÇÃO DE VARIÁVEIS
  //====================================================================================
  nx = NX / 2 ;
  ny = NY / 2 ;
  nt = NT / 2 ;
  som1 = 0.0 ;
  som2 = 0.0 ;
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA
  //====================================================================================
  AMR3D ( &integrando , NX + 1 , NY + 1 , NT + 1 ) ;
  AVR ( &fyt , NY + 1 ) ;
  AVR ( &ft , NT + 1 ) ;
  //====================================================================================
  // CÁLCULO DO INTEGRANDO
  //====================================================================================
  if ( id == 0 )
    {
      for ( k = 0 ; k <= NT ; k ++ )
	{
	  for ( i = 0 ; i <= NX - 2 ; i ++ )
	    {
	      x = i * DX ;
	      for ( j = 0 ; j <= NY - 2 ; j ++ )
		{
		  y = j * DY ;
		  int1 = ALFA ( npa , ca , ( x + DX ) , y ) *
		    ( psi [ i + 1 ] [ j ] [ k ] - psi [ i ] [ j ] [ k ] ) +
		    ALFA ( npa , ca , x , y ) *
		    ( psi [ i + 2 ] [ j ] [ k ] + 2 * psi [ i ] [ j ] [ k ]
		      - 3 * psi [ i + 1 ] [ j ] [ k ] ) ;
		  int1 *= W ;
		  int2 = ALFA ( npa , ca , x , ( y + DY ) ) *
		    ( psi [ i ] [ j + 1 ] [ k ] - psi [ i ] [ j ] [ k ] ) +
		    ALFA ( npa , ca , x , y ) *
		    ( psi [ i ] [ j + 2 ] [ k ] + 2 * psi [ i ] [ j ] [ k ]
		      - 3 * psi [ i ] [ j + 1 ] [ k ] ) ;
		  int2 *= W ;
		  integrando [ i ] [ j ] [ k ] = ( int1 + int2 ) *
		    phi [ i ] [ j ] [ k ]  ;				  
		}
	      for ( j = NY - 1 ; j <= NY ; j ++ )
		{
		  y = j * DY ;
		  int1 = ALFA ( npa , ca , ( x + DX ) , y ) *
		    ( psi [ i + 1 ] [ j ] [ k ] - psi [ i ] [ j ] [ k ] ) +
		    ALFA ( npa , ca , x , y ) *
		    ( psi [ i + 2 ] [ j ] [ k ] + 2 * psi [ i ] [ j ] [ k ]
		      - 3 * psi [ i + 1 ] [ j ] [ k ] ) ;
		  int1 *= W ;
		  int2 = ALFA ( npa , ca , x , ( y - DY ) ) *
		    ( psi [ i ] [ j - 1 ] [ k ] - psi [ i ] [ j ] [ k ] ) +
		    ALFA ( npa , ca , x , y ) *
		    ( psi [ i ] [ j - 2 ] [ k ] - 3 * psi [ i ] [ j - 1 ] [ k ] +
		      2 * psi [ i ] [ j ] [ k ] ) ;
		  int2 *= W ;
		  integrando [ i ] [ j ] [ k ] = ( int1 + int2 ) *
		    phi [ i ] [ j ] [ k ] ;
		}
	    }
	  for ( i = NX - 1 ; i <= NX ; i ++ )
	    {
	      x = i * DX ;
	      for ( j = 0 ; j <= NY - 2 ; j ++ )
		{
		  y = j * DY ;
		  int1 = ALFA ( npa , ca , ( x - DX ) , y ) *
		    ( psi [ i - 1 ] [ j ] [ k ] - psi [ i ] [ j ] [ k ] ) +
		    ALFA ( npa , ca , x , y ) *
		    ( psi [ i - 2 ] [ j ] [ k ] - 3 * psi [ i - 1 ] [ j ] [ k ] +
		      2 * psi [ i ] [ j ] [ k ] ) ;
		  int1 *= W ;
		  int2 = ALFA ( npa , ca , x , ( y + DY ) ) *
		    ( psi [ i ] [ j + 1 ] [ k ] - psi [ i ] [ j ] [ k ] ) +
		    ALFA ( npa , ca , x , y ) *
		    ( psi [ i ] [ j + 2 ] [ k ] + 2 * psi [ i ] [ j ] [ k ]
		      - 3 * psi [ i ] [ j + 1 ] [ k ] ) ;
		  int2 *= W ;
		  integrando [ i ] [ j ] [ k ] = ( int1 + int2 ) *
		    phi [ i ] [ j ] [ k ]  ;
		}
	      for ( j = NY - 1 ; j <= NY ; j ++ )
		{
		  y = j * DY ;
		  int1 = ALFA ( npa , ca , ( x - DX ) , y ) *
		    ( psi [ i - 1 ] [ j ] [ k ] - psi [ i ] [ j ] [ k ] ) +
		    ALFA ( npa , ca , x , y ) *
		    ( psi [ i - 2 ] [ j ] [ k ] - 3 * psi [ i - 1 ] [ j ] [ k ] +
		      2 * psi [ i ] [ j ] [ k ] ) ;
		  int1 *= W ;
		  int2 = ALFA ( npa , ca , x , ( y - DY ) ) *
		    ( psi [ i ] [ j - 1 ] [ k ] - psi [ i ] [ j ] [ k ] ) +
		    ALFA ( npa , ca , x , y ) *
		    ( psi [ i ] [ j - 2 ] [ k ] - 3 * psi [ i ] [ j - 1 ] [ k ] +
		      2 * psi [ i ] [ j ] [ k ] ) ;
		  int2 *= W ;
		  integrando [ i ] [ j ] [ k ] = ( int1 + int2 ) *
		    phi [ i ] [ j ] [ k ] ;
		}
	    }
	}
    }
  else
    {
      for ( k = 0 ; k <= NT ; k ++ )
	{
	  t = k * DT ;
	  for ( i = 0 ; i <= NX ; i ++ )
	    {
	      x = i * DX ;
	      for ( j = 0 ; j <= NY ; j ++ )
		{
		  y = j * DY ;
		  integrando [ i ] [ j ] [ k ] = psi [ i ] [ j ] [ k ] *
		    FXYT ( npf , cf , x , y , t ) ;
		}
	    }
	}
    }
  //====================================================================================
  // CÁLCULO DAS INTEGRAIS EM RELAÇÃO A X E Y
  //====================================================================================
  for ( k = 0 ; k <= NT ; k ++ )
    {
      for ( j = 0 ; j <= NY ; j ++ )
	{
	  som1 = 0.0 ;
	  som2 = 0.0 ;
	  for ( i = 0 ; i <= ( nx - 1 ) ; i ++ )
	    {
	      som1 += integrando [ 2 * i + 1 ] [ j ] [ k ] ;
	    }
	  for ( i = 1 ; i <= ( nx - 1 ) ; i ++ )
	    {
	      som2 += integrando [ 2 * i ] [ j ] [ k ] ;
	    }
	  fyt [ j ] = ( DX / 3.0 ) * ( integrando [ 0 ] [ j ] [ k ] +
				       integrando [ NX ] [ j ] [ k ] +
				       4 * som1 + 2 * som2 ) ;
	}
      som1 = 0.0 ;
      som2 = 0.0 ;
      for ( j = 0 ; j <= ( ny - 1 ) ; j ++ ) 
	{
	  som1 += fyt [ 2 * j + 1 ] ;
	}
      for ( j = 1 ; j <= ( ny - 1 ) ; j ++ )
	{
	  som2 += fyt [ 2 * j ] ;
	}
      ft [ k ] = ( DY / 3.0 ) * ( fyt [ 0 ] + fyt [ NY ] + 4 * som1 + 2 * som2 ) ;
    }
  //====================================================================================
  // CÁLCULO DA INTEGRAL EM RELAÇÃO A T
  //====================================================================================
  som1 = 0.0 ;
  som2 = 0.0 ;
  for ( k = 0 ; k <= ( nt - 1 ) ; k ++ )
    {
      som1 += ft [ 2 * k + 1 ] ;
    }
  for ( k = 1 ; k <= ( nt - 1 ) ; k ++ )
    {
      som2 += ft [ 2 * k ] ;
    }
  integral = ( DT / 3.0 ) * ( ft [ 0 ] + ft [ NT ] + 4 * som1 + 2 * som2 ) ;
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA
  //====================================================================================
  LMR3D ( &integrando , NX + 1 , NY + 1 ) ;
  LVR ( &fyt ) ;
  LVR ( &ft ) ;
  //====================================================================================
  // RETORNO DO RESULTADO
  //====================================================================================
  return ( integral ) ;
}
