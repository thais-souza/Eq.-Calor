//======================================================================================
// SUBROTINAS PARA O CÁLCULO DO RESÍDUO DOS PROBLEMAS PRIMAL E DUAL
//======================================================================================
// SUBROTINA PARA O CÁLCULO DO RESÍDUO DO PROBLEMA PRIMAL 
//======================================================================================
void RP ( int npa , int npf , int k , double t , double *pa , double *pf ,
	  double ***phi , double **r )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i , j ;        // VARIÁVEIS PARA OS DESLOCAMENTOS NAS DIREÇÕES X E Y
  double x , y ;     // VARIÁVEIS REAIS PARA OS DESLOCAMENTOS NAS DIREÇÕES X E Y
  //====================================================================================
  // CÁLCULO DO RESÍDUO
  //====================================================================================
  for ( j = 1 ; j < NY ; j ++ )
    {
      y = j * DY ;
      for ( i = 1 ; i < NX ; i ++ )
	{
	  x = i * DX ;
	  r [ i ] [ j ] = ALFA ( npa , pa , ( x + 0.5 * DX ) , y ) *
	    ( phi [ i + 1 ] [ j ] [ k ] - phi [ i ] [ j ] [ k ] ) -
	    ALFA ( npa , pa , ( x - 0.5 * DX ) , y ) *
	    ( phi [ i ] [ j ] [ k ] - phi [ i - 1 ] [ j ] [ k ] ) +
	    ALFA ( npa , pa , x , ( y + 0.5 * DY ) ) *
	    ( phi [ i ] [ j + 1 ] [ k ] - phi [ i ] [ j ] [ k ] ) -
	    ALFA ( npa , pa , x , ( y - 0.5 * DY ) ) *
	    ( phi [ i ] [ j ] [ k ] - phi [ i ] [ j - 1 ] [ k ] ) ;
	  r [ i ] [ j ] *= W ;
	  r [ i ] [ j ] += FXYT ( npf , pf , x , y , t ) ;
	}
    }
}

//======================================================================================
// SUBROTINA PARA O CÁLCULO DO RESÍDUO DO PROBLEMA DUAL
//======================================================================================
void RD( int npa , int k , double *pa , double ***phi , double ***phir , double ***psi,
	  double **r )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i , j ;        // VARIÁVEIS PARA OS DESLOCAMENTOS NAS DIREÇÕES X E Y
  double x , y ;     // VARIÁVEIS REAIS PARA OS DESLOCAMENTOS NAS DIREÇÕES X E Y
  //====================================================================================
  // CÁLCULO DO RESÍDUO
  //====================================================================================
  for ( j = 1 ; j < NY ; j ++ )
    {
      y = j * DY ;
      for ( i = 1 ; i < NX ; i ++ )
	{
	  x = i * DX ;
	  r [ i ] [ j ] = - ALFA ( npa , pa , ( x + 0.5 * DX ) , y ) *
	    ( psi [ i + 1 ] [ j ] [ k ] - psi [ i ] [ j ] [ k ] ) +
	    ALFA ( npa , pa , ( x - 0.5 * DX ) , y ) *
	    ( psi [ i ] [ j ] [ k ] - psi [ i - 1 ] [ j ] [ k ] ) -
	    ALFA ( npa , pa , x , ( y + 0.5 * DY ) ) *
	    ( psi [ i ] [ j + 1 ] [ k ] - psi [ i ] [ j ] [ k ] ) +
	    ALFA ( npa , pa , x , ( y - 0.5 * DY ) ) *
	    ( psi [ i ] [ j ] [ k ] - psi [ i ] [ j - 1 ] [ k ] ) ;
	  r [ i ] [ j ] *= W ;
	  r [ i ] [ j ] -= phir [ i ] [ j ] [ k ] ;
	  r [ i ] [ j ] += phi [ i ] [ j ] [ k ] ;
	}
    }
}
