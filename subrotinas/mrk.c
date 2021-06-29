//======================================================================================
// SUBROTINAS PARA A APLICAÇÃO DO MÉTODO DE RUNGE-KUTTA
//======================================================================================
// SUBROTINA PARA A APLICAÇÃO DO MÉTODO DE RUNGE-KUTTA 2ª ORDEM
//======================================================================================

void MRK20P ( int npa, double *pa , int npf , double *pf , double ***phi )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i , j , k ;    // VARIÁVEIS PARA O DESLOCAMENTO NAS DIREÇÕES X, Y E T
  double **r1 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 1
  double **r2 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 2
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA
  //====================================================================================
  AMR ( &r1 , NX + 1 , NY + 1 ) ;
  AMR ( &r2 , NX + 1 , NY + 1 ) ;
  //====================================================================================
  // CÁLCULO DA SOLUÇÃO INICIAL
  //====================================================================================
  INICIALP ( phi ) ;
  FRONTEIRAP ( phi , 0 ) ;
  //====================================================================================
  // APLICAÇÃO DO MÉTODO
  //====================================================================================
  for ( k = 1 ; k <= NT ; k ++ )
    { 
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO RESÍDUO
      //================================================================================
      RP ( npa , npf , ( k - 1 ) , ( k - 1 ) * DT , pa , pf , phi , r1 ) ;
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ ) 
	    {
	      phi [ i ] [ j ] [ k ] = phi [ i ] [ j ] [ k - 1 ] +
		DT * r1 [ i ] [ j ] ;
	    }
	}
      FRONTEIRAP ( phi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO RESÍDUO
      //================================================================================
      RP ( npa , npf , k , k * DT , pa , pf , phi , r2 ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO ESTÁGIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      phi [ i ] [ j ] [ k ]  = phi [ i ] [ j ] [ k - 1 ] +
		( DT * ( r1 [ i ] [ j ] + r2 [ i ] [ j ] ) / 2.0 ) ;
	    }
	}
      FRONTEIRAP ( phi , k ) ; 
    }
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA
  //====================================================================================
  LMR ( &r1 , NX + 1 ) ;
  LMR ( &r2 , NX + 1 ) ;
}


//======================================================================================
// SUBROTINA PARA A APLICAÇÃO DO MÉTODO DE RUNGE-KUTTA 3ª ORDEM
//======================================================================================

void MRK30P ( int npa , double *pa , int npf , double *pf , double ***phi )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i , j , k ;    // VARIÁVEIS PARA OS DESLOCAMENTOS NAS DIREÇÕES X, Y E T
  double **r1 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 1
  double **r2 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 2
  double **r3 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 3
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA
  //====================================================================================
  AMR ( &r1 , NX + 1 , NY + 1 ) ;
  AMR ( &r2 , NX + 1 , NY + 1 ) ;
  AMR ( &r3 , NX + 1 , NY + 1 ) ;
  //====================================================================================
  // CÁLCULO DA SOLUÇÃO INICIAL
  //====================================================================================
  INICIALP ( phi ) ;
  FRONTEIRAP ( phi , 0 ) ;
  //====================================================================================
  // APLICAÇÃO DO MÉTODO
  //====================================================================================
  for ( k = 1 ; k <= NT ; k ++ )
    {
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO RESÍDUO
      //================================================================================
      RP ( npa , npf , ( k - 1 ) , ( k - 1 ) * DT , pa , pf , phi , r1 ) ;
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ ) 
	    {
	      phi [ i ] [ j ] [ k ] = phi [ i ] [ j ] [ k - 1 ] + 
		DT / 2.0 * r1 [ i ] [ j ] ;
	    }
	}
      FRONTEIRAP ( phi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO RESÍDUO
      //================================================================================
      RP ( npa , npf , k , ( k - ( 1 / 2 ) ) * DT , pa , pf , phi , r2 ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      phi [ i ] [ j ] [ k ]  = phi [ i ] [ j ] [ k - 1 ] +
		( DT * 2.0 * r2 [ i ] [ j ] - DT * r1 [ i ] [ j ] ) ;
	    }
	}
      FRONTEIRAP ( phi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO TERCEIRO RESÍDUO
      //================================================================================
      RP ( npa , npf , k , k * DT , pa , pf , phi , r3 ) ;
      //================================================================================
      // DETERMINAÇÃO DO TERCEIRO ESTÁGIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      phi [ i ] [ j ] [ k ]  = phi [ i ] [ j ] [ k - 1 ] + 
		( DT / 6.0 * ( r1 [ i ] [ j ] + 4 * r2 [ i ] [ j ] +
			       r3 [ i ] [ j ] ) ) ;
	    }
	}
      FRONTEIRAP ( phi , k ) ;
    }
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA
  //====================================================================================
  LMR ( &r1 , NX + 1 ) ;
  LMR ( &r2 , NX + 1 ) ;
  LMR ( &r3 , NX + 1 ) ;
}

//======================================================================================
// SUBROTINA PARA A APLICAÇÃO DO MÉTODO DE RUNGE-KUTTA 4ª ORDEM
//======================================================================================

void MRK40P ( int npa , double *pa , int npf , double *pf , double ***phi )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i , j , k ;    // VARIÁVEIS PARA OS DESLOCAMENTOS NAS DIREÇÕES X, Y E T
  double **r1 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 1
  double **r2 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 2
  double **r3 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 3
  double **r4 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 4
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA
  //====================================================================================
  AMR ( &r1 , NX + 1 , NY + 1 ) ;
  AMR ( &r2 , NX + 1 , NY + 1 ) ;
  AMR ( &r3 , NX + 1 , NY + 1 ) ;
  AMR ( &r4 , NX + 1 , NY + 1 ) ;
  //====================================================================================
  // CÁLCULO DA SOLUÇÃO INICIAL
  //====================================================================================
  INICIALP ( phi ) ;
  FRONTEIRAP ( phi , 0 ) ;
  //====================================================================================
  // APLICAÇÃO DO MÉTODO
  //====================================================================================
  for ( k = 1 ; k <= NT ; k ++ )
    {
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO RESÍDUO
      //================================================================================
      RP ( npa , npf , ( k - 1 ) , ( k - 1 ) * DT , pa , pf , phi , r1 ) ;
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ ) 
	    {
	      phi [ i ] [ j ] [ k ] = phi [ i ] [ j ] [ k - 1 ] + 
		DT / 2.0 * r1 [ i ] [ j ] ;
	    }
	}
      FRONTEIRAP ( phi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO RESÍDUO
      //================================================================================
      RP ( npa , npf , k , ( k - ( 1 / 2 ) ) * DT , pa , pf , phi , r2 ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      phi [ i ] [ j ] [ k ]  = phi [ i ] [ j ] [ k - 1 ] +
		( DT / 2.0 * r2 [ i ] [ j ] ) ;
	    }
	}
      FRONTEIRAP ( phi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO TERCEIRO RESÍDUO
      //================================================================================
      RP ( npa , npf , k , ( k - ( 1 / 2 ) ) * DT , pa , pf , phi , r3 ) ;
      //================================================================================
      // DETERMINAÇÃO DO TERCEIRO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      phi [ i ] [ j ] [ k ]  = phi [ i ] [ j ] [ k - 1 ] +
		DT * r3 [ i ] [ j ] ;
	    }
	}
      FRONTEIRAP ( phi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO QUARTO RESÍDUO
      //================================================================================
      RP ( npa , npf , k , k * DT , pa , pf , phi , r4 ) ;
      //================================================================================
      // DETERMINAÇÃO DO TERCEIRO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      phi [ i ] [ j ] [ k ]  = phi [ i ] [ j ] [ k - 1 ] +
		( DT / 6.0 * ( r1 [ i ] [ j ] + 2 * r2 [ i ] [ j ] + 2 * r3 [ i ] [ j ] +
			       r4 [ i ] [ j ] ) );
	    }
	}
      FRONTEIRAP ( phi , k ) ;
    }
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA
  //====================================================================================
  LMR ( &r1 , NX + 1 ) ;
  LMR ( &r2 , NX + 1 ) ;
  LMR ( &r3 , NX + 1 ) ;
  LMR ( &r4 , NX + 1 ) ;
}

//======================================================================================
// SUBROTINAS PARA A APLICAÇÃO DO MÉTODO DE RUNGE-KUTTA
//======================================================================================
// SUBROTINA PARA A APLICAÇÃO DO MÉTODO DE RUNGE-KUTTA 2ª ORDEM
//======================================================================================

void MRK20D ( int npa , double *pa , double ***phi , double ***phir, double ***psi)
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i , j , k ;    // VARIÁVEIS PARA O DESLOCAMENTO NAS DIREÇÕES X, Y E T
  double **r1 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 1
  double **r2 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 2
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA
  //====================================================================================
  AMR ( &r1 , NX + 1 , NY + 1 ) ;
  AMR ( &r2 , NX + 1 , NY + 1 ) ;
  //====================================================================================
  // CÁLCULO DA SOLUÇÃO INICIAL
  //====================================================================================
  INICIALD ( psi ) ;
  FRONTEIRAD ( psi , NT ) ;
  //====================================================================================
  // APLICAÇÃO DO MÉTODO
  //====================================================================================
  for ( k = NT - 1 ; k >= 0 ; k -- )
    { 
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO RESÍDUO
      //================================================================================
      RD ( npa , ( k + 1 ) , pa , phi , phir , psi , r1 ) ;
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ ) 
	    {
	      psi [ i ] [ j ] [ k ] = psi [ i ] [ j ] [ k + 1 ] -
		DT * r1 [ i ] [ j ] ;
	    }
	}
      FRONTEIRAD ( psi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO RESÍDUO
      //================================================================================
      RD ( npa , k , pa , phi , phir , psi , r2 ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO ESTÁGIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      psi [ i ] [ j ] [ k ]  = psi [ i ] [ j ] [ k + 1 ] -
		( DT * ( r1 [ i ] [ j ] + r2 [ i ] [ j ] ) / 2.0 ) ;
	    }
	}
      FRONTEIRAD ( psi , k ) ; 
    }
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA
  //====================================================================================
  LMR ( &r1 , NX + 1 ) ;
  LMR ( &r2 , NX + 1 ) ;
}


//======================================================================================
// SUBROTINA PARA A APLICAÇÃO DO MÉTODO DE RUNGE-KUTTA 3ª ORDEM
//======================================================================================

void MRK30D ( int npa , double *pa , double ***phi , double ***phir , double ***psi )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i , j , k ;    // VARIÁVEIS PARA OS DESLOCAMENTOS NAS DIREÇÕES X, Y E T
  double **r1 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 1
  double **r2 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 2
  double **r3 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 3
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA
  //====================================================================================
  AMR ( &r1 , NX + 1 , NY + 1 ) ;
  AMR ( &r2 , NX + 1 , NY + 1 ) ;
  AMR ( &r3 , NX + 1 , NY + 1 ) ;
  //====================================================================================
  // CÁLCULO DA SOLUÇÃO INICIAL
  //====================================================================================
  INICIALD ( psi ) ;
  FRONTEIRAD ( psi , NT ) ;
  //====================================================================================
  // APLICAÇÃO DO MÉTODO
  //====================================================================================
  for ( k = NT - 1 ; k >= 0 ; k -- )
    {
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO RESÍDUO
      //================================================================================
      RD ( npa , ( k + 1 ) , pa , phi , phir , psi , r1 ) ;
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ ) 
	    {
	      psi [ i ] [ j ] [ k ] = psi [ i ] [ j ] [ k + 1 ] -
		DT / 2.0 * r1 [ i ] [ j ] ;
	    }
	}
      FRONTEIRAD ( psi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO RESÍDUO
      //================================================================================
      RD ( npa , k , pa , phi , phir , psi , r2 ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      psi [ i ] [ j ] [ k ]  = psi [ i ] [ j ] [ k + 1 ] -
		DT * 2.0 * r2 [ i ] [ j ] + DT * r1 [ i ] [ j ] ;
	    }
	}
      FRONTEIRAD ( psi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO TERCEIRO RESÍDUO
      //================================================================================
      RD ( npa , k , pa , phi , phir , psi , r3 ) ;
      //================================================================================
      // DETERMINAÇÃO DO TERCEIRO ESTÁGIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      psi [ i ] [ j ] [ k ]  = psi [ i ] [ j ] [ k + 1 ] -
		( DT / 6.0 * ( r1 [ i ] [ j ] + 4 * r2 [ i ] [ j ] +
			       r3 [ i ] [ j ] ) ) ;
	    }
	}
      FRONTEIRAD ( psi , k ) ;
    }
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA
  //====================================================================================
  LMR ( &r1 , NX + 1 ) ;
  LMR ( &r2 , NX + 1 ) ;
  LMR ( &r3 , NX + 1 ) ;
}

//======================================================================================
// SUBROTINA PARA A APLICAÇÃO DO MÉTODO DE RUNGE-KUTTA 4ª ORDEM
//======================================================================================

void MRK40D ( int npa , double *pa , double ***phi , double ***phir , double ***psi )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i , j , k ;    // VARIÁVEIS PARA OS DESLOCAMENTOS NAS DIREÇÕES X, Y E T
  double **r1 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 1
  double **r2 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 2
  double **r3 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 3
  double **r4 ;      // PONTEIRO PARA O ARMAZENAMENTO DO RESÍDUO 4
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA
  //====================================================================================
  AMR ( &r1 , NX + 1 , NY + 1 ) ;
  AMR ( &r2 , NX + 1 , NY + 1 ) ;
  AMR ( &r3 , NX + 1 , NY + 1 ) ;
  AMR ( &r4 , NX + 1 , NY + 1 ) ;
  //====================================================================================
  // CÁLCULO DA SOLUÇÃO INICIAL
  //====================================================================================
  INICIALD ( psi ) ;
  FRONTEIRAD ( psi , NT ) ;
  //====================================================================================
  // APLICAÇÃO DO MÉTODO
  //====================================================================================
  for ( k = NT - 1 ; k >= 0 ; k -- )
    {
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO RESÍDUO
      //================================================================================
      RD ( npa , ( k + 1 ) , pa , phi , phir , psi , r1 ) ;
      //================================================================================
      // DETERMINAÇÃO DO PRIMEIRO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ ) 
	    {
	      psi [ i ] [ j ] [ k ] = psi [ i ] [ j ] [ k + 1 ] - 
		DT / 2.0 * r1 [ i ] [ j ] ;
	    }
	}
      FRONTEIRAD ( psi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO RESÍDUO
      //================================================================================
      RD ( npa , k , pa , phi , phir , psi , r2 ) ;
      //================================================================================
      // DETERMINAÇÃO DO SEGUNDO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      psi [ i ] [ j ] [ k ]  = psi [ i ] [ j ] [ k + 1 ] -
		( DT / 2.0 * r2 [ i ] [ j ] ) ;
	    }
	}
      FRONTEIRAD ( psi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO TERCEIRO RESÍDUO
      //================================================================================
      RD ( npa , k , pa , phi , phir , psi , r3 ) ;
      //================================================================================
      // DETERMINAÇÃO DO TERCEIRO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      psi [ i ] [ j ] [ k ]  = psi [ i ] [ j ] [ k + 1 ] -
		DT * r3 [ i ] [ j ] ;
	    }
	}
      FRONTEIRAD ( psi , k ) ;
      //================================================================================
      // DETERMINAÇÃO DO QUARTO RESÍDUO
      //================================================================================
      RD ( npa , k , pa , phi , phir , psi , r4 ) ;
      //================================================================================
      // DETERMINAÇÃO DO TERCEIRO ESTÁGIO INTERMEDIÁRIO
      //================================================================================
      for ( j = 1 ; j < NY ; j ++ ) 
	{
	  for ( i = 1 ; i < NX ; i ++ )
	    {
	      psi [ i ] [ j ] [ k ]  = psi [ i ] [ j ] [ k + 1 ] -
		( DT / 6.0 * ( r1 [ i ] [ j ] + 2 * r2 [ i ] [ j ] + 2 * r3 [ i ] [ j ] +
			       r4 [ i ] [ j ] ) );
	    }
	}
      FRONTEIRAD ( psi , k ) ;
    }
  //====================================================================================
  // LIBERAÇÃO DE MEMÓRIA
  //====================================================================================
  LMR ( &r1 , NX + 1 ) ;
  LMR ( &r2 , NX + 1 ) ;
  LMR ( &r3 , NX + 1 ) ;
  LMR ( &r4 , NX + 1 ) ;
}