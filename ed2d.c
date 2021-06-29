//======================================================================================
// SUBROTINAS PARA A RESOLUÇÃO DOS PROBLEMAS PRIMAL E DUAL
//======================================================================================
// SUBROTINA PARA A RESOLUÇÃO DO PROBLEMA PRIMAL
//======================================================================================

void ED2DP ( int npa , double *pa , int npf , double *pf , double ***phi )
{
  //====================================================================================
  // APLICAÇÃO DO MÉTODO DE RUNGE-KUTTA
  //====================================================================================
  if ( ORDEM == 2 ) MRK20P ( npa , pa , npf , pf , phi ) ;
  else if ( ORDEM == 3 ) MRK30P ( npa , pa , npf , pf , phi ) ;
  else if ( ORDEM == 4 ) MRK40P ( npa , pa , npf , pf , phi ) ;
}

//======================================================================================
// SUBROTINA PARA A RESOLUÇÃO DO PROBLEMA DUAL
//======================================================================================

void ED2DD ( int npa , double *pa , double ***phi , double ***phir, double ***psi )
{
  //====================================================================================
  // APLICAÇÃO DO MÉTODO DE RUNGE-KUTTA
  //====================================================================================
  if ( ORDEM == 2 ) MRK20D ( npa , pa , phi , phir , psi ) ;
  else if ( ORDEM == 3 ) MRK30D ( npa , pa , phi , phir , psi ) ;
  else if ( ORDEM == 4 ) MRK40D ( npa , pa , phi , phir , psi ) ;
}
