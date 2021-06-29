//======================================================================================
// SUBROTINA PARA A LEITURA EXTERNA
//======================================================================================
void LEITURA ( char *argv [ ] , int *npa , double **pra , double **pa ,
	       int *npf , double **prf , double **pf )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;                    // VARIÁVEL AUXILIAR
  FILE *parametros ;         // PONTEIRO PARA A LEITURA DE PARÂMETROS
  //====================================================================================
  // LEITURA DA ENTRADA
  //====================================================================================
  LX = atof ( argv [ 1 ] ) ;
  LY = atof ( argv [ 2 ] ) ;
  LT = atof ( argv [ 3 ] ) ;
  NX = atoi ( argv [ 4 ] ) ;
  NT = atoi ( argv [ 5 ] ) ;
  ORDEM = atoi ( argv [ 6 ] ) ;
  GA = atoi ( argv [ 7 ] ) ;
  IGA = atoi ( argv [ 8 ] ) ;
  GRAD = atoi ( argv [ 9 ] ) ;
  MG = atoi ( argv [ 10 ] ) ;
  F1 = atoi ( argv [ 11 ] ) ;
  F2 = atoi ( argv [ 12 ] ) ;
  F3 = atoi ( argv [ 13 ] ) ;
  F4 = atoi ( argv [ 14 ] ) ;
  //====================================================================================
  // LEITURA DA QUANTIDADE DE PARÂMEROS RELACIONADOS A ALFA
  //====================================================================================
  parametros = fopen ( "parametros.inf" , "r" ) ;
  if ( parametros == NULL )
    {
      printf ( "\nAtencao: o arquivo de parametros nao existe\n" ) ;
    }
  fscanf ( parametros , "%d" , npa ) ;
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA PARA OS VETORES RELACIONADOS A ALFA
  //====================================================================================
  AVR ( pra , ( *npa ) ) ;
  AVR ( pa , ( *npa ) ) ;
  //====================================================================================
  // LEITURA DE PARÂMETROS RELACIONADOS A ALFA
  //====================================================================================
  for ( i = 0 ; i < ( *npa ) ; i ++ )
    {
      fscanf ( parametros , "%lf", & ( ( *pra ) [ i ] ) ) ;
    }
  //====================================================================================
  // LEITURA DA QUANTIDADE DE PARÂMEROS RELACIONADOS AO TERMO FONTE
  //====================================================================================
  fscanf ( parametros , "%d" , npf ) ;
  //====================================================================================
  // ALOCAÇÃO DE MEMÓRIA PARA OS VETORES RELACIONADOS AO TERMO FONTE
  //====================================================================================
  AVR ( prf , ( *npf ) ) ;
  AVR ( pf , ( *npf ) ) ;
  //====================================================================================
  // LEITURA DE PARÂMETROS RELACIONADOS AO TERMO FONTE
  //====================================================================================
  for ( i = 0 ; i < ( *npf ) ; i ++ )
    {
      fscanf ( parametros , "%lf", & ( ( *prf ) [ i ] ) ) ;
    }
  fclose ( parametros ) ;
}
