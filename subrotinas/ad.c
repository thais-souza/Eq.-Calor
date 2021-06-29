//======================================================================================
// SUBROTINAS PARA AS ALOCAÇÕES DINÂMICAS
//======================================================================================
// SUBROTINA PARA A ALOCAÇÃO DE VETORES REAIS
//======================================================================================

void AVR ( double **vetor , int colunas )
{
  ( *vetor ) = calloc ( colunas , sizeof ( double ) ) ;
  if ( ( * vetor ) == NULL )
    {
      printf ( "\nAlocacao de Vetor: Memoria Insuficiente\n" );
    }
}

//======================================================================================
// SUBROTINA PARA A LIBERAÇÃO DE VETORES REAIS
//======================================================================================

void LVR ( double **vetor )
{
  if ( ( *vetor ) == NULL )
    {
      return ;
    }
  free ( *vetor ) ;
}

//======================================================================================
// SUBROTINA PARA A ALOCAÇÃO DE MATRIZES 2D REAIS
//======================================================================================

void AMR ( double ***matriz , int linhas , int colunas )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;            // VARIÁVEL PARA O DESLOCAMENTO ENTRE AS LINHAS
  //====================================================================================
  // ALOCAÇÃO DAS LINHAS
  //====================================================================================
  ( *matriz ) = calloc ( linhas , sizeof ( double * ) ) ;
  if ( ( *matriz ) == NULL )
    {
      printf ( "\nAlocacao de Matriz 2D: Memoria Insuficiente - Linhas\n" ) ;
      return ;
    }
  //====================================================================================
  // ALOCAÇÃO DAS COLUNAS
  //====================================================================================
  for ( i = 0 ; i < linhas ; i ++ )
    {
      ( *matriz ) [ i ] = calloc ( colunas , sizeof ( double ) ) ;
      if ( ( *matriz ) [ i ] == NULL )
	{
	  printf ( "\nAlocacao de Matriz 2D: Memoria Insuficiente - Colunas\n" ) ;
	  return ;
	}
    }
}

//======================================================================================
// SUBROTINA PARA A LIBERAÇÃO DE MATRIZES 2D REAIS
//======================================================================================

void LMR ( double ***matriz , int linhas )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;            // VARIÁVEL PARA O DESLOCAMENTO ENTRE AS LINHAS
  //====================================================================================
  // VERIFICAÇÃO DA PRESENÇA DE ESPAÇOS ALOCADOS
  //====================================================================================
  if ( ( *matriz ) == NULL )
    {
      return ;
    }
  //====================================================================================
  // LIBERAÇÃO DAS COLUNAS
  //====================================================================================
  for ( i = 0 ; i < linhas ; i ++ )
    {
      free ( ( *matriz ) [ i ] ) ;
    }
  //====================================================================================
  // LIBERAÇÃO DAS LINHAS
  //====================================================================================
  free ( *matriz ) ;
}

//======================================================================================
// SUBROTINA PARA A ALOCAÇÃO DE MATRIZES 3D REAIS
//======================================================================================

void AMR3D ( double ****matriz , int linhas , int colunas , int profundidades )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;            // VARIÁVEL PARA O DESLOCAMENTO ENTRE AS LINHAS
  int j ;            // VARIÁVEL PARA O DESLOCAMENTO ENTRE AS COLUNAS
  //====================================================================================
  // ALOCAÇÃO DAS LINHAS
  //====================================================================================
  ( *matriz ) = calloc ( linhas , sizeof ( double * ) ) ;
  if ( *matriz == NULL )
    {
      printf ( "\nAlocacao de Matriz 3D: Memoria Insuficiente - Linhas\n" ) ;
      return ;
    }
  //====================================================================================
  // ALOCAÇÃO DAS COLUNAS
  //====================================================================================
  for ( i = 0 ; i < linhas ; i ++ )
    {
      ( *matriz ) [ i ] = calloc ( colunas , sizeof ( double ) ) ;
      if ( ( *matriz ) [ i ] == NULL )
	{
	  printf ("\nAlocacao de Matriz 3D: Memoria Insuficiente - Colunas\n" ) ;
	  return ;
	}
    }
  //====================================================================================
  // ALOCAÇÃO DAS PROFUNDIDADES
  //====================================================================================
  for ( i = 0 ; i < linhas ; i ++ )
    {
      for ( j = 0 ; j < colunas ; j ++ )
	{
	  ( *matriz ) [ i ] [ j ] = calloc ( profundidades , sizeof ( double ) ) ;
	  if ( ( *matriz ) [ i ] [ j ] == NULL )
	    {
	      printf ( "\nAlocacao de Matriz 3D:" ) ;
	      printf ( " Memoria Insuficiente - Produndidades\n" ) ;
	      return ;
	    }
	}
    }
}

//======================================================================================
// SUBROTINA PARA A LIBERAÇÃO DE MATRIZES 3D REAIS
//======================================================================================

void LMR3D ( double ****matriz , int linhas , int colunas )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;            // VARIÁVEL PARA O DESLOCAMENTO ENTRE AS LINHAS
  int j ;            // VARIÁVEL PARA O DESLOCAMENTO ENTRE AS COLUNAS
  //====================================================================================
  // VERIFICAÇÃO DA PRESENÇA DE ESPAÇOS ALOCADOS
  //====================================================================================
  if ( ( *matriz ) == NULL )
    {
      return ;
    }
  //====================================================================================
  // LIBERAÇÃO DAS PROFUNDIDADES E DAS COLUNAS
  //====================================================================================
  for ( i = 0 ; i < linhas ; i ++ )
    {
      for ( j = 0 ; j < colunas ; j ++ )
	{
	  free ( ( *matriz ) [ i ] [ j ] ) ;
	}
    }
  //====================================================================================
  // LIBERAÇÃO DAS LINHAS
  //====================================================================================
  free ( *matriz ) ;
}