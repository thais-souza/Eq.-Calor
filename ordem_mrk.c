//======================================================================================
// SUBROTINA PARA A VERIFICAÇÃO DA ORDEM DOS MRKs
//======================================================================================

void ORDEM_MRK( double ***phi )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;                // VARIÁVEL PARA O DESLOCAMENTO NA DIREÇÃO X
  int j ;                // VARIÁVEL PARA O DESLOCAMENTO NA DIREÇÃO Y
  int k ;                // VARÍAVEL AUXILIAR PARA O DESLOCAMENTO DO PONTEIRO
  double x ;             // VARIÁVEL REAL PARA O DESLOCAMENTO NA DIREÇÃO X
  double y ;             // VARIÁVEL REAL PARA O DESLOCAMENTO NA DIREÇÃO Y
  double a , b ;         // VARIÁVEIS AUXILIARES PARA O ARMAZENAMENTO DOS VALORES DE PHI
  double som ;           // VARIÁVEL AUXILIAR PARA O CÁLCULO DA NORMA
  char nome [ 100 ] ;    // NOME DO ARQUIVO GERADO PARA OS VALORES DE PHI
  char nomen [ 100 ] ;   // NOME DO ARQUIVO ESCRITO PARA OS VALORES DAS NORMAS
  char comando [ 100 ] ; // COMANDO PARA DELETAR ARQUIVOS
  FILE *arquivo ;        // PONTEIRO PARA O ARMAZENAMENTO DE VALORES EM ARQUIVO
  FILE *referencia ;     // PONTEIRO PARA A ABERTURA DO ARQUIVO DE REFERENCIA
  FILE *norma ;          // PONTEIRO PARA O ARMAZENAMENTO DAS NORMAS EM ARQUIVO
  FILE *grafico ;        // PONTEIRO PARA A GERAÇÃO DO GRÁFICO E DA ORDEM DO MÉTODO
  //====================================================================================
  // GERAÇÃO DO ARQUIVO COM VALORES DE PHI DO ÚLTIMO NÍVEL
  //====================================================================================
  sprintf ( nome , "%d.dat" , NT ) ;
  sprintf ( comando , "rm -f %s", nome ) ;
  system ( comando ) ;
  arquivo = fopen ( nome , "a" ) ;
  for ( j = 0 ; j <= NY ; j ++ )
    {
      y = j * DY ;
      for ( i = 0 ; i <= NX ; i++ )
	{
	  x = i * DX ;
	  fprintf ( arquivo , "% lf\t% lf\t" , x , y ) ;
	  fprintf ( arquivo , "% .25lf\n" , phi [ i ] [ j ] [ NT ] ) ;
	}
    }
  fclose ( arquivo ) ;
  //====================================================================================
  // GERAÇÃO DO ARQUIVO DE NORMAS
  //====================================================================================
  printf ( "Geracao norma %d\n" , NT ) ;
  if ( NT != 100000 )
    {
      //================================================================================
      // INICIALIZAÇÃO DE VALORES
      //================================================================================
      som = 0.0 ;
      //================================================================================
      // ABERTURA DOS ARQUIVOS
      //================================================================================
      arquivo = fopen ( nome , "r" ) ;
      referencia = fopen ( "100000.dat" , "r" ) ;
      if ( referencia == NULL )
	{
	  printf ( "\nATENCAO: O arquivo de referencia nao foi gerado.\n" ) ;
	  return ;
	}
      //================================================================================
      // CÁLCULO DO SOMATÓRIO DA NORMA
      //================================================================================
      for ( i = 0 ; i <= NX ; i ++ )
	{
	  for ( j = 0 ; j <= NY ; j ++ )
	    {
	      for ( k = 0 ; k < 3 ; k ++ )
		{
		  fscanf ( referencia , "%lf" , &a ) ;
		  fscanf ( arquivo , "%lf" , &b ) ;
		}
	      som += pow ( ( b - a ) , 2 ) ;
	    }
	}
      //================================================================================
      // CÁLCULO DA NORMA
      //================================================================================
      som = sqrt ( som * DX * DY ) ;
      //================================================================================
      // GERAÇÃO DO ARQUIVO PARA O ARMAZENAMENTO DA NORMA
      //================================================================================
      // ABERTURA DO ARQUIVO DA NORMA
      //================================================================================
      sprintf ( nomen , "mrk_%dordem.dat" , ORDEM ) ;
      norma = fopen ( nomen , "a" ) ;
      fprintf ( norma , "%d\t% .25lf\t% .25lf\t" , NT , DT , som ) ;
      fprintf ( norma , "% .25lf\t" , log10 ( DT ) ) ;
      fprintf ( norma , "% .25lf\t% .25lf\n" , log10 ( NT ) , log10 ( som ) ) ;
      //================================================================================
      // FECHAMENTO DOS ARQUIVOS ABERTOS
      //================================================================================
      fclose ( norma ) ;
      fclose ( arquivo ) ;
      fclose ( referencia ) ;
      //================================================================================
      // EXCLUSÃO DO ARQUIVO DE PHIs GERADO
      //================================================================================
      sprintf ( comando , "rm -f %s", nome ) ;
      system ( comando ) ;
      //================================================================================
      // GERAÇÃO DO GRÁFICO E CÁLCULO DA ORDEM DO MÉTODO
      //================================================================================
      if ( NT == 800 )
	{
	  //============================================================================
	  // GERAÇÃO DO COMANDO PARA O CÁLCULO DA ORDEM
	  //============================================================================
	   grafico = fopen ( "grafico_normas.gnu" , "w" ) ;
	   fprintf ( grafico , "plot '%s' u 5:6 w p pt 6 ps 1.5 lw 2\n" , nomen ) ;
	   fprintf ( grafico , "f(x) = a*x+b\n" ) ;
	   fprintf ( grafico , "fit f(x) '%s' u 5:6 via a,b" , nomen ) ;
	   fclose ( grafico ) ;
	   system ( "gnuplot grafico_normas.gnu" ) ;
	}
    }
}
