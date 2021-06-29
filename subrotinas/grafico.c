//======================================================================================
// SUBROTINA PARA A GERAÇÃO DE GRÁFICOS
//======================================================================================

void GRAFICO ( double ***phi )
{
  //====================================================================================
  // DECLARAÇÃO DE VARIÁVEIS
  //====================================================================================
  int i ;            // VARIÁVEL PARA O DESLOCAMENTO NA DIREÇÃO X
  int j ;            // VARIÁVEL PARA O DESLOCAMENTO NA DIREÇÃO Y
  int k ;            // VARIÁVEL PARA O DESLOCAMENTO TEMPORAL
  double x ;         // VARIÁVEL REAL PARA O DESLOCAMENTO NA DIREÇÃO X
  double y ;         // VARIÁVEL REAL PARA O DESLOCAMENTO NA DIREÇÃO Y
  double t ;         // VARIÁVEL REAL PARA O DESLOCAMENTO TEMPORAL
  FILE *grafico ;    // PONTEIRO PARA O ARMAZENAMENTO DE VALORES EM ARQUIVO PARA GRÁFICO
  FILE *arquivo ;    // PONTEIRO PARA O ARMAZENAMENTO DE VALORES EM ARQUIVO
  //====================================================================================
  // GERAÇÃO DOS GRÁFICOS
  //====================================================================================
  for ( k = 0 ; k <= NT ; k += IGA )
    {
      t = k * DT ;
      //============================================================================
      // DECLARAÇÃO DE VARIÁVEIS
      //============================================================================
      char nome [ 100 ] ;               // NOME DO GRÁFICO GERADO
      //============================================================================
      // GERAÇÃO DO ARQUIVO COM VALORES DE PHI
      //============================================================================
      grafico = fopen ( "graficos/phi.dat" , "w" ) ;
      arquivo = fopen ( "referencia.dat" , "a" ) ;
      if ( grafico == NULL )
	{
	  system ( "mkdir graficos" ) ;
	  grafico = fopen ( "graficos/phi.dat" , "w" ) ;
	}
      for ( j = 0 ; j <= NY ; j++ )
	{
	  y = j * DY ;
	  for ( i = 0 ; i <= NX ; i++ )
	    {
	      x = i * DX ;
	      fprintf ( grafico , "% lf\t% lf\t" , x , y ) ;
	      fprintf ( grafico , "% lf\n" , phi [ i ] [ j ] [ k ] ) ;
	      fprintf ( arquivo , "% .25lf\t% .25lf\t" , x , y ) ;
	      fprintf ( arquivo , "% .25lf\n" , phi [ i ] [ j ] [ k ] ) ;
	    }
	  fprintf ( grafico , "\n" ) ;
	  fprintf ( arquivo , "\n" ) ;
	}
      fclose ( grafico ) ;
      fclose ( arquivo ) ;
      
      grafico = fopen ( "graficos/animacao_plana.gnu" , "w" ) ;

      fprintf ( grafico , "set size ratio -1\n" ) ;
      fprintf ( grafico , "set xlabel \"x\"\n" ) ;
      fprintf ( grafico , "set ylabel \"y\"\n" ) ;
      fprintf ( grafico , "set zlabel \"$\phi(x,y,t)$\"\n" ) ;
      fprintf ( grafico , "set zrange [0:100]\n") ;
      fprintf ( grafico , "set cbrange [0:100]\n") ;  
      fprintf ( grafico , "set format x \"%s\"\n" , "% .3f" ) ;
      fprintf ( grafico , "set format y \"%s\"\n" , "% .3f" ) ;
      fprintf ( grafico , "unset key\n" ) ;
      fprintf ( grafico , "set hidden3d\n" ) ; 
      fprintf ( grafico , "set view map\n" ) ;
      fprintf ( grafico , "set pm3d interpolate 8,8\n" ) ;
      fprintf ( grafico , "set palette rgbformulae 22,13,-31\n" ) ;  
      fprintf ( grafico , "set term png\n" ) ;
      sprintf ( nome , "graficos/%d" , 10000000 + k ) ; 
      strcat ( nome , ".png" ) ;
      fprintf ( grafico , "set output \"%s\"\n" , nome ) ;
      fprintf ( grafico , "set title \"Distribuição de Temperatura em t=% .3f s\"\n" , t ) ; 
      fprintf ( grafico , "splot \"graficos/phi.dat\" w pm3d\n" ) ; 

      fclose ( grafico ) ;

      system ( "gnuplot graficos/animacao_plana.gnu" ) ;

      grafico = fopen ( "graficos/animacao3d.gnu" , "w" ) ;

      fprintf ( grafico , "set size ratio -1\n" ) ;
      fprintf ( grafico , "set xlabel \"x\"\n" ) ;
      fprintf ( grafico , "set ylabel \"y\"\n" ) ;
      fprintf ( grafico , "set zlabel \"Φ\"\n" ) ;
      fprintf ( grafico , "set zrange [0:100]\n") ;
      fprintf ( grafico , "set cbrange [0:100]\n") ;
      fprintf ( grafico , "set format x \"%s\"\n" , "% .3f" ) ;
      fprintf ( grafico , "set format y \"%s\"\n" , "% .3f" ) ;
      fprintf ( grafico , "unset key\n" ) ;
      fprintf ( grafico , "set pm3d interpolate 8,8\n" ) ;
      fprintf ( grafico , "set palette rgbformulae 22,13,-31\n" ) ;  
      fprintf ( grafico , "set term png\n" ) ;
      sprintf ( nome , "graficos/%d" , 20000000 + k ) ; 
      strcat ( nome , ".png" ) ;
      fprintf ( grafico , "set output \"%s\"\n" , nome ) ;
      fprintf ( grafico , "set title \"Distribuição de Temperatura em t=% .3f s\"\n" , t ) ; 
      fprintf ( grafico , "splot \"graficos/phi.dat\" w pm3d\n" ) ; 

      fclose ( grafico ) ;

      system ( "gnuplot graficos/animacao3d.gnu" ) ;  
    }
    
}
