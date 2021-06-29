#!/bin/sh

clear

echo "=================================================================================================================="
echo "Gerando as animacoes ..."
echo "=================================================================================================================="

cd graficos
rm -f *.gif *.avi
#convert -delay 25 1*.png animacao1.gif
#convert -delay 25 2*.png animacao2.gif
mencoder mf://1*.png -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o animacao1.avi
mencoder mf://2*.png -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o animacao2.avi
#mencoder mf://1*.png -mf w=800:h=600:fps=25:type=png -ovc raw -oac copy -o animacao3.avi
#mencoder mf://2*.png -mf w=800:h=600:fps=25:type=png -ovc raw -oac copy -o animacao4.avi
rm -f 1*.png 2*.png
cd ..

echo "=================================================================================================================="
echo "Animacoes geradas com sucesso !!!"
echo "=================================================================================================================="



