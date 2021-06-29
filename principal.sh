#======================================================================================
# SHELL-SCRIPT PARA EXECUÇÃO DO CÓDIGO-FONTE PRINCIPAL
#======================================================================================
#!/bin/sh
#======================================================================================
# LIMPA A TELA DO TERMINAL
#======================================================================================
clear
#======================================================================================
# ENTRADA DE DADOS PARA RESOLUÇÃO DA EDP
#======================================================================================
# DIMENSIONAMENTO ESPACIAL ( LX NO EIXO X E LY NO EIXO Y ) E TEMPORAL
#======================================================================================
LX=1
LY=1
LT=1
#======================================================================================
# PARTICIONAMENTO DO INTERVALO [0,LX]
#======================================================================================
#NX=60
NX=20
#======================================================================================
# ORDEM DO MÉTODO DE RUNGE-KUTTA PARA SER UTILIZADO - ORDEM PERMITIDAS: 2, 3 E 4
#======================================================================================
ORDEM=4
#======================================================================================
# FLAG PARA GERAÇÃO DE GRÁFICOS PARA ANIMAÇÃO. SE GA=0 NÃO GERAR IRÁ GERAR GRÁFICOS. SE
# GA=1 GRÁFICOS PARA ANIMAÇÃO SERÃO GERADOS.
#======================================================================================
GA=0
#======================================================================================
# CASO GA=1, INFORME ABAIXO O INTERVALO PARA GERAÇÃO DOS GRÁFICOS PARA ANIMAÇÃO.
# OBS: SE IGA=10, A CADA ITERAÇÃO DO MÉTODO DE RUNGE-KUTTA SERÁ GERADO UM GRÁFICO PNG.
#======================================================================================
IGA=10
#======================================================================================
# ESCOLHA DO TIPO DE MÉTODO PARA CÁLCULO DO VETOR GRADIENTE. SE GRAD=0, CÁLCULO VIA
# DIFERENÇAS FINITAS E SE GRAD=1, CÁLCULO VIA MÉTODO ADJUNTO.
#======================================================================================
GRAD=0
#======================================================================================
# ESCOLHA DO TIPO DE MÉTODO PARA CÁLCULO DA DIREÇÃO DE BUSCA NO MGC. SE MG=0, CÁLCULO
# VIA FLETCHER E REEVES; SE MG=1, CÁLCULO VIA POLAK E RIBIERE
#======================================================================================
MG=1
#======================================================================================
# CONDIÇÕES DE FRONTEIRA. SE O ÍNDICE FOR 0 A RESPECTIVA FRONTEIRA ESTÁ SOB CONDIÇÃO
# DE DIRICHLET E SE FOR 1 ESTÁ SOB CONDIÇÃO DE NEUMANN, SENDO:
# - F1 A FRONTEIRA ONDE 0<=X<=LX E Y=0
# - F2 A FRONTEIRA ONDE X=LX E 0<=Y<=LY
# - F3 A FRONTEIRA ONDE 0<=X<=LX E Y=LY
# - F4 A FRONTEIRA ONDE X=0 E 0<=Y<=LY
#======================================================================================
F1=1
F2=1
F3=1
F4=1
#======================================================================================
# EXECUÇÃO DO CÓDIGO-FONTE.
#======================================================================================
#NT=100000
#./principal $LX $LY $LT $NX $NT $ORDEM $GA $IGA $GRAD $MG $F1 $F2 $F3 $F4
#======================================================================================
# EXECUÇÃO DO CÓDIGO-FONTE PARA A VERIFICAÇÃO DA ORDEM.
#======================================================================================
ORDEM=2
NT=100
./principal $LX $LY $LT $NX $NT $ORDEM $GA $IGA $GRAD $MG $F1 $F2 $F3 $F4
NT=200
./principal $LX $LY $LT $NX $NT $ORDEM $GA $IGA $GRAD $MG $F1 $F2 $F3 $F4
NT=400
./principal $LX $LY $LT $NX $NT $ORDEM $GA $IGA $GRAD $MG $F1 $F2 $F3 $F4
NT=800
./principal $LX $LY $LT $NX $NT $ORDEM $GA $IGA $GRAD $MG $F1 $F2 $F3 $F4
