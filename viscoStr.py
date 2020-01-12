import microstructure as mStr
from numpy import linspace, pi, sin, concatenate, zeros_like, zeros
from copy import deepcopy

#yfrom tpp import *
#definujeme parametre na zaciatok
speed=50
linespeed=80
lineLength=10

sphereR=1.0
sphereSpeed=20

xyres=0.25

#make a slab - tu sa robia podstavce:
# parametre slabStr(xyz_sizes - velkosti podstavca vo vsetkych osiach
#, scanorder, stepsizes, speed) - zbytok je vysvetleny v dokumentacii
viscoStr=mStr.slabStr([5,5,10],[0,1,2],[xyres,0.5],speed)
viscoStr.shift([-2.5,-2.5,0])  # shift the structure

#podstavec sa okopiruje, a okopirovany sa posunie v osi x
slab2=deepcopy(viscoStr)
slab2.shift([(lineLength+2.0)*2,0,0])
#druhy podstavec vzniknuty okopirovanim sa prida do struktury ku
#prvemu - celu strukturu budeme skladat do premennej viscoStr
viscoStr.addStr(slab2)

#linspace(start, stop, num) - vrati num pocet rovnomerne rozdelenych
#hodnot z intervalu <start, stop>
x1=linspace(-pi/2,pi/2,num=25)
x2=linspace(pi/2,3*pi/2,num=50)
x3=linspace(3*pi/2,5*pi/2,num=25)

y1=sin(x1)*0.5+0.5
y2=sin(x2)
y3=sin(x3)*0.5-0.5

#[0:-1] znamena prvy riadok, posledny stlpec 
#[1:] znamena druhy riadok, vsetky stlpce
#pozn. lebo indexovanie je od nuly, concatenate spaja polia
x=concatenate((x1[0:-1],x2, x3[1:]))
y=concatenate((y1[0:-1],y2, y3[1:]))
#vytvori pole s deviatkami o rovnakej dlzke ako ma x
z=zeros_like(x)+9

#set the length:
#ptp vracia (maximum-minimum) hodnotu, cize x bude v rozsahu <0,1>
x=(x-min(x))/x.ptp()
#vynasobi vsetky cleny pola x dva krat hodnotou lineLength
x *= 2 * lineLength
#vynasobi vsetky cleny pola y 2.5 - zmeni amplitudu sinusovky:
y *= 2.5

#zeros vytvara maticu nul s prislusnymi rozmermi: riadkov bude rovnako
#ako je dlzka polia x a stlpcov 5 - to je ta matica ktora sa musi dat ako
#vstupny parameter ked vytvarame lubovolnu strukturu
sinLine=zeros((len(x),5))
#prvy stlpec matice budu hodnoty x:
#pozn. : znamena vsetky (v tomto pripade riadky)
sinLine[:,0]=x
#druhy stlpec bude y
sinLine[:,1]=y
#treti stlpec bude z (v tomto pripade vsade bude hodnota 9)
sinLine[:,2]=z
#stvrty stlpec budu same jednotky (stplec prisluchajuci shutteru)
sinLine[:,3]=1
#nasledujuci prikaz znamena, ze do posledneho riadku sa do
#stvrteho stlpca prida 0 (uzatvorenie shuttera) a do piateho stlpca
#hodnota linespeed - rychlost polymerizacie (inak ostali v celom 
#5tom stlpci nulky)
sinLine[-1,3:5]=[0,linespeed]

#tu sa vytvori objekt mikrostruktury na zaklade
#matice sinLine a prida sa do sktruktury ktora uz obsahuje podstavce
line=mStr.MicroStr(sinLine)
line.shift([2.0,0,0]) 
viscoStr.addStr(line)

#nakoniec sa vytvori a prida do struktury sfera, parametre: prve tri su 
#suradnice stredu, polomer, rychlost, xyres a dalsie dva parametre
#neviem co su zac,nema to vysvetlene ani v dokumentacii, ale shell spacing
#je vzdialenost medzi jednotlivymi vrstavami gule
sph=mStr.sphereStr(lineLength+2.0,0,9,sphereR,sphereSpeed,xyres,1.0,1,shellspacing=0.5)
viscoStr.addStr(sph)

#vykresli sa vysledna struktura
viscoStr.plot(1,markerscalef=0.1)
c=0