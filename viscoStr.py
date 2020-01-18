import microstructure as mStr
from numpy import linspace, pi, sin, cos, concatenate, zeros_like, zeros
from copy import deepcopy

# yfrom tpp import *
# definujeme parametre na zaciatok
speed = 50
linespeed = 80
lineLength = 10

sphereR = 1.0
sphereSpeed = 20

xyres = 0.25
slabXY = 5 
slabZ = 10
lineFromGround = 7

# _______________________________________________________________________
def sinusLine():
    x = linspace(0, 2 * lineLength, num=20)
    y = zeros_like(x)
    z = zeros_like(x)+9
    for i in range(len(x)):
        y[i] = 0
        z[i] = 9

    xCos = linspace(0, 2 * pi, num=10)
    yCos = cos(xCos) - 1
    for i in range(len(xCos)):       
        y[i + lineLength] = yCos[i] * 5
        z[i + lineLength] = 9

    sinLine = zeros((len(x), 5))
    sinLine[:, 0] = x
    sinLine[:, 1] = y
    sinLine[:, 2] = z
    sinLine[:, 3] = 1
    sinLine[-1, 3:5] = [0, linespeed]

    return sinLine


viscoStr = mStr.slabStr([slabXY, slabXY, slabZ], [0, 1, 2], [xyres, 0.5], speed)
viscoStr.shift([-slabXY, -slabXY/2, 0])
slab2 = deepcopy(viscoStr)
slab2.shift([(lineLength)*2 + slabXY, 0, 0])
viscoStr.addStr(slab2)

sinLine = sinusLine()
line = mStr.MicroStr(sinLine)
#line.shift([2.0, 0, 0])
viscoStr.addStr(line)

# nakoniec sa vytvori a prida do struktury sfera, parametre: prve tri su
# suradnice stredu, polomer, rychlost, xyres a dalsie dva parametre
# neviem co su zac,nema to vysvetlene ani v dokumentacii, ale shell spacing
# je vzdialenost medzi jednotlivymi vrstavami gule
sph = mStr.sphereStr(lineLength, 0, 9, sphereR,
                     sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
viscoStr.addStr(sph)

# vykresli sa vysledna struktura
viscoStr.plot(1, markerscalef=0.1)
c = 0


