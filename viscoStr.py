import microstructure as mStr
from numpy import linspace, pi, sin, cos, concatenate, zeros_like, zeros
from copy import deepcopy

# slab
slabSpeed = 50
xyres = 0.25
slabSz = [6, 5, 8] 
slabExtSz = [8, 5, 1]
slabSzDifX = slabExtSz[0] - slabSz[0]

# line
linespeed = 80
lineLength = 10
linePoints = 20
lineDrop = 1
lineElev = slabSz[2] + slabExtSz[2]/2

# sphere
sphereR = 1.0
sphereSpeed = 20

# _______________________________________________________________________
def sinusLine():
    x = linspace(0, 2 * lineLength, linePoints)
    y = zeros_like(x)
    z = zeros_like(x) + lineElev
   
    xCos = linspace(0, 2 * pi, linePoints/2)
    yCos = cos(xCos) - 1
    for i in range(len(xCos)):       
        z[i + lineLength] = lineElev + lineDrop * yCos[i] / 2

    sinLine = zeros((len(x), 5))
    sinLine[:, 0] = x
    sinLine[:, 1] = y
    sinLine[:, 2] = z
    sinLine[:, 3] = 1
    sinLine[-1, 3:5] = [0, linespeed]

    return sinLine

# vytvori stlpec a posunieho tak, aby (0,0) bolo v strede celnej steny
slab1 = mStr.slabStr(slabSz, [0, 1, 2], [xyres, 0.5], slabSpeed)
slab1.shift([-slabSz[0], -slabSz[1]/2, 0])
viscoStruct = slab1

# na stlpec ulozi rozsirenie - manzetu
slabExt = mStr.slabStr(slabExtSz, [0, 1, 2], [xyres, 0.5], slabSpeed)
slabExt.shift([-(slabSz[0] + slabSzDifX/2), -slabExtSz[1]/2, slabSz[2]])
viscoStruct.addStr(slabExt)

# skopiruje stlpec aj s rozsirenim do druhehe a posinie ho o 2 * lineLength
slab2 = deepcopy(viscoStruct)
slab2.shift([(lineLength)*2 + slabSz[0], 0, 0])
viscoStruct.addStr(slab2)

# vyutvori ciaru ktorej prva polovica bude usecka a druha sinusovka
sinLine = sinusLine()
line = mStr.MicroStr(sinLine)
viscoStruct.addStr(line)

# v polovici ciary vytvori gulu
sph = mStr.sphereStr(lineLength, 0, lineElev, sphereR,
                     sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
viscoStruct.addStr(sph)
viscoStruct.plot(1, markerscalef=0.1)


