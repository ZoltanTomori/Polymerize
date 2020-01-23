import microstructure as mStr
from numpy import linspace, pi, sin, cos, concatenate, zeros_like, zeros
from copy import deepcopy

# slab
slabSpeed = 50
xyres = 0.25
slabSz = [5, 5, 8]
slabExtSz = [7, 5, 1]
slabSzDifX = slabExtSz[0] - slabSz[0]

# line
linespeed = 80
lineLength = 10
linePoints = 50
lineDrop = 1
lineElev = slabSz[2] + slabExtSz[2]/2

# sphere
sphereR = 1.0
sphereSpeed = 20

# _______________________________________________________________________


def sinusLine():
    x1 = linspace(-pi/2, pi/2, linePoints/2)
    x2 = linspace(pi/2, 3*pi/2, linePoints)
    x3 = linspace(3*pi/2, 5*pi/2, linePoints/2)
    y1 = sin(x1)*0.5+0.5
    y2 = sin(x2)
    y3 = sin(x3)*0.5-0.5
    x = concatenate((x1[0:-1], x2, x3[1:]))
    y = concatenate((y1[0:-1], y2, y3[1:]))
    z = zeros_like(x) + lineElev

    x = (x-min(x))/x.ptp()
    x *= 2 * lineLength
    y *= lineDrop

    sinLine = zeros((len(x), 5))
    sinLine[:, 0] = x
    sinLine[:, 1] = y
    sinLine[:, 2] = z
    sinLine[:, 3] = 1
    sinLine[-1, 3:5] = [0, linespeed]

    return sinLine

def halfSinLine():
    x1 = linspace(0, lineLength, linePoints/2)
    y1 = zeros_like(x1)   
    xSin = linspace(0, 2*pi, linePoints/2)
    ySin= sin(xSin)
    x2 = lineLength + lineLength * xSin/(2*pi)
    y2 = ySin * lineDrop

    x = concatenate((x1[0:-1], x2[1:]))
    y = concatenate((y1[0:-1], y2[1:]))
    z = zeros_like(x) + lineElev


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
sinLine = halfSinLine() #sinusLine()
line = mStr.MicroStr(sinLine)
viscoStruct.addStr(line)

# v polovici ciary vytvori gulu
sph = mStr.sphereStr(lineLength, 0, lineElev, sphereR,
                     sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)

viscoStruct.addStr(sph)
viscoStruct.plot(1, markerscalef=0.1)
