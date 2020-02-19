import microstructure as mStr
from numpy import linspace, pi, sin, cos, concatenate, zeros_like, zeros
from copy import deepcopy

# slab
slabSpeed = 50
xyres = 0.25
slabSz = [5, 5, 8]
slabExtSz = [8, 5, 1]
slabSzDifX = slabExtSz[0] - slabSz[0]

# line
linespeed = 80
lineLength = 10
linePoints = 50
lineDrop = 1
lineElev = slabSz[2] + slabExtSz[2]/2

# wall 
circSpeed = 80          # rychlost polymerizacie kruhu
circR = 30              # polomer kruhu
circPoints = 50         # pocet bodov kruznice
numCircles = [4, 70]    # pocet kruznic tvoriacich polymerizovany kruh v smeroch [X, Z]
circDist = 0.2          # vzdialenost kruznic pri polymerizovani kruhu

# sphere
sphereR = 1.0
sphereSpeed = 20

# cross
crossSpeed = 20

# _______________________________________________________________________

def cross(sz):
    arr = zeros((2,5))
    arr[:, 0] = [-sz, sz]
    arr[:, 1] = [-sz, sz]
    arr[:, 2] = [0, 0.1]
    arr[:, 3] = 1
    arr[-1, 3:5] = [0, crossSpeed]
    return arr

def circle(radius, elev):
    phivec = linspace(0, 2*pi, circPoints)
    arr = zeros((len(phivec), 5))
    arr[:, 0] = cos(phivec) * radius
    arr[:, 1] = sin(phivec) * radius
    arr[:, 2] = elev
    arr[:, 3] = 1
    arr[-1, 3:5] = [0, circSpeed]
    return arr

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

# prazdna struktura
viscoStruct = mStr.MicroStr(zeros((0,5)))

# origin
origin = mStr.MicroStr(circle(0.1, 0))
origin.shift([lineLength, 0, 0])
viscoStruct.addStr(origin)

# vytvori stlpec + hlavu a posunieho tak, aby (0,0) bolo v strede celnej steny
slab1 = mStr.slabStr(slabSz, [0, 1, 2], [xyres, 0.5], slabSpeed)
slab1Ext = mStr.slabStr(slabExtSz, [0, 1, 2], [xyres, 0.5], slabSpeed)
slab1Ext.shift([-slabSzDifX/2, 0, slabSz[2]])
slab1.addStr(slab1Ext)
slab1.shift([-slabSz[0], -slabSz[1]/2, 0])
viscoStruct.addStr(slab1)

# skopiruje stlpec aj s rozsirenim do druheho a posinie ho o 2 * lineLength
slab2 = deepcopy(slab1)
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

# telo kruhovej ohrady
wall = mStr.MicroStr(zeros((0,5)))
for i in range(numCircles[1]):
    circElevation = i*circDist
    for j in range(numCircles[0]):
        circ = mStr.MicroStr(circle(circR + j*circDist, circElevation))
        wall.addStr(circ)
wall.shift([lineLength, 0, 0])
viscoStruct.addStr(wall)

viscoStruct.plot(1, markerscalef=0.1)
