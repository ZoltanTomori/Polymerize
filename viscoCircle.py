import microstructure as mStr
from numpy import linspace, pi, sin, cos, concatenate, zeros_like, zeros
from copy import deepcopy

# slab
slabSpeed = 50
xyres = 0.25
slabSz = [1, 1, 3]

# circle
circSpeed = 80
circR = 8
circPoints = 50
circElev = slabSz[2]
numCircles = 3

# sphere
sphereR = 1.0
sphereSpeed = 20

# line
lineLengthCoeff = 2  # kolko krat je rameno dlhsie ako polomer kruhu
lineSpeed = 50

# _______________________________________________________________________


def circle(radius, elev):
    phivec = linspace(0, 2*pi, circPoints)
    arr = zeros((len(phivec), 5))
    arr[:, 0] = cos(phivec) * radius
    arr[:, 1] = sin(phivec) * radius
    arr[:, 2] = elev
    arr[:, 3] = 1
    arr[-1, 3:5] = [0, circSpeed]
    return arr

verts = linspace(0, 2*pi, 4)
slabPosX = cos(verts) * circR
slabPosY = sin(verts) * circR
slab0 = mStr.slabStr([0.2, 0.2, 0.2], [0, 1, 2], [xyres, 0.2], slabSpeed)
viscoStruct = slab0
for i in range(3):
    #stlpik pod kruhom ukonceny gulickou
    slab1 = mStr.slabStr([slabSz[0], slabSz[1], slabSz[2] + sphereR/3],[0, 1, 2], [xyres, 0.2], slabSpeed)
    slab1.shift([slabPosX[i]-slabSz[0]/2, slabPosY[i]-slabSz[1]/2, 0])
    viscoStruct.addStr(slab1)
    sph1 = mStr.sphereStr(slabPosX[i], slabPosY[i], slabSz[2]+sphereR, sphereR, sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
    viscoStruct.addStr(sph1)

    # stlpik pod koncom luca ukonceny gulickou
    slab2 = mStr.slabStr([slabSz[0], slabSz[1], slabSz[2] + sphereR/3],[0, 1, 2], [xyres, 0.2], slabSpeed) 
    slab2.shift([lineLengthCoeff*slabPosX[i]-slabSz[0]/2, lineLengthCoeff*slabPosY[i]-slabSz[1]/2, 0])
    viscoStruct.addStr(slab2)
    sph2 = mStr.sphereStr(lineLengthCoeff*slabPosX[i], lineLengthCoeff*slabPosY[i], slabSz[2]+sphereR, sphereR, sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
    viscoStruct.addStr(sph2)

    # luc
    lineArr = zeros((2, 5))
    lineArr[:, 0] = [slabPosX[i], lineLengthCoeff * slabPosX[i]]
    lineArr[:, 1] = [slabPosY[i], lineLengthCoeff * slabPosY[i]]
    lineArr[:, 2] = slabSz[2] + sphereR
    lineArr[:, 3] = 1
    lineArr[-1, 3:5] = [0, lineSpeed]
    viscoStruct.addStr(mStr.MicroStr(lineArr))

# telo kruhu
for i in range(-numCircles, numCircles):
    elevation = slabSz[2] + sphereR + i*0.1
    for j in range(-numCircles, numCircles):
        circ = mStr.MicroStr(circle(circR + j*0.1, elevation))
        viscoStruct.addStr(circ)

viscoStruct.plot(1, markerscalef=0.1)
