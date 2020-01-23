import microstructure as mStr
from numpy import linspace, pi, sin, cos, concatenate, zeros_like, zeros
from copy import deepcopy

# slab
slabSpeed = 50
xyres = 0.25
slabSz = [0.5, 0.5, 5]

# circle
circSpeed = 80
circR = 10
circPoints = 50
circElev = slabSz[2]
numCircles = 5

# sphere
sphereR = 1.0
sphereSpeed = 20

# line
lineLengthCoeff = 2  # kolko krat je rameno dlhsie ako polomer kruhu
lineSpeed = 50

# _______________________________________________________________________


def circle(radius):

    phivec = linspace(0, 2*pi, circPoints)
    arr = zeros((len(phivec), 5))
    arr[:, 0] = cos(phivec) * radius
    arr[:, 1] = sin(phivec) * radius
    arr[:, 2] = slabSz[2] + sphereR
    arr[:, 3] = 1
    arr[-1, 3:5] = [0, circSpeed]

    return arr


verts = linspace(0, 2*pi, 4)
slabPosX = cos(verts) * circR
slabPosY = sin(verts) * circR
viscoStruct = mStr.slabStr([0, 0, 0], [0, 1, 2], [xyres, 0.2], slabSpeed)
for i in range(3):
    slab = mStr.slabStr(slabSz, [0, 1, 2], [xyres, 0.2], slabSpeed)
    slab2 = mStr.slabStr([slabSz[0], slabSz[1], slabSz[2]+sphereR], [0, 1, 2], [xyres, 0.2], slabSpeed)
    slab.shift([slabPosX[i], slabPosY[i], 0])
    slab2.shift([lineLengthCoeff*slabPosX[i], lineLengthCoeff*slabPosY[i], 0])
    viscoStruct.addStr(slab)
    viscoStruct.addStr(slab2)

    sph = mStr.sphereStr(slabPosX[i], slabPosY[i], slabSz[2]+sphereR, sphereR,
                         sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
    viscoStruct.addStr(sph)

    lineArr = zeros((2, 5))
    lineArr[:, 0] = [slabPosX[i], lineLengthCoeff * slabPosX[i]]
    lineArr[:, 1] = [slabPosY[i], lineLengthCoeff * slabPosY[i]]
    lineArr[:, 2] = slabSz[2] + sphereR
    lineArr[:, 3] = 1
    lineArr[-1, 3:5] = [0, lineSpeed]
    viscoStruct.addStr(mStr.MicroStr(lineArr))


for i in range(numCircles):
    circ = mStr.MicroStr(circle(circR + i*0.1))
    viscoStruct.addStr(circ)

# viscoStruct.addStr(sph)
viscoStruct.plot(1, markerscalef=0.1)
