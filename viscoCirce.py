import microstructure as mStr
from numpy import linspace, pi, sin, cos, concatenate, zeros_like, zeros
from copy import deepcopy

# slab
slabSpeed = 50
xyres = 0.25
slabSz = [1, 1, 5]

# circle
circSpeed = 80
circDiam = 20
circPoints = 50
circElev = slabSz[2]

# sphere
sphereR = 1.0
sphereSpeed = 20

# _______________________________________________________________________


def circle():
    
    phivec = linspace(0, 2*pi, circPoints)
    arr = zeros((len(phivec), 5))
    arr[:, 0] = cos(phivec)
    arr[:, 1] = sin(phivec)
    arr[:, 2] = 9
    arr[:, 3] = 1
    arr[-1, 3:5] = [0, circSpeed]
    
    return arr

slab0 = mStr.slabStr(slabSz, [0, 1, 2], [xyres, 0.5], slabSpeed)
slab0.shift([2, 0, 0])
viscoStruct = slab0

slab1 = deepcopy(slab0)
slab1.shift([4, 0, 0])
viscoStruct.addStr(slab1)

slab2 = deepcopy(slab0)
slab2.shift([6, 0, 0])
viscoStruct.addStr(slab2)


# v polovici ciary vytvori gulu
#sph = mStr.sphereStr(lineLength, 0, lineElev, sphereR,  sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
circ = mStr.MicroStr(circle())
viscoStruct.addStr(circ)

#viscoStruct.addStr(sph)
viscoStruct.plot(1, markerscalef=0.1)
