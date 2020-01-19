import microstructure as mStr
from numpy import linspace, pi, sin, cos, concatenate, zeros_like, zeros
from copy import deepcopy

# slab
slabSpeed = 50
xyres = 0.25
slabSz = [6, 5, 8] 
slabExtSz = [8, 5, 1.3]
slabSzDifX = slabExtSz[0] - slabSz[0]

# line
linespeed = 80
lineLength = 10
linePoints = 20
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
        z[i + lineLength] = lineElev + 2 * yCos[i]

    sinLine = zeros((len(x), 5))
    sinLine[:, 0] = x
    sinLine[:, 1] = y
    sinLine[:, 2] = z
    sinLine[:, 3] = 1
    sinLine[-1, 3:5] = [0, linespeed]

    return sinLine


slab1 = mStr.slabStr(slabSz, [0, 1, 2], [xyres, 0.5], slabSpeed)
slab1.shift([-slabSz[0], -slabSz[1]/2, 0])
viscoStruct = slab1

slabExt = mStr.slabStr(slabExtSz, [0, 1, 2], [xyres, 0.5], slabSpeed)
slabExt.shift([-(slabSz[0] + slabSzDifX/2), -slabExtSz[1]/2, slabSz[2]])
viscoStruct.addStr(slabExt)

slab2 = deepcopy(viscoStruct)
slab2.shift([(lineLength)*2 + slabSz[0], 0, 0])
viscoStruct.addStr(slab2)

sinLine = sinusLine()
line = mStr.MicroStr(sinLine)
viscoStruct.addStr(line)

# nakoniec sa vytvori a prida do struktury sfera, parametre: prve tri su
# suradnice stredu, polomer, rychlost, xyres a dalsie dva parametre
# neviem co su zac,nema to vysvetlene ani v dokumentacii, ale shell spacing
# je vzdialenost medzi jednotlivymi vrstavami gule
sph = mStr.sphereStr(lineLength, 0, lineElev, sphereR,
                     sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
viscoStruct.addStr(sph)
viscoStruct.plot(1, markerscalef=0.1)


