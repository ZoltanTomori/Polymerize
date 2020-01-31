import microstructure as mStr
from numpy import linspace, pi, sin, cos, concatenate, zeros_like, zeros
from copy import deepcopy

# _______________________________________________________________________

# slab
numSlabs = 4            # pocet stlpikov podopierajucich gulovite uchytky
numRays = 2             # pocet lucov
slabSpeed = 50          # rychlost polymerizovania stlpika
xyres = 0.25            # vzdialenost medzi ciarami pri budovani stlpika
slabSz = [1, 1, 3]      # rozmer stlpika [x,y,z]

# circle
circSpeed = 80          # rychlost polymerizacie kruhu
circR = 6               # polomer kruhu
circPoints = 50         # pocet bodov kruznice
numCircles = 2          # pocet kruznic tvoriacich polymerizovany kruh v oboch smeroch
circDist = 0.1          # vzdialenost kruznic pri polymerizovani kruhu

# sphere
sphereR = 1.0           # polomer gulicky na stlpikoch
sphereRayR = 0.5        # polomer gulicky "zdureniny" z ktrej ide luc z kruhu
sphereSpeed = 20        # rychlost polymerizovania gulicky

# arm (usecky alebo aj krivky ako napr. sinus)
armEndDist = 18         # vzdialenost konca ramena od stredu (dlzka ramena bude armDist-circR)
armBallDist = 12        # vzdialenost gulicky na ramene od stredu (pre sinus musi byt v strede krivky!)
armSpeed = 50           # rychlost polymerizacie luca
armPoints = 50          # pocet bodov - iba v pripade kriviek
armAmplit = 1.5         # amplituda Y (v pripade krivky)

# _______________________________________________________________________

lineLengthCoeff = armEndDist / circR

def DoLine(P1, P2):
    lineArr = zeros((2, 5))
    lineArr[:, 0] = [P1[0], P2[0]]
    lineArr[:, 1] = [P1[1], P2[1]]
    lineArr[:, 2] = slabSz[2] + sphereR
    lineArr[:, 3] = 1
    lineArr[-1, 3:5] = [0, armSpeed]
    return lineArr

def circle(radius, elev):
    phivec = linspace(0, 2*pi, circPoints)
    arr = zeros((len(phivec), 5))
    arr[:, 0] = cos(phivec) * radius
    arr[:, 1] = sin(phivec) * radius
    arr[:, 2] = elev
    arr[:, 3] = 1
    arr[-1, 3:5] = [0, circSpeed]
    return arr

def SinusLine(fromPt, toPt):
    x1 = linspace(-pi/2, pi/2, armPoints/2)
    x2 = linspace(pi/2, 3*pi/2, armPoints)
    x3 = linspace(3*pi/2, 5*pi/2, armPoints/2)
    y1 = sin(x1)*0.5+0.5
    y2 = sin(x2)
    y3 = sin(x3)*0.5-0.5
    x = concatenate((x1[0:-1], x2, x3[1:]))
    y = concatenate((y1[0:-1], y2, y3[1:]))
    z = P1[2]

    x = (x-min(x))/x.ptp()
    x = x * (toPt[0] - fromPt[0]) + fromPt[0]
    y = y * armAmplit

    sinLine = zeros((len(x), 5))
    sinLine[:, 0] = x
    sinLine[:, 1] = y
    sinLine[:, 2] = z
    sinLine[:, 3] = 1
    sinLine[-1, 3:5] = [0, armSpeed]

    return sinLine

# rozdeli 0-360 stupnov na 3 casti (4 okrajove body)
verts = linspace(pi/4, 9*pi/4, numSlabs+1)
slabPosX = cos(verts) * circR
slabPosY = sin(verts) * circR
slab0 = mStr.slabStr([0.2, 0.2, 0.2], [0, 1, 2], [xyres, 0.2], slabSpeed)
viscoStruct = slab0
for i in range(numSlabs):
    #stlpik pod kruhom ukonceny gulickou
    slab1 = mStr.slabStr([slabSz[0], slabSz[1], slabSz[2] + sphereR/3],[0, 1, 2], [xyres, 0.2], slabSpeed)
    slab1.shift([slabPosX[i]-slabSz[0]/2, slabPosY[i]-slabSz[1]/2, 0])
    viscoStruct.addStr(slab1)
    sph1 = mStr.sphereStr(slabPosX[i], slabPosY[i], slabSz[2]+sphereR, sphereR, sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
    viscoStruct.addStr(sph1)

# telo kruhu
for i in range(-numCircles, numCircles):
    circElevation = slabSz[2] + sphereR + i*circDist
    for j in range(-numCircles, numCircles):
        circ = mStr.MicroStr(circle(circR + j*circDist, circElevation))
        viscoStruct.addStr(circ)

# rozdeli 45-405 stupnov na 3 casti (4 okrajove body)
vertsRay = linspace(0, 2*pi, numRays+1)
slabRayPosX = cos(vertsRay) * circR
slabRayPosY = sin(vertsRay) * circR

 # stlpiky pod zaciatkom a koncom luca ukoncene gulickou
for i in range(numRays):
    sphRay1 = mStr.sphereStr(slabRayPosX[i], slabRayPosY[i], slabSz[2]+sphereR, sphereRayR, sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
    viscoStruct.addStr(sphRay1)

    slab2 = mStr.slabStr([slabSz[0], slabSz[1], slabSz[2] + sphereR/3],[0, 1, 2], [xyres, 0.2], slabSpeed) 
    slab2.shift([lineLengthCoeff*slabRayPosX[i]-slabSz[0]/2, lineLengthCoeff*slabRayPosY[i]-slabSz[1]/2, 0])
    viscoStruct.addStr(slab2)
    sphRay2 = mStr.sphereStr(lineLengthCoeff*slabRayPosX[i], lineLengthCoeff*slabRayPosY[i], slabSz[2]+sphereR, sphereR, sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
    viscoStruct.addStr(sphRay2)
   
# luce az nakoniec lebo su tenke
for i in range(numRays):
    P1 = [slabRayPosX[i], slabRayPosY[i], slabSz[2] + sphereR]
    P2 = [slabRayPosX[i] * lineLengthCoeff, slabRayPosY[i] * lineLengthCoeff, slabSz[2] + sphereR]
    lineArr = SinusLine(P1, P2)
    #lineArr = DoLine(P1, P2)
    viscoStruct.addStr(mStr.MicroStr(lineArr))

# gulicka do stredu luca

armCoef = armBallDist / circR
for i in range(numRays):
    sphRay2 = mStr.sphereStr(armCoef*slabRayPosX[i], armCoef*slabRayPosY[i], slabSz[2]+sphereR, sphereRayR, sphereSpeed, xyres, 1.0, 1, shellspacing=0.5)
    viscoStruct.addStr(sphRay2)
    
viscoStruct.plot(1, markerscalef=0.1)
