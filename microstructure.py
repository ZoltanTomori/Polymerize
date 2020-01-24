# CREATED BY GASZTON VIZSNYICZAI, Biological Research Centre, Szeged Hungary, gaszton@brc.hu
from matplotlib import cm
from mayavi import mlab
from copy import deepcopy
from numpy import *
from matplotlib import pyplot
import colorsys

DEFAULT_IP_STEP = 0.05


class Stroke():
    def __init__(self, coors, speed, ipstep, holo=None):
        self.coors = coors.copy()  # node coordinates
        self.speed = speed
        self.ipstep = ipstep
        self.icoors()
        self.holo = holo

    def icoors(self):
        """returns the array of interpolated coordinates and the time step dt"""

        dr = diff(self.coors, 1, 0)
        self.length = sum(sqrt(sum(dr**2, axis=1)))

        # use lists for the interpolated coordinates (numpy array append very slow)
        iX = [self.coors[0, 0]]
        iY = [self.coors[0, 1]]
        iZ = [self.coors[0, 2]]

        eps = 1./(2.**15)

        for i in range(len(dr)):
            L = sqrt(sum(dr[i, :]**2))
            # recalculate the stepsize for a round number of steps:
            stepN = ceil(L/self.ipstep-eps)
            stepsize = L/stepN
            # vector for interpolation:
            rs = dr[i, :]/L*stepsize
            # calculate the interpolated coordinates:
            ii = arange(1, stepN+1)
            iX += list(self.coors[i, 0]+ii*rs[0])
            iY += list(self.coors[i, 1]+ii*rs[1])
            iZ += list(self.coors[i, 2]+ii*rs[2])

        icoors = array([iX, iY, iZ]).transpose()
        self.dt = self.length/(len(icoors)-1)/self.speed
        return icoors, self.dt


class MicroStr():
    def __init__(self, strdata, ipstep=None, holo=None):
        assert size(strdata, 1) == 5, 'strdata must be an N-by-5 array'
        assert strdata.ndim == 2, 'strdata must be an array with 2 dimensions'

        strdata = strdata.copy()
        if ipstep is None:
            self.ipstep = DEFAULT_IP_STEP
        else:
            self.ipstep = ipstep
        self.buildStrokes(strdata, holo)

    def addStr(self, a):
        for s in a.Strokes:
            self.Strokes += [deepcopy(s)]

    def buildStrokes(self, strdata, holo=None):

        shutter = strdata[:, 3]  # where==0 marks the last element of a stroke
        strokeEnds = nonzero((shutter == 0))[0]
        strokeStart = 0
        Nstrokes = len(strokeEnds)
        self.Strokes = []

        for si in strokeEnds:
            speed = strdata[si, 4]
            coors = strdata[strokeStart:si+1, 0:3]
            self.Strokes += [Stroke(coors, speed, self.ipstep, holo)]
            strokeStart = si+1

    def set_IPstepSize(self, ipstep):
        self.ipstep = ipstep
        for s in self.Strokes:
            s.ipstep = ipstep
            s.icoors()

    def scanTime(self):
        t = 0
        for s in self.Strokes:
            t += s.length/s.speed
        return t

    def scanLength(self):
        d = 0
        for s in self.Strokes:
            d += s.length
        return d

    def shift(self, xyzShift):
        assert size(xyzShift) == 3, 'xyzShift must have 3 elements'
        for s in self.Strokes:
            s.coors[:, 0:3] += array(xyzShift)

    def scale(self, xyzScale):
        assert size(xyzScale) == 3, 'xyzScale must have 3 elements'
        for s in self.Strokes:
            s.coors[:, 0:3] *= array(xyzScale)

    def rotate(self, angle, axisstring):
        """rotates the structure coordinates with angle degrees
           around the axis specified by axisstring:'x' 'y' 'z'"""
        angle = angle/180.0*pi
        if axisstring == 'x':
            rotmat = array([[1, 0, 0],
                            [0, cos(angle), -sin(angle)],
                            [0, sin(angle), cos(angle)]])
        elif axisstring == 'y':
            rotmat = array([[cos(angle), 0, sin(angle)],
                            [0, 1, 0],
                            [-sin(angle), 0, cos(angle)]])
        elif axisstring == 'z':
            rotmat = array([[cos(angle), -sin(angle), 0],
                            [sin(angle), cos(angle), 0],
                            [0, 0, 1]])
        for s in self.Strokes:
            npoints = size(s.coors, 0)
            for i in range(npoints):
                tmp = s.coors[i, 0:3].copy()
                tmp = dot(rotmat, tmp)
                s.coors[i, 0:3] = tmp.copy()

    def rotate2(self, angle, axis):
        """rotates the structure coordinates with angle degrees
           around the axis specified by axis vector"""
        assert size(axis) == 3, 'Axis vector must have 3 elements'
        norm = sqrt(axis[0]**2+axis[1]**2+axis[2]**2)
        axis /= norm
        angle = angle/180.0*pi
        el_11 = cos(angle)+axis[0]**2*(1.-cos(angle))
        el_12 = axis[0]*axis[1]*(1.-cos(angle))-axis[2]*sin(angle)
        el_13 = axis[0]*axis[2]*(1.-cos(angle))+axis[1]*sin(angle)
        el_21 = axis[0]*axis[1]*(1.-cos(angle))+axis[2]*sin(angle)
        el_22 = cos(angle)+axis[1]**2*(1.-cos(angle))
        el_23 = axis[2]*axis[1]*(1.-cos(angle))-axis[0]*sin(angle)
        el_31 = axis[0]*axis[2]*(1.-cos(angle))-axis[1]*sin(angle)
        el_32 = axis[2]*axis[1]*(1.-cos(angle))+axis[0]*sin(angle)
        el_33 = cos(angle)+axis[2]**2*(1.-cos(angle))
        rotmat = array([[el_11, el_12, el_13],
                        [el_21, el_22, el_23],
                        [el_31, el_32, el_33]])

        for s in self.Strokes:
            npoints = size(s.coors, 0)
            for i in range(npoints):
                tmp = s.coors[i, 0:3].copy()
                tmp = dot(rotmat, tmp)
                s.coors[i, 0:3] = tmp.copy()

    def plot(self, plotmode=1, markerscalef=0.1):
        """plots the structure coordinates if mode=0 (default) 
           if plotmode=1 strokes are colored with different colors"""

        assert any(plotmode == array([0, 1])), 'plot mode must be 0 or 1'

        # get the min and max of XYZ coordinates:
        clims = array([[inf, inf, inf], [-inf, -inf, -inf]])
        for s in self.Strokes:
            mins = s.coors.min(axis=0)
            maxs = s.coors.max(axis=0)
            clims[0, :] = where(clims[0, :] > mins, mins, clims[0, :])
            clims[1, :] = where(clims[1, :] < maxs, maxs, clims[1, :])

        scf = markerscalef
        fh = mlab.figure(size=(750, 750))

        # fh.scene.disable_render=True
        cmap_lines = cm.hsv(linspace(0, 1, len(self.Strokes)+1))[1:, :]
        cmap_markers = deepcopy(cmap_lines)
        # alter the marker colors:
        for c in cmap_markers:
            hlsc = array(colorsys.rgb_to_hls(c[0], c[1], c[2]))
            # change the hls lightness:
            if hlsc[1] > 0.4:
                hlsc[1] -= 0.25
            else:
                hlsc[1] += 0.25
            # convert back to rgb:
            rgbc = colorsys.hls_to_rgb(hlsc[0], hlsc[1], hlsc[2])
            c[0:3] = rgbc

        ci = 0
        for s in self.Strokes:
            if plotmode == 1:
                linecol = tuple(cmap_lines[ci, 0:3])
                markercol = tuple(cmap_markers[ci, 0:3])
            else:
                linecol = (0, 0, 1)
                markercol = (0, 1, 0)

            mlab.points3d(s.coors[:, 0], s.coors[:, 1], s.coors[:, 2],
                          color=markercol, scale_factor=scf, mode='cube')
            mlab.plot3d(s.coors[:, 0], s.coors[:, 1], s.coors[:, 2],
                        color=linecol, tube_radius=None, line_width=2.)
            ci += 1

        # eng=mlab.get_engine()
        # for li in eng.current_scene.scene.light_manager.lights:
         #   li.intensity=1.0
        # this is needed otherwise the axes will only cover the stroke plotted last
        mlab.points3d(clims[:, 0], clims[:, 1], clims[:, 2],
                      color=(0, 0, 0), scale_factor=scf)
        mlab.axes(nb_labels=5, ranges=clims.T.flatten())
        #mlab.view(42, 73, 104, [79,  75,  76])
        mlab.show()

        # fh.scene.disable_render=False

    def plot2D(self, ax=[0, 1], col='b', plotmode=0):
        """plots the microstructure in 2D with axes ax[0] and ax[1]           
           if plotmode=1 strokes are colored with different colors"""
        assert any(plotmode == array([0, 1])), 'plot mode must be 0 or 1'

        cmap = cm.prism(linspace(0, 1, len(self.Strokes)))
        ci = 0

        for s in self.Strokes:
            col = tuple(cmap[ci, 0:3]) if plotmode == 1 else col
            pyplot.plot(s.coors[:, ax[0]], s.coors[:, ax[1]], '.-', color=col)
            # pyplot.hold(1)
            ci += 1
        pyplot.axis("equal")
        pyplot.grid()
        pyplot.show()


def slabStr(xyz_sizes, scanorder, stepsizes, speed):
    """Creates a cuboid slab with sizes along X-Y-Z defined in xyz_sizes.
       A single line is scanned back and forth along the axis defined 
       by the first value of scanorder (e.g. for [2,0,1] 2==Z-axis),
       this scanline is repeated and shifted along the axes determined by 
       the 2nd and 3rd values of scanorder in steps defined in stepsizes.
       stepsizes are recalculated to get an integer number of scanlines
       within the desired size of the structure"""

    # create a line in the first scan direction:
    lineE = zeros((2, 3))
    lineE[1, scanorder[0]] = xyz_sizes[scanorder[0]]

    # recalculate the stepsizes to match xyz_sizes:
    Nsteps = zeros(2, int)
    for i in arange(2):
        s = xyz_sizes[scanorder[i+1]]
        step = stepsizes[i]
        Nsteps[i] = round(s/step)
        stepsizes[i] = s/Nsteps[i]

    Nsteps += 1
    # replicate the line along the 2nd and 3rd scan direction:
    slab = lineReplicate3D(lineE, Nsteps, stepsizes, scanorder[1:3])

    slab[-1, 3:5] = [0, speed]

    return MicroStr(slab)


def sphereStr(cx, cy, cz, sphereR, speed, xyres, sphZratio, Nscandir, shellspacing=0.5):

    circres = 0.25  # distance of points along the circles
    circangres = 30.  # used if r is small, thus circres is too big
    Rstep = shellspacing  # distance between the shells of the sphere

    # the starting Radius:
    rstart = remainder(sphereR, Rstep)
    if (rstart == 0):
        rstart = Rstep

    sphx = array([])
    sphy = array([])
    sphz = array([])
    sh = array([])

    dir = 1
    # loop of shell radius:
    for r in linspace(rstart, sphereR, (sphereR-rstart)/Rstep+1):

        phires = xyres/r  # angular step to make the circles of shell with r
        phistart = remainder(pi, phires)/2

        phirange = arange(phistart, pi, phires)
        if (dir == -1):
            phirange = flipud(phirange)
        dir *= -1
        # make the shell with r
        for phi in phirange:

            # radius of circle at angle phi:
            R = sin(phi)*r

            if(R < 0.150):
                continue

            z = r-(cos(phi)*r)-r
            # number of points of the current circle:
            pnumber = round((2*pi*R)/circres)

            if (360/pnumber > circangres):  # refine resolution
                pnumber = 360/circangres

            pnumber = int(pnumber)
            # for i in arange(1,pnumber+1):
            angl = 360/pnumber*arange(1, pnumber+1)
            sphx = append(sphx, R*cos(angl/180*pi))
            sphy = append(sphy, R*sin(angl/180*pi))
            sphz = append(sphz, zeros(pnumber)+z)
            sh = append(sh, ones(pnumber))
        sh[-1] = 0

    sphx = expand_dims(sphx, 1)
    sphy = expand_dims(sphy, 1)
    sphz = expand_dims(sphz, 1)
    sh = expand_dims(sh, 1)
    spv = zeros_like(sh)
    spv[sh == 0] = speed

    sph = concatenate((sphx, sphy, sphz, sh, spv), axis=1)

    # swap x-z coordinates
    tmp = sph[:, 2].copy()
    sph[:, 2] = sph[:, 0]
    sph[:, 0] = tmp

    sphStr = MicroStr(sph)

    # Nscandir directions
    phi = 180./Nscandir
    t = deepcopy(sphStr)
    for i in range(1, Nscandir):
        t.rotate(phi, 'z')
        sphStr.addStr(t)

    sphStr.scale([1, 1, sphZratio])
    sphStr.shift([cx, cy, cz])

    return sphStr


def diskStr(cx, cy, cz, R, speed, xyres, zres, Nlayers, phires=None, Rmin=None, ipstep=None, scanmode=0):
    """use scanmode==1 to avoid shape deformation with large phires"""

    # the starting Radius:
    if Rmin is None:
        rstart = remainder(R, 0.5)
    else:
        rstart = Rmin

    if (rstart == 0):
        rstart = 2*xyres

    if phires is None:
        pr_mode = 0  # dynamic angular resolution
    else:
        pr_mode = 1

    disk = array([[], [], []]).transpose()
    for r in linspace(rstart, R, (R-rstart)/xyres+1):

        if (pr_mode == 1):
            tres = 2*pi*r/(360./phires)
        else:
            tres = xyres
        if tres < xyres:
            tres = xyres

        phiN = round(2*pi/(tres/r))+1
        if (scanmode == 0):
            phi = linspace(0, 2*pi, phiN)[:-1]
        else:
            phi = linspace(0, 2*pi, phiN)

        x = cos(phi)*r
        y = sin(phi)*r
        z = x*0

        if (scanmode == 1):
            # make the segment between the last 2 points xyres shorter:
            dr = array([diff(x[-2:]), diff(y[-2:])])
            rabs = sqrt(sum(dr*dr))
            rn = dr/rabs
            dr_new = rn*(rabs-xyres)
            x[-1] = x[-2]+dr_new[0]
            y[-1] = y[-2]+dr_new[1]

        temp = array([x, y, z]).transpose()
        disk = append(disk, temp, 0)

    temp = disk.copy()
    for i in arange(Nlayers-1):
        temp[:, 2] += zres
        temp = flipud(temp)
        disk = append(disk, temp, 0)

    # add the 4th and 5th columns for shutter and speed:
    disk = append(disk, ones((len(disk), 1)), 1)
    disk = append(disk, zeros((len(disk), 1)), 1)
    disk[-1, 3:5] = [0, speed]

    disk[:, 0] += cx
    disk[:, 1] += cy
    disk[:, 2] += cz

    return MicroStr(disk, ipstep=ipstep)


def helixStr(cx, cy, cz, R, Nturns, pitch, speed, xyres, zres, Nlayers, phires=None, Rmin=None, ipstep=None, scanmode=0):
    """use scanmode==1 to avoid shape deformation with large phires"""

    # the starting Radius:
    if Rmin is None:
        rstart = remainder(R, 0.5)
    else:
        rstart = Rmin

    if (rstart == 0):
        rstart = 2*xyres

    if phires is None:
        pr_mode = 0  # dynamic angular resolution
    else:
        pr_mode = 1

    helix = array([[], [], []]).transpose()
    j = 0
    for r in linspace(rstart, R, (R-rstart)/xyres+1):

        if (pr_mode == 1):
            tres = 2*pi*r/(360./phires)
        else:
            tres = xyres
        if tres < xyres:
            tres = xyres

        phiN = round((2*pi*Nturns)/(tres/r))+1
        if (scanmode == 0):
            phi = linspace(0, (2*pi*Nturns), phiN)[:-1]
            z = linspace(0, (Nturns*pitch), phiN)[:-1]
        else:
            phi = linspace(0, (2*pi*Nturns), phiN)
            z = linspace(0, (Nturns*pitch), phiN)

        x = cos(phi)*r
        y = sin(phi)*r

        if (scanmode == 1):
            # make the segment between the last 2 points xyres shorter:
            dr = array([diff(x[-2:]), diff(y[-2:])])
            rabs = sqrt(sum(dr*dr))
            rn = dr/rabs
            dr_new = rn*(rabs-xyres)
            x[-1] = x[-2]+dr_new[0]
            y[-1] = y[-2]+dr_new[1]

        temp = array([x, y, z]).transpose()
        if mod(j, 2) == 1:
            temp = flipud(temp)
        helix = append(helix, temp, 0)
        j += 1

    temp = helix.copy()
    for i in arange(Nlayers-1):
        temp[:, 2] += zres
        temp = flipud(temp)
        helix = append(helix, temp, 0)

    # add the 4th and 5th columns for shutter and speed:
    helix = append(helix, ones((len(helix), 1)), 1)
    helix = append(helix, zeros((len(helix), 1)), 1)
    helix[-1, 3:5] = [0, speed]

    helix[:, 0] += cx
    helix[:, 1] += cy
    helix[:, 2] += cz

    return MicroStr(helix, ipstep=ipstep)


def lineReplicate2D(lineE, N, d, direction):
    """replicates a polygon line defined by lineE 
       N times with d shifts along direction
       lineE must be at least 2x3 dimension"""
    lineE = lineE.copy()
    repline = lineE.copy()
    for i in arange(1, N):
        lineE[:, direction] += d
        lineE = flipud(lineE)
        repline = append(repline, lineE, axis=0)

    # add the 4th and 5th columns for shutter and speed:
    repline = append(repline, ones((len(repline), 1)), 1)
    repline = append(repline, zeros((len(repline), 1)), 1)

    return repline


def lineReplicate3D(lineE, Ns, ds, dirs):
    """replicates a polygon line defined by lineE 
       in the 2 directions defined by dirs
       with 2 number of replications defined by Ns
       with the 2 spacing defined in ds
       lineE must be at least 2x3 dimension"""

    lineE = lineE.copy()
    levelE = lineE.copy()

    for i in arange(1, Ns[0]):
        lineE[:, dirs[0]] += ds[0]
        lineE = flipud(lineE)
        levelE = append(levelE, lineE, axis=0)

    repline = levelE.copy()

    for i in arange(1, Ns[1]):
        levelE[:, dirs[1]] += ds[1]
        levelE = flipud(levelE)
        repline = append(repline, levelE, axis=0)

    # add the 4th and 5th columns for shutter and speed:
    repline = append(repline, ones((len(repline), 1)), 1)
    repline = append(repline, zeros((len(repline), 1)), 1)

    return repline


def letter0(lsize, cx, cy, speed):

    R = lsize/2

    degres = 15
    phistart = 0
    phiend = 360
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N
    zeroframe = zeros((N+1, 3))
    phivec = arange(phistart, phiend+degres, degres)/180.*pi
    zeroframe[:, 0] = cos(phivec)*R*0.66
    zeroframe[:, 1] = sin(phivec)*R

    # add the 4th and 5th columns for shutter and speed:
    zeroframe = append(zeroframe, ones((len(zeroframe), 1)), 1)
    zeroframe = append(zeroframe, zeros((len(zeroframe), 1)), 1)
    zeroframe[-1, 3:5] = array([0, speed])

    zero = zeroframe
    zero[:, 0:2] = zero[:, 0:2]-zero[:, 0:2].min(axis=0)+array([cx, cy])

    return MicroStr(zero)


def letter1(lsize, cx, cy, speed):

    height = lsize
    width = lsize
    one = zeros((3, 5))

    one[0, 1] = height/3.*2
    one[1, 0] = width/3.
    one[1, 1] = height
    one[2, 0] = width/3.

    one[:, 3] = 1
    one[-1, 3:5] = array([0, speed])

    one[:, 0:2] = one[:, 0:2]+array([cx, cy])

    return MicroStr(one)


def letter2(lsize, cx, cy, speed):

    width = lsize*0.66
    height = lsize
    R = width/2.

    degres = 15
    phistart = 165
    phiend = 15
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N

    two = zeros((N+3, 5))
    phivec = arange(phistart, phiend+degres, degres)/180.*pi
    two[:N+1, 0] = cos(phivec)*R+R
    two[:N+1, 1] = sin(phivec)*R+height-R

    two[N+2, 0] = width
    two[:, 3] = 1
    two[-1, 3:5] = array([0, speed])

    two[:, 0:2] += array([cx, cy])

    return MicroStr(two)


def letter3(lsize, cx, cy, speed):

    width = lsize
    height = lsize
    R = width/4.

    degres = 15.
    phistart = 165.
    phiend = -70.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N
    three = zeros((N+1, 5))
    phivec = arange(phistart, phiend+degres, degres)/180.*pi
    three[:, 0] = cos(phivec)*R+R
    three[:, 1] = sin(phivec)*R+height-R

    three[:, 3] = 1

    np = array([[three[-1, 0]-0.150, height/2., 0, 1, 0],
                [width/4., height/2., 0, 0, speed]])

    three = append(three, np, axis=0)

    m = size(three, 0)

    np = three[:m+1, :].copy()
    np[:, 1] = abs(three[0:m+1, 1]-height)
    np[:, 0:2] = flipud(np[:, 0:2])

    three = append(three, np, axis=0)
    three[-1, 3:5] = array([0, speed])
    three[:, 0:2] += array([cx, cy])

    return MicroStr(three)


def letter4(lsize, cx, cy, speed):

    height = lsize
    width = lsize
    four = zeros((4, 5))

    four[0, 0] = width/2.+width/10.
    four[0, 1] = height/5.*2

    four[1, 1] = height/5.*2

    four[2, 0] = width/2.
    four[2, 1] = height
    four[3, 0] = width/2.

    four[:, 3] = 1
    four[-1, 3:5] = array([0, speed])

    four[:, 0:2] += array([cx, cy])

    return MicroStr(four)


def letter5(lsize, cx, cy, speed):

    width = lsize
    height = lsize
    R = width*0.33

    degres = 15
    phistart = 110
    phiend = -135
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N
    phivec = arange(phistart, phiend+degres, degres)/180.*pi

    five = zeros((N+2, 5))
    five[1:, 0] = cos(phivec)*R+R
    five[1:, 1] = sin(phivec)*R+R

    five[0, 0] = 2*R
    five[0:2, 1] = height
    five[1, 0] = five[2, 0]

    five[:, 3] = 1
    five[-1, 3:5] = array([0, speed])

    five[:, 0:2] += array([cx, cy])

    return MicroStr(five)


def letter6(lsize, cx, cy, speed):

    width = lsize
    height = lsize
    R = width*0.33

    degres = 15.
    phistart = 0.
    phiend = 360.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N

    six = zeros((N+1, 5))
    phivec = (arange(phistart, phiend+degres, degres)+180)/180.*pi
    six[:, 0] = cos(phivec)*R+R
    six[:, 1] = sin(phivec)*R+R

    six[:, 3] = 1
    six[-1, 3:5] = array([0, speed])

    degres = 15.
    phistart = 175.
    phiend = 80.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N
    phivec = (arange(phistart, phiend+degres, degres))/180.*pi
    np = zeros((N+1, 5))
    np[:, 0] = cos(phivec)*1.5*R+1.5*R
    np[:, 1] = sin(phivec)*2*R+1.2*R

    six = append(six, np, axis=0)
    six[:, 3] = 1
    six[-1, 3:5] = array([0, speed])

    six[:, 0:2] += array([cx, cy])

    return MicroStr(six)


def letter7(lsize, cx, cy, speed):

    height = lsize
    width = lsize/2.
    seven = zeros((3, 5))

    seven[0, 1] = height
    seven[1, 0] = width
    seven[1, 1] = height

    seven[:, 3] = 1
    seven[-1, 3:5] = array([0, speed])

    seven[:, 0:2] += array([cx, cy])

    return MicroStr(seven)


def letter8(lsize, cx, cy, speed):

    width = lsize
    height = lsize
    R = width/4.

    degres = 15.
    phistart = -70.
    phiend = 250.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N
    eight = zeros((N+1, 5))
    phivec = (arange(phistart, phiend+degres, degres))/180.*pi
    eight[:, 0] = cos(phivec)*R+R
    eight[:, 1] = sin(phivec)*R+height-R

    eight[:, 3] = 1
    eight[-1, 3:5] = array([0, speed])

    degres = 15.
    phistart = 90.
    phiend = -270.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N

    np = zeros((N+1, 5))
    phivec = (arange(phistart, phiend+degres, degres))/180.*pi
    np[:, 0] = cos(phivec)*R+R
    np[:, 1] = sin(phivec)*R+R
    eight = append(eight, np, axis=0)

    eight[:, 3] = 1
    eight[-1, 3:5] = array([0, speed])

    eight[:, 0:2] += array([cx, cy])

    return MicroStr(eight)


def letter9(lsize, cx, cy, speed):

    width = lsize
    height = lsize
    R = width*0.33

    degres = 15.
    phistart = 0.
    phiend = 360.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N

    nine = zeros((N+1, 5))
    phivec = (arange(phistart, phiend+degres, degres)+180)/180.*pi
    nine[:, 0] = cos(phivec)*R+R
    nine[:, 1] = sin(phivec)*R+R

    nine[:, 3] = 1
    nine[-1, 3:5] = array([0, speed])

    degres = 15.
    phistart = 175.
    phiend = 80.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N
    phivec = (arange(phistart, phiend+degres, degres))/180.*pi
    np = zeros((N+1, 5))
    np[:, 0] = cos(phivec)*1.5*R+1.5*R
    np[:, 1] = sin(phivec)*2*R+1.2*R

    nine = append(nine, np, axis=0)
    nine[:, 3] = 1
    nine[-1, 3:5] = array([0, speed])

    nine[:, 0:2] = -nine[:, 0:2]
    nine[:, 0:2] -= nine[:, 0:2].min(axis=0)

    nine[:, 0:2] += array([cx, cy])

    return MicroStr(nine)


def letterA(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize
    A = zeros((5, 5))
    A[:, 3] = 1

    A[1, 0] = width/2.
    A[1, 1] = height
    A[2, 0] = width

    A[2, 3:5] = array([0, speed])

    A[3, 1] = height/3.
    A[3, 0] = (width/2.)/height*A[3, 1]+0.150
    A[4, 1] = height/3.
    A[4, 0] = width-(width/2.)/height*A[3, 1]-0.150

    A[-1, 3:5] = array([0, speed])
    A[:, 0:2] += array([cx, cy])

    return MicroStr(A)


def half_circle(R, N, cx, cy):

    circ = zeros((N, 5))

    angres = pi/(N+1)
    an = arange(N, 0, -1)
    circ[:, 0] = cos(an*angres-pi/2)*R
    circ[:, 1] = sin(an*angres-pi/2)*R

    circ[:, 0:2] += array([cx, cy])
    circ[:, 3] = 1
    return circ


def letterB(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize
    B = zeros((3, 5))

    B[1, 1] = height
    B[2, 0] = width/2.
    B[2, 1] = height
    B[:, 3] = 1
    circN = 10
    circ = half_circle(height/4., circN, 0.5*width, 0.75*height)
    B = append(B, circ, axis=0)

    np = array([[width/2., height/2., 0, 1, 0],
                [0.150, height/2., 0, 1, 0]])
    B = append(B, np, axis=0)
    B[-1, 3:5] = array([0, speed])

    circ = half_circle(height/4., circN, 0.5*width, 0.25*height)

    B = append(B, circ, axis=0)

    np = array([[width/2., 0, 0, 1, 0],
                [0.150, 0, 0, 1, 0]])
    B = append(B, np, axis=0)

    B[-1, 3:5] = array([0, speed])

    B[:, 0:2] += array([cx, cy])

    return MicroStr(B)


def letterC(lsize, cx, cy, speed):

    R = lsize/2.
    degres = 2.*arcsin(0.3/R)/pi*180

    phistart = -50.
    phiend = -310.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N
    C = zeros((N+1, 5))

    phivec = (arange(phistart, phiend+degres, degres))/180.*pi
    C[:, 0] = cos(phivec)*R
    C[:, 1] = sin(phivec)*R

    C[:, 3] = 1
    C[-1, 3:5] = array([0, speed])
    C[:, 0:2] -= C[:, 0:2].min(axis=0)
    C[:, 0:2] += array([cx, cy])

    return MicroStr(C)


def letterD(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize
    D = zeros((3, 5))

    D[1, 1] = height
    D[2, 0] = width/2.
    D[2, 1] = height
    circN = 10
    circ = half_circle(height/2., circN, 0.5*width, 0.5*height)
    D = append(D, circ, axis=0)

    np = array([[width/2., 0, 0, 1, 0],
                [0.150, 0, 0, 1, 0]])
    D = append(D, np, axis=0)

    D[:, 3] = 1
    D[-1, 3:5] = array([0, speed])
    D[:, 0:2] += array([cx, cy])

    return MicroStr(D)


def letterE(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize
    E = zeros((6, 5))
    E[:, 3] = 1

    E[0, 0] = width
    E[2, 1] = height
    E[3, 0] = width
    E[3, 1] = height
    E[3, 3:5] = array([0, speed])

    E[4, 0] = 0.150
    E[4, 1] = height/2.
    E[5, 0] = width*0.75
    E[5, 1] = height/2.

    E[-1, 3:5] = array([0, speed])
    E[:, 0:2] += array([cx, cy])

    return MicroStr(E)


def letterF(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize
    F = zeros((5, 5))
    F[:, 3] = 1

    F[1, 1] = height
    F[2, 0] = width
    F[2, 1] = height
    F[2, 3:5] = array([0, speed])

    F[3, 0] = 0.150
    F[3, 1] = height/2.
    F[4, 0] = width*0.75
    F[4, 1] = height/2.

    F[-1, 3:5] = array([0, speed])
    F[:, 0:2] += array([cx, cy])

    return MicroStr(F)


def letterG(lsize, cx, cy, speed):

    R = lsize/2.
    degres = 2.*arcsin(0.3/R)/pi*180

    phistart = -degres
    phiend = -300.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N

    G = zeros((N+3, 5))
    G[1, 0] = R

    phivec = (arange(phistart, phiend+degres, degres))/180.*pi
    G[2:N+3, 0] = cos(phivec)*R
    G[2:N+3, 1] = sin(phivec)*R

    G[:, 3] = 1
    G[-1, 3:5] = array([0, speed])
    G[:, 0:2] = G[:, 0:2]-G[:, 0:2].min(axis=0)+array([cx, cy])

    return MicroStr(G)


def letterH(lsize, cx, cy, speed):
    width = lsize*0.75
    height = lsize
    H = zeros((6, 5))
    H[:, 3] = 1

    H[1, 1] = height
    H[1, 3:5] = array([0, speed])

    H[2:4, 0] = width
    H[3, 1] = height
    H[3, 3:5] = array([0, speed])

    H[4, 0] = 0.150
    H[4, 1] = height/2.
    H[5, 0] = width-0.150
    H[5, 1] = height/2.

    H[-1, 3:5] = array([0, speed])
    H[:, 0:2] += array([cx, cy])

    return MicroStr(H)


def letterI(lsize, cx, cy, speed):

    I = zeros((2, 5))
    I[0, 3] = 1
    I[1, 1] = lsize
    I[1, 3:5] = array([0, speed])
    I[:, 0:2] += array([cx, cy])

    return MicroStr(I)


def letterJ(lsize, cx, cy, speed):

    width = lsize
    height = lsize

    J = zeros((2, 5))
    J[0, 0] = width
    J[0, 1] = height
    J[1, 0] = width
    J[1, 1] = height/3.

    circN = 10
    circ = half_circle(width/3., circN, width/3.*2., height/3.)
    circ = circ[-4:, :]
    J = append(J, circ, axis=0)
    J = append(J, array([J[-1, :]]), axis=0)
    J[-1, 0] -= width/6.
    J[:, 3] = 1
    J[-1, 3:5] = array([0, speed])

    J[:, 0] = J[:, 0]-min(J[:, 0])
    J[:, 0] = J[:, 0]-max(J[:, 0])/2.
    J[:, 0:2] -= J[:, 0:2].min(axis=0)
    J[:, 0:2] += array([cx, cy])

    return MicroStr(J)


def letterK(lsize, cx, cy, speed):

    width = lsize
    height = lsize
    K = zeros((6, 5))
    K[:, 3] = 1
    K[1, 1] = height
    K[1, 3:5] = array([0, speed])

    K[2, 0] = 0.150
    K[2, 1] = height/2.+0.150
    K[3, 0] = width/4.*3.
    K[3, 1] = height

    K[3, 3:5] = array([0, speed])

    K[4, 0] = 0.150
    K[4, 1] = height/2.-0.150
    K[5, 0] = width/4.*3.

    K[-1, 3:5] = array([0, speed])

    K[:, 0:2] += array([cx, cy])

    return MicroStr(K)


def letterL(lsize, cx, cy, speed):

    width = lsize
    height = lsize
    L = zeros((3, 5))
    L[:, 3] = 1
    L[0, 0] = width*0.75
    L[2, 1] = height
    L[-1, 3:5] = array([0, speed])
    L[:, 0:2] += array([cx, cy])

    return MicroStr(L)


def letterM(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize
    M = zeros((5, 5))

    M[1, 1] = height
    M[2, 0] = width/2.
    M[2, 1] = height/3.
    M[3, 0] = width
    M[3, 1] = height
    M[4, 0] = width

    M[:, 3] = 1
    M[-1, 3:5] = array([0, speed])
    M[:, 0:2] += array([cx, cy])

    return MicroStr(M)


def letterN(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize
    N = zeros((4, 5))

    N[1, 1] = height
    N[2, 0] = width
    N[3, 0] = width
    N[3, 1] = height

    N[:, 3] = 1
    N[-1, 3:5] = array([0, speed])
    N[:, 0:2] += array([cx, cy])

    return MicroStr(N)


def letterO(lsize, cx, cy, speed):

    R = lsize/2.
    degres = 20.
    phistart = 0.
    phiend = 360.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    degres = (phiend-phistart)/N
    O = zeros((N+1, 5))

    phivec = (arange(phistart, phiend+degres, degres))/180.*pi
    O[:, 0] = cos(phivec)*R*0.75
    O[:, 1] = sin(phivec)*R

    O[:, 3] = 1
    O[-1, 3:5] = array([0, speed])
    O[:, 0:2] = O[:, 0:2]-O[:, 0:2].min(axis=0)+array([cx, cy])

    return MicroStr(O)


def letterP(lsize, cx, cy, speed):

    width = lsize/2.
    height = lsize
    P = zeros((3, 5))

    P[1, 1] = height
    P[2, 0] = width/2.
    P[2, 1] = height
    circN = 10
    circ = half_circle(height/4., circN, 0.5*width, 0.75*height)
    P = append(P, circ, axis=0)

    np = array([[width/2., height/2., 0, 1, 0],
                [0.150, height/2., 0, 1, 0]])
    P = append(P, np, axis=0)

    P[:, 3] = 1
    P[-1, 3:5] = array([0, speed])
    P[:, 0:2] += array([cx, cy])

    return MicroStr(P)


def letterQ(lsize, cx, cy, speed):

    R = lsize/2.
    degres = 2.*arcsin(0.3/R)/pi*180.

    phistart = 0.
    phiend = -360.
    N = abs((phiend-phistart)/degres)
    N = round(N)

    Q = zeros((N, 5))
    phivec = (linspace(phistart, phiend, N))/180.*pi
    Q[:, 0] = cos(phivec)*R*0.75
    Q[:, 1] = sin(phivec)*R

    Q[:, 3] = 1
    Q[-1, 3:5] = array([0, speed])

    Qline = zeros((2, 5))
    Qline[0, 0] = lsize/8.
    Qline[0, 1] = -lsize/8.
    Qline[1, 0] = lsize/2.
    Qline[1, 1] = -lsize/2.

    Qline[:, 3] = 1
    Qline[-1, 3:5] = array([0, speed])

    Q = append(Q, Qline, axis=0)
    Q[:, 0:2] -= Q[:, 0:2].min(axis=0)
    Q[:, 0:2] += array([cx, cy])

    return MicroStr(Q)


def letterR(lsize, cx, cy, speed):

    width = lsize/2
    height = lsize
    R = zeros((3, 5))

    R[1, 1] = height
    R[2, 0] = width/2.
    R[2, 1] = height

    circN = 10
    circ = half_circle(height/4., circN, 0.5*width, 0.75*height)
    R = append(R, circ, axis=0)
    R[:, 3] = 1
    np = array([[width/2., height/2., 0, 1, 0],
                [0.150, height/2., 0, 0, speed],
                [width/2., height/2., 0, 1, 0],
                [width, 0, 0, 0, speed]])
    R = append(R, np, axis=0)

    R[:, 0:2] += array([cx, cy])

    return MicroStr(R)


def letterS(lsize, cx, cy, speed):

    R = lsize/4.

    degres = 20.
    phistart = 50.
    phiend = 270.
    N = abs((phiend-phistart)/degres)
    N = round(N)
    phivec = (linspace(phistart, phiend, N))/180.*pi
    S = zeros((len(phivec), 5))
    S[:, 0] = cos(phivec)*R
    S[:, 1] = sin(phivec)*R+lsize*0.75

    phistart = phiend-180.
    phiend = -140.
    N = abs((phiend-phistart)/degres)
    N = round(N)

    phivec = (linspace(phistart, phiend, N))/180.*pi
    phivec = phivec[1:]
    np = zeros((len(phivec), 5))
    np[:, 0] = cos(phivec)*R
    np[:, 1] = sin(phivec)*R+0.25*lsize

    S = append(S, np, axis=0)

    S[:, 3] = 1
    S[-1, 3:5] = array([0, speed])
    S[:, 0:2] -= S[:, 0:2].min(axis=0)
    S[:, 0:2] += array([cx, cy])

    return MicroStr(S)


def letterT(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize

    T = array([[width/2., 0, 0, 1, 0],
               [width/2., height-0.150, 0, 0, speed],
               [0, height, 0, 1, 0],
               [width, height, 0, 0, speed]])
    T[:, 0:2] += array([cx, cy])

    return MicroStr(T)


def letterU(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize

    U = zeros((2, 5))
    U[0, 0] = width
    U[0, 1] = height
    U[1, 0] = width
    U[1, 1] = height/2.

    circN = 10
    R = width/2.
    circ = zeros((circN, 5))
    angres = pi/(circN+1)
    an = arange(circN, 0, -1)
    circ[:, 0] = cos(an*angres-pi)*R
    circ[:, 1] = sin(an*angres-pi)*R
    circ[:, 0:2] += array([width/2., width/2.])
    circ[:, 3] = 1

    U = append(U, circ, axis=0)

    np = array([[0, height/2., 0, 1, 0],
                [0, height, 0, 1, 0]])

    U = append(U, np, axis=0)
    U[:, 3] = 1
    U[-1, 3:5] = array([0, speed])
    U[:, 0:2] += array([cx, cy])

    return MicroStr(U)


def letterV(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize
    V = zeros((3, 5))

    V[0, 1] = height
    V[1, 0] = width/2.
    V[2, 0] = width
    V[2, 1] = height

    V[:, 3] = 1
    V[-1, 3:5] = array([0, speed])
    V[:, 0:2] += array([cx, cy])

    return MicroStr(V)


def letterW(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize
    W = zeros((5, 5))

    W[0, 1] = height
    W[1, 0] = width/4.
    W[2, 0] = width/2.
    W[2, 1] = height/3.*2
    W[3, 0] = width/4.*3
    W[4, 0] = width
    W[4, 1] = height

    W[:, 3] = 1
    W[-1, 3:5] = array([0, speed])
    W[:, 0:2] += array([cx, cy])

    return MicroStr(W)


def letterX(lsize, cx, cy, speed):

    width = lsize/4.*3
    height = lsize

    X = array([[0, 0, 0, 1, 0],
               [width, height, 0, 0, speed],
               [0, height, 0, 1, 0],
               [width, 0, 0, 0, speed]])

    X[:, 0:2] += array([cx, cy])

    return MicroStr(X)


def letterY(lsize, cx, cy, speed):

    width = lsize/4.*3
    height = lsize
    Y = zeros((5, 5))
    Y[:, 3] = 1

    Y[0, 0] = width/2.
    Y[1, 0] = width/2.
    Y[1, 1] = height/2.
    Y[2, 1] = height

    Y[2, 3:5] = array([0, speed])

    Y[3, 0] = width/2.
    Y[3, 1] = height/2.
    Y[4, 0] = width
    Y[4, 1] = height
    Y[4, 3:5] = array([0, speed])
    Y[:, 0:2] += array([cx, cy])

    return MicroStr(Y)


def letterZ(lsize, cx, cy, speed):

    width = lsize*0.75
    height = lsize
    Z = array([[0, height, 0, 1, 0],
               [width, height, 0, 1, 0],
               [0, 0, 0, 1, 0],
               [width, 0, 0, 0, speed]])
    Z[:, 0:2] += array([cx, cy])

    return MicroStr(Z)


def letterPoint(lsize, cx, cy, speed):

    Point = diskStr(cx, cy, 0, lsize*0.1, speed, lsize*0.1, lsize*0.1, 1)

    return Point


def textStr(txt, textsize, speed):

    ln = len(txt)
    space = 0.5*textsize
    ypos = 0
    xpos = 0
    for i in range(ln):

        l = txt[i].upper()
        if l == '.':
            l = 'Point'
        if (l == '"'):
            ypos = ypos-textsize-space
            xpos = 0
            continue
        if (l == ' '):
            xpos += space
            continue
        fstring = "lstr=letter" + l + "("
        fstring += str(textsize) + "," + str(xpos) + "," + \
            str(ypos) + "," + str(speed) + ")"
        exec(fstring)

        if i == 0:
            txstr = deepcopy(lstr)
        else:
            txstr.addStr(lstr)

        # get the min and max of XYZ coordinates:
        clims = array([[inf, inf, inf], [-inf, -inf, -inf]])
        for s in lstr.Strokes:
            mins = s.coors.min(axis=0)
            maxs = s.coors.max(axis=0)
            clims[0, :] = where(clims[0, :] > mins, mins, clims[0, :])
            clims[1, :] = where(clims[1, :] < maxs, maxs, clims[1, :])
        xpos += space+clims[1, 0]-clims[0, 0]

    txstr.rotate(180, 'z')
    for s in txstr.Strokes:
        s.coors[:, 0] *= -1

#    txstr.strdata[:,0]*=-1
#    txstr.strdata[:,0]-=min(txstr.strdata[:,0])

    return txstr
