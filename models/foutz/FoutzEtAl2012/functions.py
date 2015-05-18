def find_section_coordinates(sec):
    """ Determine xyz coordinates for a given section """
    from numpy import array
    from neuron import h
    x=[];y=[];z=[];
    n=h.n3d(sec=sec).__int__()
    for ii in xrange(n):
        x.append(h.x3d(ii,sec=sec))
        y.append(h.y3d(ii,sec=sec))
        z.append(h.z3d(ii,sec=sec))
    return array([array(x),array(y),array(z)])
def find_mean_section_coordinates(sec):
    """ Determine average coordinate for a given section """
    from neuron import h
    from numpy import array, mean
    x=[];y=[];z=[];
    n=h.n3d(sec=sec).__int__()
    for ii in xrange(n):
        x.append(h.x3d(ii,sec=sec))
        y.append(h.y3d(ii,sec=sec))
        z.append(h.z3d(ii,sec=sec))
    return array([mean(x),mean(y),mean(z)])
def rotate_coordinates(X,Y,Z,v,axis='z'):
    #print "NEW RUN: ",v
    if axis=='z':
        from numpy import array, c_, arctan,sin,cos,matrix,ones,pi,sqrt,sign
        X = array(X)
        Y = array(Y)
        Z = array(Z)
        assert X.shape==Y.shape==Z.shape
        v = matrix(v)
        dx = v[1,0]-v[0,0]
        dy = v[1,1]-v[0,1]
        dz = v[1,2]-v[0,2]
        #print dx,dy,dz
        length = sqrt(dx**2 + dy**2 + dz**2)
        xyz = matrix(c_[X.flatten(),
                        Y.flatten(),
                        Z.flatten()])
        # this is where we translate/rotate
        T = matrix([[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [-v[1,0], -v[1,1], -v[1,2], 1]])
        xyz = (c_[xyz,ones(X.flatten().shape)] * T)[:,:3]
        v   = (c_[v,  ones(v[:,0].shape)] * T)[:,:3]

        # Rotate about X axis
        dx = v[1,0]-v[0,0]
        dy = v[1,1]-v[0,1]
        dz = v[1,2]-v[0,2]
        #print dx,dy,dz
        #assert length == sqrt(dx**2 + dy**2 + dz**2)
        #print length, sqrt(dx**2 + dy**2 + dz**2)
        if dz==0:
            theta_x = -1 * sign(dy) * pi/2.0
        else:
            theta_x = -arctan(dy/dz)
        #print 'THETA X: ',theta_x
        #print v

        Rx = matrix([[1,0,0],
                     [0,cos(theta_x),-sin(theta_x)],
                     [0,sin(theta_x),cos(theta_x)]])
        xyz *= Rx
        v   *= Rx

        # Rotate about Y axis
        dx = v[1,0]-v[0,0]
        dy = v[1,1]-v[0,1]
        dz = v[1,2]-v[0,2]
        #print dx,dy,dz
        #assert length == sqrt(dx**2 + dy**2 + dz**2)
        #print length, sqrt(dx**2 + dy**2 + dz**2)
        if dz==0:
            theta_y = sign(dx) * pi/2.0
        else:
            theta_y = arctan(dx/dz)
        Ry = matrix([[cos(theta_y),0,sin(theta_y)],
                     [0,1,0],
                     [-sin(theta_y),0,cos(theta_y)]])
        xyz *= Ry
        v   *= Ry
        if sign((v[1,2]-v[0,2]))==-1: # aligned in the wrong direction
            flip = matrix([[1,0,0],
                           [0,-1,0],
                           [0,0,-1]])
            xyz *= flip
            v   *= flip
        #assert length == sqrt(dx**2 + dy**2 + dz**2)
        #print length, sqrt(dx**2 + dy**2 + dz**2)
        #print "v", v

        # Rotate about Z axis
##        dx = v[1,0]-v[0,0]
##        dy = v[1,1]-v[0,1]
##        if dy != 0:
##            theta_z = -arctan(dx/dy)
##            print "THETA Z: ",theta_z
##            print "v", v
##            Rz = matrix([[cos(theta_z),-sin(theta_z),0],
##                         [sin(theta_z),cos(theta_z),0],
##                         [0,0,1]])
##            xyz *= Rz
##            v   *= Rz

        xyz = array(xyz)
        v   = array(v)
        return xyz[:,0].reshape(X.shape),xyz[:,1].reshape(Y.shape),xyz[:,2].reshape(Z.shape),v
    else:
        raise StandardError,'Not implemented'
def find_cylindrical_coords(X,Y,Z,v):
    " Align v with positive z axis "
    from numpy import sqrt
    X,Y,Z,v = rotate_coordinates(X,Y,Z,v)
    R = sqrt(X**2 + Y**2)
    return R,Z
def rotate_point_around_vector(xyz, uvw,angle,angle_units='degree'):
    """ Rotate the point (x,y,z) around the vector (u,v,w) """
    from numpy import cos, sin, array, pi, deg2rad
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    u = uvw[0]
    v = uvw[1]
    w = uvw[2]
    if angle_units=='degree':
        angle = deg2rad(angle)
    ux = u * x
    uy = u * y
    uz = u * z
    vx = v * x
    vy = v * y
    vz = v * z
    wx = w * x
    wy = w * y
    wz = w * z
    sa = sin(angle)
    ca = cos(angle)
    new_x = u * (ux + vy + wz) + (x * (v * v + w * w) - u * (vy + wz)) * ca + (-wy + vz) * sa
    new_y = v * (ux + vy + wz) + (y * (u * u + w * w) - v * (ux + wz)) * ca + (wx - uz) * sa
    new_z = w * (ux + vy + wz) + (z * (u * u + v * v) - w * (ux + vy)) * ca + (-vx + uy) * sa
    return array([new_x, new_y, new_z])
def point_along_vector(origin, direction, length):
    from numpy import array
    from numpy.linalg import norm
    origin = array(origin).astype(float)
    direction = array(direction).astype(float)
    unit_direction = direction/norm(direction)
    return origin + length * unit_direction
def _float_approx_equal(x, y, tol=1e-18, rel=1e-7):
    ## {{{ http://code.activestate.com/recipes/577124/ (r1)

    if tol is rel is None:
        raise TypeError('cannot specify both absolute and relative errors are None')
    tests = []
    if tol is not None: tests.append(tol)
    if rel is not None: tests.append(rel*abs(x))
    assert tests
    return abs(x - y) <= max(tests)
def approx_equal(x, y, *args, **kwargs):
    """approx_equal(float1, float2[, tol=1e-18, rel=1e-7]) -> True|False
    approx_equal(obj1, obj2[, *args, **kwargs]) -> True|False

    Return True if x and y are approximately equal, otherwise False.

    If x and y are floats, return True if y is within either absolute error
    tol or relative error rel of x. You can disable either the absolute or
    relative check by passing None as tol or rel (but not both).

    For any other objects, x and y are checked in that order for a method
    __approx_equal__, and the result of that is returned as a bool. Any
    optional arguments are passed to the __approx_equal__ method.

    __approx_equal__ can return NotImplemented to signal that it doesn't know
    how to perform that specific comparison, in which case the other object is
    checked instead. If neither object have the method, or both defer by
    returning NotImplemented, approx_equal falls back on the same numeric
    comparison used for floats.

    >>> almost_equal(1.2345678, 1.2345677)
    True
    >>> almost_equal(1.234, 1.235)
    False

    """
    ## {{{ http://code.activestate.com/recipes/577124/ (r1)
    if not (type(x) is type(y) is float):
        # Skip checking for __approx_equal__ in the common case of two floats.
        methodname = '__approx_equal__'
        # Allow the objects to specify what they consider "approximately equal",
        # giving precedence to x. If either object has the appropriate method, we
        # pass on any optional arguments untouched.
        for a,b in ((x, y), (y, x)):
            try:
                method = getattr(a, methodname)
            except AttributeError:
                continue
            else:
                result = method(b, *args, **kwargs)
                if result is NotImplemented:
                    continue
                return bool(result)
    # If we get here without returning, then neither x nor y knows how to do an
    # approximate equal comparison (or are both floats). Fall back to a numeric
    # comparison.
    return _float_approx_equal(x, y, *args, **kwargs)
def make_legend(**extra_args):
    from matplotlib import pyplot
    leg = pyplot.legend(numpoints=1, labelspacing=0, scatterpoints=1,
                        **extra_args)
    for t in leg.get_texts():
        t.set_fontsize('small')    # the legend text fontsize
def find_rho(radius=0.2,n=1.36,NA=0.37):
    """ Rho is a measure of light spread
    radius: (mm) In Gradinaru paper, they report the diameter as 400 um
    n: index of refraction of gray matter
    NA: numerical aperture of the optical fiber
    """
    from numpy import sqrt
    return radius * sqrt((n/(NA**2)) - 1)
def spreading(dist,radius=0.2,n=1.36,NA=0.37):
    rho = find_rho(radius,n,NA)
    return (rho ** 2) / ((rho + dist) ** 2)
