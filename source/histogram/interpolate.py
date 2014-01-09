
__author__ = "Andreas Weiler <andreas.weiler@desy.de>"
__version__ = "0.1"

######################################################################
# Helpers
######################################################################

def _myDet(p, q, r):
    """Calc. determinant of a special matrix with three 2D points.

    The sign, "-" or "+", determines the side, right or left,
    respectivly, on which the point r lies, when measured against
    a directed vector from p to q.
    """

    # We use Sarrus' Rule to calculate the determinant.
    # (could also use the Numeric package...)
    sum1 = q[0]*r[1] + p[0]*q[1] + r[0]*p[1]
    sum2 = q[0]*p[1] + r[0]*q[1] + p[0]*r[1]

    return sum1 - sum2


def _isRightTurn((p, q, r)):
    "Do the vectors pq:qr form a right turn, or not?"

    assert p != q and q != r and p != r

    if _myDet(p, q, r) < 0:
        return 1
    else:
        return 0


def convexHull(P):
    "Calculate the convex hull of a set of points."

    # Get a local list copy of the points and sort them lexically.
    points = map(None, P)
    points.sort()

    # Build upper half of the hull.
    upper = [points[0], points[1]]
    for p in points[2:]:
        upper.append(p)
    while len(upper) > 2 and not _isRightTurn(upper[-3:]):
        del upper[-2]

    # Build lower half of the hull.
    points.reverse()
    lower = [points[0], points[1]]
    for p in points[2:]:
        lower.append(p)
    while len(lower) > 2 and not _isRightTurn(lower[-3:]):
        del lower[-2]

    # Remove duplicates.
    del lower[0]
    del lower[-1]

    # Concatenate both halfs and return.
    return tuple(upper + lower)


def isPointInPolygon(r, P):
    "Is point r inside a given polygon P?"

    # We assume the polygon is a list of points, listed clockwise!
    for i in xrange(len(P[:-1])):
        p, q = P[i], P[i+1]
        # if not _isRightTurn((p, q, r)):
        #     return 0 # Out!

    return 1 # It's within!


######################################################################
# Interpolation
######################################################################

def Interpolate1D(listIn, x):
    '''linear interpolated function y = f(x)

    Expects a (not necessarily sorted) list (xi ,yi)

    listIn = [ [1.,2.] , [2,-3.], [4.,3.] ]

    and calculates linearly interpolated y at x.
    '''

    listIn = [[float(xi), float(yi)] for xi, yi in listIn]
    #print listIn
    listIn.sort()   # sort in place
    minX = listIn[0][0]
    maxX = listIn[-1][0]    # last element
    if minX > maxX:
        errmsg = "No valid list given: " + listIn + " !"
        raise Exception(errmsg)
    elif not minX <= x <= maxX:
        errmsg = "x-value (" + str(x) + ")outside boundaries (" + \
        str(minX) + "," + str(maxX) + ")! \n"
        print errmsg
        raise Exception(errmsg)

    distlist = [[abs(xel - x), listIn.index([xel, yel])] for xel, yel in listIn]
    distlist.sort()     # list of distance and index number
    x1Index = distlist[0][1]   # index of closest element
    x2Index = distlist[1][1]

    x1, y1 = listIn[x1Index]
    x2, y2 = listIn[x2Index]

    if x1 == x2:
        errmsg = "two entries for same x-value : " + str(x1)
        raise Exception(errmsg)

    yinterpolate = (y1 - y2) / (x1 - x2) * x - (x2 * y1 - x1 * y2) / (x1 - x2)

    return yinterpolate


def Interpolate2D(listIn, x):
    '''linear interpolated function z = f(x,y)

    Expects a (not necessarily sorted) list (xi ,yi, zi)

    listIn = [ [1., 2. ,2.] , [2, 1,-3.], [4.,2, 3.] ,, [2.,2, 3.] ]

    and calculates bilinearly interpolated z at (x,y). Checks if interpolated
    point is in convex hull of supplied list.
    '''

    listCoordinates = [[float(xi), float(yi)] for xi, yi, zi in listIn]
    listIn = [[float(xi), float(yi), float(zi)] for xi, yi, zi in listIn]

    p = (float(x[0]), float(x[1]))

    # check if points is part of grid

    for pis in listCoordinates:
        if p[0] == pis[0] and p[1] == pis[1]:
            print "no interpolation needed"
            return listIn[listCoordinates.index(pis)]

    if len(listIn) < 4:
        raise Exception("Need at least 4 points")

    boundaryPolygon = convexHull(listCoordinates)

    if not isPointInPolygon(p, boundaryPolygon):
        #errmsg = "Point" + str(x) + "is not inside interpolated region!"
        #print errmsg
        #raise Exception(errmsg)
        return -100

    #  Find four closest points and lineraly interpolate

    distlist = [[abs(xel - p[0]) ** 2 + abs(yel - p[1]) ** 2, \
                    listIn.index([xel, yel, zel])] for xel, yel, zel in listIn]
    distlist.sort()     # list of distance and index number
    x1Index = distlist[0][1]   # index of closest element
    x2Index = distlist[1][1]
    x3Index = distlist[2][1]
    x4Index = distlist[3][1]

    x1, y1, z1 = listIn[x1Index]
    x2, y2, z2 = listIn[x2Index]
    x3, y3, z3 = listIn[x3Index]

    # check if first three points are in a line
    # if so go ahead in list of alternatives

    determinant = -x2 * y1 + x3 *y1 + x1 *y2 - x3 *y2 - x1 *y3 + x2* y3

    epsilon = 10**(-16)
    ind = 1
    while abs(determinant) < epsilon:
        #print "det123"
        x3, y3, z3 = listIn[distlist[2+ind][1]]
        ind += 1
        determinant = -x2 * y1 + x3 *y1 + x1 *y2 - x3 *y2 - x1 *y3 + x2* y3

    x4, y4, z4 = listIn[2 + ind]
    determinant124 = -x2 * y1 + x4 *y1 + x1 *y2 - x4 *y2 - x1 *y4 + x2* y4
    determinant234 = -x2 * y4 + x3 *y4 + x4 *y2 - x3 *y2 - x4 *y3 + x2* y3
    determinant134 = -x4 * y1 + x3 *y1 + x1 *y4 - x3 *y4 - x1 *y3 + x4* y3

    # check the same for the 4th point

    while abs(determinant124) < epsilon or abs(determinant234) < epsilon or abs(determinant134) < epsilon:
        #print "det 4th point"
        ind += 1
        x4, y4, z4 = listIn[2 + ind]
        determinant124 = -x2 * y1 + x4 *y1 + x1 *y2 - x4 *y2 - x1 *y4 + x2* y4
        determinant234 = -x2 * y4 + x3 *y4 + x4 *y2 - x3 *y2 - x4 *y3 + x2* y3
        determinant134 = -x4 * y1 + x3 *y1 + x1 *y4 - x3 *y4 - x1 *y3 + x4* y3

# We fit a bivariate polynomial p[x_, y_] := a x y + b y + c x + d


    a = (x2*y3*z1 - x2*y4*z1 - x1*y3*z2 + x1*y4*z2 - x2*y1*z3 + x1*y2*z3 - x1*y4*z3 + \
        x2*y4*z3 + x4*(y2*z1 - y3*z1 - y1*z2 + y3*z2 + y1*z3 - y2*z3) + x2*y1*z4 - \
        x1*y2*z4 + x1*y3*z4 - x2*y3*z4 + x3*(y4*z1 + y1*z2 - y4*z2 - y1*z4 + y2*(-z1 + z4)))/ \
        (x2*(x3*(y2 - y3)*(y1 - y4) - x4*(y1 - y3)*(y2 - y4)) + x1*(x4*(y2 - y3)*(y1 - y4) -\
        x3*(y1 - y3)*(y2 - y4) + x2*(y1 - y2)*(y3 - y4)) + x3*x4*(y1 - y2)*(y3 - y4))

    b = (x1*x4*(y1 - y4)*(z2 - z3) + x3*(x4*(y3 - y4)*(z1 - z2) - x1*(y1 - y3)*(z2 - z4)) + \
        x2*(-(x4*(y2 - y4)*(z1 - z3)) + x3*(y2 - y3)*(z1 - z4) + x1*(y1 - y2)*(z3 - z4)))/ \
        (x2*(x3*(y2 - y3)*(y1 - y4) - x4*(y1 - y3)*(y2 - y4)) + x1*(x4*(y2 - y3)*(y1 - y4) - \
        x3*(y1 - y3)*(y2 - y4) + x2*(y1 - y2)*(y3 - y4)) + x3*x4*(y1 - y2)*(y3 - y4))

    c = (-(x4*y2*y4*z1) + x4*y3*y4*z1 + x1*y1*y3*z2 - x1*y1*y4*z2 + x4*y1*y4*z2 - x4*y3*y4*z2 -\
        x1*y1*y2*z3 + x1*y1*y4*z3 - x4*y1*y4*z3 + x4*y2*y4*z3 + x1*y1*(y2 - y3)*z4 + \
         x3*y3*(y2*z1 - y4*z1 - y1*z2 + y4*z2 + y1*z4 - y2*z4) + x2*y2*(y4*z1 + y1*z3 - y4*z3 -\
        y1*z4 + y3*(-z1 + z4)))/   (x2*(x3*(y2 - y3)*(y1 - y4) - x4*(y1 - y3)*(y2 - y4)) + \
       x1*(x4*(y2 - y3)*(y1 - y4) - x3*(y1 - y3)*(y2 - y4) + x2*(y1 - y2)*(y3 - y4)) + x3*x4*(y1 - y2)*(y3 - y4))

    d = (-(x1*x4*(y1 - y4)*(y3*z2 - y2*z3)) + x3* \
        (-(x4*(y3 - y4)*(y2*z1 - y1*z2)) + x1*(y1 - y3)*(y4*z2 - y2*z4)) +    \
        x2*(x4*(y2 - y4)*(y3*z1 - y1*z3) - x3*(y2 - y3)*(y4*z1 - y1*z4) - x1*(y1 - y2)*(y4*z3 - y3*z4)))/\
        (x2*(x3*(y2 - y3)*(y1 - y4) - x4*(y1 - y3)*(y2 - y4)) + \
         x1*(x4*(y2 - y3)*(y1 - y4) - x3*(y1 - y3)*(y2 - y4) + x2*(y1 - y2)*(y3 - y4)) + \
         x3*x4*(y1 - y2)*(y3 - y4))

    zinterpolate = a * p[0] * p[1] + b * p[1] + c * p[0] + d

    # except ZeroDivisionError:
    #     print 'points :', listIn[x1Index], listIn[x2Index], listIn[x3Index], listIn[x4Index]
    #     print 'xint: ', p
    #     import csv
    #     f = open('debug.m', 'w')
    #     f.write(str(listIn))
    #     f.close()
    #     print '----------'

    return zinterpolate
