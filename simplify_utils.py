import numpy as np


try:
    import numba
    numba_available = True
    #import simplify_utils_getSquareSegmentDistance_numba
    print("NOTE: NUMBA is available on this system. NUMBA loaded.")
    from lib.simplify.simplify_utils_getSquareSegmentDistance_numba import getSquareSegmentDistance

except ImportError:
    numba_available = False
    print("!!! WARNING !!! NUMBA is NOT available on this system. Installation is strongly advised to speed up simplification (10-100 fold)")
    #import simplify_utils_getSquareSegmentDistance_no_numba
    from lib.simplify.simplify_utils_getSquareSegmentDistance_no_numba import getSquareSegmentDistance

if numba_available:
    from numba import jit

#
# @jit
# def getSquareSegmentDistance(p, p1, p2):
#     """
# Square distance between point and a segment
# """
#     x = p1[0]
#     y = p1[1]
#
#     dx = p2[0] - x
#     dy = p2[1] - y
#
#     if dx != 0 or dy != 0:
#         t = ((p[0] - x) * dx + (p[1] - y) * dy) / (dx * dx + dy * dy)
#
#         if t > 1:
#             x = p2[0]
#             y = p2[1]
#         elif t > 0:
#             x += dx * t
#             y += dy * t
#
#     dx = p[0] - x
#     dy = p[1] - y
#
#     return dx * dx + dy * dy
#

def getSquareDistanceOriginal(p1, p2):
    """
Square distance between two points
"""
    dx = p1['x'] - p2['x']
    dy = p1['y'] - p2['y']

    return dx * dx + dy * dy

def getSquareDistance_slow(p1, p2):
    """
Square distance between two points
"""
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]

    return dx * dx + dy * dy

def getSquareSegmentDistanceOriginal(p, p1, p2):
    """
Square distance between point and a segment
"""
    x = p1['x']
    y = p1['y']

    dx = p2['x'] - x
    dy = p2['y'] - y

    if dx != 0 or dy != 0:
        t = ((p['x'] - x) * dx + (p['y'] - y) * dy) / (dx * dx + dy * dy)

        if t > 1:
            x = p2['x']
            y = p2['y']
        elif t > 0:
            x += dx * t
            y += dy * t

    dx = p['x'] - x
    dy = p['y'] - y

    return dx * dx + dy * dy


def getSquareSegmentDistance(p, p1, p2):
    """
    Square distance between point and a segment
    """
    x = p1[0]
    y = p1[1]

    dx = p2[0] - x
    dy = p2[1] - y

    if dx != 0 or dy != 0:
        t = ((p[0] - x) * dx + (p[1] - y) * dy) / (dx * dx + dy * dy)

        if t > 1:
            x = p2[0]
            y = p2[1]
        elif t > 0:
            x += dx * t
            y += dy * t

    dx = p[0] - x
    dy = p[1] - y

    return dx * dx + dy * dy
#

# NUMBA SYNTAX: newFct = jit (returnParam (inputParam1, inputParam2)) (oldFct)
# getSquareSegmentDistance = jit(double (double[2], double[2],double[2]))(getSquareSegmentDistance_slow)


def simplifyRadialDistance_slow(points, tolerance):
    length = len(points)
    prev_point = points[0]
    new_points = [prev_point]

    for i in range(length):
        point = points[i]

        if getSquareDistance_slow(point, prev_point) > tolerance:
            new_points.append(point)
            prev_point = point

    if prev_point.all != point.all:
        new_points.append(point)

    return new_points

# NUMBA SYNTAX: newFct = jit (returnParam (inputParam1, inputParam2)) (oldFct)
#simplifyRadialDistance = jit(double[:,:] (double[:,:], double))(simplifyRadialDistance_slow)


def simplifyDouglasPeucker_slow(points, tolerance):
    length = len(points)
    markers = [0] * length # Maybe not the most efficent way?

    first = 0
    last = length - 1

    first_stack = []
    last_stack = []

    new_points = []

    markers[first] = 1
    markers[last] = 1

    while last:
        max_sqdist = 0

        for i in range(first, last):
            sqdist = getSquareSegmentDistance(points[i], points[first], points[last])

            if sqdist > max_sqdist:
                index = i
                max_sqdist = sqdist

        if max_sqdist > tolerance:
            markers[index] = 1

            first_stack.append(first)
            last_stack.append(index)

            first_stack.append(index)
            last_stack.append(last)

        # Can pop an empty array in Javascript, but not Python, so check
        # the length of the list first
        if len(first_stack) == 0:
            first = None
        else:
            first = first_stack.pop()

        if len(last_stack) == 0:
            last = None
        else:
            last = last_stack.pop()

    for i in range(length):
        if markers[i]:
            new_points.append(points[i])

    return new_points

# NUMBA SYNTAX: newFct = jit (returnParam (inputParam1, inputParam2)) (oldFct)
# simplifyDouglasPeucker = jit(double[:,:] (double[:,:], double))(simplifyDouglasPeucker_slow)


def simplify(points, tolerance=0.1, highestQuality=True):
    sqtolerance = tolerance * tolerance

    if not highestQuality:
        points = simplifyRadialDistance_slow(points, sqtolerance)

    points = simplifyDouglasPeucker_slow(points, sqtolerance)

    return points


def simplify_plot_arrays(dep_array, time_array, tolerance=0.1):
    """

    :param dep_arr: dependent variable array, must be the same shape as time_array
    :param time_array: array of time values, must be the same shape as dep_array
    :param tolerance: value to tell the simplification algorithm the distance from the original curve will be tolerated.
    As this is an absolute tolerance, this value should be reduce
    :return: simp_dep, simp_t: simplified arrays
    """

    pt_array = np.vstack((dep_array, time_array)).transpose()

    simple_pts = np.array(simplify(pt_array, tolerance=tolerance))

    simp_dep = simple_pts[:, 0]
    simp_t = simple_pts[:, 1]

    return simp_dep, simp_t