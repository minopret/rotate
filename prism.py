#!/usr/bin/env python3
# pylint: disable=missing-docstring
# pylint: disable=unsubscriptable-object
import doctest
import math
import linear
import quaternion

class Point(linear.Vector):
    __slots__ = tuple()

    def ccw(self, pt0):
        '''
        Returns whether a path in two dimensions turns counterclockwise if,
        upon completing a line segment from origin to ``pt0``, it proceeds
        by a line segment to ``self``.
        '''
        pt1 = self
        return (pt0.y * pt1.x < 0 if pt0.x == 0 else
             pt0.x * pt1.y > 0 if pt1.x == 0 else
             pt0.y/pt0.x < pt1.y/pt1.x if pt0.x * pt1.x > 0 else
             pt0.y/pt0.x > pt1.y/pt1.x)

    def __sub__(self, pt1):
        pt0 = self
        return Point(pt0.x - pt1.x, pt0.y - pt1.y)

POLYGON32 = (lambda power_of_angle, angle: [
    Point(*(lambda z: (2 * z.real, 2 * z.imag))(power_of_angle := power_of_angle * angle))
    for i in range(32)])(1, math.cos(math.tau/32) + math.sin(math.tau/32) * 1j)
POLYGON16 = [Point(*p) for p in zip(
    [0,  1,    2, 3,   4, 4.125, 4, 3,   2, 1,    0, -1,  -2, -2.125,-2,-1,  ],
    [0, -0.25, 0, 0.5, 2, 3,     4, 5.5, 6, 6.25, 6,  5.5, 4,  3,     2, 0.5,])]
OCTAGON = POLYGON16[::2]
DART_SW = [Point(*p)
        for p in zip([0, 2, 1, 1], [0, 0, 1, 2])]
DART_NW = [Point(*p)
        for p in zip([0, 2, 2, -1], [0, 0, 3, -1])]

class ConvexPolygon:
    def __init__(self, path, validate=True):
        self.path = path
        if validate and not self.is_valid():
            raise ValueError('path is not both convex and closed')
    def is_valid(self):
        '''
        >>> ConvexPolygon(DART_SW, False).is_valid()
        False
        >>> ConvexPolygon(DART_NW, False).is_valid()
        False
        >>> ConvexPolygon(OCTAGON, False).is_valid()
        True
        >>> ConvexPolygon(POLYGON16, False).is_valid()
        True
        '''
        path = self.path
        iterator = iter(path + path[0:2])
        origin = next(iterator)
        start = next(iterator)
        for end in iterator:
            relstart = start - origin
            relend = end - origin
            if not relend.ccw(relstart):
                return False
            origin = start
            start = end
        return True
    @classmethod
    def get(cls, path, default):
        '''
        >>> len(ConvexPolygon.get(POLYGON16, None).path)
        16
        >>> ConvexPolygon.get(DART_NW, None) is None
        True
        '''
        polygon = ConvexPolygon(path, validate=False)
        if polygon.is_valid():
            return polygon
        return default

class Prism:
    '''
    A solid formed by sliding a polygon from an
    initial position to a final position.
    It can be placed at any coordinate in a 3D space,
    rotated in any way, and mirror-reversed.
    As the polygon slides, its size can remain the
    same, linearly increase, or linearly decrease.

    The polygon is defined as a path that returns to
    its starting point. The polygon's slide is defined
    by the vector named ``sheer``. It is the vector
    traveled by that starting point, as the polygon
    moves from its initial position to its final position.
    '''
    def __init__(self,
            polygon: ConvexPolygon,
            rotreflscale: quaternion.Quaternion,
            start: linear.Vector3,
            sheer: linear.Vector3):
        self.polygon = polygon
        self.rotreflscale = rotreflscale
        self.start = start
        self.sheer = sheer
    def vertices(self):
        v_i, v_j, _v_k = quaternion.Versor(
                self.rotreflscale.vector, scalar=self.rotreflscale.scalar
                ).rotation()
        start = self.start
        path1 = [v_i * p_x + v_j * p_y + start
                for p_x, p_y in self.polygon.path]
        v_i, v_j, _v_k = self.rotreflscale.rotation()
        sheer = start + self.sheer
        path2 = [v_i * p_x + v_j * p_y + sheer
                for p_x, p_y in self.polygon.path]
        return [path1, path2]

if __name__ == '__main__':
    doctest.testmod()
    prism = Prism(ConvexPolygon(OCTAGON),
            quaternion.Quaternion(2, linear.Vector3(1.0, 1.0, 1.0)),
            linear.Vector3(2.0, 3.0, 5.0),
            linear.Vector3(4.0, -4.0, 4.0))
    print(prism.vertices())
