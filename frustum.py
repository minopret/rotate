#!/usr/bin/env python3
# Aaron Mansheim <minopret@gmail.com> 2022
# pylint: disable=missing-docstring
# pylint: disable=unsubscriptable-object
import doctest
import math
import linear
import ply
import quaternion

class Point(linear.Vector2):
    __slots__ = tuple()

POLYGON32 = (lambda power_of_angle, angle: [
    Point(*(lambda z: (2 * z.real, 2 * z.imag))(power_of_angle := power_of_angle * angle))
    for i in range(32)])(1, math.cos(math.tau/32) + math.sin(math.tau/32) * 1j)
POLYGON16 = [Point(*p) for p in zip(
    [ 7,  6,  5,  3,  0, -3, -5, -6, -7, -6, -5, -3,  0,  3,  5,  6, ],
    [ 0,  3,  5,  6,  7,  6,  5,  3,  0, -3, -5, -6, -7, -6, -5, -3, ])]
OCTAGON = POLYGON16[::2]
QUADRILATERAL = POLYGON16[::4]
DART_SW = [Point(*p)
        for p in zip([0, 2, 1,  1], [0, 0, 1,  2])]
DART_NW = [Point(*p)
        for p in zip([0, 2, 2, -1], [0, 0, 3, -1])]

def fan(size):
    '''
    Specify a sequence for visiting vertices of
    a convex polygon, so that the polygon is the
    disjoint union of the triangles that are
    defined by each three vertices in the sequence.
    '''

    # Is it helpful that it's a flat structure?
    focus = size - 2
    start = size - 1
    for end in range(focus):
        yield focus
        yield start
        yield end
        start = end

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
        if len(path) < 3:
            return False
        origin = path[-2]
        start = path[-1]
        for end in path:
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

class Frustum:
    '''
    Define the boundary of a solid by a shell interpolated linearly
    between a starting path that lies in some plane and an ending path
    that lies in some parallel plane.

    The starting path must be a convex polygon.

    The ending path must be "the same polygon" in the sense
    that it is parallel to the starting path and geometrically
    similar to the starting path.

    The ending path must be related to the starting path by a
    linear transformation that is either a translation or else a
    composition of a translation with a scaling.

    The starting path must be expressed as ``polygon``, which must
    be a counter-clockwise path in two-dimensional coordinates.

    Its placement in space is defined by vector ``start`` and quaternion
    ``conformal``. That quaternion is interpreted as a conformal linear
    transformation. The placement is as follows.

        - Rotate and possibly reflect ``polygon``, by the conjugate
          operation of the unit quaternion that is the result of
          normalizing ``conformal``.
        - Add the vector ``start``.

    The placement of the ending path in space is defined as follows:

        - Rotate and possibly reflect ``polygon``, by the conjugate
          operation of ``conformal``.
        - Add the vector ``start``.
        - Add the vector ``shear``.
    '''
    def __init__(self,
            polygon: ConvexPolygon,
            conformal: quaternion.Quaternion,
            start: linear.Vector3,
            shear: linear.Vector3):
        self.polygon = polygon
        self.conformal = conformal
        self.start = start
        self.shear = shear
        self.pose = None
    def apex(self):
        '''
        Return the point ``apex``, defined as follows.

        If ``conformal`` is not a unit quaternion and ``path`` contains
        at least one edge not parallel to ``shear``, then it has a
        unique ``apex``, a point contained simultaneously in every line
        that passes through a vertex of the starting path and the
        corresponding vertex of the ending path.
        The apex is in the direction ``-shear`` from ``start``.

        Otherwise, ``conformal`` is a unit quaternion or ``path`` contains
        only edges parallel to ``shear``. Define that the distance
        from ``start`` to ``apex`` is (positive) infinite.

        In the following example, I produce a starting path
        ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)) by defining ``polygon`` as
        ((0, 0   ), (0, -1  ), (1, -1  ), (1, 0   )), then defining ``conformal`` as
        a half-turn on the x-axis.

        >>> frustum = Frustum(ConvexPolygon([Point(0, 0), Point(0, -1), Point(1, -1), Point(1, 0)]),
        ...         quaternion.Quaternion(0.0, linear.Vector3(1/math.sqrt(2.0), 0.0, 0.0)),
        ...         linear.Vector3(), linear.Vector3(0.25, 0.25, 1.0))
        >>> print(frustum.apex())
        Vector3(0.4999999999999999, 0.4999999999999999, 1.9999999999999996)
        '''

        pose = self.vertices()
        vertex1 = pose[0][0] # == self.start
        vertex2 = pose[0][1]
        dir1 = pose[1][0] - vertex1 # == self.shear
        dir2 = pose[1][1] - vertex2
        return vertex1 + dir1 * vertex1.line_approaches_line(dir1, vertex2, dir2)
    def vertices(self) -> tuple:
        r'''
        >>> frustum = Frustum(ConvexPolygon(OCTAGON),
        ...         quaternion.Quaternion(2, linear.Vector3(1.0, 1.0, 1.0)),
        ...         linear.Vector3(2.0,  3.0, 5.0),
        ...         linear.Vector3(4.0, -4.0, 4.0))
        >>> print(";\n".join([",\n".join([str(vertex) for vertex in path])
        ... for path in frustum.vertices()]))
        Vector3(5.0, 1.0, 11.0),
        Vector3(8.428571428571427, 3.7142857142857144, 7.857142857142856),
        Vector3(7.999999999999999, 6.0, 3.0),
        Vector3(4.142857142857142, 6.571428571428571, -0.7142857142857135),
        Vector3(-1.0, 5.0, -0.9999999999999991),
        Vector3(-4.428571428571427, 2.2857142857142856, 2.1428571428571437),
        Vector3(-3.999999999999999, 0.0, 7.0),
        Vector3(-0.1428571428571419, -0.5714285714285712, 10.714285714285714);
        Vector3(27.0, -15.0, 51.0),
        Vector3(51.0, 4.0, 29.0),
        Vector3(48.0, 20.0, -5.0),
        Vector3(21.0, 24.0, -31.0),
        Vector3(-15.0, 13.0, -33.0),
        Vector3(-39.0, -6.0, -11.0),
        Vector3(-36.0, -22.0, 23.0),
        Vector3(-9.0, -26.0, 49.0)
        '''
        if self.pose is None:
            v_i, v_j, _v_k = quaternion.Versor(
                    self.conformal.vector,
                    scalar=self.conformal.scalar
                    ).rotation()
            start = self.start
            path1 = [v_i * p_x + v_j * p_y + start
                    for p_x, p_y in self.polygon.path]
            v_i, v_j, _v_k = self.conformal.rotation()
            shear = start + self.shear
            path2 = [v_i * p_x + v_j * p_y + shear
                    for p_x, p_y in self.polygon.path]
            self.pose = (path1, path2)
        return self.pose

    def triangulation(self):
        r'''
        >>> frustum = Frustum(ConvexPolygon(QUADRILATERAL),
        ...         quaternion.Quaternion(2, linear.Vector3(1.0, 1.0, 1.0)),
        ...         linear.Vector3(2.0,  3.0, 5.0),
        ...         linear.Vector3(4.0, -4.0, 4.0))
        >>> print("\n".join(str(point) for point in frustum.triangulation()))
        Vector3(-1.0, 5.0, -0.9999999999999991)
        Vector3(-3.999999999999999, 0.0, 7.0)
        Vector3(5.0, 1.0, 11.0)
        Vector3(-1.0, 5.0, -0.9999999999999991)
        Vector3(5.0, 1.0, 11.0)
        Vector3(7.999999999999999, 6.0, 3.0)
        Vector3(-15.0, 13.0, -33.0)
        Vector3(-36.0, -22.0, 23.0)
        Vector3(27.0, -15.0, 51.0)
        Vector3(-15.0, 13.0, -33.0)
        Vector3(27.0, -15.0, 51.0)
        Vector3(48.0, 20.0, -5.0)
        Vector3(48.0, 20.0, -5.0)
        Vector3(7.999999999999999, 6.0, 3.0)
        Vector3(5.0, 1.0, 11.0)
        Vector3(48.0, 20.0, -5.0)
        Vector3(5.0, 1.0, 11.0)
        Vector3(27.0, -15.0, 51.0)
        Vector3(-15.0, 13.0, -33.0)
        Vector3(-1.0, 5.0, -0.9999999999999991)
        Vector3(7.999999999999999, 6.0, 3.0)
        Vector3(-15.0, 13.0, -33.0)
        Vector3(7.999999999999999, 6.0, 3.0)
        Vector3(48.0, 20.0, -5.0)
        Vector3(-36.0, -22.0, 23.0)
        Vector3(-3.999999999999999, 0.0, 7.0)
        Vector3(-1.0, 5.0, -0.9999999999999991)
        Vector3(-36.0, -22.0, 23.0)
        Vector3(-1.0, 5.0, -0.9999999999999991)
        Vector3(-15.0, 13.0, -33.0)
        Vector3(27.0, -15.0, 51.0)
        Vector3(5.0, 1.0, 11.0)
        Vector3(-3.999999999999999, 0.0, 7.0)
        Vector3(27.0, -15.0, 51.0)
        Vector3(-3.999999999999999, 0.0, 7.0)
        Vector3(-36.0, -22.0, 23.0)
        '''

        # Correct, but can be simplified.
        (path1, path2) = self.vertices()
        size = len(path1)
        result = [path1[i] for i in fan(size)]
        result += [path2[i] for i in fan(size)]
        result += self._parallelograms(path1, path2)
        return result

    @staticmethod
    def _parallelograms(path1, path2):
        tri4 = list(fan(4))
        indices = list(range(len(path1)))
        parallelograms = [(path1[a], path2[a], path2[b], path1[b])
                for (a, b) in zip(indices, indices[1:] + [indices[0]])]
        return [parallelogram[i]
                for parallelogram in parallelograms
                for i in tri4]

    def vertex_list(self):
        return tuple((point for path in self.vertices() for point in path))

    def face_list(self):
        size = len(self.vertex_list()) // 2
        result = list(fan(size))
        result += [i + size for i in result]
        path1 = list(range(size))
        path2 = [x + size for x in path1]
        result += self._parallelograms(path1, path2)
        return [result[i:i + 3] for i in range(0, len(result), 3)]

    def ply(self):
        """
        >>> frustum = Frustum(ConvexPolygon(OCTAGON),
        ...         quaternion.Quaternion(0.8, linear.Vector3(0.7, 0.7, 0.7)),
        ...         linear.Vector3(6.0, -3.0, -7.0),
        ...         linear.Vector3(b=-18.0))
        >>> print(frustum.ply())
        ply
        format ascii 1.0
        element vertex 16
        property float x
        property float y
        property float z
        element face 28
        property list uchar int vertex_index
        end_header
        6.49763033175356 -3.46445497630332 -0.03317535545024
        11.33175355450237 -2.97630331753555 -2.35545023696683
        12.96682464454976 -2.50236966824644 -7.46445497630332
        10.62085308056872 -2.31279620853080 -12.30805687203792
        5.50236966824644 -2.53554502369668 -13.96682464454976
        0.66824644549763 -3.02369668246445 -11.64454976303318
        -0.96682464454976 -3.49763033175356 -6.53554502369668
        1.37914691943128 -3.68720379146920 -1.69194312796208
        7.05000000000000 -21.98000000000000 7.70000000000000
        17.25000000000000 -20.95000000000000 2.80000000000000
        20.70000000000000 -19.95000000000000 -7.98000000000000
        15.75000000000000 -19.55000000000000 -18.20000000000000
        4.95000000000000 -20.02000000000000 -21.70000000000000
        -5.25000000000000 -21.05000000000000 -16.80000000000000
        -8.70000000000000 -22.05000000000000 -6.02000000000000
        -3.75000000000000 -22.45000000000000 4.20000000000000
        3 6 7 0
        3 6 0 1
        3 6 1 2
        3 6 2 3
        3 6 3 4
        3 6 4 5
        3 14 15 8
        3 14 8 9
        3 14 9 10
        3 14 10 11
        3 14 11 12
        3 14 12 13
        3 9 1 0
        3 9 0 8
        3 10 2 1
        3 10 1 9
        3 11 3 2
        3 11 2 10
        3 12 4 3
        3 12 3 11
        3 13 5 4
        3 13 4 12
        3 14 6 5
        3 14 5 13
        3 15 7 6
        3 15 6 14
        3 8 0 7
        3 8 7 15
        """
        return ply.PLY('ascii',
                list(tuple((f'{value:.14f}' for value in vertex)) for vertex in self.vertex_list()),
                list(tuple([str(len(face))] + [str(value) for value in face])
                    for face in self.face_list()))

if __name__ == '__main__':
    doctest.testmod()
