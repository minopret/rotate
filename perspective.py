#!/usr/bin/env python3

import doctest
import prism

def line_meets_plane(online1, online2, inplane, normal):
    numerator = normal * (inplane - online1)
    denominator = normal * (online2 - online1)
    if denominator == 0:
        # No unique solution.
        return online1 if numerator == 0 else None
    return online1 + (online2 - online1) * numerator / denominator

class FlattenedPrism(prism.Prism):
    def __init__(self, *args, pose=None, **kwargs):
        if pose is None:
            pose = tuple()
        self.fixed = pose
        super().__init__(*args, **kwargs)
    def vertices(self) -> tuple:
        return self.fixed

class Camera(prism.Prism):
    def view(self, prism2):
        r"""
        >>> import linear
        >>> import math
        >>> import quaternion
        >>> base = prism.ConvexPolygon([
        ...         prism.Point(-0.5, -0.5), prism.Point( 0.5, -0.5),
        ...         prism.Point( 0.5,  0.5), prism.Point(-0.5,  0.5)])
        >>> conformal = quaternion.Quaternion(0.0, linear.Vector3(a=1/math.sqrt(2.0)))
        >>> shear = linear.Vector3(0, 0, 1.0)
        >>> camera = Camera(base, conformal,linear.Vector3(), shear)
        >>> prism2 = prism.Prism(base, conformal, linear.Vector3(c=-5), shear)
        >>> camera.view(prism2).vertices()
        ([Vector3(-0.14285714285714285, 0.14285714285714285, 0.0), Vector3(0.14285714285714285, 0.14285714285714285, 0.0), Vector3(0.14285714285714285, -0.14285714285714285, 0.0), Vector3(-0.14285714285714285, -0.14285714285714285, 0.0)], [Vector3(-0.0833333333333333, 0.0833333333333333, 0.0), Vector3(0.0833333333333333, 0.0833333333333333, 0.0), Vector3(0.0833333333333333, -0.0833333333333333, 0.0), Vector3(-0.0833333333333333, -0.0833333333333333, 0.0)])
        >>> print(camera.view(prism2).ply())
        ply
        format ascii 1.0
        element vertex 8
        property float x
        property float y
        property float z
        element face 12
        property list uchar int vertex_index
        end_header
        -0.14285714285714 0.14285714285714 0.00000000000000
        0.14285714285714 0.14285714285714 0.00000000000000
        0.14285714285714 -0.14285714285714 0.00000000000000
        -0.14285714285714 -0.14285714285714 0.00000000000000
        -0.08333333333333 0.08333333333333 0.00000000000000
        0.08333333333333 0.08333333333333 0.00000000000000
        0.08333333333333 -0.08333333333333 0.00000000000000
        -0.08333333333333 -0.08333333333333 0.00000000000000
        3 2 3 0
        3 2 0 1
        3 6 7 4
        3 6 4 5
        3 5 1 0
        3 5 0 4
        3 6 2 1
        3 6 1 5
        3 7 3 2
        3 7 2 6
        3 4 0 3
        3 4 3 7
        """
        def projection(eye, inplane, normal):
            def project(point):
                return line_meets_plane(eye, point, inplane, normal)
            return project
        eye = self.apex()
        inplane = self.vertices()[0][0]
        normal = self.conformal.rotation()[2]
        project = projection(eye, inplane, normal)
        pose = prism2.vertices()
        pose = ([project(point) for point in pose[0]],
                [project(point) for point in pose[1]])
        return FlattenedPrism(prism2.polygon, prism2.conformal,
                prism2.start, prism2.shear, pose=pose)

if __name__ == '__main__':
    doctest.testmod()
