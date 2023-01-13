#!env python3
import math
from linear import Vector3
from perspective import Camera
import quaternion
# from frustum import ConvexPolygon, OCTAGON, Frustum
import frustum

base = frustum.ConvexPolygon([
        frustum.Point(-0.5, -0.5), frustum.Point( 0.5, -0.5),
        frustum.Point( 0.5,  0.5), frustum.Point(-0.5,  0.5)])
conformal = quaternion.Quaternion(0.0, Vector3(a=1/math.sqrt(2.0)))
shear = Vector3(0, 0, 1.0)
camera = Camera(base, conformal, Vector3(), shear)
frustum2 = frustum.Frustum(base, conformal, Vector3(c=-5), shear)

# frustum2 = Frustum(ConvexPolygon(OCTAGON),
#         quaternion.Quaternion(0.8, Vector3(0.7, 0.7, 0.7)),
#         Vector3(6.0, -3.0, -7.0),
#         Vector3(b=-18.0))
print(camera.view(frustum2).ply())
