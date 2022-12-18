#!env python3
import math
from linear import Vector3
from perspective import Camera
import quaternion
# from prism import ConvexPolygon, OCTAGON, Prism
import prism

base = prism.ConvexPolygon([
        prism.Point(-0.5, -0.5), prism.Point( 0.5, -0.5),
        prism.Point( 0.5,  0.5), prism.Point(-0.5,  0.5)])
conformal = quaternion.Quaternion(0.0, Vector3(a=1/math.sqrt(2.0)))
shear = Vector3(0, 0, 1.0)
camera = Camera(base, conformal, Vector3(), shear)
prism2 = prism.Prism(base, conformal, Vector3(c=-5), shear)

# prism2 = Prism(ConvexPolygon(OCTAGON),
#         quaternion.Quaternion(0.8, Vector3(0.7, 0.7, 0.7)),
#         Vector3(6.0, -3.0, -7.0),
#         Vector3(b=-18.0))
print(camera.view(prism2).ply())
