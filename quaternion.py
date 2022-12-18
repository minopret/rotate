#!/usr/bin/env python3
# Aaron Mansheim <minopret@gmail.com> 2022
import abc
import doctest
from math import atan2, hypot, sin, sqrt
import numbers
import linear

class Vector4(linear.Vector):
    __slots__ = ('x', 'y', 'z', 'w')
    # pylint: disable=invalid-name
    def __init__(self, /, x=0.0, y=0.0, z=0.0, w=0.0):
        self.x, self.y, self.z, self.w = x, y, z, w
    @staticmethod
    def indices():
        return (0, 1, 2, 3)
    def __len__(self):
        return 4
    def _getitem(self, key):
        indices = self.__class__.indices()
        return (self.x, self.y, self.z, self.w)[indices[key]]

class DivRingElement(linear.AddableElement):
    def __init_subclass__(cls, divring=None, **kwargs):
        super().__init_subclass__(addgroup=divring, **kwargs)
        cls.ONE = cls(1)
    @abc.abstractmethod
    def __mul__(self, other):
        pass
    @abc.abstractmethod
    def __truediv__(self, other):
        pass

class Quaternion(DivRingElement):
    '''
    Approximate an element of H ("Hamilton's"), the incomputable set of all
    quaternions in mathematics.
    '''
    def __init__(self, scalar=0, vector=None):
        '''
        >>> q1 = Quaternion(31/22.0, linear.Vector3(13/22.0, 15/22.0, 19/22.0))
        >>> q1.scalar
        1.4090909090909092
        >>> q1.vector
        Vector3(0.5909090909090909, 0.6818181818181818, 0.8636363636363636)
        '''
        self.scalar, self.vector = scalar, linear.Vector3(*(vector or tuple()))
    def __abs__(self):
        'Returns the real distance from 0. Called for abs(self).'
        return hypot(self.scalar, *self.vector)
    def __eq__(self, other):
        '''
        self == other

        >>> q1 = Quaternion(31/22.0, linear.Vector3(13/22.0, 15/22.0, 19/22.0))
        >>> q2 = Quaternion(31/22.0, linear.Vector3(13/22.0, 15/22.0, 19/22.0))
        >>> q3 = Quaternion(vector=linear.Vector3())
        >>> q1 == q2
        True
        >>> q1 == q3
        False
        '''
        return (self.scalar, self.vector) == (other.scalar, other.vector)
    def __hash__(self):
        return hash((self.scalar, self.vector))
    def __mul__(self, other):
        'self * other'
        if isinstance(other, (Quaternion, Versor)):
            return Quaternion(scalar=self.scalar * other.scalar - self.vector * other.vector,
                    vector=(other.vector * self.scalar + self.vector * other.scalar +
                        self.vector.cross(other.vector)))
        if isinstance(other, numbers.Number):
            return Quaternion(scalar=self.scalar * other, vector=self.vector * other)
        return NotImplemented
    def __repr__(self):
        '''
        >>> q1 = Quaternion(31/22.0, linear.Vector3(13/22.0, 15/22.0, 19/22.0))
        >>> print(q1)
        Quaternion(1.40909090909091, 0.59090909090909, 0.68181818181818, 0.86363636363636)
        '''
        scalar, vector1, vector2, vector3 = (
                round(x, 14) for x in [self.scalar] + list(self.vector))
        return f'{self.__class__.__name__}({scalar}, {vector1}, {vector2}, {vector3})'
    def __sub__(self, other):
        'self - other'
        return Quaternion(scalar=self.scalar - other.scalar,
                vector=(self.vector[i] - other.vector[i] for i in range(3)))
    def __truediv__(self, other):
        'self times multiplicative inverse of other'
        if isinstance(other, Versor):
            return self * other.conjugate()
        if isinstance(other, Quaternion):
            square = other.scalar**2 + other.vector**2
            return Quaternion(scalar=self.scalar / square, vector=self.vector / square
                    ) * other.conjugate()
        if isinstance(other, numbers.Number):
            return Quaternion(scalar=self.scalar / other, vector=self.vector / other)
        return NotImplemented
    def conjugate(self):
        'Quaternion(scalar=a, vector=v).conjugate() == Quaternion(scalar=a, vector=-v).'
        return Quaternion(scalar=self.scalar, vector=-self.vector)
    def apply(self, point):
        extent = point.proj(self.vector)
        deviation = point - extent
        absdev = abs(deviation)
        return (self.conjugate() * Quaternion(0, deviation * (1 / absdev)) * self
                ).vector * absdev + extent
    def rotate(self, point):
        return (self.conjugate() * Quaternion(0, point) * self).vector
    def rotation(self):
        '''
        Compute vectors ri, rj, rk such that, for each vector (a, b, c),
        this quaternion rotates it to ``a ri + b rj + c rk``.
        '''
        i, j, k = linear.Vector3(a=1), linear.Vector3(b=1), linear.Vector3(c=1)
        return [self.rotate(i), self.rotate(j), self.rotate(k)]

Quaternion.__neg__.__doc__ = \
        '''
        -self

        >>> q1 = -Quaternion(31/22.0, linear.Vector3(13/22.0, 15/22.0, 19/22.0))
        >>> q1.scalar
        -1.4090909090909092
        >>> q1.vector
        Vector3(-0.5909090909090909, -0.6818181818181818, -0.8636363636363636)
        '''

Quaternion.__pos__.__doc__ = \
        '''
        +self

        >>> q1 = +Quaternion(31/22.0, linear.Vector3(13/22.0, 15/22.0, 19/22.0))
        >>> q2 =  Quaternion(31/22.0, linear.Vector3(13/22.0, 15/22.0, 19/22.0))
        >>> q1 == q2
        True
        '''

class Versor(Quaternion, divring=Quaternion):
    '''
    Approximate a subset V of H, where H is the incomputable set of all
    quaternions in mathematics and V = {q for q in H if abs(q) == 1}.
    The set V is closed under products, quotients, and additive inverses,
    but not addition.
    '''

    CLASS_VERSION='1.0.0'

    def __init__(self, vector=None, positive=True, *, scalar=None):
        '''
        >>> ( Versor(vector=(1/11, 2/11, 4/11)) ).scalar * 11
        10.0
        >>> ( Versor(vector=(1/11, 2/11, 4/11), positive=False) ).scalar * 11
        -10.0
        >>> (Versor.ONE.scalar, Versor.ONE.vector)
        (1.0, Vector3(0.0, 0.0, 0.0))
        '''
        if scalar is None:
            if vector == 1 or vector is None:
                scalar, vector = 1.0, linear.Vector3()
            else:
                vector = linear.Vector3(*vector)
                scalar = sqrt(1.0 - abs(vector)**2)
        scalar *= 1.0 if positive else -1.0
        scalar, *vector = Vector4(scalar, *vector).normalized()
        super().__init__(scalar, vector)
    def __abs__(self):
        '''
        Returns the real distance from 0. Called for abs(self).

        >>> abs( Versor(vector=(1/11, 2/11, 4/11)) )
        1.0
        '''
        return 1.0
    def __add__(self, other):
        '''
        self + other

        >>> print(Versor(vector=(1/11, 2/11, 4/11)) + Versor(vector=(.5, .5, .5)))
        Quaternion(1.40909090909091, 0.59090909090909, 0.68181818181818, 0.86363636363636)
        '''
        # Is a versor if, and only if, one of the following is the case:
        # *) ``other`` equals zero;
        # *) TODO the perpendicular plane that bisects the line segment
        #    (self, other) contains the origin.
        return super().__add__(other) if other else self
    def __bool__(self):
        '''
        True if self != 0. Called for bool(self).

        >>> bool( Versor(vector=(1/11, 2/11, 4/11)) )
        True
        '''
        return True
    def __mul__(self, other):
        '''
        self * other

        >>> from math import cos, isclose, pi, sqrt
        >>> (Versor() * Versor(vector=(1/11, 2/11, 4/11))).scalar * 11
        10.0
        >>> ipos, jpos, kpos = (Versor(vector=v) for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
        >>> ineg, jneg, kneg = (Versor(vector=v) for v in ((-1, 0, 0), (0, -1, 0), (0, 0, -1)))
        >>> assert ipos * jpos == kpos, f'{ipos * jpos} <> {kpos}'
        >>> assert jpos * kpos == ipos, f'{jpos * kpos} <> {ipos}'
        >>> assert kpos * ipos == jpos, f'{kpos * ipos} <> {jpos}'
        >>> assert ipos * kpos == jneg, f'{ipos * kpos} <> {jneg}'
        >>> assert jpos * ipos == kneg, f'{jpos * ipos} <> {kneg}'
        >>> assert kpos * jpos == ineg, f'{kpos * jpos} <> {ineg}'
        >>> assert abs( Versor(vector=(0.5, 0.5, 0.5)) * Versor(vector=(0, sin(pi/6), 0))
        ...     - Versor(vector=(sin(pi/12)/sqrt(2), cos(pi/12)/sqrt(2),
        ...     cos(pi/12)/sqrt(2))) ) < 1e-14
        '''
        if isinstance(other, Versor):
            return Versor(vector=(other.vector * self.scalar + self.vector * other.scalar +
                        self.vector.cross(other.vector)))
        if isinstance(other, Quaternion):
            return super().__mul__(other)
        return NotImplemented
    def __neg__(self):
        '''
        -self

        >>> print(-Versor(vector=(1/11, 2/11, 4/11)))
        Versor(-0.90909090909091, -0.09090909090909, -0.18181818181818, -0.36363636363636)
        '''
        return Versor(vector=-self.vector, positive=(self.scalar < 0))
    def __rtruediv__(self, other):
        'other / self'
        if isinstance(other, numbers.Number):
            return self.conjugate() * other
        return NotImplemented
    def __repr__(self):
        scalar, vector1, vector2, vector3 = (
                round(x, 14) for x in [self.scalar] + list(self.vector))
        return f'{self.__class__.__name__}({scalar}, {vector1}, {vector2}, {vector3})'
    #def __sub__(self, other):
    #    'self - other'
        # Is a versor if, and only if, one of the following is the case:
        # *) ``other`` equals zero;
        # *) TODO the perpendicular plane that bisects the line segment
        #    (self, -other) contains the origin.
    #    return super().__sub__(other)
    def __truediv__(self, other):
        '''
        self times multiplicative inverse of other

        >>> print( Versor.ONE / Versor(vector=(1/11, 2/11, 4/11)) )
        Versor(0.90909090909091, -0.09090909090909, -0.18181818181818, -0.36363636363636)
        >>> print( Versor(vector=(1/11, 2/11, 4/11)) / Versor(vector=(1/11, 2/11, 4/11)) )
        Versor(1.0, 0.0, 0.0, 0.0)
        '''
        if self == other == Versor.ONE or isinstance(other, numbers.Number) and other == 1:
            return self
        if isinstance(other, Versor):
            return self * other.conjugate()
        return super().__truediv__(other)
    def conjugate(self):
        '''
        Versor(vector=v).conjugate() == -Versor(vector=v, positive=(v.scalar < 0))

        >>> print(Versor(vector=(1/11, 2/11, 4/11)).conjugate())
        Versor(0.90909090909091, -0.09090909090909, -0.18181818181818, -0.36363636363636)
        >>> print(( -Versor((1/11, 2/11, 4/11)) ).conjugate())
        Versor(-0.90909090909091, 0.09090909090909, 0.18181818181818, 0.36363636363636)
        '''
        return Versor(vector=-self.vector, positive=(self.scalar >= 0))
    def angle(self):
        '''
        Express axis-angle rotation as a versor:
        Q =(cos θ/2, k sin θ/2) where k is the axis, a unit; θ is the angle.
        k = 0 if (θ = 0) else x / sin θ/2

        >>> from math import cos, isclose, pi
        >>> theta = pi/2
        >>> assert isclose(theta,
        ...     Versor(vector=(cos(pi/4)*sin(theta/2), sin(pi/4)*sin(theta/2), 0)).angle())
        '''
        return 2 * atan2(abs(self.vector), self.scalar)
    def axis(self):
        '''
        >>> print(tuple(round(x**2, 14) for x in Versor(vector=(.5, .5, .5)).axis()))
        (0.33333333333333, 0.33333333333333, 0.33333333333333)
        '''
        theta = self.angle()
        return 0 if theta == 0 else self.vector / sin(theta/2)
    @classmethod
    def fromangleaxis(cls, angle, axis):
        return Versor(vector=axis.normalized() * sin(angle/2))

Versor.apply.__doc__ = \
        '''
        >>> from math import isclose, pi
        >>> p1 = linear.Vector3(1, 1, 1)
        >>> r = linear.Vector3(1/11, 2/11, 4/11)
        >>> u1 = Versor.fromangleaxis(pi / 2, r.normalized())
        >>> p2 = u1.apply(p1)
        >>> s1 = p1 - p1.proj(u1.vector)
        >>> s2 = p2 - p1.proj(u1.vector)
        >>> assert isclose(u1.vector.anglebetween(p1), u1.vector.anglebetween(p2))
        >>> assert isclose(abs(p1), abs(p2))
        >>> assert isclose(u1.vector * s2, 0, abs_tol=1e-9)
        >>> assert isclose(u1.vector * s1, 0, abs_tol=1e-9)
        >>> assert isclose(s2 * p1, 0, abs_tol=1e-9), f'{s2 * p1} >> 0'
        >>> assert isclose(s2 * s1, 0, abs_tol=1e-9), f'{s2 * s1} >> 0'
        >>> assert isclose(p2 * s1, 0, abs_tol=1e-9), f'{p2 * s1} >> 0'
        >>> u2 = Versor.fromangleaxis(pi / 3, r.normalized())
        >>> p3 = u2.apply(p1)
        >>> s3 = p3 - p1.proj(u2.vector)
        >>> assert isclose(u2.vector.anglebetween(p1), u2.vector.anglebetween(p3)
        ...     ), f'{u2.vector.anglebetween(p1)} <> {u2.vector.anglebetween(p3)}'
        >>> assert isclose(s1.anglebetween(s3), pi / 3), f'{s1.anglebetween(s3)} <> {pi / 3}'
        >>> assert isclose(abs(p1), abs(p3)), f'{abs(p1)} <> {abs(p3)}'
        >>> assert isclose(u2.vector * s3, 1e-9, abs_tol=1e-9), f'{u2.vector * s3} <> 0'
        >>> assert isclose(u2.vector * s1, 1e-9, abs_tol=1e-9), f'{u2.vector * s1} >> 0'
        '''

if __name__ == '__main__':
    doctest.testmod()
