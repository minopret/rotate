#!/usr/bin/env python3
import abc
import doctest
import math
import numbers

class Vector(tuple):
    def __new__(cls, *args):
        return super().__new__(cls, args)
    def __add__(self, *args):
        return self.__class__(*( sum(values) for values in zip(*([self] + list(args))) ))
    def __mul__(self, other):
        if isinstance(other, Vector):
            return sum(ours * theirs for ours, theirs in zip(self, other))
        return self.__class__(*(ours * other for ours in self))
    def __neg__(self):
        return self.__class__(*(ours * -1 for ours in self))
    def __pow__(self, exp):
        if exp == 2:
            # The likely case.
            return self * self
        if exp >= 0 and isinstance(exp, int):
            return (1 if exp == 0 else
                    self if exp == 1 else
                    self * (self * self)**((exp - 1) / 2) if exp % 2 else
                    (self * self)**(exp / 2)
        raise ValueError
    def __repr__(self):
        return f'{self.__class__.__name__}' + super().__repr__()
    def __sub__(self, other):
        return self + other * -1
    def __truediv__(self, other):
        return self * (1 / other)
    def proj(self, other):
        '''
        Projection of self ``a`` in direction of another vector ``b``: ``proj_b a``
        '''
        return other * (self * other) / (other * other)
    def norm(self):
        return math.sqrt(self * self)
    def normalized(self):
        return self / self.norm()
    def azimuth(self):
        basisi = Vector3(1, 0, 0)
        projij = Vector3(*self[:2], 0)
        return basisi.anglebetween(projij) * (+1 if self[1] > 0 else -1)
    def elevation(self):
        projij = Vector3(*self[:2], 0)
        return self.anglebetween(projij) * (+1 if self[2] > 0 else -1)
    def anglebetween(self, other):
        return math.acos(abs(other.normalized() * self.normalized()))

class Vector3(Vector):
    def __new__(cls, /, a, b, c):
        return super().__new__(cls, a, b, c)
    def cross(self, other):
        return Vector3(self[1] * other[2] - other[1] * self[2],
                self[2] * other[0] - other[2] * self[0],
                self[0] * other[1] - other[0] * self[1])

class MulGroup(abc.ABC):
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.ONE = cls(1)
    @abc.abstractmethod
    def __mul__(self, other):
        pass
    def __pos__(self):
        '+self'
        return self
    @abc.abstractmethod
    def __truediv__(self, other):
        pass

class DivRingElement(abc.ABC):
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.ZERO = cls()
        cls.ONE = cls(1)
    def __add__(self, other):
        'self + other'
        return self - (self.__class__.ZERO - other)
    def __bool__(self):
        'True if self != 0. Called for bool(self).'
        return self != self.__class__.ZERO
    @abc.abstractmethod
    def __eq__(self, other):
        pass
    @abc.abstractmethod
    def __mul__(self, other):
        pass
    def __neg__(self):
        '-self'
        return self.__class__.ZERO - self
    def __pos__(self):
        '+self'
        return self
    @abc.abstractmethod
    def __sub__(self, other):
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
        >>> q1 = Quaternion(0.5 + 10/11.0, Vector3(0.5 + 1/11.0, 0.5 + 2/11.0, 0.5 + 4/11.0))
        >>> q1.scalar
        1.4090909090909092
        >>> q1.vector
        Vector3(0.5909090909090909, 0.6818181818181819, 0.8636363636363636)
        '''
        self.scalar, self.vector = scalar, Vector3(*vector) if vector else Vector3(0.0, 0.0, 0.0)
    def quaternion(self):
        return self
    def __abs__(self):
        'Returns the real distance from 0. Called for abs(self).'
        return math.sqrt(self.scalar**2 + self.vector**2)
    def __eq__(self, other):
        '''
        self == other

        >>> q1 = Quaternion(0.5 + 10/11.0, Vector3(0.5 + 1/11.0, 0.5 + 2/11.0, 0.5 + 4/11.0))
        >>> q2 = Quaternion(0.5 + 10/11.0, Vector3(0.5 + 1/11.0, 0.5 + 2/11.0, 0.5 + 4/11.0))
        >>> q3 = Quaternion(vector=Vector3(0, 0, 0))
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
        >>> q1 = Quaternion(0.5 + 10/11.0, Vector3(0.5 + 1/11.0, 0.5 + 2/11.0, 0.5 + 4/11.0))
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

Quaternion.__neg__.__doc__ = \
        '''
        -self

        >>> q1 = -Quaternion(0.5 + 10/11.0, Vector3(0.5 + 1/11.0, 0.5 + 2/11.0, 0.5 + 4/11.0))
        >>> q1.scalar
        -1.4090909090909092
        >>> q1.vector
        Vector3(-0.5909090909090909, -0.6818181818181819, -0.8636363636363636)
        '''

Quaternion.__pos__.__doc__ = \
        '''
        +self

        >>> q1 = +Quaternion(0.5 + 10/11.0, Vector3(0.5 + 1/11.0, 0.5 + 2/11.0, 0.5 + 4/11.0))
        >>> q2 =  Quaternion(0.5 + 10/11.0, Vector3(0.5 + 1/11.0, 0.5 + 2/11.0, 0.5 + 4/11.0))
        >>> q1 == q2
        True
        '''

class Versor(MulGroup):
    '''
    Approximate a subset V of H, where H is the incomputable set of all
    quaternions in mathematics and V = {q for q in H if abs(q) == 1}.
    The set V is closed under products, quotients, and additive inverses,
    but not addition.
    '''

    CLASS_VERSION='1.0.0'

    def __init__(self, vector=None, positive=True):
        '''
        >>> ( Versor((1/11, 2/11, 4/11)) ).scalar * 11
        10.0
        >>> ( Versor((1/11, 2/11, 4/11), positive=False) ).scalar * 11
        -10.0
        >>> (Versor.ONE.scalar, Versor.ONE.vector)
        (1.0, Vector3(0.0, 0.0, 0.0))
        '''
        self.vector = Vector3(*vector) if vector and vector != 1 else Vector3(0.0, 0.0, 0.0)
        self.scalar = (1.0 if positive else -1.0) * math.sqrt(1.0 - self.vector**2)
    def quaternion(self):
        return Quaternion(self.scalar, self.vector)
    def __abs__(self):
        '''
        Returns the real distance from 0. Called for abs(self).

        >>> abs( Versor((1/11, 2/11, 4/11)) )
        1.0
        '''
        return 1.0
    def __add__(self, other):
        '''
        self + other

        >>> print(Versor((1/11, 2/11, 4/11)) + Versor((.5, .5, .5)))
        Quaternion(1.40909090909091, 0.59090909090909, 0.68181818181818, 0.86363636363636)
        '''
        # Is a versor if, and only if, one of the following is the case:
        # *) ``other`` equals zero;
        # *) TODO the perpendicular plane that bisects the line segment
        #    (self, other) contains the origin.
        return self.quaternion() + other if other else self
    def __bool__(self):
        '''
        True if self != 0. Called for bool(self).

        >>> bool( Versor((1/11, 2/11, 4/11)) )
        True
        '''
        return True
    def __eq__(self, other):
        '''
        self == other

        >>> Versor() == Versor.ONE
        True
        '''
        return self.quaternion() == other.quaternion()
    def __hash__(self):
        return hash(self.quaternion())
    def __mul__(self, other):
        '''
        self * other

        >>> (Versor() * Versor((1/11, 2/11, 4/11))).scalar * 11
        10.0
        >>> ipos, jpos, kpos = Versor((1, 0, 0)), Versor((0, 1, 0)), Versor((0, 0, 1))
        >>> ineg, jneg, kneg = Versor((-1, 0, 0)), Versor((0, -1, 0)), Versor((0, 0, -1))
        >>> assert ipos * jpos == kpos, f'{ipos * jpos} {kpos}'
        >>> assert jpos * kpos == ipos, f'{jpos * kpos} {ipos}'
        >>> assert kpos * ipos == jpos, f'{kpos * ipos} {jpos}'
        >>> assert ipos * kpos == jneg, f'{ipos * kpos} {jneg}'
        >>> assert jpos * ipos == kneg, f'{jpos * ipos} {kneg}'
        >>> assert kpos * jpos == ineg, f'{kpos * jpos} {ineg}'
        >>> assert abs( Versor((0.5, 0.5, 0.5)) * Versor((0, math.sin(math.pi/6), 0))
        ...     - Versor((math.sin(math.pi/12)/math.sqrt(2), math.cos(math.pi/12)/math.sqrt(2),
        ...     math.cos(math.pi/12)/math.sqrt(2))) ) < 1e-14
        '''
        if isinstance(other, Versor):
            return Versor(vector=(other.vector * self.scalar + self.vector * other.scalar +
                        self.vector.cross(other.vector)))
        if isinstance(other, Quaternion):
            return self.quaternion() * other
        return NotImplemented
    def __neg__(self):
        '''
        -self

        >>> print(-Versor((1/11, 2/11, 4/11)))
        Versor(-0.90909090909091, -0.09090909090909, -0.18181818181818, -0.36363636363636)
        '''
        return Versor(-self.vector, positive=(self.scalar < 0))
    def __rtruediv__(self, other):
        'other / self'
        if isinstance(other, numbers.Number):
            return self.conjugate() * other
        return NotImplemented
    def __repr__(self):
        scalar, vector1, vector2, vector3 = (
                round(x, 14) for x in [self.scalar] + list(self.vector))
        return f'{self.__class__.__name__}({scalar}, {vector1}, {vector2}, {vector3})'
    def __sub__(self, other):
        'self - other'
        # Is a versor if, and only if, one of the following is the case:
        # *) ``other`` equals zero;
        # *) TODO the perpendicular plane that bisects the line segment
        #    (self, -other) contains the origin.
        return self.quaternion() - other
    def __truediv__(self, other):
        '''
        self times multiplicative inverse of other

        >>> print( Versor.ONE / Versor((1/11, 2/11, 4/11)) )
        Versor(0.90909090909091, -0.09090909090909, -0.18181818181818, -0.36363636363636)
        >>> print( Versor((1/11, 2/11, 4/11)) / Versor((1/11, 2/11, 4/11)) )
        Versor(1.0, 0.0, 0.0, 0.0)
        '''
        if self == other == Versor.ONE or isinstance(other, numbers.Number) and other == 1:
            return self
        if isinstance(other, Versor):
            return self * other.conjugate()
        return self.quaternion() / other
    def conjugate(self):
        '''
        Versor(vector=v).conjugate() == -Versor(vector=v, positive=(v.scalar < 0))

        >>> print(Versor((1/11, 2/11, 4/11)).conjugate())
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

        >>> theta = math.pi/2
        >>> abs(theta - Versor((math.cos(math.pi/4)*math.sin(theta/2),
        ...     math.sin(math.pi/4)*math.sin(theta/2), 0)).angle()) < 1e-14
        True
        '''
        return 2 * math.atan2(math.sqrt(self.vector * self.vector), self.scalar)
    def axis(self):
        '''
        >>> print(tuple(round(x**2, 14) for x in Versor((.5, .5, .5)).axis()))
        (0.33333333333333, 0.33333333333333, 0.33333333333333)
        '''
        theta = self.angle()
        return 0 if theta == 0 else self.vector / math.sin(theta/2)
    @classmethod
    def fromangleaxis(cls, angle, axis):
        return Versor(axis.normalized() * math.sin(angle/2))
    def apply(self, point):
        '''
        >>> p1 = Vector3(1, 1, 1)
        >>> r = Vector3(1/11, 2/11, 4/11)
        >>> u1 = Versor.fromangleaxis(math.pi / 2, r.normalized())
        >>> p2 = u1.apply(p1)
        >>> s1 = p1 - p1.proj(u1.vector)
        >>> s2 = p2 - p1.proj(u1.vector)
        >>> abs( u1.vector.anglebetween(p1) - u1.vector.anglebetween(p2) ) < 1e-7
        True
        >>> abs( p1.norm() - p2.norm() ) < 1e-7
        True
        >>> abs( u1.vector * s2 ) < 1e-7
        True
        >>> abs( u1.vector * s1 ) < 1e-7
        True
        >>> abs( s2 * p1 ) < 1e-7
        True
        >>> abs( s2 * s1 ) < 1e-7
        True
        >>> abs( p2 * s1 ) < 1e-7
        True
        >>> u2 = Versor.fromangleaxis(math.pi / 3, r.normalized())
        >>> p3 = u2.apply(p1)
        >>> s3 = p3 - p1.proj(u2.vector)
        >>> abs( u2.vector.anglebetween(p1) - u2.vector.anglebetween(p3) ) < 1e-7
        True
        >>> abs( s1.anglebetween(s3) - math.pi / 3 ) < 1e-7
        True
        >>> abs( p1.norm() - p3.norm() ) < 1e-7
        True
        >>> abs( u2.vector * s3 ) < 1e-7
        True
        >>> abs( u2.vector * s1 ) < 1e-7
        True
        '''
        extent = point.proj(self.vector)
        deviation = point - extent
        lendev = deviation.norm()
        return (self.conjugate() * Quaternion(0, deviation * (1 / lendev)) * self
                ).vector * lendev + extent
    def rotate(self, point):
        return (self.conjugate() * Quaternion(0, point) * self).vector
    def rotation(self):
        i, j, k = Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1)
        return [self.rotate(i), self.rotate(j), self.rotate(k)]

if __name__ == '__main__':
    doctest.testmod()
