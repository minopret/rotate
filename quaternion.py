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
            return self * self
        # These other cases are conventional, but not expected to see use.
        if exp < 0:
            raise ValueError
        if exp == 0:
            return 1
        if exp == 1:
            return self
        even_exp = (exp - 1) / 2 if exp % 2 else exp / 2
        return (self if exp % 2 else 1) * (self * self)**even_exp
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
    def __new__(cls, *args):
        assert len(args) == 3
        return super().__new__(cls, *args)
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

class QuaQua(abc.ABC):
    '''
    "quasi-Quaternion", the interface of Quaternion.

    In particular, Versor implements QuaQua and does not subclass Quaternion.

    * Versor does not have a constructor __class__(scalar, vector).
    * Versor does not extend DivRingElement, the superclass of Quaternion.
    * No Versor is a "zero", an additive identity for instances of Versor.
    * Versor does not have a class variable Versor.ZERO.
    '''
    @abc.abstractmethod
    def __abs__(self):
        pass
    @abc.abstractmethod
    def __add__(self, other):
        pass
    @abc.abstractmethod
    def __bool__(self):
        pass
    @abc.abstractmethod
    def __eq__(self, other):
        pass
    @abc.abstractmethod
    def __hash__(self):
        pass
    @abc.abstractmethod
    def __mul__(self, other):
        pass
    @abc.abstractmethod
    def __neg__(self):
        pass
    @abc.abstractmethod
    def __pos__(self):
        pass
    @abc.abstractmethod
    def __sub__(self, other):
        pass
    @abc.abstractmethod
    def __truediv__(self, other):
        pass
    @abc.abstractmethod
    def conjugate(self):
        pass

class Quaternion(DivRingElement, QuaQua):
    '''
    Approximate an element of H, the incomputable set of all ("Hamilton's")
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
        self.scalar = scalar
        self.vector = Vector3(*vector) if vector else Vector3(0.0, 0.0, 0.0)
    def quaternion(self):
        return self
    def __abs__(self):
        'Returns the real distance from 0. Called for abs(self).'
        return math.sqrt(self.scalar**2 + self.vector**2)
        # return math.hypot(self.scalar, *self.vector) # Python 3.9
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
        if isinstance(other, QuaQua):
            return Quaternion(
                    scalar=self.scalar * other.scalar - self.vector * other.vector,
                    vector=(
                        other.vector * self.scalar +
                        self.vector * other.scalar +
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
        return f'Quaternion({scalar}, {vector1}, {vector2}, {vector3})'
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
            return Quaternion(scalar=self.scalar / square,
                    vector=self.vector / square) * other.conjugate()
        if isinstance(other, numbers.Number):
            return Quaternion(scalar=self.scalar / other,
                    vector=self.vector / other)
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

class Versor(MulGroup, QuaQua):
    '''
    Approximate a subset V of H, where H is the incomputable set of all
    quaternions in mathematics and V = {q for q in H if abs(q) == 1}.

    Having the invariant abs(q) == 1, the set V is closed under products,
    quotients, and additive inverses, but not under addition. V does
    not contain the zero quaternion, which is the zero element of
    quaternion addition.

    It is inevitable that some mathematical equalities for versors will
    not be computational equalities for class Versor. No representation
    for quaternions, or versors, or real numbers, can provide a
    representation for every such object that is both exact and small
    enough to use. We can only choose which such objects will have
    exact, small-enough representations. Consequently, all operations
    will have inputs for which the computed results are not the
    mathematically exact results. The computed results will never
    satisfy all of the properties of the mathematically exact results.

    This class represents versors as 4-tuples of binary64 values. An
    example of a versor that has no exact representation in that form
    is (0, 4/21, 5/21, 20/21). In binary64 math, that versor does not
    satisfy the invariant abs(q) == 1, if abs would be implemented as
    the sum of squared components.

    If abs would be implemented as math.hypot in Python 3.10, then
    that versor actually does satisfy abs(q) == 1. But the inevitable
    remains inevitable.
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
        self.scalar = (1.0 if positive else -1.0) * math.sqrt(1.0 - self.vector * self.vector)
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
            return Versor(vector=(
                        other.vector * self.scalar +
                        self.vector * other.scalar +
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
        return f'Versor({scalar}, {vector1}, {vector2}, {vector3})'
    def __sub__(self, other):
        'self - other'
        # Is a versor if, and only if, one of the following is the case:
        # *) ``other`` equals zero;
        # *) the perpendicular plane that bisects the line segment
        #    (self, -other) contains the origin.
        return self.quaternion() - other
        #return Quaternion(scalar=self.scalar - other.scalar,
        #        vector=(self.vector[i] - other.vector[i] for i in range(3)))
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
        Q =(cos θ/2, k sin θ/2)
        where
        * k is the axis, a unit;
        * θ is the angle.

        From Q=(s, x) a versor, obtain axis and angle:

        θ = 2 atan2(|x|, s)

        k = 0 if (θ = 0) else x / sin θ/2

        >>> theta = math.pi/2
        >>> abs(theta - Versor((math.cos(math.pi/4)*math.sin(theta/2),
        ...     math.sin(math.pi/4)*math.sin(theta/2), 0)).angle()) < 1e-14
        True
        '''
        return 2 * math.atan2(math.sqrt(self.vector * self.vector), self.scalar)
        # return 2 * math.atan2(math.hypot(*self.vector), self.scalar) # Python 3.9
    def axis(self):
        '''
        >>> print(tuple(round(x**2, 14) for x in Versor((.5, .5, .5)).axis()))
        (0.33333333333333, 0.33333333333333, 0.33333333333333)
        '''
        theta = self.angle()
        return 0 if theta == 0 else self.vector / math.sin(theta/2)
    @classmethod
    def fromangleaxis(cls, angle, axis):
        return Versor(axis * math.sin(angle/2))
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

if __name__ == '__main__':
    doctest.testmod()
