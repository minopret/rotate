#!/usr/bin/env python3
import abc
import collections
import doctest
import math

class AddableInitializer:
    def __init__(self, addgroup):
        self.addgroup = addgroup
    def __get__(self, instance, owner=None):
        owner.ZERO = self.addgroup()
        return owner.ZERO

class AddableElement(abc.ABC):
    def __init_subclass__(cls, addgroup=None, **kwargs):
        super().__init_subclass__(**kwargs)
        if addgroup is None:
            addgroup = cls
        cls.addgroup = addgroup
        cls.ZERO = AddableInitializer(addgroup)
    @abc.abstractmethod
    def __sub__(self, other):
        pass
    def __eq__(self, other):
        return self - other == self.__class__.ZERO
    def __bool__(self):
        'True if self != 0. Called for bool(self).'
        return self != self.__class__.ZERO
    def __neg__(self):
        '-self'
        return self.__class__.ZERO - self
    def __pos__(self):
        '+self'
        return self
    def __add__(self, other):
        'self + other'
        return self - (-other)

class AddableSequence(collections.abc.Sequence, AddableElement):
    def __sub__(self, other):
        return self.__class__( *(
                    addend - subtrahend
                    for addend, subtrahend in zip(self, other)
                ) )
    def __eq__(self, other):
        return isinstance(other, collections.abc.Sequence) and tuple(self) == tuple(other)

class Vector(AddableSequence):
    @staticmethod
    @abc.abstractmethod
    def indices():
        'Return indices of components'
    @abc.abstractmethod
    def __len__(self):
        'Number of components'
    @abc.abstractmethod
    def _getitem(self, key):
        'Retrieve component by index number'
    def __getitem__(self, key):
        if isinstance(key, slice):
            indices = self.__class__.indices()
            return tuple(self._getitem(key) for index in indices[key])
        if isinstance(key, int):
            return self._getitem(key)
        syntax_warning = TypeError(
                f'Vector indices must be integers or slices, not {type(key).__name__}')
        raise syntax_warning
    def __hash__(self):
        return hash(tuple(self))
    def __mul__(self, other):
        if isinstance(other, Vector):
            return sum(ours * theirs for ours, theirs in zip(self, other))
        return self.__class__(*(ours * other for ours in self))
    def __pow__(self, exp):
        if exp == 2:
            # The likely case.
            return self * self
        if exp >= 0 and isinstance(exp, int):
            return (1 if exp == 0 else
                    self if exp == 1 else
                    self * (self * self)**((exp - 1) / 2) if exp % 2 else
                    (self * self)**(exp / 2))
        raise ValueError()
    def __repr__(self):
        return f'{self.__class__.__name__}(' + ', '.join(repr(elem) for elem in self) + ')'
    def __truediv__(self, other):
        return self * (1 / other)
    def proj(self, other):
        '''
        Projection of self ``a`` in direction of another vector ``b``: ``proj_b a``
        If ``a`` is not zero, then the projection is zero if, and only if,
        the two are orthogonal.
        '''
        return other * (self * other) / (other * other)
    def __abs__(self):
        return math.hypot(*self)
    def normalized(self):
        return self / abs(self)
    def anglebetween(self, other):
        return math.acos(abs(other.normalized() * self.normalized()))
    def line_approaches_line(self, dir1, vector2, dir2):
        vector1 = self
        '''
        What scalar, multiplied by dir1 and added to vector1, gives the
        closest point to the line through vector2 that has direction dir2?
        '''
        normal = dir1.cross(dir2)
        if normal == dir1.__class__():
            # The two are parallel lines.
            return math.inf

        # Define proj2 as vector2 projected onto vector1's plane.
        proj2 = vector2 + (vector1.proj(normal) - vector2.proj(normal))

        # Now find the intersection of the lines through vector1 and proj2
        # with respective directions dir1 and dir2, as the parameter t
        # for the formula vector1 + t dir1.

        # https://math.stackexchange.com/a/271366
        return abs(dir2.cross(proj2 - vector1)) / abs(normal)

class Vector2(Vector):
    __slots__ = ('x', 'y', )
    def __init__(self, a=0.0, b=0.0):
        self.x, self.y = a, b
    @staticmethod
    def indices():
        return (0, 1)
    def __len__(self):
        return 2
    def _getitem(self, key):
        indices = self.__class__.indices()
        return (self.x, self.y)[indices[key]]
    def ccw(self, pt0):
        '''
        Returns whether the plane angle from the origin to ``pt0`` to
        ``self`` is counterclockwise from above.
        '''
        pt1 = self
        return pt1.x * pt0.y < pt0.x * pt1.y

class Vector3(Vector):
    __slots__ = ('x', 'y', 'z', )
    def __init__(self, /, a=0.0, b=0.0, c=0.0):
        self.x, self.y, self.z = a, b, c
    @staticmethod
    def indices():
        return (0, 1, 2)
    def __len__(self):
        return 3
    def _getitem(self, key):
        indices = self.__class__.indices()
        return (self.x, self.y, self.z)[indices[key]]
    def cross(self, other):
        return Vector3(self[1] * other[2] - other[1] * self[2],
                self[2] * other[0] - other[2] * self[0],
                self[0] * other[1] - other[0] * self[1])
    def azimuth(self):
        basisi = Vector3(1)
        projij = Vector3(*self[:2])
        return basisi.anglebetween(projij) * (+1 if self[1] > 0 else -1)
    def elevation(self):
        projij = Vector3(*self[:2])
        return self.anglebetween(projij) * (+1 if self[2] > 0 else -1)

if __name__ == '__main__':
    doctest.testmod()
