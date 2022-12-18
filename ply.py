# References:
# http://paulbourke.net/dataformats/ply/
# https://web.archive.org/web/20161204152348/http://www.dcs.ed.ac.uk/teaching/cs4/www/graphics/Web/ply.html
# https://www.loc.gov/preservation/digital/formats/fdd/fdd000501.shtml

import abc
import textwrap

class PLY(abc.ABC):
    HEADER_1_0 = textwrap.dedent('''\
    ply
    format {fmt} 1.0
    element vertex {vertex_count}
    property float x
    property float y
    property float z
    element face {face_count}
    property list uchar int vertex_index
    end_header
    ''')
    def __init__(self, fmt, vertices: list, faces: list):
        self.fmt = fmt
        self.vertices = vertices
        self.faces = faces
    def header(self):
        cls = self.__class__
        return cls.HEADER_1_0.format(fmt=self.fmt,
                vertex_count=len(self.vertices), face_count=len(self.faces))
    def vertex_list(self):
        return '\n'.join(' '.join(vertex) for vertex in self.vertices)
    def face_list(self):
        return '\n'.join(' '.join(face) for face in self.faces)
    def __str__(self):
        return self.header() + '\n'.join([self.vertex_list(), self.face_list()])
