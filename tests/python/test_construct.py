import os
import unittest

import numpy as np
from nose2.tools import params

from tinymesh import Mesh

CWD = os.getcwd()

model_dir = 'data/models'
filenames = [
    'box.obj',
    'torus.obj',
    'fandisk.ply',
    'bunny.ply',
]


class TestConstruct(unittest.TestCase):
    @params(*filenames)
    def test_mesh_io(self, filename):
        filename = os.path.join(CWD, model_dir, filename)
        base, ext = os.path.splitext(filename)
        try:
            mesh = Mesh(filename)
        except Exception:
            self.fail('Failed to load mesh!')

        self.assertGreater(mesh.num_vertices(), 0)
        self.assertGreater(mesh.num_faces(), 0)
        self.assertGreater(mesh.num_edges(), 0)
        self.assertGreater(mesh.num_halfedges(), 0)

        try:
            mesh.save(base + '_test' + ext)
        except Exception:
            self.fail('Failed to save mesh!')

    def test_tetrahedron(self):
        # Construct a simple tetrahedron
        a = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='double')
        b = np.asarray([0, 1, 2, 0, 2, 3, 0, 3, 1, 3, 2, 1], dtype='uint32')
        try:
            _ = Mesh(a, b)
        except Exception:
            self.fail('Mesh construction from vertices/indices was failed!')
