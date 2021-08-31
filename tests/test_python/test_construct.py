try:
    import unittest2 as unittest
except ImportError:
    import unittest

import os

import numpy as np

from tinymesh import Mesh

CWD = os.getcwd()


class TestConstruct(unittest.TestCase):
    def test_io(self):
        filename = os.path.join(CWD, 'data/models/bunny.ply')
        base, ext = os.path.splitext(filename)
        mesh = Mesh(filename)
        mesh.save(base + '_test' + ext)

    def test_tetrahedron(self):
        # Construct a simple tetrahedron
        a = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='double')
        b = np.asarray([0, 1, 2, 0, 2, 3, 0, 3, 1, 3, 2, 1], dtype='uint32')
        try:
            _ = Mesh(a, b)
        except Exception:
            self.fail('Mesh construction from vertices/indices was failed!')
