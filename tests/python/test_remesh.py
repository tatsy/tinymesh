import os
import unittest

from nose2.tools import params

import tinymesh as tm

CWD = os.getcwd()

model_dir = 'data/models'
filenames = [
    'box.obj',
    'torus.obj',
    'fandisk.ply',
    'bunny_mini.ply',
]


class TestRemesh(unittest.TestCase):

    @params(*filenames)
    def test_hole_fill(self, filename):
        filename = os.path.join(CWD, model_dir, filename)
        mesh = tm.Mesh(filename)

        try:
            mesh.fill_holes()
        except Exception:
            self.fail('Failed!')

        self.assertTrue(mesh.verify())

    @params(*filenames)
    def test_remesh_triangular(self, filename):
        filename = os.path.join(CWD, model_dir, filename)
        mesh = tm.Mesh(filename)
        mesh.fill_holes()

        try:
            tm.remesh_triangular(mesh)
        except Exception:
            self.fail('Failed!')

    @params(*filenames)
    def test_simplify_qem(self, filename):
        filename = os.path.join(CWD, model_dir, filename)
        mesh = tm.Mesh(filename)
        mesh.fill_holes()

        try:
            tm.simplify_qem(mesh, mesh.num_faces() // 10)
        except Exception:
            self.fail('Failed!')
