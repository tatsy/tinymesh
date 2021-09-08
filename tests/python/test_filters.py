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


class TestFilters(unittest.TestCase):
    @params(*filenames)
    def test_smooth_laplacian(self, filename):
        filename = os.path.join(CWD, model_dir, filename)
        mesh = tm.Mesh(filename)
        tm.hole_fill(mesh)

        try:
            tm.smooth_laplacian(mesh)
        except Exception:
            self.fail('Failed!')

    @params(*filenames)
    def test_smooth_taubin(self, filename):
        filename = os.path.join(CWD, model_dir, filename)
        mesh = tm.Mesh(filename)
        tm.hole_fill(mesh)

        try:
            tm.smooth_taubin(mesh)
        except Exception:
            self.fail('Failed!')

    @params(*filenames)
    def test_implicit_fairing(self, filename):
        filename = os.path.join(CWD, model_dir, filename)
        mesh = tm.Mesh(filename)
        tm.hole_fill(mesh)

        try:
            tm.implicit_fairing(mesh)
        except Exception:
            self.fail('Failed!')

    @params(*filenames)
    def test_denoise_normal_gaussian(self, filename):
        filename = os.path.join(CWD, model_dir, filename)
        mesh = tm.Mesh(filename)
        tm.hole_fill(mesh)

        try:
            tm.denoise_normal_gaussian(mesh)
        except Exception:
            self.fail('Failed!')

    @params(*filenames)
    def test_denoise_normal_bilateral(self, filename):
        filename = os.path.join(CWD, model_dir, filename)
        mesh = tm.Mesh(filename)
        tm.hole_fill(mesh)

        try:
            tm.denoise_normal_bilateral(mesh)
        except Exception:
            self.fail('Failed!')

    @params(*filenames)
    def test_denoise_l0_smooth(self, filename):
        filename = os.path.join(CWD, model_dir, filename)
        mesh = tm.Mesh(filename)
        tm.hole_fill(mesh)

        try:
            tm.denoise_l0_smooth(mesh)
        except Exception:
            self.fail('Failed!')
