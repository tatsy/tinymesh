import os
from itertools import product

import pytest

import tinymesh as tms

CWD = os.getcwd()

model_dir = 'data/models'
filenames = [
    'torus.obj',
    'fandisk.ply',
    'bunny_mini.ply',
]

methods = [
    tms.smooth_laplacian,
    tms.smooth_taubin,
    tms.implicit_fairing,
    tms.denoise_normal_gaussian,
    tms.denoise_normal_bilateral,
    tms.denoise_l0_smooth,
]


@pytest.mark.parametrize("method, filename", product(methods, filenames))
def test_smooth_method(method, filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tms.Mesh(filename)
    mesh.fill_holes()

    try:
        method(mesh)
    except Exception:
        pytest.fail('Failed!')
