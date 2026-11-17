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
    tms.remesh_triangular,
    lambda m: tms.simplify_qem(m, m.num_faces() // 2),
    lambda m: tms.simplify_qem(m, m.num_faces() // 10),
]


@pytest.mark.parametrize('filename', filenames)
def test_hole_fill(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tms.Mesh(filename)

    try:
        mesh.fill_holes()
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), 'Mesh verification failed!'


@pytest.mark.parametrize('method, filename', product(methods, filenames))
def test_remesh_method(method, filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tms.Mesh(filename)
    mesh.fill_holes()

    try:
        method(mesh)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), 'Mesh verification failed!'
