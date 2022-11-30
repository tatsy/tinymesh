import os

import pytest

import tinymesh as tm

CWD = os.getcwd()

model_dir = 'data/models'
filenames = [
    'box.obj',
    'torus.obj',
    'fandisk.ply',
    'bunny_mini.ply',
]


@pytest.mark.parametrize("filename", filenames)
def test_hole_fill(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tm.Mesh(filename)

    try:
        mesh.fill_holes()
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"


@pytest.mark.parametrize("filename", filenames)
def test_remesh_triangular(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tm.Mesh(filename)
    mesh.fill_holes()

    try:
        tm.remesh_triangular(mesh)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"


@pytest.mark.parametrize("filename", filenames)
def test_simplify_qem(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tm.Mesh(filename)
    mesh.fill_holes()

    try:
        tm.simplify_qem(mesh, mesh.num_faces() // 10)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"
