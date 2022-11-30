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
def test_smooth_laplacian(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tm.Mesh(filename)
    mesh.fill_holes()

    try:
        tm.smooth_laplacian(mesh)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"


@pytest.mark.parametrize("filename", filenames)
def test_smooth_taubin(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tm.Mesh(filename)
    mesh.fill_holes()

    try:
        tm.smooth_taubin(mesh)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"


@pytest.mark.parametrize("filename", filenames)
def test_implicit_fairing(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tm.Mesh(filename)
    mesh.fill_holes()

    try:
        tm.implicit_fairing(mesh)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"


@pytest.mark.parametrize("filename", filenames)
def test_denoise_normal_gaussian(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tm.Mesh(filename)
    mesh.fill_holes()

    try:
        tm.denoise_normal_gaussian(mesh)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"


@pytest.mark.parametrize("filename", filenames)
def test_denoise_normal_bilateral(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tm.Mesh(filename)
    mesh.fill_holes()

    try:
        tm.denoise_normal_bilateral(mesh)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"


@pytest.mark.parametrize("filename", filenames)
def test_denoise_l0_smooth(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tm.Mesh(filename)
    mesh.fill_holes()

    try:
        tm.denoise_l0_smooth(mesh)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"
