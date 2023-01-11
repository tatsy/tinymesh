import os
from functools import partial
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
    partial(tms.get_mesh_laplacian, type=tms.MeshLaplace.ADJACENT),
    partial(tms.get_mesh_laplacian, type=tms.MeshLaplace.COTANGENT),
    partial(tms.get_mesh_laplacian, type=tms.MeshLaplace.BELKIN08),
    tms.get_principal_curvatures,
    tms.get_principal_curvatures_with_derivatives,
    tms.get_feature_line_field,
]


@pytest.mark.parametrize("filename", filenames)
def test_hks(filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tms.Mesh(filename)
    mesh.fill_holes()

    L = tms.get_mesh_laplacian(mesh, tms.MeshLaplace.COTANGENT)
    K = min(L.shape[0], 200)

    try:
        tms.get_heat_kernel_signatures(L, K)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"


@pytest.mark.parametrize("method, filename", product(methods, filenames))
def test_ops_method(method, filename):
    filename = os.path.join(CWD, model_dir, filename)
    mesh = tms.Mesh(filename)
    mesh.fill_holes()

    try:
        method(mesh)
    except Exception:
        pytest.fail('Failed!')

    assert mesh.verify(), "Mesh verification failed!"
