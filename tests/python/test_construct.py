import os

import numpy as np
import pytest

from tinymesh import Mesh

CWD = os.getcwd()

model_dir = 'data/models'
filenames = [
    'torus.obj',
    'fandisk.ply',
    'bunny.ply',
]


@pytest.mark.parametrize("filename", filenames)
def test_mesh_io(filename):
    filename = os.path.join(CWD, model_dir, filename)
    base, ext = os.path.splitext(filename)

    try:
        mesh = Mesh(filename)
    except Exception:
        pytest.fail('Failed to load mesh!')

    assert mesh.num_vertices() > 0
    assert mesh.num_faces() > 0
    assert mesh.num_edges() > 0
    assert mesh.num_halfedges() > 0

    try:
        mesh.save(base + '_test' + ext)
    except Exception:
        pytest.fail('Failed to save mesh!')


@pytest.mark.parametrize("filename", ["not_exist.ply", "invalid_ext.obj"])
def test_mesh_invalid_io(filename):
    with pytest.raises(RuntimeError):
        _ = Mesh(filename)


def test_tetrahedron():
    # Construct a simple tetrahedron
    a = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='double')
    b = np.asarray([0, 1, 2, 0, 2, 3, 0, 3, 1, 3, 2, 1], dtype='uint32')
    try:
        _ = Mesh(a, b)
    except Exception:
        pytest.fail('Mesh construction from vertices/indices was failed!')
