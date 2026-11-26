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


@pytest.mark.parametrize('filename', filenames)
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


@pytest.mark.parametrize('filename', ['not_exist.ply', 'invalid_ext.obj'])
def test_mesh_invalid_io(filename):
    with pytest.raises(RuntimeError):
        _ = Mesh(filename)


@pytest.mark.parametrize('filename', ['unsupported.stl', 'missing.ply'])
def test_mesh_load_failure_preserves_existing_mesh(tmp_path, filename):
    mesh = Mesh(
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
        [0, 1, 2],
    )
    mesh.vertex(0).lock()
    vertices = np.array(mesh.get_vertices())
    indices = mesh.get_vertex_indices()

    with pytest.raises(RuntimeError):
        mesh.load(str(tmp_path / filename))

    np.testing.assert_array_equal(mesh.get_vertices(), vertices)
    assert mesh.get_vertex_indices() == indices
    assert mesh.vertex(0).is_locked()
    assert mesh.verify()


def test_tetrahedron():
    # Construct a simple tetrahedron
    a = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='double')
    b = np.asarray([0, 1, 2, 0, 2, 3, 0, 3, 1, 3, 2, 1], dtype='uint32')
    try:
        _ = Mesh(a, b)
    except Exception:
        pytest.fail('Mesh construction from vertices/indices was failed!')


def test_mesh_clone_preserves_topology_and_lock_state():
    mesh = Mesh(
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
        [0, 1, 2],
    )
    mesh.vertex(0).lock()
    mesh.vertex(1).set_pos(mesh.vertex(0).pos())

    copy = mesh.clone()

    assert copy.num_vertices() == 3
    assert copy.get_vertex_indices() == mesh.get_vertex_indices()
    assert copy.verify()
    assert copy.vertex(0).is_locked()

    copy.vertex(0).unlock()
    assert mesh.vertex(0).is_locked()
