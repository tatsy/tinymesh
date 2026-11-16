import importlib.util
from pathlib import Path

import numpy as np
import pytest


matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg', force=True)

import matplotlib.pyplot as plt  # noqa: E402
from tinymesh.visualization import Viewer  # noqa: E402

from tinymesh import Mesh  # noqa: E402


def test_viewer_renders_and_saves_mesh(tmp_path, monkeypatch):
    vertices = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    indices = np.array([0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3], dtype=np.uint32)
    mesh = Mesh(vertices, indices)
    viewer = Viewer(400, 300)
    viewer.lookat(3, 4, 5, 0, 0, 0, 0, 1, 0)
    viewer.perspective(45, viewer.aspect, 0.1, 100)
    monkeypatch.setattr(plt, 'show', lambda: None)
    output = tmp_path / 'mesh.png'

    try:
        viewer.mesh_visualization(mesh, save=True, filename=str(output), wireframe=True)
        assert output.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
        collection = plt.gcf().axes[0].collections[0]
        assert len(collection.get_paths()) == 4
    finally:
        plt.close('all')


def test_example_viewer_import_remains_compatible():
    path = Path(__file__).resolve().parents[2] / 'examples/python/viewer.py'
    spec = importlib.util.spec_from_file_location('example_viewer', path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.Viewer is Viewer
