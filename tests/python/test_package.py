import sys
import subprocess
from textwrap import dedent
from importlib import metadata, resources

import pytest

import tinymesh
from tinymesh import _core


def test_existing_api_is_reexported():
    for name in dir(_core):
        if not name.startswith('_'):
            assert getattr(tinymesh, name) is getattr(_core, name)
    assert tinymesh.__version__ == metadata.version('tinymesh')


def test_type_information_is_installed():
    package = resources.files('tinymesh')
    assert package.joinpath('py.typed').is_file()
    assert package.joinpath('_core.pyi').is_file()


@pytest.mark.parametrize('import_error', ['ModuleNotFoundError', 'ImportError'])
def test_core_import_without_visualization_dependencies(tmp_path, import_error):
    result = subprocess.run(
        [
            sys.executable,
            '-c',
            dedent("""\
                import sys
                import builtins
                from importlib.abc import MetaPathFinder

                error_type = getattr(builtins, sys.argv[1])

                class BlockVisualization(MetaPathFinder):
                    def find_spec(self, fullname, path=None, target=None):
                        if fullname.split('.')[0] in {'matplotlib', 'open3d'}:
                            raise error_type(fullname)

                sys.meta_path.insert(0, BlockVisualization())
                import tinymesh
                assert tinymesh.Mesh().num_vertices() == 0
                assert 'tinymesh.visualization' not in sys.modules

                try:
                    import tinymesh.visualization
                except ImportError as exc:
                    assert 'uv sync --extra viz' in str(exc)
                    assert 'pip install "tinymesh[viz]"' in str(exc)
                    assert isinstance(exc.__cause__, error_type)
                else:
                    raise AssertionError('Expected a missing optional dependency')
            """),
            import_error,
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
