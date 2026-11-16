"""Optional mesh visualization. Install tinymesh[viz] to use it."""

try:
    from .viewer import Viewer
except ImportError as exc:
    raise ImportError(
        'Visualization could not be imported. Install the optional dependencies with '
        '`uv sync --extra viz` or `pip install "tinymesh[viz]"`.'
    ) from exc


__all__ = ['Viewer']
