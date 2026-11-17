"""Optional mesh visualization. Install tinymesh[viz] to use it."""

try:
    from .viewer import Viewer
except ImportError as exc:
    raise ImportError(
        'Visualization could not be imported. Install the optional dependencies with `pip install tinymesh`.'
    ) from exc


__all__ = ['Viewer']
