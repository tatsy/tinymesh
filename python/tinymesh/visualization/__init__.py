"""Optional mesh visualization. Install tinymesh[viz] to use it."""

try:
    from .viewer import Viewer
except ModuleNotFoundError as exc:
    if exc.name != 'matplotlib':
        raise
    raise ModuleNotFoundError(
        'Visualization requires Matplotlib. Install tinymesh[viz] or run uv sync --extra viz.',
        name='matplotlib',
    ) from exc


__all__ = ['Viewer']
