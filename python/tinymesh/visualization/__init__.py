try:
    from .viewer import Viewer
except ImportError as exc:
    raise ImportError('Visualization could not be imported. Install the module with `pip install tinymesh`.') from exc


__all__ = ['Viewer']
