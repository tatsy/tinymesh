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
