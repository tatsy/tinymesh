import re
import pathlib
import platform

from setuptools.errors import CCompilerError, PackageDiscoveryError
from setuptools.command.build_ext import build_ext

try:
    from pybind11.setup_helpers import Pybind11Extension
except ImportError:
    from setuptools import Extension as Pybind11Extension


class TinyMeshBuildExt(build_ext):
    def run(self):
        try:
            super(TinyMeshBuildExt, self).run()
        except FileNotFoundError:
            raise Exception("File not found. Could not compile C extension")

    def build_extension(self, ext):
        # common settings
        for e in self.extensions:
            e.define_macros.extend(
                [
                    ("TINYMESH_PYTHON_MODULE", 1),
                ]
            )

        # OS specific settings
        if platform.system() == "Darwin":
            for e in self.extensions:
                e.extra_compile_args.extend(
                    [
                        "-mmacosx-version-min=10.15",
                    ]
                )

        elif platform.system() == "Linux":
            for e in self.extensions:
                e.extra_compile_args.extend(
                    [
                        "-fopenmp",
                    ]
                )
                e.extra_link_args.extend(
                    [
                        "-fopenmp",
                        "-lstdc++fs",
                    ]
                )

        # compiler specific settings
        if self.compiler.compiler_type == "unix":
            for e in self.extensions:
                e.extra_compile_args.extend(
                    [
                        "-std=c++17",
                    ]
                )

        elif self.compiler.compiler_type == "msvc":
            for e in self.extensions:
                e.extra_compile_args.extend(["/utf-8", "/openmp"])
                e.define_macros.extend(
                    [
                        ("_CRT_SECURE_NO_WARNINGS", 1),
                        ("_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING", 1),
                    ]
                )

        # building
        try:
            super(TinyMeshBuildExt, self).build_extension(ext)
        except (CCompilerError, PackageDiscoveryError, ValueError):
            raise Exception("Could not compile C extension")


exclude = ["src/tinymesh/ext"]
sources = pathlib.Path().glob("src/**/*.cpp")
sources = [str(path).replace("\\", "/") for path in sources]
sources = [path for path in sources if all([not re.search(e, path) for e in exclude])]

include = [
    "src/tinymesh",
    "src/tinymesh/ext/tinyobjloader",
    "src/tinymesh/ext/tinyply/source",
    "src/tinymesh/ext/eigen",
    "src/tinymesh/ext/spectra/include",
]

ext_modules = [
    Pybind11Extension(
        "tinymesh",
        sources,
        language="c++",
        include_dirs=include,
    )
]


def build(setup_kwargs):
    setup_kwargs.update(
        {
            "ext_modules": ext_modules,
            "cmdclass": {"build_ext": TinyMeshBuildExt},
        }
    )
