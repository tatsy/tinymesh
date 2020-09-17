import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CommandMixin(object):
    user_options = [
        ('cmake-build-args=', None, 'CMake additonal build arguments')
    ]

    def initialize_options(self):
        super().initialize_options()
        self.cmake_build_args = None


class InstallCommand(CommandMixin, install):
    description = "Install modules"
    user_options = getattr(install, 'user_options', []) + CommandMixin.user_options


class CMakeBuild(CommandMixin, build_ext):
    description = "Build modules using CMake"
    user_options = getattr(build_ext, 'user_options', []) + CommandMixin.user_options

    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_args = ['-DPYTHON_MODULE_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-DBUILD_PYTHON_MODULE=ON']

        if self.cmake_build_args is not None:
            args = self.cmake_build_args.split('-D')
            args = [arg.strip() for arg in args]
            for arg in args:
                if arg != '':
                    cmake_args.append('-D' + arg)

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--config', 'Release']
            build_args += ['--parallel', '2']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--parallel', '2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)


setup(
    name='tinymesh',
    version='0.1.0',
    author='Tatsuya Yatagawa',
    author_email='tatsy.mail@gmail.com',
    description='TinyMesh is a light-weight mesh processing library',
    long_description='',
    license='MPL-v2',
    ext_modules=[CMakeExtension('tinymesh')],
    cmdclass={
        'install': InstallCommand,
        'build_ext': CMakeBuild
    }
)