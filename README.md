TinyMesh
===

[![Codacy Badge](https://app.codacy.com/project/badge/Grade/2fd9a7f621e44654ad8b81bc38138662)](https://www.codacy.com/manual/tatsy/tinymesh?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tatsy/tinymesh&amp;utm_campaign=Badge_Grade)
[![Build Status](https://travis-ci.com/tatsy/tinymesh.svg?branch=master)](https://travis-ci.com/tatsy/tinymesh)
![Windows CI](https://github.com/tatsy/tinymesh/workflows/Windows%20CI/badge.svg)
![MacOS CI](https://github.com/tatsy/tinymesh/workflows/MacOS%20CI/badge.svg)
![Ubuntu CI](https://github.com/tatsy/tinymesh/workflows/Ubuntu%20CI/badge.svg)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)

> TinyMesh is a light-weight mesh processing library in C/C++

Modules
---

Here is the list of modules and reference papers for that.

*   **Smoothing**
  * Laplacian smoothing
  * Taubin smoothing [[Taubin 1995]](https://dl.acm.org/doi/10.1145/218380.218473)
  * Implicit fairing [[Desbrun 1999]](https://dl.acm.org/doi/10.1145/311535.311576)
*   **Remeshing**
  * Uniform triangulation [[Hoppe 1996]](https://dl.acm.org/doi/10.1145/237170.237216)
*   **Simplification**
  * Quadric error metrics (QEM) [[Garland and Heckbert 1997]](https://dl.acm.org/doi/10.1145/258734.258849)
*   **Hole filling**
  * Min-area hole filling [[Barequet and Sharir 1995]](https://www.sciencedirect.com/science/article/pii/016783969400011G?via%3Dihub)

Install
---

### Dependency

The module is tested its compilation using the following compilers.

*   **Windows** - Visual Studio 2019
*   **MacOS** - Apple Clang 11.0 (MacOS 10.15)
*   **Linux** - LLVM Clang 10.0, GNU C Compiler 9.3.0

Also, the module depends on Eigen for the linear algebra.

*   [Eigen 3.x](http://eigen.tuxfamily.org/index.php)

### C++ Library and Examples

You can build a shared library and all the examples by `CMake` with the following commands.

```shell
$ git clone https://github.com/tatsy/tinymesh.git --depth 1
$ cd tinymesh
$ git submodule update --init --recursive
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON ..
$ cmake --build . --config Release --parallel 2
```

#### Run examples (C++)

```shell
$ ./build/bin/example_simplify data/models/bunny.ply
```

### Python module

By the following shell command, you can build the module, which will be placed on the root of this repo.

```shell
$ python setup.py build_ext -i
```

Also, you can pass the CMake build arguments to `setup.py`.

```shell
$ python setup.py build_ext -i --cmake-build-args "-DEIGEN3_INCLUDE_DIR=[Your Eigen3 Location]"
```

Otherwise, you can build the module also using `CMake`.

```shell
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_MODULE=ON ..
$ cmake --build . --config Release --parallel 2
```

#### Run examples (Python)

```shell
$ python examples/python/fill_and_fair.py data/models/bunny.ply
```

Screen shots
---

| Input | Before | Simplified |
|:--------:|:--------:|:--------:|
|![](figures/dragon.jpg)|![](figures/dragon_before_simplify.jpg)|![](figures/dragon_after_simplify.jpg)|

Notice
---

The modules are not perfect and their process will fail for problematic meshes, e.g., with non-manifold faces. In such cases, you can fix the problem by repairing the mesh using [MeshFix](https://github.com/MarcoAttene/MeshFix-V2.1).

License
---

Mozilla Public License v2 (c) Tatsuya Yatagawa 2020
