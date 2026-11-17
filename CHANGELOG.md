# CHANGELOG

## v0.2.5

### Added

* Add an opt-in C++ mesh construction benchmark and a Python runner with generated/file inputs, validation, memory reporting, and CSV/JSON output.
* Ship Python type stubs and `py.typed`, and expose the packaged viewer through `tinymesh.visualization.Viewer`.

### Changed

* Replace native dependency submodules with pinned, SHA256-verified CMake FetchContent dependencies.
* Require CMake 3.24 or later and a C++20 compiler with `std::format` support. Preserve standalone C++ library and example builds without Python.
* Move Python builds to scikit-build-core and development environment management to uv. Install the native extension as `tinymesh._core` while preserving the top-level Python API.
* Support Python 3.10–3.13, with NumPy, SciPy, and Matplotlib as runtime dependencies. Keep notebook tooling in the `examples` dependency group.
* Move the example-local viewer into the Python package and update notebooks to import it.
* Use reserved hash tables for halfedge construction lookups.
* Consolidate CI for C++ builds, examples, and tests, and Python distribution builds and installed-wheel tests on Linux, macOS, and Windows.

### Fixed

* Validate triangle index counts, index ranges, duplicated edges, inconsistent winding, non-manifold edges, and disconnected vertex fans during mesh construction.
* Clear stored indices when loading or reconstructing a mesh, and explicitly identify boundary faces.
* Correct boundary handling in vertex normals, Voronoi areas, cotangent weights, and mean face area. Return zero curvature at boundary vertices and correct the mean curvature scaling.
* Read PLY float32/float64 coordinates and uint8/uint16/int32/uint32 triangle indices, validate their ranges, and report unsupported data types and file errors.
* Correct vector-scalar subtraction, component-wise bounds operations, triangle closest-point queries, and rotation between antiparallel vectors.
* Build a shared library when `TINYMESH_BUILD_SHARED=ON`, with the appropriate Windows export/import definitions.
* Handle absent source extensions and filenames containing whitespace in the formatting script.

### C++ API compatibility

* Rename `Vertex::volonoiArea()` to `Vertex::voronoiArea()` (including the face overload).
* Disable `Mesh` copy construction and copy assignment. Use `clone()` for an independent copy or move the mesh.

## v0.2.4

* Update CI scripts to build wheels for 5 latest Python versions.

## v0.2.3

* Update some method names (e.g., `isStatic` -> `isLocked`)
* Delete `InHalfedgeIterator`, and rename `OutHalfedgeIterator` as `HalfedgeIterator`.
* Add const version of iterators.
* Advancing front method [Zhao et al. 2007] is implemented.
* Change dev-env manager to Poetry.
* Add a method to compute principal curvatures (and thier directions) [Rusinkiewicz 2004].
* Add a method of context-based coherent surface completion [Harary et al. 2016].
* Add own print function to support both C++/Python standard outputs.

## v0.2.2

* Change Python env manager to Pipenv.
* Add methods computing approximated Laplacian-Beltrami operators (adjacent, cotangent, Belkin+2008).
* Add IPython notebooks as Pythons examples rather than simple scripts.
* Change traversal order of elements around a vertex to counter-clockwise order.
* Add heat kernel signature [Sun et al. 2009].

## v0.2.1

* Minor bug fix for GCC v8.

## v0.2.0

* Add unit test for each functions (just to check run or not).
* Add feature preservation in remeshing.
* Include Eigen to repo as a submodule.
* Update Python setup to use Poetry.
* Add mesh denoising via bilateral filter [Zheng et al. 2011] and L0 smoothing [He and Schaefer 2013].
* Update Python module build to use Python-side `pybind11` installation.
* Add mesh denoising via normal Gaussian filter [Ohtake et al. 2001].
* Implicit fairing `implicit_fair` is updated to its normalized version.
* Hole filling `hole_fill` is updated to take the upper bound of allowed dihedral angle.
* Update `setup.py` to support the operation `install`.

## v0.1.0

* First release
