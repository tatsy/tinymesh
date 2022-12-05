v0.2.4
---
*   Delete `InHalfedgeIterator', and rename `OutHalfedgeIterator` as `HalfedgeIterator`.
*   Add const version of iterators.
*   Advancing front method [Zhao et al. 2007] is implemented.
*   Change dev-env manager to Poetry.

v0.2.3
---
*   Update some method names (e.g., `isStatic` -> `isLocked`)

v0.2.2
---
*   Change Python env manager to Pipenv.
*   Add methods computing approximated Laplacian-Beltrami operators (adjacent, cotangent, Belkin+2008).
*   Add IPython notebooks as Pythons examples rather than simple scripts.
*   Change traversal order of elements around a vertex to counter-clockwise order.
*   Add heat kernel signature [Sun et al. 2009].

v0.2.1
---
*   Minor bug fix for GCC v8.

v0.2.0
---
*   Add unit test for each functions (just to check run or not).
*   Add feature preservation in remeshing.
*   Include Eigen to repo as a submodule.
*   Update Python setup to use Poetry.
*   Add mesh denoising via bilateral filter [Zheng et al. 2011] and L0 smoothing [He and Schaefer 2013].
*   Update Python module build to use Python-side `pybind11` installation.
*   Add mesh denoising via normal Gaussian filter [Ohtake et al. 2001].
*   Implicit fairing `implicit_fair` is updated to its normalized version.
*   Hole filling `hole_fill` is updated to take the upper bound of allowed dihedral angle.
*   Update `setup.py` to support the operation `install`.

v0.1.0
---
*   First release
