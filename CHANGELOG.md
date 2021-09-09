v 0.2.1
---
*  Minor bug fix for GCC v8.

v.0.2.0
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

* First release
