v.0.1.1
---

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
