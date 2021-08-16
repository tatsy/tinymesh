#include "tinymesh/tinymesh.h"
using namespace tinymesh;

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

// type caster: Matrix <-> NumPy-array
namespace pybind11 {
namespace detail {

template <typename Float, int Dims>
struct type_caster<Vec<Float, Dims>> {
public:
    using Type = Vec<Float, Dims>;
    PYBIND11_TYPE_CASTER(Type, _("Vec<Float, Dims>"));

    // Conversion part 1 (Python -> C++)
    bool load(py::handle src, bool convert) {
        if (!convert && !py::array_t<Float>::check_(src)) return false;

        auto buf = py::array_t<Float, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf) return false;

        auto dims = buf.ndim();
        if (dims != 1) return false;
        auto elems = buf.shape()[0];
        if (elems != Dims) return false;

        value = Vec<Float, Dims>();
        const Float *ptr = buf.data();
        for (int d = 0; d < Dims; d++) {
            value[d] = ptr[d];
        }

        return true;
    }

    //Conversion part 2 (C++ -> Python)
    static py::handle cast(const Vec<Float, Dims> &src, py::return_value_policy policy, py::handle parent) {
        std::vector<size_t> shape  (1, Dims);
        std::vector<size_t> strides(1, sizeof(Float));
        py::array a(std::move(shape), std::move(strides), (Float*)&src);
        return a.release();
    }
};

}  // namespace detail
}  // namespace pybind11

PYBIND11_MODULE(tinymesh, m) {
    py::class_<Mesh>(m, "Mesh")
        .def(py::init<>())
        .def(py::init<const std::string &>())
        .def(py::init<const std::vector<Vec3>&, const std::vector<uint32_t>&>())
        .def("load", &Mesh::load)
        .def("save", &Mesh::save)
        .def("vertex", &Mesh::vertex, py::return_value_policy::reference)
        .def("face", &Mesh::face, py::return_value_policy::reference)
        .def("get_vertices", &Mesh::getVertices)
        .def("get_vertex_indices", &Mesh::getVertexIndices)
        .def("num_vertices", &Mesh::num_vertices)
        .def("num_faces", &Mesh::num_faces);

    py::class_<Vertex, std::shared_ptr<Vertex>>(m, "Vertex")
        .def(py::init<>())
        .def("is_boundary", &Vertex::isBoundary)
        .def("is_static", &Vertex::isStatic);

    py::class_<Face, std::shared_ptr<Face>>(m, "Face")
        .def(py::init<>())
        .def("is_boundary", &Face::isBoundary)
        .def("is_static", &Face::isStatic)
        .def("set_is_static", &Face::setIsStatic);


    /*** Smoothing ***/

    m.def("smooth_laplacian", &smoothLaplacian,
          "Laplacian smoothing",
          py::arg("mesh"),
          py::arg("epsilon") = 1.0,
          py::arg("cotangent_weight") = false,
          py::arg("iterations") = 3);

    m.def("smooth_taubin", &smoothTaubin,
          "Taubin smoothing",
          py::arg("mesh"),
          py::arg("shrink") = 1.0,
          py::arg("inflate") = 1.0,
          py::arg("iterations") = 3);

    m.def("implicit_fairing", &implicitFairing,
          "Implicit fairing",
          py::arg("mesh"),
          py::arg("epsilon") = 1.0,
          py::arg("iterations") = 1);

    /*** Remesh ***/

    m.def("remesh_triangular", &remeshTriangular,
          "Triangle remeshing",
          py::arg("mesh"),
          py::arg("short_length") = 0.8,
          py::arg("long_length") = 1.333,
          py::arg("angle_thresh") = 0.2,
          py::arg("iterations") = 5);

    /*** Simplification ***/

    m.def("simplify_qem", &simplifyQEM,
          "QEM-based simplification",
          py::arg("mesh"),
          py::arg("n_triangles"));

    /*** Hole filling ***/

    m.def("hole_fill", &holeFill,
          "Max-area hole filling",
          py::arg("mesh"),
          py::arg("dihedral_bound") = Pi);
}
