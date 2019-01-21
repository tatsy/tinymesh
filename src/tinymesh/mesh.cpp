#define TINYMESH_API_EXPORT
#include "mesh.h"

#include <functional>
#include <unordered_map>

#include "vector.h"
#include "point.h"
#include "halfedge.h"
#include "face.h"
#include "tiny_obj_loader.h"

namespace std {

template <>
struct hash<tinymesh::Vector> {
    size_t operator()(const tinymesh::Vector &v) const {
        size_t h = 0;
        h = std::hash<double>()(v.x) ^ (h << 1);
        h = std::hash<double>()(v.y) ^ (h << 1);
        h = std::hash<double>()(v.z) ^ (h << 1);
        return h;
    }
};

}

namespace tinymesh {

Mesh::Mesh() {}

Mesh::Mesh(const std::string &filename) {
    load(filename);
}

void Mesh::load(const std::string &filename) {
    // Open
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> mats;
    std::string errMsg;
    bool success = tinyobj::LoadObj(&attrib, &shapes, &mats, &errMsg, filename.c_str());
    if (!errMsg.empty()) {
        Warning("%s", errMsg.c_str());
    }

    if (!success) {
        FatalError("Failed to load *.obj file: %s", filename.c_str());
    }

    // Clear
    m_pts.clear();
    m_hes.clear();
    m_faces.clear();
    m_indices.clear();

    // Traverse triangles
    std::unordered_map<Vector, uint32_t> uniqueVertices;
    for (const auto &shape : shapes) {
        for (const auto &index : shape.mesh.indices) {
            Vector v;
            if (index.vertex_index >= 0) {
                v = Vector(attrib.vertices[index.vertex_index * 3 + 0],
                           attrib.vertices[index.vertex_index * 3 + 1],
                           attrib.vertices[index.vertex_index * 3 + 2]);
            }

            if (uniqueVertices.count(v) == 0) {
                uniqueVertices[v] = static_cast<uint32_t>(m_pts.size());
                m_pts.emplace_back(v);
            }

            m_indices.push_back(uniqueVertices[v]);
        }
    }

    // Setup halfedge structure
    for (int i = 0; i < m_indices.size(); i += 3) {
        auto he0 = std::make_shared<Halfedge>();
        auto he1 = std::make_shared<Halfedge>();
        auto he2 = std::make_shared<Halfedge>();

        he0->m_start = m_pts[m_indices[i + 0]];
        he0->m_end = m_pts[m_indices[i + 1]];

        he1->m_start = m_pts[m_indices[i + 1]];
        he1->m_end = m_pts[m_indices[i + 2]];

        he2->m_start = m_pts[m_indices[i + 2]];
        he2->m_end = m_pts[m_indices[i + 0]];

        he0->m_next = he1;
        he1->m_next = he2;
        he2->m_next = he0;

        m_pts[m_indices[i + 0]]->m_he = he0;
        m_pts[m_indices[i + 1]]->m_he = he1;
        m_pts[m_indices[i + 2]]->m_he = he2;

        auto face = std::make_shared<Face>();
        face->m_he = he0;

        he0->m_face = face;
        he1->m_face = face;
        he2->m_face = face;

        m_hes.push_back(he0);
        m_hes.push_back(he1);
        m_hes.push_back(he2);
        m_faces.push_back(face);
    }

    // Set opposite half edges
    for (const auto &he0 : m_hes) {
        for (const auto &he1 : he0->m_end->out_he_iterator()) {
            
        }
    }
}

}  // namespace tinymesh
