#define TINYMESH_API_EXPORT
#include "mesh.h"

#include <iostream>
#include <fstream>
#include <functional>
#include <unordered_map>

#include "vector.h"
#include "vertex.h"
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

// ----------
// Mesh
// ----------

Mesh::Mesh() {}

Mesh::Mesh(const std::string &filename) {
    load(filename);
}

void Mesh::load(const std::string &filename) {
    // Open
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> mats;
    std::string warnMsg;
    std::string errMsg;
    bool success = tinyobj::LoadObj(&attrib, &shapes, &mats, &warnMsg, &errMsg, filename.c_str());
    if (!errMsg.empty()) {
        Warning("%s", errMsg.c_str());
    }

    if (!success) {
        FatalError("Failed to load *.obj file: %s", filename.c_str());
    }

    // Clear
    m_verts.clear();
    m_hes.clear();
    //m_faces.clear();
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
                uniqueVertices[v] = static_cast<uint32_t>(m_verts.size());
                m_verts.push_back(std::make_shared<Vertex>(v));
            }

            m_indices.push_back(uniqueVertices[v]);
        }
    }

    // Setup halfedge structure
    std::unordered_map<Vertex*, std::vector<Halfedge*>> perVertexHEs;
    for (int i = 0; i < m_indices.size(); i += 3) {
        auto he0 = std::make_shared<Halfedge>();
        auto he1 = std::make_shared<Halfedge>();
        auto he2 = std::make_shared<Halfedge>();

        he0->m_start = m_verts[m_indices[i + 0]].get();
        he0->m_end = m_verts[m_indices[i + 1]].get();

        he1->m_start = m_verts[m_indices[i + 1]].get();
        he1->m_end = m_verts[m_indices[i + 2]].get();

        he2->m_start = m_verts[m_indices[i + 2]].get();
        he2->m_end = m_verts[m_indices[i + 0]].get();

        he0->m_next = he1.get();
        he1->m_next = he2.get();
        he2->m_next = he0.get();

        he0->m_prev = he2.get();
        he1->m_prev = he0.get();
        he2->m_prev = he1.get();

        m_verts[m_indices[i + 0]]->m_he = he0.get();
        m_verts[m_indices[i + 1]]->m_he = he1.get();
        m_verts[m_indices[i + 2]]->m_he = he2.get();

        auto face = std::make_shared<Face>();
        face->m_he = he0.get();

        he0->m_face = face.get();
        he1->m_face = face.get();
        he2->m_face = face.get();

        m_verts[m_indices[i + 0]]->m_face = face.get();
        m_verts[m_indices[i + 1]]->m_face = face.get();
        m_verts[m_indices[i + 2]]->m_face = face.get();

        perVertexHEs[m_verts[m_indices[i + 0]].get()].push_back(he0.get());
        perVertexHEs[m_verts[m_indices[i + 1]].get()].push_back(he1.get());
        perVertexHEs[m_verts[m_indices[i + 2]].get()].push_back(he2.get());

        m_hes.push_back(he0);
        m_hes.push_back(he1);
        m_hes.push_back(he2);
        m_faces.push_back(face);
    }

    // Set opposite half edges
    for (const auto &he0 : m_hes) {
        for (Halfedge *he1 : perVertexHEs[he0->m_end]) {
            if (he0->m_start == he1->m_end && he0->m_end == he1->m_start) {
                he0->m_opposite = he1;
                break;
            }
        }
    }
}

void Mesh::save(const std::string &filename) {
    std::ofstream writer(filename.c_str(), std::ios::out);
    if (writer.fail()) {
        FatalError("Failed to open file: %s", filename.c_str());
    }

    for (const auto &p : m_verts) {
        writer << "v " << p->pt().x << " " << p->pt().y << " " << p->pt().z << std::endl;
    }

    for (int i = 0; i < m_indices.size(); i += 3) {
        writer << "f " << (m_indices[i + 0] + 1) << " " << (m_indices[i + 1] + 1) << " " << (m_indices[i + 2] + 1) << std::endl;
    }

    writer.close();
}

Mesh::VertexIterator Mesh::v_begin() {
    return Mesh::VertexIterator(m_verts);
}

Mesh::VertexIterator Mesh::v_end() {
    return Mesh::VertexIterator(m_verts, m_verts.size());
}

// ----------
// VertexIterator
// ----------

Mesh::VertexIterator::VertexIterator(std::vector<std::shared_ptr<Vertex>> &vertices, int index)
    : m_verts{ vertices }
    , m_index{ index } {
}

bool Mesh::VertexIterator::operator!=(const Mesh::VertexIterator &it) const {
    return m_index != it.m_index;
}

Vertex &Mesh::VertexIterator::operator*() {
    return *m_verts[m_index];
}

Vertex *Mesh::VertexIterator::operator->() const {
    return m_verts[m_index].get();
}

Mesh::VertexIterator &Mesh::VertexIterator::operator++() {
    ++m_index;
    return *this;
}

Mesh::VertexIterator Mesh::VertexIterator::operator++(int) {
    int prev = m_index;
    ++m_index;
    return Mesh::VertexIterator(m_verts, prev);
}

}  // namespace tinymesh
