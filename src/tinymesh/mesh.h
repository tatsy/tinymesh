#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_MESH_H
#define TINYMESH_MESH_H

#include <string>
#include <vector>
#include <memory>

#include "common.h"

namespace tinymesh {

class TINYMESH_EXPORTS Mesh {
public:
    // Forward declaration
    class VertexIterator;

public:
    Mesh();
    Mesh(const std::string &filename);
    void load(const std::string &filename);
    void save(const std::string &filename);

    VertexIterator v_begin();
    VertexIterator v_end();

    size_t num_vertices() { return m_verts.size(); }

private:
    std::vector<std::shared_ptr<Vertex>> m_verts;
    std::vector<std::shared_ptr<Halfedge>> m_hes;
    std::vector<std::shared_ptr<Face>> m_faces;
    std::vector<uint32_t> m_indices;

    friend class Mesh::VertexIterator;
};

class TINYMESH_EXPORTS Mesh::VertexIterator {
public:
    explicit VertexIterator(std::vector<std::shared_ptr<Vertex>> &m_vtx, int index = 0);
    bool operator!=(const VertexIterator &it) const;
    Vertex &operator*();
    Vertex *operator->() const;
    VertexIterator &operator++();
    VertexIterator operator++(int);

private:
    std::vector<std::shared_ptr<Vertex>> &m_verts;
    int m_index;
};

}  // namespace tinymesh

#endif  // TINYMESH_MESH_H
