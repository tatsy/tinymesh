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
    class HalfedgeIterator;

public:
    Mesh();
    Mesh(const std::string &filename);
    void load(const std::string &filename);
    void save(const std::string &filename);
    void flipHE(Halfedge *he);
    bool splitHE(Halfedge *he);
    bool collapseHE(Halfedge *he);

    void verify() const;

    VertexIterator v_begin();
    VertexIterator v_end();
    HalfedgeIterator he_begin();
    HalfedgeIterator he_end();

    size_t num_vertices() { return m_verts.size(); }
    size_t num_halfedges() { return m_hes.size(); }
    size_t num_faces() { return m_faces.size(); }

private:
    void addVertex(Vertex *v);
    void addHalfedge(Halfedge *he);
    void addFace(Face *f);
    void removeVertex(Vertex *v);
    void removeHalfedge(Halfedge *he);
    void removeFace(Face *f);

    void verifyVertex(Vertex *v) const;

    std::vector<std::shared_ptr<Vertex>> m_verts;
    std::vector<std::shared_ptr<Halfedge>> m_hes;
    std::vector<std::shared_ptr<Face>> m_faces;
    std::vector<uint32_t> m_indices;

    friend class Mesh::VertexIterator;
};

class TINYMESH_EXPORTS Mesh::VertexIterator {
public:
    explicit VertexIterator(std::vector<std::shared_ptr<Vertex>> &verts, int index = 0);
    bool operator!=(const VertexIterator &it) const;
    Vertex &operator*();
    Vertex *operator->() const;
    Vertex *ptr() const;
    VertexIterator &operator++();
    VertexIterator operator++(int);
    VertexIterator &operator--();
    VertexIterator operator--(int);

private:
    std::vector<std::shared_ptr<Vertex>> &m_verts;
    int m_index;
};

class TINYMESH_EXPORTS Mesh::HalfedgeIterator {
public:
    explicit HalfedgeIterator(std::vector<std::shared_ptr<Halfedge>> &hes, int index = 0);
    bool operator!=(const HalfedgeIterator &it) const;
    Halfedge &operator*();
    Halfedge *operator->() const;
    Halfedge *ptr() const;
    HalfedgeIterator &operator++();
    HalfedgeIterator operator++(int);
    HalfedgeIterator &operator--();
    HalfedgeIterator operator--(int);

private:
    std::vector<std::shared_ptr<Halfedge>> &m_hes;
    int m_index;
};

}  // namespace tinymesh

#endif  // TINYMESH_MESH_H
