#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_MESH_H
#define TINYMESH_MESH_H

#include <string>
#include <vector>
#include <memory>

#include "core/common.h"

namespace tinymesh {

class TINYMESH_API Mesh {
public:
    // Forward declaration
    class VertexIterator;
    class HalfedgeIterator;
    class FaceIterator;

public:
    Mesh();
    Mesh(const std::string &filename);
    void load(const std::string &filename);
    void save(const std::string &filename);

    bool flipHE(Halfedge *he);
    bool splitHE(Halfedge *he);
    bool collapseHE(Halfedge *he);

    Vertex *vertex(int i) { return vertices_[i].get(); }

    bool verify() const;

    VertexIterator v_begin();
    VertexIterator v_end();
    HalfedgeIterator he_begin();
    HalfedgeIterator he_end();
    FaceIterator f_begin();
    FaceIterator f_end();

    size_t num_vertices() { return vertices_.size(); }
    size_t num_halfedges() { return halfedges_.size(); }
    size_t num_faces() { return faces_.size(); }

private:
    void loadOBJ(const std::string &filename);
    void loadPLY(const std::string &filename);

    void addVertex(Vertex *v);
    void addHalfedge(Halfedge *he);
    void addFace(Face *f);
    void removeVertex(Vertex *v);
    void removeHalfedge(Halfedge *he);
    void removeFace(Face *f);

    bool verifyVertex(Vertex *v) const;

    std::vector<std::shared_ptr<Vertex>> vertices_;
    std::vector<std::shared_ptr<Halfedge>> halfedges_;
    std::vector<std::shared_ptr<Face>> faces_;
    std::vector<uint32_t> indices_;

    friend class Mesh::VertexIterator;
};

class TINYMESH_API Mesh::VertexIterator {
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
    std::vector<std::shared_ptr<Vertex>> &vertices_;
    int index_;
};

class TINYMESH_API Mesh::HalfedgeIterator {
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
    std::vector<std::shared_ptr<Halfedge>> &halfedges_;
    int index_;
};

class TINYMESH_API Mesh::FaceIterator {
public:
    explicit FaceIterator(std::vector<std::shared_ptr<Face>> &faces, int index = 0);
    bool operator!=(const FaceIterator &it) const;
    Face &operator*();
    Face *operator->() const;
    Face *ptr() const;
    FaceIterator &operator++();
    FaceIterator operator++(int);
    FaceIterator &operator--();
    FaceIterator operator--(int);

private:
    std::vector<std::shared_ptr<Face>> &faces_;
    int index_;
};

}  // namespace tinymesh

#endif  // TINYMESH_MESH_H
