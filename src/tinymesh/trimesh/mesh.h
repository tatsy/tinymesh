#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_MESH_H
#define TINYMESH_MESH_H

#include <string>
#include <vector>
#include <set>
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
    void save(const std::string &filename) const;

    bool flipHE(Halfedge *he);
    bool splitHE(Halfedge *he);
    bool collapseHE(Halfedge *he);
    bool triangulate(Face *f);

    Vertex *vertex(int i) { return vertices_[i].get(); }

    bool verify() const;

    Vertex *vertex(int index) const {
        Assertion(index >= 0 && index < vertices_.size(), "Vertex index out of bounds!");
        return vertices_[index].get();
    }

    Halfedge* halfedge(int index) const {
        Assertion(index >= 0 && index < halfedges_.size(), "Halfedge index out of bounds!");
        return halfedges_[index].get();
    }

    Face* face(int index) const {
        Assertion(index >= 0 && index < faces_.size(), "Face index out of bounds!");
        return faces_[index].get();
    }

    size_t num_vertices() { return vertices_.size(); }
    size_t num_halfedges() { return halfedges_.size(); }
    size_t num_faces() { return faces_.size(); }

private:
    void loadOBJ(const std::string &filename);
    void loadPLY(const std::string &filename);
    void saveOBJ(const std::string &filename) const;
    void savePLY(const std::string &filename) const;

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
};

}  // namespace tinymesh

#endif  // TINYMESH_MESH_H
