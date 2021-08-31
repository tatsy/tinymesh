#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_MESH_H
#define TINYMESH_MESH_H

#include <string>
#include <vector>
#include <set>
#include <memory>

#include "core/api.h"
#include "core/vec.h"

namespace tinymesh {

/**
 * Halfedge data structure for polygonal mesh.
 */
class TINYMESH_API Mesh {
public:
    // Public methods
    Mesh();
    explicit Mesh(const std::string &filename);
    Mesh(const std::vector<Vec3> &vertices, const std::vector<uint32_t> &indices);

    void load(const std::string &filename);
    void save(const std::string &filename) const;
    void construct(const std::vector<Vec3> &vertices, const std::vector<uint32_t> &indices);

    std::vector<Vec3> getVertices() const;
    std::vector<uint32_t> getVertexIndices() const;

    bool flipHE(Halfedge *he);
    bool splitHE(Halfedge *he);
    bool collapseHE(Halfedge *he);
    bool triangulate(Face *f, double dihedralBound = Pi);

    bool verify() const;

    Vertex* vertex(int index) const {
        Assertion(index >= 0 && index < (int)vertices_.size(), "Vertex index out of bounds!");
        return vertices_[index].get();
    }

    Edge* edge(int index) const {
        Assertion(index >= 0 && index < (int)edges_.size(), "Edge index out of bounds!");
        return edges_[index].get();
    }

    Halfedge* halfedge(int index) const {
        Assertion(index >= 0 && index < (int)halfedges_.size(), "Halfedge index out of bounds!");
        return halfedges_[index].get();
    }

    Face* face(int index) const {
        Assertion(index >= 0 && index < (int)faces_.size(), "Face index out of bounds!");
        return faces_[index].get();
    }

    //! Compute Gaussian curvature at a vertex
    double K(Vertex *v) const;
    //! Compute mean curvature at a vertex
    double H(Vertex *v) const;

    size_t numVertices() { return vertices_.size(); }
    size_t numEdges() { return edges_.size(); }
    size_t numHalfedges() { return halfedges_.size(); }
    size_t numFaces() { return faces_.size(); }

private:
    // Private methods
    void loadOBJ(const std::string &filename);
    void loadPLY(const std::string &filename);
    void saveOBJ(const std::string &filename) const;
    void savePLY(const std::string &filename) const;

    void construct();

    void addVertex(Vertex *v);
    void addEdge(Edge *e);
    void addHalfedge(Halfedge *he);
    void addFace(Face *f);
    void removeVertex(Vertex *v);
    void removeEdge(Edge *e);
    void removeHalfedge(Halfedge *he);
    void removeFace(Face *f);

    bool verifyVertex(Vertex *v) const;

    // Private parameters
    std::vector<std::shared_ptr<Vertex>> vertices_;
    std::vector<std::shared_ptr<Edge>> edges_;
    std::vector<std::shared_ptr<Halfedge>> halfedges_;
    std::vector<std::shared_ptr<Face>> faces_;
    std::vector<uint32_t> indices_;
};

}  // namespace tinymesh

#endif  // TINYMESH_MESH_H
