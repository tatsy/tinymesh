#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_MESH_H
#define TINYMESH_MESH_H

#include <string>
#include <vector>
#include <set>
#include <memory>
#include <unordered_map>

#include "core/api.h"
#include "core/vec.h"
#include "core/types.h"
#include "core/debug.h"

namespace tinymesh {

/**
 * "Mesh" is a class for halfedge data structure for polygonal mesh.
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
    Mesh clone();

    std::vector<Vec3> getVertices() const;
    std::vector<uint32_t> getVertexIndices() const;

    bool flipHE(Halfedge *he);
    bool splitHE(Halfedge *he);
    bool collapseHE(Halfedge *he);
    void fillHoles();

    double getAvgEdgeLength() const;

    bool verify() const;

    Vertex *vertex(int index) const {
        Assertion(index >= 0 && index < (int)vertices_.size(), "Vertex index out of bounds!");
        return vertices_[index].get();
    }

    Edge *edge(int index) const {
        Assertion(index >= 0 && index < (int)edges_.size(), "Edge index out of bounds!");
        return edges_[index].get();
    }

    Halfedge *halfedge(int index) const {
        Assertion(index >= 0 && index < (int)halfedges_.size(), "Halfedge index out of bounds!");
        return halfedges_[index].get();
    }

    Face *face(int index) const {
        Assertion(index >= 0 && index < (int)faces_.size(), "Face index out of bounds!");
        return faces_[index].get();
    }

    size_t numVertices() {
        return vertices_.size();
    }
    size_t numEdges() {
        return edges_.size();
    }
    size_t numHalfedges() {
        return halfedges_.size();
    }
    size_t numFaces() {
        return faces_.size();
    }

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

    void triangulate(Face *f);
    bool verifyVertex(Vertex *v) const;

    // ----- Mesh completion methods (in completion.cpp) -----

    //! This method adds a new triangle using vertices specified with "boundary" and indices in the list.
    Face *addNewTriangle(const std::vector<Vertex *> &boundary, const std::tuple<uint32_t, uint32_t, uint32_t> &tri,
                         std::unordered_map<IndexPair, Halfedge *> &pair2he,
                         std::unordered_map<Halfedge *, IndexPair> &he2pair);
    void holeFillMinDihedral_(Face *face, double dihedralBound = Pi);
    void holeFillAdvancingFront_(Face *face);

    // Private parameters
    std::vector<std::shared_ptr<Vertex>> vertices_;
    std::vector<std::shared_ptr<Edge>> edges_;
    std::vector<std::shared_ptr<Halfedge>> halfedges_;
    std::vector<std::shared_ptr<Face>> faces_;
    std::vector<uint32_t> indices_;

    // Friend methods
    friend TINYMESH_API void holeFillMinDihedral(Mesh &mesh, Face *face, double dihedralBound);
    friend TINYMESH_API void holeFillAdvancingFront(Mesh &mesh, Face *face);
};

}  // namespace tinymesh

#endif  // TINYMESH_MESH_H
