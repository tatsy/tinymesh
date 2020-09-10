#define TINYMESH_API_EXPORT
#include "mesh.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <queue>
#include <tuple>
#include <functional>
#include <unordered_map>

#include <Eigen/Core>

#include "core/vec.h"
#include "core/filesystem.h"
#include "polymesh/vertex.h"
#include "polymesh/edge.h"
#include "polymesh/halfedge.h"
#include "polymesh/face.h"
#include "tiny_obj_loader.h"
#include "tinyply.h"

namespace fs = std::filesystem;

using IndexPair = std::pair<uint32_t, uint32_t>;

struct IndexPairHash : public std::function<size_t(IndexPair)> {
    std::size_t operator()(const IndexPair& k) const {
        return std::get<0>(k) ^ std::get<1>(k);
    }
};

namespace tinymesh {

// ----------
// Mesh
// ----------

Mesh::Mesh() {}

Mesh::Mesh(const std::string &filename) { load(filename); }

Mesh::Mesh(const std::vector<Vec3> &vertices, const std::vector<uint32_t> &indices) {
    construct(vertices, indices);
}

void Mesh::load(const std::string &filename) {
    // Clear existing elements
    vertices_.clear();
    edges_.clear();
    halfedges_.clear();
    faces_.clear();

    // Load new mesh
    const std::string ext = fs::path(filename.c_str()).extension().string();
    if (ext == ".obj") {
        loadOBJ(filename);
    } else if (ext == ".ply") {
        loadPLY(filename);
    } else {
        FatalError("Unsupported file extension: %s", ext.c_str());
    }

    construct();
}

void Mesh::construct(const std::vector<Vec3> &vertices, const std::vector<uint32_t> &indices) {
    // Clear existing elements
    vertices_.clear();
    edges_.clear();
    halfedges_.clear();
    faces_.clear();

    std::unordered_map<Vec3, uint32_t> uniqueVertices;
    vertices_.clear();
    for (uint32_t i : indices) {
        const Vec3 &pos = vertices[i];

        if (uniqueVertices.count(pos) == 0) {
            uniqueVertices[pos] = static_cast<uint32_t>(vertices_.size());
            vertices_.push_back(std::make_shared<Vertex>(pos));
        }
        indices_.push_back(uniqueVertices[pos]);
    }

    construct();
}

std::vector<Vec3> Mesh::getVertices() const {
    std::vector<Vec3> ret;
    for (const auto &v : vertices_) {
        ret.push_back(v->pos());
    }
    return ret;
}

std::vector<uint32_t> Mesh::getVertexIndices() const {
    std::vector<uint32_t> ret;
    for (const auto &f : faces_) {
        if (f->isBoundary()) {
            continue;
        }

        const int i0 = f->halfedge_->src()->index();
        const int i1 = f->halfedge_->next()->src()->index();
        const int i2 = f->halfedge_->prev()->src()->index();
        ret.push_back(i0);
        ret.push_back(i1);
        ret.push_back(i2);
    }
    return ret;
}

void Mesh::construct() {
    for (int i = 0; i < indices_.size(); i += 3) {
        if (indices_[i + 0] == indices_[i + 1] || indices_[i + 1] == indices_[i + 2] ||
            indices_[i + 2] == indices_[i + 0]) {
            int size = (int)indices_.size();
            std::swap(indices_[i + 0], indices_[size - 3]);
            std::swap(indices_[i + 1], indices_[size - 2]);
            std::swap(indices_[i + 2], indices_[size - 1]);
            indices_.resize(size - 3);
        }
    }

    // Put vertex indices
    for (int i = 0; i < vertices_.size(); i++) {
        vertices_[i]->index_ = i;
    }

    // Check if all triangles consist of three distinct vertices
    for (int i = 0; i < indices_.size(); i += 3) {
        if (indices_[i + 0] == indices_[i + 1] || indices_[i + 1] == indices_[i + 2] ||
            indices_[i + 2] == indices_[i + 0]) {
            FatalError("Each triangle must have three distinct vertices!");
        }
    }

    // Setup halfedge structure
    static const int degree = 3;
    std::map<IndexPair, Halfedge *> pairToHalfedge;
    std::map<uint32_t, uint32_t> vertexDegree;
    for (int i = 0; i < indices_.size(); i += degree) {
        // Check polygon duplication
        bool isDuplicated = false;
        for (int j = 0; j < degree; j++) {
            const uint32_t a = indices_[i + j];
            const uint32_t b = indices_[i + (j + 1) % degree];
            IndexPair ab(a, b);
            if (pairToHalfedge.find(ab) != pairToHalfedge.end()) {
                Warn("An edge with vertices #%d and #%d is detected! Skip the face that includes this edge.", a, b);
                isDuplicated = true;
                break;
            }
        }

        if (isDuplicated) {
            continue;
        }

        // Traverse face vertices
        auto face = new Face();
        std::vector<Halfedge *> faceHalfedges;
        for (int j = 0; j < degree; j++) {
            // Count up vertex degree
            const uint32_t a = indices_[i + j];
            if (vertexDegree.find(a) == vertexDegree.end()) {
                vertexDegree[a] = 1;
            } else {
                vertexDegree[a] += 1;
            }

            // Check duplicated halfedges
            const uint32_t b = indices_[i + (j + 1) % degree];
            IndexPair ab(a, b);
            if (pairToHalfedge.find(ab) != pairToHalfedge.end()) {
                FatalError("An edge with vertices #%d and #%d is duplicated\n", a, b);
            }

            auto he = new Halfedge();
            addHalfedge(he);
            pairToHalfedge[ab] = he;

            vertices_[a]->halfedge_ = he;
            he->src_ = vertices_[a].get();
            he->face_ = face;
            face->halfedge_ = he;

            faceHalfedges.push_back(he);

            // Set opposite halfedge if it exists.
            // Besides, add an edge if opposite halfedge is successfully paired.
            IndexPair ba(b, a);
            auto iba = pairToHalfedge.find(ba);
            if (iba != pairToHalfedge.end()) {
                Halfedge *rev = iba->second;
                he->rev_ = rev;
                rev->rev_ = he;

                auto edge = new Edge();
                addEdge(edge);
                he->edge_ = edge;
                rev->edge_ = edge;
                edge->halfedge_ = he;
            } else {
                he->rev_ = nullptr;
            }
        }
        addFace(face);

        // Set next halfedges
        for (int j = 0; j < degree; j++) {
            const int k = (j + 1) % degree;
            faceHalfedges[j]->next_ = faceHalfedges[k];
        }
    }

    // For convenience, advance halfedge of each boundary vertex to one that is on the boundary
    for (auto v : vertices_) {
        Halfedge *he = v->halfedge_;
        do {
            if (he->rev_ == nullptr) {
                v->halfedge_ = he;
                break;
            }

            he = he->rev_->next_;
        } while (he != v->halfedge_);
    }

    // Construct new faces (possibly non-triangle) for each boundary components
    const int numHalfedges = halfedges_.size();
    for (int i = 0; i < numHalfedges; i++) {
        auto &he = halfedges_[i];
        Assertion(he.get() != nullptr, "Null halfedge detected!");
        if (he->rev_ == nullptr) {
            auto face = new Face();
            addFace(face);
            face->isBoundary_ = true;

            std::vector<Halfedge *> boundaryHalfedges;
            Halfedge *it = he.get();
            do {
                auto rev = new Halfedge();
                addHalfedge(rev);
                boundaryHalfedges.push_back(rev);

                auto edge = new Edge();
                addEdge(edge);

                it->rev_ = rev;
                rev->rev_ = it;
                edge->halfedge_ = it;
                it->edge_ = edge;
                rev->edge_ = edge;

                rev->face_ = face;
                rev->src_ = it->next_->src_;

                // Advance it to the next halfedge along the current boundary loop
                it = it->next_;
                while (it != he.get() && it->rev_ != nullptr) {
                    it = it->rev_;
                    it = it->next_;
                }
            } while (it != he.get());

            face->halfedge_ = boundaryHalfedges.front();

            const int degree = boundaryHalfedges.size();
            for (int j = 0; j < degree; j++) {
                const int k = (j - 1 + degree) % degree;
                boundaryHalfedges[j]->next_ = boundaryHalfedges[k];
            }
        }
    }

    // To make the later traversal easier, update the vertex halfedges to be not on the boundary.
    for (auto v : vertices_) {
        v->halfedge_ = v->halfedge_->rev_->next_;
    }

    // Check if the mesh is non-manifold
    for (auto v : vertices_) {
        if (v->halfedge_ == nullptr) {
            FatalError("Some vertices are not referenced by any polygon!");
        }

        uint32_t count = 0;
        Halfedge *he = v->halfedge_;
        do {
            if (!he->face()->isBoundary()) {
                count += 1;
            }

            he = he->rev()->next();
        } while (he != v->halfedge_);

        if (count != vertexDegree[v->index()]) {
            FatalError("At least one of the vertices is non-manifold: %d vs %d\n", count, vertexDegree[v->index()]);
        }
    }
}

void Mesh::loadOBJ(const std::string &filename) {
    // Open
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> mats;
    std::string warnMsg;
    std::string errMsg;
    bool success = tinyobj::LoadObj(&attrib, &shapes, &mats, &warnMsg, &errMsg, filename.c_str());
    if (!errMsg.empty()) {
        Warn("%s", errMsg.c_str());
    }

    if (!success) {
        FatalError("Failed to load *.obj file: %s", filename.c_str());
    }

    // Traverse triangles
    vertices_.clear();
    indices_.clear();
    std::unordered_map<Vec3, uint32_t> uniqueVertices;
    for (const auto &shape : shapes) {
        for (const auto &index : shape.mesh.indices) {
            Vec3 v;
            if (index.vertex_index >= 0) {
                v = Vec3(attrib.vertices[index.vertex_index * 3 + 0], attrib.vertices[index.vertex_index * 3 + 1],
                         attrib.vertices[index.vertex_index * 3 + 2]);
            }

            if (uniqueVertices.count(v) == 0) {
                uniqueVertices[v] = static_cast<uint32_t>(vertices_.size());
                vertices_.push_back(std::make_shared<Vertex>(v));
            }
            indices_.push_back(uniqueVertices[v]);
        }
    }

    Assertion(indices_.size() % 3 == 0, "Non-triangle face might be contained!");
}

void Mesh::loadPLY(const std::string &filename) {
    using tinyply::PlyFile;
    using tinyply::PlyData;

    try {
        // Open
        std::ifstream reader(filename.c_str(), std::ios::binary);
        if (reader.fail()) {
            FatalError("Failed to open file: %s", filename.c_str());
        }

        // Read header
        PlyFile file;
        file.parse_header(reader);

        // Request vertex data
        std::shared_ptr<PlyData> vert_data, norm_data, uv_data, face_data;
        try {
            vert_data = file.request_properties_from_element("vertex", { "x", "y", "z" });
        } catch (std::exception &e) {
            std::cerr << "tinyply exception: " << e.what() << std::endl;
        }

        try {
            norm_data = file.request_properties_from_element("vertex", { "nx", "ny", "nz" });
        } catch (std::exception &e) {
            std::cerr << "tinyply exception: " << e.what() << std::endl;
        }

        try {
            uv_data = file.request_properties_from_element("vertex", { "u", "v" });
        } catch (std::exception &e) {
            std::cerr << "tinyply exception: " << e.what() << std::endl;
        }

        try {
            face_data = file.request_properties_from_element("face", { "vertex_indices" }, 3);
        } catch (std::exception &e) {
            std::cerr << "tinyply exception: " << e.what() << std::endl;
        }

        // Read vertex data
        file.read(reader);

        // Copy vertex data
        const size_t numVerts = vert_data->count;
        std::vector<float> raw_vertices(numVerts * 3);
        std::memcpy(raw_vertices.data(), vert_data->buffer.get(), sizeof(float) * numVerts * 3);
        
        const size_t numFaces = face_data->count;
        std::vector<uint32_t> raw_indices(numFaces * 3);
        std::memcpy(raw_indices.data(), face_data->buffer.get(), sizeof(uint32_t) * numFaces * 3);

        std::unordered_map<Vec3, uint32_t> uniqueVertices;
        vertices_.clear();
        for (uint32_t i : raw_indices) {
            const Vec3 pos = Vec3(raw_vertices[i * 3 + 0], raw_vertices[i * 3 + 1], raw_vertices[i * 3 + 2]);
            if (uniqueVertices.count(pos) == 0) {
                uniqueVertices[pos] = static_cast<uint32_t>(vertices_.size());
                vertices_.push_back(std::make_shared<Vertex>(pos));
            }
            indices_.push_back(uniqueVertices[pos]);
        }
    } catch (const std::exception &e) {
        std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
    }
}

void Mesh::save(const std::string &filename) const {
    // Load new mesh
    const std::string ext = fs::path(filename.c_str()).extension().string();
    if (ext == ".obj") {
        saveOBJ(filename);
    } else if (ext == ".ply") {
        savePLY(filename);
    } else {
        FatalError("Unsupported file extension: %s", ext.c_str());
    }
}

void Mesh::saveOBJ(const std::string &filename) const {
    std::ofstream writer(filename.c_str(), std::ios::out);
    if (writer.fail()) {
        FatalError("Failed to open file: %s", filename.c_str());
    }

    for (const auto &v : vertices_) {
        writer << "v " << v->pos().x() << " " << v->pos().y() << " " << v->pos().z() << std::endl;
    }

    for (const auto &f : faces_) {
        if (f->isBoundary()) {
            continue;
        }

        const int i0 = f->halfedge_->src()->index();
        const int i1 = f->halfedge_->next()->src()->index();
        const int i2 = f->halfedge_->prev()->src()->index();
        if (i0 < vertices_.size() && i1 <= vertices_.size() && i2 <= vertices_.size()) {
            writer << "f " << (i0 + 1) << " " << (i1 + 1) << " " << (i2 + 1) << std::endl;
        }
    }

    writer.close();
}

void Mesh::savePLY(const std::string &filename) const {
    using namespace tinyply;

    std::filebuf buffer;
    buffer.open(filename.c_str(), std::ios::out | std::ios::binary);

    std::ostream outstream(&buffer);
    if (outstream.fail()) {
        FatalError("Failed to open file: %s", filename.c_str());
    }

    std::vector<float> vertexData(vertices_.size() * 3);
    for (int i = 0; i < vertices_.size(); i++) {
        const Vec3 v = vertices_[i]->pos();
        vertexData[i * 3 + 0] = v.x();
        vertexData[i * 3 + 1] = v.y();
        vertexData[i * 3 + 2] = v.z();
    }

    std::vector<uint32_t> indexData;
    indexData.reserve(faces_.size() * 3);
    for (const auto &f : faces_) {
        if (f->isBoundary()) {
            continue;
        }

        if (!f->halfedge_) {
            Warn("Face without halfedge is detected!");
            continue;
        }

        const int i0 = f->halfedge_->src()->index();
        const int i1 = f->halfedge_->next()->src()->index();
        const int i2 = f->halfedge_->prev()->src()->index();

        indexData.push_back(i0);
        indexData.push_back(i1);
        indexData.push_back(i2);
    }

    PlyFile plyfile;

    plyfile.add_properties_to_element("vertex", { "x", "y", "z" }, Type::FLOAT32, vertexData.size() / 3,
                                      reinterpret_cast<uint8_t *>(vertexData.data()), Type::INVALID, 0);

    plyfile.add_properties_to_element("face", { "vertex_indices" }, Type::UINT32, indexData.size() / 3,
                                      reinterpret_cast<uint8_t *>(indexData.data()), Type::UINT8, 3);

    plyfile.get_comments().push_back("generated by tinyply 2.2");

    plyfile.write(outstream, true);
}

bool Mesh::splitHE(Halfedge *he) {
    Halfedge *rev = he->rev();
    Vertex *v0 = he->src();
    Vertex *v1 = rev->src();
    if (v0->degree() < v1->degree()) {
        std::swap(he, rev);
        std::swap(v0, v1);
    }

    if (v0->degree() < 3 || v1->degree() < 3) {
        return false;
    }

    // Prepare halfedge / face
    auto he_new = new Halfedge();
    auto he01 = new Halfedge();
    auto he02 = new Halfedge();
    auto rev_new = new Halfedge();
    auto he11 = new Halfedge();
    auto he12 = new Halfedge();

    auto e_new = new Edge();
    auto e0 = new Edge();
    auto e1 = new Edge();

    auto he0 = he->next_;
    auto he1 = he0->next_;
    auto he2 = rev->next_;
    auto he3 = he2->next_;

    auto f0 = he->face();
    auto f1 = rev->face();
    if (f0->isBoundary() || f1->isBoundary()) {
        return false;
    }

    auto f0_new = new Face();
    auto f1_new = new Face();

    auto v2 = he1->src_;
    auto v3 = he3->src_;

    // Prepare new vertex
    const Vec3 newPos = 0.5 * (he->src()->pos() + he->dst()->pos());
    auto v_new = new Vertex(newPos);

    // Update next half-edges
    he->next_ = he0;
    he0->next_ = he01;
    he01->next_ = he;

    he_new->next_ = he02;
    he02->next_ = he1;
    he1->next_ = he_new;

    rev->next_ = he11;
    he11->next_ = he3;
    he3->next_ = rev;

    rev_new->next_ = he2;
    he2->next_ = he12;
    he12->next_ = rev_new;

    // Set edge-halfedge links
    e_new->halfedge_ = he_new;
    he_new->edge_ = e_new;
    rev_new->edge_ = e_new;

    e0->halfedge_ = he01;
    he01->edge_ = e0;
    he02->edge_ = e0;

    e1->halfedge_ = he11;
    he11->edge_ = e1;
    he12->edge_ = e1;

    // Update rev half-edges
    he->rev_ = rev;
    rev->rev_ = he;

    he_new->rev_ = rev_new;
    rev_new->rev_ = he_new;

    he01->rev_ = he02;
    he02->rev_ = he01;

    he11->rev_ = he12;
    he12->rev_ = he11;

    // Update faces
    he->face_ = f0;
    he0->face_ = f0;
    he01->face_ = f0;

    he_new->face_ = f0_new;
    he02->face_ = f0_new;
    he1->face_ = f0_new;

    rev->face_ = f1;
    he3->face_ = f1;
    he11->face_ = f1;

    rev_new->face_ = f1_new;
    he12->face_ = f1_new;
    he2->face_ = f1_new;

    // Update halfedge origins
    he->src_ = v_new;
    he0->src_ = v1;
    he01->src_ = v2;

    he_new->src_ = v0;
    he02->src_ = v_new;
    he1->src_ = v2;

    rev->src_ = v1;
    he11->src_ = v_new;
    he3->src_ = v3;

    rev_new->src_ = v_new;
    he2->src_ = v0;
    he12->src_ = v3;

    // Update vertex halfedges
    v0->halfedge_ = he2;
    v1->halfedge_ = he0;
    v_new->halfedge_ = he;

    // Update face halfedges
    f0->halfedge_ = he;
    f0_new->halfedge_ = he_new;
    f1->halfedge_ = rev;
    f1_new->halfedge_ = rev_new;

    // Add new components to the list
    addFace(f0_new);
    addFace(f1_new);
    addVertex(v_new);
    addEdge(e_new);
    addEdge(e0);
    addEdge(e1);
    addHalfedge(he_new);
    addHalfedge(rev_new);
    addHalfedge(he01);
    addHalfedge(he02);
    addHalfedge(he11);
    addHalfedge(he12);

    return true;
}

bool Mesh::collapseHE(Halfedge *he) {
    // Collapse halfedge v0->v2.
    //       v1
    //     /  \
    //    /    \
    //  v0 ---- v2
    //    \    /
    //     \  /
    //      v3

    Edge *e = he->edge_;
    Halfedge *rev = he->rev_;
    Assertion(he->index_ < halfedges_.size(), "Daemon halfedge detected!");
    Assertion(rev->index_ < halfedges_.size(), "Daemon halfedge detected!");

    // Merge vertex to the one with smaller degree
    const int d0 = he->src()->degree();
    const int d1 = rev->src()->degree();
    if (d0 > d1) {
        std::swap(he, rev);
    }

    // Skip halfedge at the boundary
    if (he->face()->isBoundary() || he->rev()->face()->isBoundary()) {
        return false;
    }

    // Check degrees of end points
    Vertex *v0 = he->src();
    Vertex *v1 = he->prev()->src();
    Vertex *v2 = rev->src();
    Vertex *v3 = rev->prev()->src();

    Assertion(v0->index_ < vertices_.size(), "Daemon vertex detected!");
    Assertion(v1->index_ < vertices_.size(), "Daemon vertex detected!");
    Assertion(v2->index_ < vertices_.size(), "Daemon vertex detected!");
    Assertion(v3->index_ < vertices_.size(), "Daemon vertex detected!");

    if (v0->degree() <= 3 || v1->degree() <= 4 || v2->degree() <= 3 || v3->degree() <= 4) {
        return false;
    }

    // Collect neighboring vertices of the one to be removed
    std::vector<Vertex *> neighbors;
    for (auto it = v2->ohe_begin(); it != v2->ohe_end(); ++it) {
        Assertion(it->dst() != nullptr, "Null vertex is found!");
        neighbors.push_back(it->dst());
    }

    Assertion(!neighbors.empty(), "No neighboring vertex!");

    // Check unsafe collapse (# of shared vertices)
    int numUnion = 0;
    std::set<Vertex *> neighborSet(neighbors.begin(), neighbors.end());
    for (auto it = v0->ohe_begin(); it != v0->ohe_end(); ++it) {
        Assertion(it->dst() != nullptr, "Null vertex is found!");
        if (neighborSet.find(it->dst()) != neighborSet.end()) {
            numUnion += 1;
        }
    }

    if (numUnion != 2) {
        return false;
    }

    // Check unsafe collapse (face flip)
    Vec3 norm(0.0);
    for (int i = 0; i < neighbors.size(); i++) {
        const int j = (i + 1) % neighbors.size();
        const Vec3 e1 = neighbors[i]->pos() - v2->pos();
        const Vec3 e2 = neighbors[j]->pos() - v2->pos();
        norm += cross(e1, e2);
    }

    // Degenerated face if norm is zero-length
    if (length(norm) == 0.0) {
        return false;
    }
    norm = normalize(norm);

    double sign = 0.0;
    for (int i = 0; i < neighbors.size(); i++) {
        const int j = (i + 1) % neighbors.size();
        const Vec3 e1 = neighbors[i]->pos() - v0->pos();
        const Vec3 e2 = neighbors[j]->pos() - v0->pos();
        const double s = dot(norm, cross(e1, e2));
        if (sign * s < 0.0) {
            return false;
        }
        sign = s;
    }

    // Update origins of all outward halfedges
    Vertex *v_remove = v2;
    Vertex *v_remain = v0;
    Assertion(v_remain->index() < vertices_.size(), "Something is wrong1!");
    Assertion(v_remove->index() < vertices_.size(), "Something is wrong2!");
    Assertion(v_remove != v_remain, "Something is wrong3!");
    for (auto it = v_remove->ohe_begin(); it != v_remove->ohe_end(); ++it) {
        Assertion(it->src() == v_remove, "Invalid halfedge origin detected!");
        it->src_ = v_remain;
    }

    for (auto it = v_remain->ohe_begin(); it != v_remain->ohe_end(); ++it) {
        Assertion(it->src() == v_remain, "Invalid halfedge origin detected!");
    }

    // Update edges and opposite halfedges
    Edge *e0 = he->next()->next()->edge();
    Edge *e1 = he->next()->edge();
    Halfedge *he0 = he->next()->next()->rev();
    Halfedge *he1 = he->next()->rev();
    e0->halfedge_ = he0;
    he0->rev_ = he1;
    he0->edge_ = e0;
    he1->rev_ = he0;
    he1->edge_ = e0;

    Edge *e2 = rev->next()->next()->edge();
    Edge *e3 = rev->next()->edge();
    Halfedge *he2 = rev->next()->next()->rev();
    Halfedge *he3 = rev->next()->rev();
    e2->halfedge_ = he2;
    he2->rev_ = he3;
    he2->edge_ = e2;
    he3->rev_ = he2;
    he3->edge_ = e2;

    // Update halfedge of the origin vertex (v0)
    v_remain->halfedge_ = he0;

    // Update halfedge for wing-vertices
    v0->halfedge_ = he0;
    v1->halfedge_ = he1;
    v2->halfedge_ = he2;
    v3->halfedge_ = he3;

    // Remove edges
    removeEdge(e);
    removeEdge(e1);
    removeEdge(e3);

    // Remove halfedges
    Face *f0 = he->face_;
    Face *f1 = rev->face_;
    Assertion(!f0->isBoundary(), "Invalid face removal!");
    Assertion(!f1->isBoundary(), "Invalid face removal!");

    removeHalfedge(he->next_->next_);
    removeHalfedge(he->next_);
    removeHalfedge(he);
    removeHalfedge(rev->next_->next_);
    removeHalfedge(rev->next_);
    removeHalfedge(rev);

    // Remove faces
    removeFace(f0);
    removeFace(f1);

    // Remove vertices
    removeVertex(v_remove);

    return true;
}

bool Mesh::flipHE(Halfedge *he) {
    Halfedge *rev = he->rev_;
    if (!rev) {
        FatalError("Flip is called boundary halfedge!");
    }

    Vertex *u0 = he->src();
    Vertex *u1 = he->prev()->src();
    Vertex *u2 = rev->src();
    Vertex *u3 = rev->prev()->src();
    if (u0->degree() <= 4 || u2->degree() <= 4) {
        return false;
    }

    // Cannot flip halfedge if ends vertices of flipped
    // halfedge are already the ends of other edges.
    for (auto it = u1->v_begin(); it != u1->v_end(); ++it) {
        if (it.ptr() == u3) {
            return false;
        }
    }

    // Get surrounding vertices, halfedges and faces
    Halfedge *he0 = he->next_;
    Halfedge *he1 = he->next_->next_;
    Halfedge *he2 = rev->next_;
    Halfedge *he3 = rev->next_->next_;

    Vertex *v0 = he0->src_;
    Vertex *v1 = he1->src_;
    Vertex *v2 = he2->src_;
    Vertex *v3 = he3->src_;

    Face *f0 = he->face_;
    Face *f1 = rev->face_;

    // Update halfedges of start/end vertices
    v0->halfedge_ = he0;
    v2->halfedge_ = he2;

    // Uate_ halfedge's start and end
    he->src_ = v3;
    rev->src_ = v1;

    // Update face circulation
    he->next_ = he1;
    he1->next_ = he2;
    he2->next_ = he;
    rev->next_ = he3;
    he3->next_ = he0;
    he0->next_ = rev;

    // Update faces
    f0->halfedge_ = he;
    he->face_ = f0;
    he1->face_ = f0;
    he2->face_ = f0;

    f1->halfedge_ = rev;
    rev->face_ = f1;
    he0->face_ = f1;
    he3->face_ = f1;

    return true;
}

bool Mesh::triangulate(Face *face) {
    /*
     * This is an triangulation algorithm described in the following paper.
     * Barequet and Sharir, "Filling Gaps in the Boundary of a Polyhedron", 1995.
     */

    using E = IndexPair;
    using T = std::tuple<uint32_t, uint32_t, uint32_t>;

    std::vector<Vertex *> vs;
    std::unordered_map<E, Halfedge *, IndexPairHash> pair2he;
    std::unordered_map<Halfedge *, E> he2pair;

    Halfedge *it = face->halfedge_;
    do {
        vs.push_back(it->src());
        it = it->next_;
    } while (it != face->halfedge_);
    const int n_verts = (int)vs.size();

    int count = 0;
    it = face->halfedge_;
    do {
        const int next = (count + 1) % n_verts;
        pair2he.insert(std::make_pair(E(count, next), it));
        he2pair.insert(std::make_pair(it, E(count, next)));
        it = it->next_;
        count += 1;
    } while (it != face->halfedge_);

    Eigen::MatrixXd W(n_verts, n_verts);
    Eigen::MatrixXi O(n_verts, n_verts);
    W.setConstant(1.0e12);
    O.setConstant(-1);
    for (int i = 0; i < n_verts - 1; i++) {
        W(i, i + 1) = 0.0;
        O(i, i + 1) = -1;
        if (i < n_verts - 2) {
            const Vec3 v0 = vs[i]->pos();
            const Vec3 v1 = vs[i + 1]->pos();
            const Vec3 v2 = vs[i + 2]->pos();
            W(i, i + 2) = 0.5 * length(cross(v1 - v0, v2 - v0));
            O(i, i + 2) = i + 1;
        }
    }

    for (int j = 3; j < n_verts; j++) {
        for (int i = 0; i < n_verts - j; i++) {
            const int k = i + j;
            for (int m = i + 1; m < k; m++) {
                const Vec3 v0 = vs[i]->pos();
                const Vec3 v1 = vs[m]->pos();
                const Vec3 v2 = vs[k]->pos();
                const double F = 0.5 * length(cross(v1 - v0, v2 - v0));
                const double Wnew = W(i, m) + W(m, k) + F;
                if (Wnew < W(i, k)) {
                    W(i, k) = Wnew;
                    O(i, k) = m;
                }
            }
        }
    }

    std::queue<E> que;
    std::vector<T> tris;

    que.push(E(0, n_verts - 1));
    while (!que.empty()) {
        const E ik = que.front();
        que.pop();

        const int i = std::get<0>(ik);
        const int k = std::get<1>(ik);
        const int m = O(i, k);

        if (i + 2 == k) {
            tris.emplace_back(i, m, k);
        } else {
            if (i + 1 != m) {
                que.push(E(i, m));
            }

            tris.emplace_back(i, m, k);

            if (m != k - 1) {
                que.push(E(m, k));
            }
        }
    }

    std::vector<Face *> new_faces;
    for (const auto &t : tris) {
        const int i0 = std::get<0>(t);
        const int i1 = std::get<1>(t);
        const int i2 = std::get<2>(t);

        Halfedge *new_he0 = nullptr;
        if (pair2he.count(E(i0, i1)) != 0) {
            new_he0 = pair2he[E(i0, i1)];
        } else {
            new_he0 = new Halfedge();
            pair2he.insert(std::make_pair(E(i0, i1), new_he0));
            he2pair.insert(std::make_pair(new_he0, E(i0, i1)));
            addHalfedge(new_he0);
        }

        Halfedge *new_he1 = nullptr;
        if (pair2he.count(E(i1, i2)) != 0) {
            new_he1 = pair2he[E(i1, i2)];
        } else {
            new_he1 = new Halfedge();
            pair2he.insert(std::make_pair(E(i1, i2), new_he1));
            he2pair.insert(std::make_pair(new_he1, E(i1, i2)));
            addHalfedge(new_he1);
        }

        Halfedge *new_he2 = nullptr;
        if (pair2he.count(E(i2, i0)) != 0) {
            new_he2 = pair2he[E(i2, i0)];
        } else {
            new_he2 = new Halfedge();
            pair2he.insert(std::make_pair(E(i2, i0), new_he2));
            he2pair.insert(std::make_pair(new_he2, E(i2, i0)));
            addHalfedge(new_he2);
        }

        new_he0->src_ = vs[i0];
        new_he1->src_ = vs[i1];
        new_he2->src_ = vs[i2];

        new_he0->next_ = new_he1;
        new_he1->next_ = new_he2;
        new_he2->next_ = new_he0;

        auto new_face = new Face();
        addFace(new_face);
        new_faces.push_back(new_face);

        new_he0->face_ = new_face;
        new_he1->face_ = new_face;
        new_he2->face_ = new_face;
        new_face->halfedge_ = new_he0;
    }

    for (auto f : new_faces) {
        Halfedge *he = f->halfedge_;
        do {
            if (he->rev_ == NULL) {
                auto e = he2pair[he];
                int i = std::get<0>(e);
                int j = std::get<1>(e);

                auto r = E(j, i);
                if (pair2he.count(r) != 0) {
                    Halfedge *rev = pair2he[r];
                    he->rev_ = rev;
                    rev->rev_ = he;

                    Edge *new_edge = new Edge();
                    he->edge_ = new_edge;
                    rev->edge_ = new_edge;
                    new_edge->halfedge_ = he;
                    addEdge(new_edge);
                } else {
                    printf("%d %d %p\n", i, j, he);
                }
            }

            he = he->next_;
        } while (he != f->halfedge_);
    }

    removeFace(face);

    return true;
}

double Mesh::K(Vertex *v) const {
    std::vector<Vertex *> neighbors;
    for (auto it = v->v_begin(); it != v->v_end(); ++it) {
        neighbors.push_back(it.ptr());
    }

    const int N = static_cast<int>(neighbors.size());
    double sumAngles = 0.0;
    double sumAreas = 0.0;
    for (int i = 0; i < N; i++) {
        const int j = (i + 1) % N;
        const Vec3 e1 = neighbors[i]->pos() - v->pos();
        const Vec3 e2 = neighbors[j]->pos() - v->pos();
        sumAngles += std::atan2(length(cross(e1, e2)), dot(e1, e2));
        sumAreas += length(cross(e1, e2)) / 6.0;
    }
    return (2.0 * Pi - sumAngles) / sumAreas;
}

double Mesh::H(Vertex *v) const {
    std::vector<Vertex *> neighbors;
    for (auto it = v->v_begin(); it != v->v_end(); ++it) {
        neighbors.push_back(it.ptr());
    }

    const int N = static_cast<int>(neighbors.size());
    Vec3 laplace = Vec3(0.0);
    double sumAreas = 0.0;
    for (int i = 0; i < N; i++) {
        const int prev = (i - 1 + N) % N;
        const int post = (i + 1) % N;

        const Vec3 ea1 = v->pos() - neighbors[prev]->pos();
        const Vec3 ea2 = neighbors[i]->pos() - neighbors[prev]->pos();
        const double cota = length(ea1) * length(ea2) / length(cross(ea1, ea2));
        const Vec3 eb1 = v->pos() - neighbors[post]->pos();
        const Vec3 eb2 = v->pos() - neighbors[post]->pos();
        const double cotb = length(eb1) * length(eb2) / length(cross(eb1, eb2));
        laplace += (cota + cotb) * (neighbors[i]->pos() - v->pos());

        const Vec3 e1 = neighbors[i]->pos() - v->pos();
        const Vec3 e2 = neighbors[post]->pos() - v->pos();
        sumAreas += length(cross(e1, e2)) / 6.0;
    }
    return length(laplace) / (2.0 * sumAreas);
}

bool Mesh::verify() const {
    bool success = true;
    for (int i = 0; i < vertices_.size(); i++) {
        Vertex *v = vertices_[i].get();
        if (v->index() != i) {
            Warn("Vertex index does not match array index: v[%d].index = %d", i, v->index());
            success = false;
        }
        success &= verifyVertex(v);
    }

    for (int i = 0; i < edges_.size(); i++) {
        Edge *e = edges_[i].get();
        if (e->index() != i) {
            Warn("Edge index does not match array index: e[%d].index = %d", i, e->index());
            success = false;
        }

        if (e->halfedge_ == nullptr) {
            Warn("Edge[%d] is not with non-null halfedge!", e->index());
            success = false;
        }
    }

    for (int i = 0; i < halfedges_.size(); i++) {
        auto he = halfedges_[i];
        if (he->index() != i) {
            Warn("Halfedge index does not match array index: he[%d].index = %d", i, he->index());
            success = false;
        }

        if (he->edge_ == nullptr) {
            Warn("Halfedge[%d] is not with non-null edge!", he->index());
            success = false;
        }
    }

    for (int i = 0; i < faces_.size(); i++) {
        auto f = faces_[i];
        if (f->index() != i) {
            Warn("Face index does not match array index: f[%d].index = %d", i, f->index());
            success = false;
        }

        if (f->halfedge_ == nullptr) {
            Warn("Face[%d] is not with non-null halfedge!", f->index());
            success = false;
        }
    }

    return success;
}

bool Mesh::verifyVertex(Vertex *v) const {
    bool success = true;
    if (v->index() >= vertices_.size()) {
        fprintf(stderr, "Vertex index out of range: %d > %d\n", v->index(), (int)vertices_.size());
        success = false;
    }

    Vec3 p = v->pos();
    if (std::isinf(p[0]) || std::isnan(p[0]) || std::isinf(p[1]) || std::isnan(p[1]) || std::isinf(p[2]) ||
        std::isnan(p[2])) {
        fprintf(stderr, "Inf of NaN found at v[%d]: (%f, %f, %f)\n", v->index(), p[0], p[1], p[2]);
        success = false;
    }

    for (auto it = v->ohe_begin(); it != v->ohe_end(); ++it) {
        if (it->src() != v) {
            fprintf(stderr, "Origin of an outward halfedge is invalid at v[%d]!\n", v->index());
            success = false;
        }
    }

    return success;
}

void Mesh::addVertex(Vertex *v) {
    v->index_ = static_cast<int>(vertices_.size());
    vertices_.emplace_back(v);
}

void Mesh::addEdge(Edge *e) {
    e->index_ = static_cast<int>(edges_.size());
    edges_.emplace_back(e);
}

void Mesh::addHalfedge(Halfedge *he) {
    he->index_ = static_cast<int>(halfedges_.size());
    halfedges_.emplace_back(he);
}

void Mesh::addFace(Face *f) {
    f->index_ = static_cast<int>(faces_.size());
    faces_.emplace_back(f);
}

void Mesh::removeVertex(Vertex *v) {
    Assertion((*v) == (*vertices_[v->index_]), "Invalid vertex indexing!");
    Assertion(v->index_ < vertices_.size(), "Vertex index out of range!");

    if (v->index_ < vertices_.size() - 1) {
        std::swap(vertices_[v->index_], vertices_[vertices_.size() - 1]);
        std::swap(vertices_[v->index_]->index_, vertices_[vertices_.size() - 1]->index_);
    }
    vertices_.resize(vertices_.size() - 1);
}

void Mesh::removeEdge(Edge *e) {
    Assertion((*e) == (*edges_[e->index_]), "Invalid edge indexing!");
    Assertion(e->index_ < edges_.size(), "Edge index out of range!");

    if (e->index_ < edges_.size() - 1) {
        std::swap(edges_[e->index_], edges_[edges_.size() - 1]);
        std::swap(edges_[e->index_]->index_, edges_[edges_.size() - 1]->index_);
    }
    edges_.resize(edges_.size() - 1);
}

void Mesh::removeHalfedge(Halfedge *he) {
    Assertion((*he) == (*halfedges_[he->index_]), "Invalid halfedge indexing!");

    if (he->index_ < halfedges_.size() - 1) {
        std::swap(halfedges_[he->index_], halfedges_[halfedges_.size() - 1]);
        std::swap(halfedges_[he->index_]->index_, halfedges_[halfedges_.size() - 1]->index_);
    }
    halfedges_.resize(halfedges_.size() - 1);
}

void Mesh::removeFace(Face *f) {
    Assertion((*f) == (*faces_[f->index_]), "Invalid face indexing!");

    if (f->index_ < faces_.size() - 1) {
        std::swap(faces_[f->index_], faces_[faces_.size() - 1]);
        std::swap(faces_[f->index_]->index_, faces_[faces_.size() - 1]->index_);
    }
    faces_.resize(faces_.size() - 1);
}

}  // namespace tinymesh
