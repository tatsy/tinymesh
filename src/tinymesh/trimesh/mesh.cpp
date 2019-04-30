#define TINYMESH_API_EXPORT
#include "mesh.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <functional>
#include <unordered_map>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

#include "core/vec.h"
#include "trimesh/vertex.h"
#include "trimesh/halfedge.h"
#include "trimesh/face.h"
#include "tiny_obj_loader.h"
#include "tinyply.h"

using IndexPair = std::pair<uint32_t, uint32_t>;

namespace tinymesh {

// ----------
// Mesh
// ----------

Mesh::Mesh() {}

Mesh::Mesh(const std::string &filename) {
    load(filename);
}

void Mesh::load(const std::string &filename) {
    // Clear existing elements
    halfedges_.clear();
    vertices_.clear();
    faces_.clear();

    // Load new mesh
    const std::string ext = fs::path(filename.c_str()).extension();
    if (ext == ".obj") {
        loadOBJ(filename);
    } else if (ext == ".ply") {
        loadPLY(filename);
    } else {
        FatalError("Unsupported file extension: %s", ext.c_str());
    }

    // Put vertex indices
    for (int i = 0; i < vertices_.size(); i++) {
        vertices_[i]->index_ = i;
    }

    // Check if all triangles consist of three distinct vertices
    for (int i = 0; i < indices_.size(); i += 3) {
        if (indices_[i] == indices_[i + 1] || indices_[i + 1] == indices_[i + 2] || indices_[i + 1] == indices_[i + 2]) {
            FatalError("Each triangle must have three distinct vertices!");
        }
    }

    // Setup halfedge structure
    static const int degree = 3;
    std::map<IndexPair, Halfedge*> pairToHalfedge;
    std::map<uint32_t, uint32_t> vertexDegree;
    for (int i = 0; i < indices_.size(); i += degree) {
        // Traverse face vertices
        auto face = std::make_shared<Face>();
        std::vector<Halfedge*> faceHalfedges;
        for (int j = 0; j < degree; j++) {
            const uint32_t a = indices_[i + j];
            if (vertexDegree.find(a) == vertexDegree.end()) {
                vertexDegree[a] = 1;
            } else {
                vertexDegree[a] += 1;
            }

            const uint32_t b = indices_[i + (j + 1) % degree];
            IndexPair ab(a, b);
            auto he = std::make_shared<Halfedge>();

            if (pairToHalfedge.find(ab) != pairToHalfedge.end()) {
                FatalError("Duplicated halfedges are found!");
            }

            halfedges_.push_back(he);
            pairToHalfedge[ab] = he.get();

            vertices_[a]->halfedge_ = he.get();
            he->src_ = vertices_[a].get();
            he->face_ = face.get();
            face->halfedge_ = he.get();

            faceHalfedges.push_back(he.get());

            // Set opposite halfedge if exists
            IndexPair ba(b, a);
            auto iba = pairToHalfedge.find(ba);
            if (iba != pairToHalfedge.end()) {
                Halfedge *rev = iba->second;

                he->rev_ = rev;
                rev->rev_ = he.get();
            } else {
                he->rev_ = nullptr;
            }
        }
        faces_.push_back(face);

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
    for (auto he : halfedges_) {
        if (he->rev_ == nullptr) {
            auto face = std::make_shared<Face>();
            faces_.push_back(face);
            face->isBoundary_ = true;

            std::vector<Halfedge*> boundaryHalfedges;
            Halfedge *it = he.get();
            do {
                auto rev = std::make_shared<Halfedge>();
                halfedges_.push_back(rev);
                boundaryHalfedges.push_back(rev.get());

                it->rev_ = rev.get();
                rev->rev_ = it;
                rev->face_ = face.get();
                rev->src_ = it->next_->src_;

                // Advance it to the next halfedge along the current boundary loop
                it = it->next_;
                while (it != he.get() && it->rev_ != nullptr) {
                    it = it->rev_->next_;
                }
            } while (it != he.get());

            face->halfedge_ = boundaryHalfedges.front();

            const int degree = boundaryHalfedges.size();
            printf("deg = %d\n", degree);
            for (int j = 0; j < degree; j++) {
                const int k = (j - 1 + degree) % degree;
                boundaryHalfedges[j]->next_ = boundaryHalfedges[k];
            }
        }
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
            FatalError("At least one of the vertices is non-manifold!");
        }
    }

    // Put halfedge indices
    for (int i = 0; i < halfedges_.size(); i++) {
        halfedges_[i]->index_ = i;
    }

    // Put face indices
    for (int i = 0; i < faces_.size(); i++) {
        faces_[i]->index_ = i;
    }
}

void Mesh::save(const std::string &filename) {
    std::ofstream writer(filename.c_str(), std::ios::out);
    if (writer.fail()) {
        FatalError("Failed to open file: %s", filename.c_str());
    }

    for (const auto &v : vertices_) {
        writer << "v " << v->pos().x << " " << v->pos().y << " " << v->pos().z << std::endl;
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
    std::unordered_map<Vec, uint32_t> uniqueVertices;
    for (const auto &shape : shapes) {
        for (const auto &index : shape.mesh.indices) {
            Vec v;
            if (index.vertex_index >= 0) {
                v = Vec(attrib.vertices[index.vertex_index * 3 + 0],
                        attrib.vertices[index.vertex_index * 3 + 1],
                        attrib.vertices[index.vertex_index * 3 + 2]);
            }

            if (uniqueVertices.count(v) == 0) {
                uniqueVertices[v] = static_cast<uint32_t>(vertices_.size());
                vertices_.push_back(std::make_shared<Vertex>(v));
            }

            indices_.push_back(uniqueVertices[v]);
        }
    }
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
        std::vector<float> raw_vertices, raw_normals, raw_uv;
        if (vert_data) {
            raw_vertices.resize(numVerts * 3);
            std::memcpy(raw_vertices.data(), vert_data->buffer.get(), sizeof(float) * numVerts * 3);
        }

        const size_t numFaces = face_data->count;
        std::vector<uint32_t> raw_indices(numFaces * 3);
        std::memcpy(raw_indices.data(), face_data->buffer.get(), sizeof(uint32_t) * numFaces * 3);

        std::unordered_map<Vec, uint32_t> uniqueVertices;
        vertices_.clear();
        for (uint32_t i : raw_indices) {
            Vec pos;

            if (vert_data) {
                pos = Vec(raw_vertices[i * 3 + 0],
                          raw_vertices[i * 3 + 1],
                          raw_vertices[i * 3 + 2]);
            }

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

bool Mesh::splitHE(Halfedge *he) {
    Halfedge *rev = he->rev();
    Vertex *v0 = he->src();
    Vertex *v1 = rev->src();
    if (v0->degree() < v1->degree()) {
        std::swap(he, rev);
        std::swap(v0, v1);
    }

    if (v0->degree() < 6) {
        return false;
    }

    // Prepare halfedge / face
    auto he_00 = new Halfedge();
    auto he_01 = new Halfedge();
    auto he_02 = new Halfedge();
    auto he_10 = new Halfedge();
    auto he_11 = new Halfedge();
    auto he_12 = new Halfedge();
    
    auto f0 = new Face();
    auto f1 = new Face();

    // Prepare new vertex
    const Vec newPos = 0.5 * (he->src()->pos() + he->dst()->pos());
    auto v_new = new Vertex(newPos);

    // Setup connection between new components
    he_01->next_ = he_02;
    he_02->next_ = he_00;
    he_00->next_ = he_01;
    he_10->next_ = he_11;
    he_11->next_ = he_12;
    he_12->next_ = he_10;

    he_00->face_ = f0;
    he_01->face_ = f0;
    he_02->face_ = f0;
    f0->halfedge_ = he_00;

    he_10->face_ = f1;
    he_11->face_ = f1;
    he_12->face_ = f1;
    f1->halfedge_ = he_10;

    // Update opposite halfedges
    Halfedge *whe0 = he->prev()->rev()->prev();
    Halfedge *wre0 = whe0->rev();
    Halfedge *whe1 = rev->next()->rev()->next();
    Halfedge *wre1 = whe1->rev();

    whe0->rev_ = he_01;
    wre0->rev_ = he_02;
    whe1->rev_ = he_12;
    wre1->rev_ = he_11;

    he_01->rev_ = whe0;
    he_02->rev_ = wre0;
    he_12->rev_ = whe1;
    he_11->rev_ = wre1;

    he_00->rev_ = he_10;
    he_10->rev_ = he_00;

    // Update halfedge origins
    he_00->src_ = v0;
    he_01->src_ = v_new;
    he_02->src_ = whe0->src_;
    he_10->src_ = v_new;
    he_11->src_ = v0;
    he_12->src_ = wre1->src_;

    he->src_ = v_new;
    he->prev()->rev()->src_ = v_new;
    rev->next()->src_ = v_new;
    whe1->src_ = v_new;

    // Update halfedges for vertices
    v0->halfedge_ = he_00;
    v_new->halfedge_ = he_10;

    // Add new components to the list
    addFace(f0);
    addFace(f1);
    addVertex(v_new);
    addHalfedge(he_00);
    addHalfedge(he_01);
    addHalfedge(he_02);
    addHalfedge(he_10);
    addHalfedge(he_11);
    addHalfedge(he_12);

    return true;
}

bool Mesh::collapseHE(Halfedge* he) {
    Halfedge *rev = he->rev_;
    const int d0 = he->src()->degree();
    const int d1 = rev->src()->degree();
    // Merge vertex to the one with smaller degree
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
    if (v0->degree() <= 3 || v1->degree() <= 4 || v2->degree() <= 3 || v3->degree() <= 4) {
        return false;
    }

    // Check unsafe collapse
    std::set<Vertex*> neighbors;
    for (auto it = v0->ohe_begin(); it != v0->ohe_end(); ++it) {
        neighbors.insert(it->dst());
    }

    int numUnion = 0;
    for (auto it = v2->ohe_begin(); it != v2->ohe_end(); ++it) {
        if (neighbors.find(it->dst()) != neighbors.end()) {
            numUnion += 1;            
        }
    }

    if (numUnion != 2) {
        return false;
    }

    // Update halfedge of origin
    he->src()->halfedge_ = he->prev()->rev();

    // Update origins of all outward halfedges
    Vertex *v_remove = rev->src_;
    Vertex *v_remain = he->src_;
    for (auto it = he->src_->ohe_begin(); it != he->src_->ohe_end(); ++it) {
        it->src_ = v_remain;
    }

    for (auto it = rev->src_->ohe_begin(); it != rev->src_->ohe_end(); ++it) {
        it->src_ = v_remain;
    }

    // Update opposite halfedges
    Halfedge *he0 = he->next_->rev_;
    Halfedge *he1 = he->next_->next_->rev_;
    he0->rev_ = he1;
    he1->rev_ = he0;

    Halfedge *he2 = rev->next_->rev_;
    Halfedge *he3 = rev->next_->next_->rev_;
    he2->rev_ = he3;
    he3->rev_ = he2;

    // Update halfedge for wing-vertices
    he->prev()->src()->halfedge_ = he->next_->rev_;
    rev->prev()->src()->halfedge_ = rev->next_->rev_;

    // Remove faces
    Face *f0 = he->face_;
    Face *f1 = rev->face_;
    removeFace(f0);
    removeFace(f1);

    // Remove vertices
    removeVertex(v_remove);

    // Remove halfedges
    removeHalfedge(he->next_->next_);
    removeHalfedge(he->next_);
    removeHalfedge(he);
    removeHalfedge(rev->next_->next_);
    removeHalfedge(rev->next_);
    removeHalfedge(rev);

    return true;
}

bool Mesh::flipHE(Halfedge* he) {    
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

bool Mesh::verify() const {
    bool success = true;
    for (int i = 0; i < vertices_.size(); i++) {
        Vertex *v = vertices_[i].get();
        if (v->index() != i) {
            fprintf(stderr, "Vertex index does not match array index: v[%d].index = %d\n", i, v->index());
            success = false;
        }
        success &= verifyVertex(v);
    }

    for (int i = 0; i < halfedges_.size(); i++) {
        auto he = halfedges_[i];
        if (he->index() != i) {
            fprintf(stderr, "Halfedge index does not match array index: he[%d].index = %d\n", i, he->index());
            success = false;
        }

        if (!he->isBorder() && he.get() != he->rev()->rev()) {
            fprintf(stderr, "Halfedge opposition conflict: %p != %p\n", he.get(), he->rev()->rev());
            success = false;
        }

        success &= verifyVertex(he->src());
        success &= verifyVertex(he->dst());
    }

//    for (int i = 0; i < faces_.size(); i++) {
//        auto f = faces_[i];
//        if (f->index_ != i) {
//            printf("v: %d %d\n", f->index_, i);
//        }
//
//        if (f->m_he->index_ >= halfedges_.size()) {
//            printf("v: %d %d\n", f->m_he->index_, halfedges_.size());
//        }
//
//        const int i0 = f->m_he->src()->index();
//        const int i1 = f->m_he->next()->src()->index();
//        const int i2 = f->m_he->prev()->src()->index();
//        if (i0 >= vertices_.size() || i1 >= vertices_.size() || i2 >= vertices_.size()) {
//            printf("%d %d %d, %d\n", i0, i1, i2, vertices_.size());
//        }
//    }

    return success;
}

bool Mesh::verifyVertex(Vertex* v) const {
    bool success = true;
    if (v->index() >= vertices_.size()) {
        fprintf(stderr, "Vertex index out of range: %d > %d\n", v->index(), (int)vertices_.size());
        success = false;
    }

    Vec p = v->pos();
    if (std::isinf(p.x) || std::isnan(p.x) || std::isinf(p.y) || std::isnan(p.y) || std::isinf(p.z) || std::isnan(p.z)) {
        fprintf(stderr, "Inf of NaN found at v[%d]: (%f, %f, %f)\n", v->index(), p.x, p.y, p.z);
        success = false;
    }
    return success;
}

Mesh::VertexIterator Mesh::v_begin() {
    return Mesh::VertexIterator(vertices_);
}

Mesh::VertexIterator Mesh::v_end() {
    return Mesh::VertexIterator(vertices_, vertices_.size());
}

Mesh::HalfedgeIterator Mesh::he_begin() {
    return Mesh::HalfedgeIterator(halfedges_);
}

Mesh::HalfedgeIterator Mesh::he_end() {
    return Mesh::HalfedgeIterator(halfedges_, halfedges_.size());
}

Mesh::FaceIterator Mesh::f_begin() {
    return Mesh::FaceIterator(faces_);
}

Mesh::FaceIterator Mesh::f_end() {
    return Mesh::FaceIterator(faces_, faces_.size());
}

void Mesh::addVertex(Vertex *v) {
    v->index_ = vertices_.size();
    vertices_.emplace_back(v);
}

void Mesh::addHalfedge(Halfedge *he) {
    he->index_ = halfedges_.size();
    halfedges_.emplace_back(he);
}

void Mesh::addFace(Face *f) {
    f->index_ = faces_.size();
    faces_.emplace_back(f);
}

void Mesh::removeVertex(Vertex* v) {
    if (v->index_ < vertices_.size() - 1) {
        std::swap(vertices_[v->index_], vertices_[vertices_.size() - 1]);
        std::swap(vertices_[v->index_]->index_ , vertices_[vertices_.size() - 1]->index_);
    } 
    vertices_.resize(vertices_.size() - 1);
}

void Mesh::removeHalfedge(Halfedge *he) {
    if (he->index_ < halfedges_.size() - 1) {
        std::swap(halfedges_[he->index_], halfedges_[halfedges_.size() - 1]);
        std::swap(halfedges_[he->index_]->index_, halfedges_[halfedges_.size() - 1]->index_);
    }
    halfedges_.resize(halfedges_.size() - 1);
}

void Mesh::removeFace(Face *f) {
    if (f->index_ < faces_.size() - 1) {
        std::swap(faces_[f->index_], faces_[faces_.size() - 1]);
        std::swap(faces_[f->index_]->index_, faces_[faces_.size() - 1]->index_);
    }
    faces_.resize(faces_.size() - 1);
}


// ----------
// VertexIterator
// ----------

Mesh::VertexIterator::VertexIterator(std::vector<std::shared_ptr<Vertex>> &vertices, int index)
    : vertices_{ vertices }
    , index_{ index } {
}

bool Mesh::VertexIterator::operator!=(const Mesh::VertexIterator &it) const {
    return index_ != it.index_;
}

Vertex &Mesh::VertexIterator::operator*() {
    return *vertices_[index_];
}

Vertex *Mesh::VertexIterator::operator->() const {
    return index_ < vertices_.size() ? vertices_[index_].get() : nullptr;
}

Vertex* Mesh::VertexIterator::ptr() const {
    return index_ < vertices_.size() ? vertices_[index_].get() : nullptr;
}

Mesh::VertexIterator &Mesh::VertexIterator::operator++() {
    ++index_;
    return *this;
}

Mesh::VertexIterator Mesh::VertexIterator::operator++(int) {
    int prev = index_;
    ++index_;
    return Mesh::VertexIterator(vertices_, prev);
}

Mesh::VertexIterator &Mesh::VertexIterator::operator--() {
    --index_;
    return *this;
}

Mesh::VertexIterator Mesh::VertexIterator::operator--(int) {
    int prev = index_;
    --index_;
    return Mesh::VertexIterator(vertices_, prev);
}

// ----------
// HalfedgeIterator
// ----------

Mesh::HalfedgeIterator::HalfedgeIterator(std::vector<std::shared_ptr<Halfedge>>& hes, int index)
    : halfedges_{ hes }
    , index_{ index } {    
}

bool Mesh::HalfedgeIterator::operator!=(const HalfedgeIterator& it) const {
    return index_ != it.index_;
}

Halfedge &Mesh::HalfedgeIterator::operator*() {
    return *halfedges_[index_];
}

Halfedge* Mesh::HalfedgeIterator::operator->() const {
    return index_ < halfedges_.size() ? halfedges_[index_].get() : nullptr;
}

Halfedge* Mesh::HalfedgeIterator::ptr() const {
    return index_ < halfedges_.size() ? halfedges_[index_].get() : nullptr;
}

Mesh::HalfedgeIterator &Mesh::HalfedgeIterator::operator++() {
    ++index_;
    return *this;
}

Mesh::HalfedgeIterator Mesh::HalfedgeIterator::operator++(int) {
    int prev = index_;
    ++index_;
    return Mesh::HalfedgeIterator(halfedges_, prev);
}

Mesh::HalfedgeIterator &Mesh::HalfedgeIterator::operator--() {
    --index_;
    return *this;
}

Mesh::HalfedgeIterator Mesh::HalfedgeIterator::operator--(int) {
    int prev = index_;
    --index_;
    return Mesh::HalfedgeIterator(halfedges_, prev);
}

// ----------
// FaceIterator
// ----------

Mesh::FaceIterator::FaceIterator(std::vector<std::shared_ptr<Face>> &faces, int index)
    : faces_{ faces }
    , index_{ index } {
}

bool Mesh::FaceIterator::operator!=(const FaceIterator &it) const {
    return index_ != it.index_;
}

Face &Mesh::FaceIterator::operator*() {
    return *faces_[index_];
}

Face *Mesh::FaceIterator::operator->() const {
    return index_ < faces_.size() ? faces_[index_].get() : nullptr;
}

Face *Mesh::FaceIterator::ptr() const {
    return index_ < faces_.size() ? faces_[index_].get() : nullptr;
}

Mesh::FaceIterator &Mesh::FaceIterator::operator++() {
    ++index_;
    return *this;
}

Mesh::FaceIterator Mesh::FaceIterator::operator++(int) {
    int prev = index_;
    ++index_;
    return Mesh::FaceIterator(faces_, prev);
}

Mesh::FaceIterator &Mesh::FaceIterator::operator--() {
    --index_;
    return *this;
}

Mesh::FaceIterator Mesh::FaceIterator::operator--(int) {
    int prev = index_;
    --index_;
    return Mesh::FaceIterator(faces_, prev);
}

}  // namespace tinymesh
