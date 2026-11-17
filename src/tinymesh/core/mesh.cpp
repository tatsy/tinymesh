#include "mesh.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <queue>
#include <tuple>
#include <functional>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>

#include <Eigen/Core>

#include "core/debug.h"
#include "core/utils.h"
#include "core/vec.h"
#include "core/vertex.h"
#include "core/face.h"
#include "core/halfedge.h"
#include "restore/restore.h"
#include "tiny_obj_loader.h"
#include "tinyply.h"

namespace fs = std::filesystem;

namespace tinymesh {

// ----------
// Mesh
// ----------

Mesh::Mesh() {
}

Mesh::Mesh(const std::string &filename) {
    load(filename);
}

Mesh::Mesh(const std::vector<Vec3> &vertices, const std::vector<uint32_t> &indices) {
    construct(vertices, indices);
}

size_t Mesh::numVertices() const {
    return vertices_.size();
}

size_t Mesh::numEdges() const {
    std::unordered_set<Halfedge *> hes;
    for (size_t i = 0; i < halfedges_.size(); i++) {
        Halfedge *he = halfedges_[i].get();
        if (hes.count(he) != 0 || hes.count(he->rev_) != 0) continue;
        hes.insert(he);
    }
    return hes.size();
}

size_t Mesh::numHalfedges() const {
    return halfedges_.size();
}

size_t Mesh::numFaces() const {
    return faces_.size();
}

void Mesh::load(const std::string &filename) {
    // Clear existing elements
    vertices_.clear();
    halfedges_.clear();
    faces_.clear();
    indices_.clear();

    // Load new mesh
    const std::string ext = fs::path(filename.c_str()).extension().string();
    if (ext == ".obj") {
        loadOBJ(filename);
    } else if (ext == ".ply") {
        loadPLY(filename);
    } else {
        Error("Unsupported file extension: {}", ext.c_str());
    }

    construct();
}

void Mesh::construct(const std::vector<Vec3> &vertices, const std::vector<uint32_t> &indices) {
    // Clear existing elements
    vertices_.clear();
    halfedges_.clear();
    faces_.clear();
    indices_.clear();

    std::unordered_map<Vec3, uint32_t> uniqueVertices;
    for (uint32_t i : indices) {
        if (i >= vertices.size()) {
            Error("Vertex index {} is out of range.", i);
        }

        const Vec3 &pos = vertices[i];

        if (uniqueVertices.count(pos) == 0) {
            uniqueVertices[pos] = static_cast<uint32_t>(vertices_.size());
            vertices_.push_back(std::make_shared<Vertex>(pos));
        }
        indices_.push_back(uniqueVertices[pos]);
    }

    construct();
}

Mesh Mesh::clone() {
    return Mesh(this->getVertices(), this->getVertexIndices());
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
        if (f->isHole()) {
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

double Mesh::getMeanEdgeLength() const {
    std::unordered_set<Halfedge *> visited;
    double mean = 0.0;
    for (int i = 0; i < halfedges_.size(); i++) {
        Halfedge *he = halfedges_[i].get();
        if (visited.count(he) != 0) continue;
        if (visited.count(he->rev_) != 0) continue;
        visited.insert(he);

        mean += halfedges_[i]->length();
    }
    return mean / (double)visited.size();
}

double Mesh::getMeanDihedralAngle() const {
    std::unordered_set<Halfedge *> visited;
    double mean = 0.0;
    for (int i = 0; i < halfedges_.size(); i++) {
        Halfedge *he = halfedges_[i].get();
        if (visited.count(he) != 0) continue;
        if (visited.count(he->rev_) != 0) continue;
        visited.insert(he);

        Halfedge *rev = he->rev_;
        Vertex *vh1 = he->src();
        Vertex *vh2 = he->next()->dst();
        Vertex *vh3 = rev->src();
        Vertex *vh4 = rev->next()->dst();

        const Vec3 p1 = vh1->pos();
        const Vec3 p2 = vh2->pos();
        const Vec3 p3 = vh3->pos();
        const Vec3 p4 = vh4->pos();
        mean += dihedral(p2, p1, p3, p4);
    }
    return mean / (double)visited.size();
}

double Mesh::getMeanFaceArea() const {
    double totalArea = 0.0;
    size_t count = 0;
    for (const auto &f : faces_) {
        if (f->isHole()) {
            continue;
        }
        totalArea += f->area();
        count++;
    }

    if (count == 0) return 0.0;

    return totalArea / static_cast<double>(count);
}

void Mesh::construct() {
    if (indices_.size() % 3 != 0) {
        Error("Triangle index count must be a multiple of three: {}", indices_.size());
    }

    // Put vertex indices
    for (int i = 0; i < (int)vertices_.size(); i++) {
        vertices_[i]->index_ = i;
    }

    // Check if all triangles consist of three distinct vertices
    for (size_t i = 0; i < indices_.size(); i += 3) {
        if (indices_[i + 0] == indices_[i + 1] || indices_[i + 1] == indices_[i + 2] ||
            indices_[i + 2] == indices_[i + 0]) {
            Error("Each triangle must have three distinct vertices!");
        }
    }

    // Setup halfedge structure
    static const int ngon = 3;
    std::unordered_map<IndexPair, Halfedge *, IndexPairHash> pairToHalfedge;
    std::unordered_map<uint32_t, uint32_t> vertexDegree;
    pairToHalfedge.reserve(indices_.size());
    vertexDegree.reserve(vertices_.size());
    for (int i = 0; i < (int)indices_.size(); i += ngon) {
        // Check polygon duplication
        bool isDuplicated = false;
        for (int j = 0; j < ngon; j++) {
            const uint32_t a = indices_[i + j];
            const uint32_t b = indices_[i + (j + 1) % ngon];

            if (pairToHalfedge.contains(IndexPair{ a, b })) {
                if (pairToHalfedge.contains(IndexPair{ b, a })) {
                    Error("Non-manifold edge #{}--#{}: more than two incident faces (face #{})", a, b, i / ngon);
                }

                Error("Inconsistent winding or duplicated face at edge #{}--#{} (face #{})", a, b, i / ngon);
            }
        }

        // Traverse face vertices
        auto face = new Face();
        std::vector<Halfedge *> faceHalfedges;
        for (int j = 0; j < ngon; j++) {
            // Count up #faces around each vertex
            const uint32_t a = indices_[i + j];
            if (vertexDegree.find(a) == vertexDegree.end()) {
                vertexDegree[a] = 1;
            } else {
                vertexDegree[a] += 1;
            }

            // Check duplicated halfedges
            const uint32_t b = indices_[i + (j + 1) % ngon];
            IndexPair ab(a, b);
            if (pairToHalfedge.find(ab) != pairToHalfedge.end()) {
                Error("An edge with vertices #{} and #{} is duplicated\n", a, b);
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
            } else {
                he->rev_ = nullptr;
            }
        }
        addFace(face);

        // Set next halfedges
        for (int j = 0; j < ngon; j++) {
            const int k = (j + 1) % ngon;
            faceHalfedges[j]->next_ = faceHalfedges[k];
        }
    }

    // Find outgoing boundary halfedges.
    std::vector<Halfedge *> boundaryStart(vertices_.size(), nullptr);

    for (const auto &he : halfedges_) {
        if (he->rev_ != nullptr) {
            continue;
        }

        const auto index = he->src_->index_;
        if (boundaryStart[index] != nullptr) {
            Error("Non-manifold vertex #{}: multiple outgoing boundary edges.", index);
        }

        boundaryStart[index] = he.get();
    }

    // Check connectivity of incident triangles around each vertex.
    for (const auto &v : vertices_) {
        if (v->halfedge_ == nullptr) {
            Error("Vertex #{} is not referenced by any triangle.", v->index_);
        }

        Halfedge *start = boundaryStart[v->index_];
        if (start == nullptr) {
            start = v->halfedge_;
        }

        const uint32_t expected = vertexDegree[v->index_];
        uint32_t count = 0;
        Halfedge *he = start;

        do {
            if (++count > expected) {
                Error("Non-manifold vertex #{}: invalid face fan.", v->index_);
            }

            he = he->next_->next_->rev_;
        } while (he != nullptr && he != start);

        if (count != expected) {
            Error("Non-manifold vertex #{}: disconnected face fans ({} of {} faces).", v->index_, count, expected);
        }

        v->halfedge_ = start;
    }

    // Construct new faces (possibly non-triangle) for each boundary components
    const auto numHalfedges = static_cast<int>(halfedges_.size());
    for (int i = 0; i < numHalfedges; i++) {
        Halfedge *const he = halfedges_[i].get();
        Assertion(he != nullptr, "Null halfedge detected!");
        if (he->rev_ == nullptr) {
            auto face = new Face();
            face->isHole_ = true;
            addFace(face);

            std::vector<Halfedge *> boundaryHalfedges;
            Halfedge *it = he;
            do {
                auto rev = new Halfedge();
                addHalfedge(rev);
                boundaryHalfedges.push_back(rev);

                it->rev_ = rev;
                rev->rev_ = it;

                rev->face_ = face;
                rev->src_ = it->next_->src_;

                it->isBoundary_ = true;
                rev->isBoundary_ = true;

                // Advance it to the next halfedge along the current boundary loop
                it = it->next_;
                while (it != he && it->rev_ != nullptr) {
                    it = it->rev_;
                    it = it->next_;
                }
            } while (it != he);

            face->halfedge_ = boundaryHalfedges.front();

            const auto degree = static_cast<int>(boundaryHalfedges.size());
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
        Warn("{}", errMsg);
    }

    if (!success) {
        Error("Failed to load *.obj file: {}", filename);
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
    using tinyply::PlyData;
    using tinyply::PlyFile;

    // Open
    std::ifstream reader(filename.c_str(), std::ios::binary);
    if (reader.fail()) {
        Error("Failed to open file: {}", filename);
    }

    // Read header
    PlyFile file;
    file.parse_header(reader);

    // Request vertex data
    std::shared_ptr<PlyData> vert_data, norm_data, uv_data, face_data;
    try {
        vert_data = file.request_properties_from_element("vertex", { "x", "y", "z" });
    } catch (std::exception &e) {
        Error("Failed to read PLY vertex data: {}", e.what());
    }

    try {
        norm_data = file.request_properties_from_element("vertex", { "nx", "ny", "nz" });
    } catch (std::exception &) {
        // std::cerr << "tinyply exception: " << e.what() << std::endl;
    }

    try {
        uv_data = file.request_properties_from_element("vertex", { "u", "v" });
    } catch (std::exception &) {
        // std::cerr << "tinyply exception: " << e.what() << std::endl;
    }

    try {
        face_data = file.request_properties_from_element("face", { "vertex_indices" }, 3);
    } catch (std::exception &e) {
        Error("Failed to read PLY face data: {}", e.what());
    }

    // Read vertex data
    file.read(reader);

    // Copy vertex data
    const size_t numVerts = vert_data->count;
    std::vector<uint8_t> raw_vertices;

    if (vert_data->t == tinyply::Type::FLOAT32) {
        raw_vertices.resize(numVerts * 3 * sizeof(float));
    } else if (vert_data->t == tinyply::Type::FLOAT64) {
        raw_vertices.resize(numVerts * 3 * sizeof(double));
    } else {
        Error("Unsupported vertex data type: {}", (int)vert_data->t);
    }
    std::memcpy(raw_vertices.data(), vert_data->buffer.get(), raw_vertices.size());

    const size_t numFaces = face_data->count;
    std::vector<uint8_t> raw_indices;
    if (face_data->t == tinyply::Type::UINT8) {
        raw_indices.resize(numFaces * 3 * sizeof(uint8_t));
    } else if (face_data->t == tinyply::Type::UINT16) {
        raw_indices.resize(numFaces * 3 * sizeof(uint16_t));
    } else if (face_data->t == tinyply::Type::INT32 || face_data->t == tinyply::Type::UINT32) {
        raw_indices.resize(numFaces * 3 * sizeof(uint32_t));
    } else {
        Error("Unsupported face data type: {}", (int)face_data->t);
    }
    std::memcpy(raw_indices.data(), face_data->buffer.get(), raw_indices.size());

    std::unordered_map<Vec3, uint32_t> uniqueVertices;
    vertices_.clear();
    for (size_t i = 0; i < numFaces * 3; i++) {
        uint32_t j = 0;
        if (face_data->t == tinyply::Type::UINT8) {
            j = reinterpret_cast<uint8_t *>(raw_indices.data())[i];
        } else if (face_data->t == tinyply::Type::UINT16) {
            j = reinterpret_cast<uint16_t *>(raw_indices.data())[i];
        } else if (face_data->t == tinyply::Type::INT32) {
            j = reinterpret_cast<int32_t *>(raw_indices.data())[i];
        } else if (face_data->t == tinyply::Type::UINT32) {
            j = reinterpret_cast<uint32_t *>(raw_indices.data())[i];
        } else {
            Error("Unsupported face data type: {}", (int)face_data->t);
        }

        if (j >= numVerts) {
            Error("Vertex index {} is out of range.", j);
        }

        Vec3 pos;
        if (vert_data->t == tinyply::Type::FLOAT32) {
            const float x = reinterpret_cast<float *>(raw_vertices.data())[j * 3 + 0];
            const float y = reinterpret_cast<float *>(raw_vertices.data())[j * 3 + 1];
            const float z = reinterpret_cast<float *>(raw_vertices.data())[j * 3 + 2];
            pos = Vec3(x, y, z);
        } else if (vert_data->t == tinyply::Type::FLOAT64) {
            const double x = reinterpret_cast<double *>(raw_vertices.data())[j * 3 + 0];
            const double y = reinterpret_cast<double *>(raw_vertices.data())[j * 3 + 1];
            const double z = reinterpret_cast<double *>(raw_vertices.data())[j * 3 + 2];
            pos = Vec3(x, y, z);
        } else {
            Error("Unsupported vertex data type: {}", (int)vert_data->t);
        }

        if (uniqueVertices.count(pos) == 0) {
            uniqueVertices[pos] = static_cast<uint32_t>(vertices_.size());
            vertices_.push_back(std::make_shared<Vertex>(pos));
        }
        indices_.push_back(uniqueVertices[pos]);
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
        Error("Unsupported file extension: {}", ext);
    }
}

void Mesh::saveOBJ(const std::string &filename) const {
    std::ofstream writer(filename.c_str(), std::ios::out);
    if (writer.fail()) {
        Error("Failed to open file: {}", filename);
    }

    for (const auto &v : vertices_) {
        writer << "v " << v->pos().x() << " " << v->pos().y() << " " << v->pos().z() << std::endl;
    }

    const int nv = (int)vertices_.size();
    for (const auto &f : faces_) {
        if (f->isHole()) {
            continue;
        }

        const int i0 = f->halfedge_->src()->index();
        const int i1 = f->halfedge_->next()->src()->index();
        const int i2 = f->halfedge_->prev()->src()->index();
        if (i0 < nv && i1 < nv && i2 < nv) {
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
        Error("Failed to open file: {}", filename);
    }

    std::vector<float> vertexData(vertices_.size() * 3);
    for (size_t i = 0; i < vertices_.size(); i++) {
        const Vec3 v = vertices_[i]->pos();
        vertexData[i * 3 + 0] = (float)v.x();
        vertexData[i * 3 + 1] = (float)v.y();
        vertexData[i * 3 + 2] = (float)v.z();
    }

    std::vector<uint32_t> indexData;
    indexData.reserve(faces_.size() * 3);
    for (const auto &f : faces_) {
        if (f->isHole()) {
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

    if (he->isBoundary() || he->rev_->isBoundary()) {
        return false;
    }

    // Prepare halfedge / face
    auto he_new = new Halfedge();
    auto he01 = new Halfedge();
    auto he02 = new Halfedge();
    auto rev_new = new Halfedge();
    auto he11 = new Halfedge();
    auto he12 = new Halfedge();

    auto he0 = he->next_;
    auto he1 = he0->next_;
    auto he2 = rev->next_;
    auto he3 = he2->next_;

    auto f0 = he->face();
    auto f1 = rev->face();

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
    addHalfedge(he_new);
    addHalfedge(rev_new);
    addHalfedge(he01);
    addHalfedge(he02);
    addHalfedge(he11);
    addHalfedge(he12);

    return true;
}

bool Mesh::collapseHE(Halfedge *he) {
    /*
     *  Collapse halfedge v0->v2.
     *         v1
     *       /  \
     *      /    \
     *    v0 ---- v2
     *      \    /
     *       \  /
     *        v3
     */

    Halfedge *rev = he->rev_;
    Assertion(he->index_ < (int)halfedges_.size(), "Daemon halfedge detected!");
    Assertion(rev->index_ < (int)halfedges_.size(), "Daemon halfedge detected!");

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

    Assertion(v0->index_ < (int)vertices_.size(), "Daemon vertex detected!");
    Assertion(v1->index_ < (int)vertices_.size(), "Daemon vertex detected!");
    Assertion(v2->index_ < (int)vertices_.size(), "Daemon vertex detected!");
    Assertion(v3->index_ < (int)vertices_.size(), "Daemon vertex detected!");

    if (v0->degree() <= 3 || v1->degree() <= 4 || v2->degree() <= 3 || v3->degree() <= 4) {
        return false;
    }

    // Collect neighboring vertices of the one to be removed
    std::vector<Vertex *> neighbors;
    for (auto it = v2->he_begin(); it != v2->he_end(); ++it) {
        Assertion(it->dst() != nullptr, "Null vertex is found!");
        neighbors.push_back(it->dst());
    }

    Assertion(!neighbors.empty(), "No neighboring vertex!");

    // Check unsafe collapse (# of shared vertices)
    int numUnion = 0;
    std::set<Vertex *> neighborSet(neighbors.begin(), neighbors.end());
    for (auto it = v0->he_begin(); it != v0->he_end(); ++it) {
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
    for (size_t i = 0; i < neighbors.size(); i++) {
        const size_t j = (i + 1) % neighbors.size();
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
    for (size_t i = 0; i < neighbors.size(); i++) {
        const size_t j = (i + 1) % neighbors.size();
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
    Assertion(v_remain->index() < (int)vertices_.size(), "Something is wrong1!");
    Assertion(v_remove->index() < (int)vertices_.size(), "Something is wrong2!");
    Assertion(v_remove != v_remain, "Something is wrong3!");
    for (auto it = v_remove->he_begin(); it != v_remove->he_end(); ++it) {
        Assertion(it->src() == v_remove, "Invalid halfedge origin detected!");
        it->src_ = v_remain;
    }

    for (auto it = v_remain->he_begin(); it != v_remain->he_end(); ++it) {
        Assertion(it->src() == v_remain, "Invalid halfedge origin detected!");
    }

    // Update edges and opposite halfedges
    Halfedge *he0 = he->next()->next()->rev();
    Halfedge *he1 = he->next()->rev();
    he0->rev_ = he1;
    he1->rev_ = he0;

    Halfedge *he2 = rev->next()->next()->rev();
    Halfedge *he3 = rev->next()->rev();
    he2->rev_ = he3;
    he3->rev_ = he2;

    // Update halfedge of the origin vertex (v0)
    v_remain->halfedge_ = he0;

    // Update halfedge for wing-vertices
    v0->halfedge_ = he0;
    v1->halfedge_ = he1;
    v2->halfedge_ = he2;
    v3->halfedge_ = he3;

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

bool Mesh::collapseFace(Face *f) {
    std::vector<Halfedge *> hes;
    for (auto it = f->he_begin(); it != f->he_end(); ++it) {
        hes.push_back(it.ptr());
    }
    Assertion(hes.size() == 3, "Non-triangle face detected!");

    Halfedge *he1 = hes[0];
    Halfedge *he2 = he1->next_->rev_;
    bool ret = collapseHE(he1);
    Assertion(he2->index() >= 0 && he2->index() < numHalfedges(), "Error!!");

    ret &= collapseHE(he2);

    return ret;

    // if (!collapseHE(he1) || !collapseHE(he2)) {
    //     return false;
    // }

    // return true;
}

bool Mesh::flipHE(Halfedge *he) {
    Halfedge *rev = he->rev_;
    if (!rev) {
        Error("Flip is called boundary halfedge!");
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

void Mesh::fillHoles() {
    std::vector<int> holeFaces;
    for (int i = 0; i < (int)this->numFaces(); i++) {
        Face *f = this->face(i);
        if (f->isHole()) {
            holeFaces.push_back(i);
        }
    }

    for (auto i : holeFaces) {
        Face *f = this->face(i);
        this->triangulate(f);
    }
}

void Mesh::triangulate(Face *face) {
    this->holeFillMinDihedral_(face);
}

bool Mesh::verify() const {
    bool success = true;
    for (int i = 0; i < (int)vertices_.size(); i++) {
        Vertex *v = vertices_[i].get();
        if (v->index() != i) {
            Warn("Vertex index does not match array index: v[{}].index = {}", i, v->index());
            success = false;
        }
        success &= verifyVertex(v);
    }

    for (int i = 0; i < (int)halfedges_.size(); i++) {
        auto he = halfedges_[i];
        if (he->index() != i) {
            Warn("Halfedge index does not match array index: he[{}].index = {}", i, he->index());
            success = false;
        }
    }

    for (int i = 0; i < (int)faces_.size(); i++) {
        auto f = faces_[i];
        if (f->index() != i) {
            Warn("Face index does not match array index: f[{}].index = {}", i, f->index());
            success = false;
        }

        if (f->halfedge_ == nullptr) {
            Warn("Face[{}] is not with non-null halfedge!", f->index());
            success = false;
        }
    }

    return success;
}

bool Mesh::verifyVertex(Vertex *v) const {
    bool success = true;
    if (v->index() >= (int)vertices_.size()) {
        Warn("Vertex index out of range: {} > {}", v->index(), (int)vertices_.size());
        success = false;
    }

    Vec3 p = v->pos();
    if (std::isinf(p[0]) || std::isnan(p[0]) || std::isinf(p[1]) || std::isnan(p[1]) || std::isinf(p[2]) ||
        std::isnan(p[2])) {
        Warn("Inf of NaN found at v[{}]: ({}, {}, {})\n", v->index(), p[0], p[1], p[2]);
        success = false;
    }

    for (auto it = v->he_begin(); it != v->he_end(); ++it) {
        if (it->src() != v) {
            Warn("Origin of an outward halfedge is invalid at v[{}]\n", v->index());
            success = false;
        }
    }

    return success;
}

void Mesh::addVertex(Vertex *v) {
    v->index_ = static_cast<int>(vertices_.size());
    vertices_.emplace_back(v);
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
    const int nVerts = (int)vertices_.size();
    Assertion((*v) == (*vertices_[v->index_]), "Invalid vertex indexing!");
    Assertion(v->index_ < nVerts, "Vertex index out of range!");

    if (v->index_ < nVerts - 1) {
        std::swap(vertices_[v->index_], vertices_[nVerts - 1]);
        std::swap(vertices_[v->index_]->index_, vertices_[nVerts - 1]->index_);
        vertices_[nVerts - 1]->index_ = -1;
    }
    vertices_.resize(nVerts - 1);
}

void Mesh::removeHalfedge(Halfedge *he) {
    const int nHEs = (int)halfedges_.size();
    Assertion((*he) == (*halfedges_[he->index_]), "Invalid halfedge indexing!");
    Assertion(he->index_ < nHEs, "Halfedge index out of range!");

    if (he->index_ < nHEs - 1) {
        std::swap(halfedges_[he->index_], halfedges_[nHEs - 1]);
        std::swap(halfedges_[he->index_]->index_, halfedges_[nHEs - 1]->index_);
        halfedges_[nHEs - 1]->index_ = -1;
    }
    halfedges_.resize(nHEs - 1);
}

void Mesh::removeFace(Face *f) {
    const int nFaces = (int)faces_.size();
    Assertion((*f) == (*faces_[f->index_]), "Invalid face indexing!");
    Assertion(f->index_ < nFaces, "Face index out of range!");

    if (f->index_ < nFaces - 1) {
        std::swap(faces_[f->index_], faces_[nFaces - 1]);
        std::swap(faces_[f->index_]->index_, faces_[nFaces - 1]->index_);
        faces_[nFaces - 1]->index_ = -1;
    }
    faces_.resize(faces_.size() - 1);
}

}  // namespace tinymesh
