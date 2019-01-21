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
    Mesh();
    Mesh(const std::string &filename);
    void load(const std::string &filename);

private:
    std::vector<std::shared_ptr<Point>> m_pts;
    std::vector<std::shared_ptr<Halfedge>> m_hes;
    std::vector<std::shared_ptr<Face>> m_faces;
    std::vector<uint32_t> m_indices;
};

}  // namespace tinymesh

#endif  // TINYMESH_MESH_H
