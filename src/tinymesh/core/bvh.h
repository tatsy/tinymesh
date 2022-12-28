#ifdef _MSC_VER
#pragma once
#endif  // _MSC_VER

#ifndef TINYMESH_BVH_H
#define TINYMESH_BVH_H

#include <vector>
#include <memory>

#include "api.h"
#include "vec.h"
#include "bounds.h"
#include "triangle.h"

namespace tinymesh {

struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() {
    }
    BVHPrimitiveInfo(int pid, const Bounds3 &b)
        : primIdx(pid)
        , centroid()
        , bounds(b) {
        centroid = (b.posMax() + b.posMin()) * 0.5;
    }

    int primIdx;
    Vec3 centroid;
    Bounds3 bounds;
};

struct BVHNode {
    void initLeaf(const Bounds3 &b, int pid) {
        this->bounds = b;
        this->splitAxis = 0;
        this->primIdx = pid;
    }

    void initFork(const Bounds3 &b, BVHNode *l, BVHNode *r, int axis) {
        this->bounds = b;
        this->left = l;
        this->right = r;
        this->splitAxis = axis;
        this->primIdx = -1;
    }

    bool isLeaf() const {
        return primIdx >= 0;
    }

    Bounds3 bounds;
    BVHNode *left = nullptr;
    BVHNode *right = nullptr;
    int splitAxis;
    int primIdx;
};

class TINYMESH_API BVH {
public:
    BVH() = default;
    BVH(const Mesh &mesh);
    BVH(const std::vector<Vec3> &vertices, const std::vector<uint32_t> &indices);
    void construct(const std::vector<Vec3> &vertices, const std::vector<uint32_t> &indices);

    double distance(const Vec3 &p) const;
    Vec3 closestPoint(const Vec3 &p) const;

private:
    BVHNode *constructRec(std::vector<BVHPrimitiveInfo> &buildData, int start, int end);

    BVHNode *root = nullptr;
    std::vector<Triangle> tris;
    std::vector<std::shared_ptr<BVHNode>> nodes;
};

}  // namespace tinymesh

#endif  // TINYMESH_BVH_H
