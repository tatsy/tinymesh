#define TINYMESH_API_EXPORT
#include "bvh.h"

#include <stack>

#include "mesh.h"

namespace tinymesh {

struct BucketInfo {
    int count;
    Bounds3 bounds;
    BucketInfo()
        : count(0)
        , bounds() {
    }
};

struct ComparePoint {
    int dim;
    explicit ComparePoint(int d)
        : dim(d) {
    }
    bool operator()(const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) const {
        return a.centroid[dim] < b.centroid[dim];
    }
};

struct CompareToBucket {
    int splitBucket, nBuckets, dim;
    const Bounds3 &centroidBounds;

    CompareToBucket(int split, int num, int d, const Bounds3 &b)
        : splitBucket(split)
        , nBuckets(num)
        , dim(d)
        , centroidBounds(b) {
    }

    bool operator()(const BVHPrimitiveInfo &p) const {
        const double cmin = centroidBounds.posMin()[dim];
        const double cmax = centroidBounds.posMax()[dim];
        const double inv = (1.0) / (std::abs(cmax - cmin) + 1.0e-12);
        const double diff = std::abs(p.centroid[dim] - cmin);
        int b = static_cast<int>(nBuckets * diff * inv);
        if (b >= nBuckets) {
            b = nBuckets - 1;
        }
        return b <= splitBucket;
    }
};

BVH::BVH(const Mesh &mesh)
    : BVH{ mesh.getVertices(), mesh.getVertexIndices() } {
}

BVH::BVH(const std::vector<Vec3> &vertices, const std::vector<uint32_t> &indices) {
    construct(vertices, indices);
}

void BVH::construct(const std::vector<Vec3> &vertices, const std::vector<uint32_t> &indices) {
    tris.clear();
    std::vector<BVHPrimitiveInfo> buildData;
    for (size_t i = 0; i < indices.size(); i += 3) {
        const uint32_t i0 = indices[i + 0];
        const uint32_t i1 = indices[i + 1];
        const uint32_t i2 = indices[i + 2];
        const size_t idx = tris.size();
        tris.emplace_back(vertices[i0], vertices[i1], vertices[i2]);
        buildData.emplace_back(idx, tris[idx].bounds());
    }

    root = constructRec(buildData, 0, tris.size());
}

double BVH::distance(const Vec3 &p) const {
    const Vec3 pc = closestPoint(p);
    return length(pc - p);
}

Vec3 BVH::closestPoint(const Vec3 &query) const {
    std::stack<BVHNode *> stk;
    stk.push(root);

    double dist = 1.0e20;
    Vec3 ret;
    while (!stk.empty()) {
        BVHNode *node = stk.top();
        stk.pop();

        if (node->isLeaf()) {
            const Vec3 pc = tris[node->primIdx].closestPoint(query);
            const double dc = length(query - pc);
            if (dc < dist) {
                dist = dc;
                ret = pc;
            }
        } else {
            // distance to furthest corner
            const bool isInside = node->bounds.inside(query);
            const double dc = node->bounds.distance(query);
            if (isInside || dc <= dist) {
                if (node->left) stk.push(node->left);
                if (node->right) stk.push(node->right);
            }
        }
    }
    return ret;
}

BVHNode *BVH::constructRec(std::vector<BVHPrimitiveInfo> &buildData, int start, int end) {
    if (start == end) return nullptr;

    BVHNode *node = new BVHNode();
    nodes.emplace_back(node);

    Bounds3 bounds;
    for (int i = start; i < end; i++) {
        bounds = Bounds3::merge(bounds, buildData[i].bounds);
    }

    int nprims = end - start;
    if (nprims == 1) {
        // Leaf node
        node->initLeaf(bounds, buildData[start].primIdx);
    } else {
        // Fork node
        Bounds3 centroidBounds;
        for (int i = start; i < end; i++) {
            centroidBounds.merge(buildData[i].centroid);
        }

        int splitAxis = centroidBounds.maxExtent();
        int mid = (start + end) / 2;
        if (nprims <= 8) {
            std::nth_element(buildData.begin() + start, buildData.begin() + mid, buildData.begin() + end,
                             ComparePoint(splitAxis));
        } else {
            // Seperate with SAH (surface area heuristics)
            const int nBuckets = 16;
            BucketInfo buckets[nBuckets];

            const double cmin = centroidBounds.posMin()[splitAxis];
            const double cmax = centroidBounds.posMax()[splitAxis];
            const double idenom = 1.0 / (std::abs(cmax - cmin) + 1.0e-12);
            for (int i = start; i < end; i++) {
                const double numer = buildData[i].centroid[splitAxis] - centroidBounds.posMin()[splitAxis];
                int b = static_cast<int>(nBuckets * std::abs(numer) * idenom);
                if (b == nBuckets) {
                    b = nBuckets - 1;
                }

                buckets[b].count++;
                buckets[b].bounds = Bounds3::merge(buckets[b].bounds, buildData[i].bounds);
            }

            double bucketCost[nBuckets - 1] = { 0 };
            for (int i = 0; i < nBuckets - 1; i++) {
                Bounds3 b0, b1;
                int cnt0 = 0, cnt1 = 0;
                for (int j = 0; j <= i; j++) {
                    b0 = Bounds3::merge(b0, buckets[j].bounds);
                    cnt0 += buckets[j].count;
                }
                for (int j = i + 1; j < nBuckets; j++) {
                    b1 = Bounds3::merge(b1, buckets[j].bounds);
                    cnt1 += buckets[j].count;
                }
                bucketCost[i] += 0.125 + (cnt0 * b0.area() + cnt1 * b1.area()) / bounds.area();
            }

            double minCost = bucketCost[0];
            int minCostSplit = 0;
            for (int i = 1; i < nBuckets - 1; i++) {
                if (minCost > bucketCost[i]) {
                    minCost = bucketCost[i];
                    minCostSplit = i;
                }
            }

            if (minCost < nprims) {
                auto it = std::partition(buildData.begin() + start, buildData.begin() + end,
                                         CompareToBucket(minCostSplit, nBuckets, splitAxis, centroidBounds));
                mid = it - buildData.begin();
            }
        }

        BVHNode *left = constructRec(buildData, start, mid);
        BVHNode *right = constructRec(buildData, mid, end);
        node->initFork(bounds, left, right, splitAxis);
    }
    return node;
}

}  // namespace tinymesh
