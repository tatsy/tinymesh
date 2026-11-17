#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include "tinymesh/tinymesh.h"

namespace {

using Clock = std::chrono::steady_clock;
using tinymesh::Mesh;

struct Input {
    std::vector<Vec3> vertices;
    std::vector<uint32_t> indices;
    size_t holes = 0;
    size_t boundaryHalfedges = 0;
};

double seconds(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double>(end - start).count();
}

uint32_t number(const std::string &text) {
    size_t used = 0;
    const auto value = std::stoull(text, &used);
    if (used != text.size() || text.empty() || text.front() == '-' ||
        value > std::numeric_limits<uint32_t>::max()) {
        throw std::invalid_argument("Expected an unsigned 32-bit integer: " + text);
    }
    return static_cast<uint32_t>(value);
}

Input generate(const std::string &shape, uint32_t n) {
    const size_t cells = static_cast<size_t>(n) * n;
    // Mesh indices and element counts currently use signed int internally.
    if (n < 3 || cells > static_cast<size_t>(std::numeric_limits<int>::max()) / 8) {
        throw std::invalid_argument("Size must be >= 3 and keep all element counts within signed int range.");
    }

    Input input;
    if (shape == "grid") {
        input.vertices.reserve(static_cast<size_t>(n + 1) * (n + 1));
        input.indices.reserve(6 * cells);
        for (uint32_t y = 0; y <= n; ++y) {
            for (uint32_t x = 0; x <= n; ++x) {
                input.vertices.emplace_back(x, y, 0.0);
            }
        }
        for (uint32_t y = 0; y < n; ++y) {
            for (uint32_t x = 0; x < n; ++x) {
                const uint32_t a = y * (n + 1) + x;
                const uint32_t b = a + 1, c = a + n + 1, d = c + 1;
                input.indices.insert(input.indices.end(), { a, b, d, a, d, c });
            }
        }
        input.holes = 1;
        input.boundaryHalfedges = 4 * n;
    } else if (shape == "torus") {
        input.vertices.reserve(cells);
        input.indices.reserve(6 * cells);
        for (uint32_t y = 0; y < n; ++y) {
            const double v = 2 * Pi * y / n;
            for (uint32_t x = 0; x < n; ++x) {
                const double u = 2 * Pi * x / n;
                const double r = 3.0 + std::cos(v);
                input.vertices.emplace_back(r * std::cos(u), r * std::sin(u), std::sin(v));
            }
        }
        for (uint32_t y = 0; y < n; ++y) {
            for (uint32_t x = 0; x < n; ++x) {
                const uint32_t a = y * n + x, b = y * n + (x + 1) % n;
                const uint32_t c = ((y + 1) % n) * n + x, d = ((y + 1) % n) * n + (x + 1) % n;
                input.indices.insert(input.indices.end(), { a, b, d, a, d, c });
            }
        }
    } else if (shape == "triangles") {
        input.vertices.reserve(3 * cells);
        input.indices.reserve(3 * cells);
        for (uint32_t y = 0; y < n; ++y) {
            for (uint32_t x = 0; x < n; ++x) {
                const auto a = static_cast<uint32_t>(input.vertices.size());
                input.vertices.emplace_back(3.0 * x, 3.0 * y, 0.0);
                input.vertices.emplace_back(3.0 * x + 1, 3.0 * y, 0.0);
                input.vertices.emplace_back(3.0 * x, 3.0 * y + 1, 0.0);
                input.indices.insert(input.indices.end(), { a, a + 1, a + 2 });
            }
        }
        input.holes = cells;
        input.boundaryHalfedges = 3 * cells;
    } else if (shape == "fan") {
        input.vertices.reserve(cells + 1);
        input.indices.reserve(3 * cells);
        input.vertices.emplace_back(0.0, 0.0, 0.0);
        for (size_t i = 0; i < cells; ++i) {
            const double theta = 2 * Pi * i / cells;
            input.vertices.emplace_back(std::cos(theta), std::sin(theta), 0.0);
        }
        for (uint32_t i = 0; i < cells; ++i) {
            input.indices.insert(input.indices.end(), { 0, i + 1, (i + 1) % static_cast<uint32_t>(cells) + 1 });
        }
        input.holes = 1;
        input.boundaryHalfedges = cells;
    } else {
        throw std::invalid_argument("Unknown case: " + shape);
    }
    return input;
}

void shuffleFaces(Input &input, uint32_t seed) {
    std::vector<size_t> order(input.indices.size() / 3);
    std::iota(order.begin(), order.end(), 0);
    std::mt19937 random(seed);
    std::shuffle(order.begin(), order.end(), random);
    std::vector<uint32_t> indices;
    indices.reserve(input.indices.size());
    for (const auto face : order) {
        indices.insert(indices.end(), input.indices.begin() + 3 * face, input.indices.begin() + 3 * face + 3);
    }
    input.indices.swap(indices);
}

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

size_t validate(const Mesh &mesh, const Input *input) {
    require(mesh.verify(), "Mesh::verify failed");
    size_t holes = 0, realFaces = 0, boundary = 0;
    for (size_t i = 0; i < mesh.numFaces(); ++i) {
        const auto *face = mesh.face(static_cast<int>(i));
        if (face->isHole()) {
            ++holes;
        } else {
            ++realFaces;
            require(face->numCorners() == 3, "A real face is not triangular");
        }
    }
    for (size_t i = 0; i < mesh.numHalfedges(); ++i) {
        const auto *he = mesh.halfedge(static_cast<int>(i));
        require(he->rev() && he->next() && he->src() && he->face(), "Null halfedge link");
        require(he->rev()->rev() == he, "Reverse pairing is not reciprocal");
        require(he->src() == he->rev()->dst() && he->dst() == he->rev()->src(), "Reverse endpoints mismatch");
        require(he->next()->face() == he->face(), "Next halfedge has a different face");
        require(he->isBoundary() == he->rev()->isBoundary(), "Boundary flags mismatch");
        boundary += he->isBoundary();
    }
    require(mesh.numHalfedges() % 2 == 0, "Odd halfedge count");
    if (input) {
        require(mesh.numVertices() == input->vertices.size(), "Unexpected vertex welding or missing vertices");
        require(realFaces == input->indices.size() / 3, "Unexpected triangle count");
        require(holes == input->holes, "Unexpected boundary loop count");
        require(boundary == 2 * input->boundaryHalfedges, "Unexpected boundary halfedge count");
        require(mesh.numHalfedges() == 3 * realFaces + input->boundaryHalfedges, "Unexpected halfedge count");
    }
    return holes;
}

void peakRss() {
#if defined(__linux__) || defined(__APPLE__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
#if defined(__APPLE__)
        std::cout << usage.ru_maxrss;
#else
        std::cout << static_cast<uint64_t>(usage.ru_maxrss) * 1024;
#endif
        return;
    }
#endif
    std::cout << "null";
}

void info() {
    std::cout << "{\"build_config\":" << std::quoted(TINYMESH_BENCH_BUILD_CONFIG);
#if defined(__clang__)
    std::cout << ",\"compiler\":\"Clang\",\"compiler_version\":" << std::quoted(__clang_version__);
#elif defined(_MSC_VER)
    std::cout << ",\"compiler\":\"MSVC\",\"compiler_version\":\"" << _MSC_FULL_VER << '"';
#elif defined(__GNUC__)
    std::cout << ",\"compiler\":\"GCC\",\"compiler_version\":" << std::quoted(__VERSION__);
#endif
    std::cout << ",\"sizeof_vertex\":" << sizeof(tinymesh::Vertex)
              << ",\"sizeof_halfedge\":" << sizeof(tinymesh::Halfedge)
              << ",\"sizeof_face\":" << sizeof(tinymesh::Face);
#if defined(TINYMESH_DISABLE_ASSERT)
    std::cout << ",\"assertions_enabled\":false";
#else
    std::cout << ",\"assertions_enabled\":true";
#endif
    std::cout << "}\n";
}

}  // namespace

int main(int argc, char **argv) {
    try {
        std::string shape = "grid", order = "ordered", file;
        uint32_t size = 128, seed = 42;
        for (int i = 1; i < argc; ++i) {
            const std::string key = argv[i];
            if (key == "--info") {
                info();
                return 0;
            }
            if (i + 1 == argc) throw std::invalid_argument("Missing value for " + key);
            const std::string value = argv[++i];
            if (key == "--case") shape = value;
            else if (key == "--size") size = number(value);
            else if (key == "--seed") seed = number(value);
            else if (key == "--order") order = value;
            else if (key == "--file") file = value;
            else throw std::invalid_argument("Unknown option: " + key);
        }
        if (order != "ordered" && order != "shuffled") throw std::invalid_argument("Unknown face order");

        const auto begin = Clock::now();
        Input input;
        if (file.empty()) {
            input = generate(shape, size);
            if (order == "shuffled") shuffleFaces(input, seed);
        }
        const auto prepared = Clock::now();
        std::unique_ptr<Mesh> mesh;
        if (file.empty()) mesh = std::make_unique<Mesh>(input.vertices, input.indices);
        else mesh = std::make_unique<Mesh>(file);
        const auto built = Clock::now();

        // Capture the construction peak before the untimed validation pass.
        std::cout << std::setprecision(12) << "{\"peak_rss_bytes\":";
        peakRss();
        const auto checkBegin = Clock::now();
        const auto holes = validate(*mesh, file.empty() ? &input : nullptr);
        const auto checked = Clock::now();
        const auto vertices = mesh->numVertices(), faces = mesh->numFaces(), halfedges = mesh->numHalfedges();
        const auto destroyBegin = Clock::now();
        mesh.reset();
        const auto destroyed = Clock::now();
        std::cout << ",\"preparation_seconds\":" << seconds(begin, prepared)
                  << ",\"build_seconds\":" << seconds(prepared, built)
                  << ",\"verification_seconds\":" << seconds(checkBegin, checked)
                  << ",\"destruction_seconds\":" << seconds(destroyBegin, destroyed)
                  << ",\"input_vertices\":" << (file.empty() ? std::to_string(input.vertices.size()) : "null")
                  << ",\"input_triangles\":" << (file.empty() ? std::to_string(input.indices.size() / 3) : "null")
                  << ",\"vertices\":" << vertices << ",\"triangles\":" << faces - holes
                  << ",\"faces_including_holes\":" << faces << ",\"halfedges\":" << halfedges
                  << ",\"boundary_loops\":" << holes << ",\"valid\":true}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
