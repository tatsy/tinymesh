#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_TEST_UTILS_H
#define TINYMESH_TEST_UTILS_H

#include "gtest/gtest.h"
#include "test_config.h"

namespace fs = std::filesystem;

class TinyMeshBaseTest : public ::testing::Test {
protected:
    TinyMeshBaseTest() {
    }
    virtual ~TinyMeshBaseTest() {
    }
};

class TinyMeshBaseTestWithParam : public TinyMeshBaseTest, public ::testing::WithParamInterface<std::string> {
public:
    TinyMeshBaseTestWithParam() {
        fs::path modelDir(MODEL_DIRECTORY);
        fs::path filePath(GetParam().c_str());
        filename = (modelDir / filePath).string();
    }

    virtual ~TinyMeshBaseTestWithParam() {
    }

protected:
    std::string filename;
    Mesh mesh;
};

#endif  // TINYMESH_TEST_UTILS_H
