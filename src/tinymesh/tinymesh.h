#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_H
#define TINYMESH_H

#include "core/api.h"
#include "core/debug.h"
#include "core/vec.h"
#include "core/openmp.h"
#include "core/filesystem.h"

#include "core/mesh.h"
#include "core/vertex.h"
#include "core/edge.h"
#include "core/halfedge.h"
#include "core/face.h"

#include "ops/ops.h"
#include "remesh/remesh.h"
#include "filters/filters.h"
#include "restore/restore.h"

#endif  // TINYMESH_H
