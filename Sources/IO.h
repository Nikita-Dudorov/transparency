// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#pragma once

#include <string>
#include <memory>

#include "Scene.h"
#include "Mesh.h"
#include "Material.h"

namespace IO {

/// Loads an OFF mesh file. See https://en.wikipedia.org/wiki/OFF_(file_format).
std::shared_ptr<Mesh> loadOFFMesh (const std::string & filename);

}
