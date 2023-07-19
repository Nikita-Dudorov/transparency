// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#define _USE_MATH_DEFINES

#include <cmath>
#include <algorithm>

#include "Mesh.h"

using namespace std;

Mesh::~Mesh () {
	clear ();
}

void Mesh::recomputePerVertexNormals (bool angleBased) {
	m_vertexNormals.clear ();
	// Change the following code to compute a proper per-vertex normal
	m_vertexNormals.resize (m_vertexPositions.size (), glm::vec3 (0.0, 0.0, 0.0));
	for (auto & t : m_triangleIndices) {
		glm::vec3 e0 (m_vertexPositions[t[1]] - m_vertexPositions[t[0]]);
		glm::vec3 e1 (m_vertexPositions[t[2]] - m_vertexPositions[t[0]]);
		glm::vec3 n = normalize (cross (e0, e1));
		for (size_t i = 0; i < 3; i++)
			m_vertexNormals[t[i]] += n;
	}
	for (auto & n : m_vertexNormals)
		n = normalize (n);
}

void Mesh::clear () {
	m_vertexPositions.clear ();
	m_vertexNormals.clear ();
	m_triangleIndices.clear ();
}