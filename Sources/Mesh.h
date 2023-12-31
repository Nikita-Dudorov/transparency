// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#pragma once

#include <vector>
#include <memory>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "Transform.h"

class Mesh : public Transform {
public:
	virtual ~Mesh ();

	inline const std::vector<glm::vec3> & vertexPositions () const { return m_vertexPositions; } 
	inline std::vector<glm::vec3> & vertexPositions () { return m_vertexPositions; }
	inline const std::vector<glm::vec3> & vertexNormals () const { return m_vertexNormals; } 
	inline std::vector<glm::vec3> & vertexNormals () { return m_vertexNormals; } 
	inline const std::vector<glm::vec2> & vertexTextCoords () const { return m_vertexTextCoords; } 
	inline std::vector<glm::vec2> & vertexTextCoords () { return m_vertexTextCoords; } 
	inline const std::vector<glm::uvec3> & triangleIndices () const { return m_triangleIndices; }
	inline std::vector<glm::uvec3> & triangleIndices () { return m_triangleIndices; }

	void recomputePerVertexNormals (bool angleBased = false);

	void clear ();

private:
	std::vector<glm::vec3> m_vertexPositions;
	std::vector<glm::vec3> m_vertexNormals;
	std::vector<glm::vec2> m_vertexTextCoords;
	std::vector<glm::uvec3> m_triangleIndices;
};
