// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

class LightSource {
public:
	LightSource(const glm::vec3& direction = glm::vec3 (0.f, 0.f, -1.f), const glm::vec3& color = glm::vec3 (1.f, 1.f, 1.f), float intensity = 1.f) : m_direction(direction), m_color(color), m_intensity(intensity) {}

	inline const glm::vec3& direction() const { return m_direction; }

	inline const glm::vec3& color() const { return m_color; }

	inline float intensity() const { return m_intensity; }

private:
	glm::vec3 m_direction;
	glm::vec3 m_color;
	float m_intensity;
};