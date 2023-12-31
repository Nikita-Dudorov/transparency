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
#include <algorithm>
#include <cmath>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

class Material {
public:
	Material (const glm::vec3& diffuseAlbedo = glm::vec3(0.6f, 0.6f, 0.6f), 
			  float roughness = 0.5f, 
			  float metallicness = 0.f,
			  float transparency = 0.f, 
			  float refraction = 1.f,
			  float reflectance = 0.f)
		: m_albedo (diffuseAlbedo), 
		  m_roughness (roughness), 
		  m_metallicness (metallicness),
		  m_transparency (transparency),
		  m_refraction (refraction),
		  m_reflectance (reflectance) {}

	virtual ~Material () {}

	inline glm::vec3 albedo () const { return m_albedo; }
	inline void setAlbedo (glm::vec3 albedo) { m_albedo = albedo; }
	inline float roughness () const { return m_roughness; }
	inline void setRoughness (float roughness) { m_roughness = roughness; }
	inline float metallicness () const { return m_metallicness; }
	inline void setMetallicness (float metallicness) { m_metallicness = metallicness; }
	inline float transparency () const { return m_transparency; }
	inline void setTransparency (float transparency) { m_transparency = transparency; }
	inline float refraction () const { return m_refraction; }
	inline void setRefraction (float refraction) { m_refraction = refraction; }
	inline float reflectance () const { return m_reflectance; }
	inline void setReflectance (float reflectance) { m_reflectance = reflectance; }

private:
	glm::vec3 m_albedo;
	float m_roughness;
	float m_metallicness;
	float m_transparency = 0; // opaque
	float m_refraction = 1; // no refraction 
	float m_reflectance = 0; // no reflection
};

