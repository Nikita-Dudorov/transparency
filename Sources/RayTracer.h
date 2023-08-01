// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#pragma once

#include <cmath>
#include <algorithm>
#include <limits>
#include <memory>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "Renderer.h"
#include "Image.h"
#include "Scene.h"
#include "BVH.h"

using namespace std;

class RayTracer : public Renderer {
public:
	
	RayTracer();
	virtual ~RayTracer();

	inline void setResolution (int width, int height) { m_imagePtr = make_shared<Image> (width, height); }
	inline std::shared_ptr<Image> image () { return m_imagePtr; }
	inline const std::shared_ptr<Image> image () const { return m_imagePtr; }
	void init (const std::shared_ptr<Scene> scenePtr);
	void buildBVH(const std::shared_ptr<Scene> scene);
	virtual void render (const std::shared_ptr<Scene> scenePtr) final;

private:
	template<typename T>
	inline T barycentricInterpolation (const T & p0, const T & p1, const T & p2, float w, float u, float v) const {	return w * p0 + u * p1 + v * p2; }

	struct Hit {
		size_t m_meshIndex;
		size_t m_triangleIndex;
		float m_uCoord;
		float m_vCoord;
		float m_distance;
	};

	bool rayTrace(const Ray& ray, const std::shared_ptr<Scene> scene, size_t originMeshIndex, size_t originTriangleIndex, Hit& hit, bool anyHit);
	inline bool rayTrace(const Ray& ray, const std::shared_ptr<Scene> scene, size_t originMeshIndex, size_t originTriangleIndex)
	{ 
		Hit hit = Hit();
		return rayTrace (ray, scene, originMeshIndex, originTriangleIndex, hit, true);
	};
	glm::vec3 lightRadiance (const std::shared_ptr<LightSource> lightPtr, const glm::vec3 & position) const;
	glm::vec3 materialReflectance (const std::shared_ptr<Scene> scenePtr, 
								   const std::shared_ptr<Material> material, 
								   const glm::vec3& wi, 
								   const glm::vec3& wo, 
								   const glm::vec3& n,
								   const glm::vec2& textCoord) const;
	glm::vec3 shade(const std::shared_ptr<Scene> scenePtr, const Ray & ray, const Hit& hit);
	glm::vec3 sample (const std::shared_ptr<Scene> scenePtr, const Ray & ray, size_t originMeshIndex, size_t originTriangleIndex);

	std::shared_ptr<Image> m_imagePtr;
	
    std::unique_ptr<BVH> m_bvh;
};