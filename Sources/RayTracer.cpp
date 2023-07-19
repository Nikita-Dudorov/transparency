// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#include "RayTracer.h"

#include<chrono>

#include "Console.h"
#include "Camera.h"
#include "PBR.h"

RayTracer::RayTracer() :
	Renderer (), 
	m_imagePtr (std::make_shared<Image>()) {}

RayTracer::~RayTracer() {}

void RayTracer::init (const std::shared_ptr<Scene> scenePtr) {
	buildBVH(scenePtr);
}

void RayTracer::buildBVH(const std::shared_ptr<Scene> scene) {

	if (m_bvh != nullptr)
		m_bvh.reset();

	std::vector<std::pair<size_t,size_t>> list; 
	for (size_t mIndex = 0; mIndex < scene->numOfMeshes (); mIndex++) {
		const auto& T = scene->mesh(mIndex)->triangleIndices();
		for (size_t tIndex = 0; tIndex < T.size(); tIndex++) 
			list.push_back (std::pair<size_t,size_t> (mIndex, tIndex));
	}
	
	m_bvh = std::make_unique<BVH>(scene,list);
}

void RayTracer::render (const std::shared_ptr<Scene> scenePtr) {
	size_t width = m_imagePtr->width();
	size_t height = m_imagePtr->height();
	m_imagePtr->clear (scenePtr->backgroundColor ());
	const auto cameraPtr = scenePtr->camera();
	for (int y = 0; y < height; y++) {
#pragma omp parallel for
		for (int x = 0; x < width; x++) {
			glm::vec3 colorResponse (0.f, 0.f, 0.f);
			Ray ray = cameraPtr->rayAt((float(x) + 0.5) / width, 1.f - (float(y) + 0.5) / height);
			colorResponse += sample (scenePtr, ray, 0, 0);
			m_imagePtr->operator()(x, y) = colorResponse;
		}
	}
}

bool RayTracer::rayTrace(const Ray& ray, const std::shared_ptr<Scene> scene, size_t originMeshIndex, size_t originTriangleIndex, Hit& hit, bool anyHit) {
	float closest = std::numeric_limits<float>::max();
	bool intersectionFound = false;
	std::vector<std::pair<size_t,size_t>> candidateMeshTrianglePairs;
	//candidateMeshTrianglePairs.reserve(16);
	std::vector<std::pair<size_t,size_t>> list; 
	for (size_t mIndex = 0; mIndex < scene->numOfMeshes (); mIndex++) {
		const auto& T = scene->mesh(mIndex)->triangleIndices();
		for (size_t tIndex = 0; tIndex < T.size(); tIndex++) 
			list.push_back (std::pair<size_t,size_t> (mIndex, tIndex));
	}
	m_bvh->intersect(ray, candidateMeshTrianglePairs);
	//if(candidateMeshTrianglePairs.size()>0)
	//	std::cout << "*****" << candidateMeshTrianglePairs.size() << "*********" << std::endl;
	// for (size_t mIndex = 0; mIndex < scene->numOfMeshes (); mIndex++) {
	// 	const auto& T = scene->mesh(mIndex)->triangleIndices();
	// 	for (size_t tIndex = 0; tIndex < T.size(); tIndex++) 
	// 		candidateMeshTrianglePairs.push_back (std::pair<size_t,size_t> (mIndex, tIndex));
	// }
	for (size_t i = 0; i < candidateMeshTrianglePairs.size(); i++) {
		size_t mIndex = candidateMeshTrianglePairs[i].first;
		size_t tIndex = candidateMeshTrianglePairs[i].second;
		if (anyHit && mIndex == originMeshIndex && tIndex == originTriangleIndex)
				continue;
		const auto& mesh = scene->mesh(mIndex);
		const auto& triangleIndices = mesh->triangleIndices();
		const glm::uvec3& triangle = triangleIndices[tIndex];
		const auto& P = mesh->vertexPositions();
		glm::mat4 modelMatrix = mesh->computeTransformMatrix ();
		float ut, vt, dt;
		if (ray.triangleIntersect(glm::vec3 (modelMatrix * glm::vec4 (P[triangle[0]], 1.0)), 
								  glm::vec3 (modelMatrix * glm::vec4 (P[triangle[1]], 1.0)), 
								  glm::vec3 (modelMatrix * glm::vec4 (P[triangle[2]], 1.0)), 
								  ut, vt, dt) == true) {
			if (dt > 0.f) {
				if (anyHit)
					return true;
				else if (dt < closest) {
					intersectionFound = true;
					closest = dt;
					hit.m_meshIndex = mIndex;
					hit.m_triangleIndex = tIndex;
					hit.m_uCoord = ut;
					hit.m_vCoord = vt;
					hit.m_distance = dt;
				}
			}
		}
	}
	return intersectionFound;
}

glm::vec3 RayTracer::lightRadiance (const std::shared_ptr<LightSource> lightPtr, const glm::vec3 & position) const {
	return lightPtr->color() * lightPtr->intensity() * glm::pi<float>(); 
}

glm::vec3 RayTracer::materialReflectance (const std::shared_ptr<Scene> scenePtr,
										  const std::shared_ptr<Material> materialPtr, 
										  const glm::vec3& wi, 
										  const glm::vec3& wo, 
										  const glm::vec3& n) const {
		return BRDF (wi, wo, n, materialPtr->albedo (), materialPtr->roughness (), materialPtr->metallicness ()); 
}

glm::vec3 RayTracer::shade(const std::shared_ptr<Scene> scenePtr, const Ray & ray, const Hit& hit) {
	const auto& mesh = scenePtr->mesh(hit.m_meshIndex);
	const std::shared_ptr<Material> materialPtr = scenePtr->material(scenePtr->mesh2material(hit.m_meshIndex));
	const auto& P = mesh->vertexPositions();
	const auto& N = mesh->vertexNormals();
	glm::mat4 modelMatrix = mesh->computeTransformMatrix ();
	const glm::uvec3 & triangle = mesh->triangleIndices()[hit.m_triangleIndex];
	float w = 1.f - hit.m_uCoord - hit.m_vCoord;
	glm::vec3 hitPosition = barycentricInterpolation(P[triangle[0]], P[triangle[1]], P[triangle[2]], w, hit.m_uCoord, hit.m_vCoord);
	hitPosition = glm::vec3 (modelMatrix * glm::vec4 (hitPosition, 1.0));
	glm::vec3 unormalizedHitNormal = barycentricInterpolation(N[triangle[0]], N[triangle[1]], N[triangle[2]], w, hit.m_uCoord, hit.m_vCoord);
	glm::mat4 normalMatrix = glm::transpose (glm::inverse (modelMatrix));
	glm::vec3 hitNormal = normalize (glm::vec3 (normalMatrix * glm::vec4 (normalize (unormalizedHitNormal), 1.0)));
	glm::vec3 wo = normalize(-ray.direction ());
	glm::vec3 colorResponse (0.f, 0.f, 0.f);
	for (size_t i = 0; i < scenePtr->numOfLightSources (); ++i) {
		const std::shared_ptr<LightSource> light = scenePtr->lightSource (i);
		glm::vec3 wi = normalize(-light->direction());
		float wiDotN = max(0.f, dot(wi, hitNormal));
		if (wiDotN <= 0.f)
			continue;
		colorResponse += lightRadiance (light, hitPosition) * materialReflectance(scenePtr, materialPtr, wi, wo, hitNormal) * wiDotN;
	}
	return colorResponse;
}

// Preparing for Monte Carlo Path Tracing...
glm::vec3 RayTracer::sample (const std::shared_ptr<Scene> scenePtr, const Ray & ray, size_t originMeshIndex, size_t originTriangleIndex) {
	Hit hit;
	bool intersectionFound = rayTrace(ray, scenePtr, originMeshIndex, originTriangleIndex, hit, false);
	if (intersectionFound && hit.m_distance > 0.f) {
		return shade(scenePtr, ray, hit);
	} else 
		return scenePtr->backgroundColor (); 
}