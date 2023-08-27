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
	list = BVH::makeIndexPairSet(scene);
	m_bvh = std::make_unique<BVH>(scene, list);
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
			colorResponse += sample (scenePtr, ray, -1, -1, false);
			m_imagePtr->operator()(x, y) = colorResponse;
		}
	}
}

bool RayTracer::rayTrace(const Ray& ray, const std::shared_ptr<Scene> scene, size_t originMeshIndex, size_t originTriangleIndex, Hit& hit, bool insideHit) {
	float closest = std::numeric_limits<float>::max();
	bool intersectionFound = false;
	std::vector<std::pair<size_t,size_t>> candidateMeshTrianglePairs;
	candidateMeshTrianglePairs.reserve(16);
	std::vector<std::pair<size_t,size_t>> list; 
	for (size_t mIndex = 0; mIndex < scene->numOfMeshes (); mIndex++) {
		const auto& T = scene->mesh(mIndex)->triangleIndices();
		for (size_t tIndex = 0; tIndex < T.size(); tIndex++) 
			list.push_back (std::pair<size_t,size_t> (mIndex, tIndex));
	}
	m_bvh->intersect(ray, candidateMeshTrianglePairs);

	for (size_t i = 0; i < candidateMeshTrianglePairs.size(); i++) {
		size_t mIndex = candidateMeshTrianglePairs[i].first;
		size_t tIndex = candidateMeshTrianglePairs[i].second;
		if (!insideHit && mIndex == originMeshIndex && tIndex == originTriangleIndex)
			continue;
		if (insideHit && mIndex == originMeshIndex) // ray is currently inside an object, 
			continue;							    // approximation: we limit the number of hits with the same mesh by 2
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
				if (dt < closest) {
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
										  const glm::vec3& n,
										  const glm::vec2& textCoord) const {
		return BRDF (wi, wo, n, textCoord, materialPtr->albedo (), materialPtr->roughness (), materialPtr->metallicness ()); 
}

glm::vec3 RayTracer::shade(const std::shared_ptr<Scene> scenePtr, const Ray & ray, const Hit& hit) {
	const auto& mesh = scenePtr->mesh(hit.m_meshIndex);
	const std::shared_ptr<Material> materialPtr = scenePtr->material(scenePtr->mesh2material(hit.m_meshIndex));
	const auto& P = mesh->vertexPositions();
	const auto& N = mesh->vertexNormals();
	const auto& T = mesh->vertexTextCoords();
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
	glm::vec2 textCoord (0.f, 0.f);
	if (!T.empty()){
		textCoord = barycentricInterpolation(T[triangle[0]], T[triangle[1]], T[triangle[2]], w, hit.m_uCoord, hit.m_vCoord);
	}  
	for (size_t i = 0; i < scenePtr->numOfLightSources (); ++i) {
		const std::shared_ptr<LightSource> light = scenePtr->lightSource (i);
		glm::vec3 wi = normalize(-light->direction());
		float wiDotN = max(0.f, dot(wi, hitNormal));
		if (wiDotN <= 0.f)
			continue;
		colorResponse += lightRadiance (light, hitPosition) * materialReflectance(scenePtr, materialPtr, wi, wo, hitNormal, textCoord) * wiDotN;
	}
	return colorResponse;
}

std::pair<glm::vec3, glm::vec3> RayTracer::hitPN (const std::shared_ptr<Scene> scenePtr, const Hit& hit) {
	const auto& mesh = scenePtr->mesh(hit.m_meshIndex);
	const auto& P = mesh->vertexPositions();
	const auto& N = mesh->vertexNormals();
	glm::mat4 modelMatrix = mesh->computeTransformMatrix ();
	const glm::uvec3 & triangle = mesh->triangleIndices()[hit.m_triangleIndex];
	float w = 1.f - hit.m_uCoord - hit.m_vCoord;
	glm::vec3 hitPosition = barycentricInterpolation(P[triangle[0]], P[triangle[1]], P[triangle[2]], w, hit.m_uCoord, hit.m_vCoord);
	glm::vec3 unormalizedHitNormal = barycentricInterpolation(N[triangle[0]], N[triangle[1]], N[triangle[2]], w, hit.m_uCoord, hit.m_vCoord);
	glm::mat4 normalMatrix = glm::transpose (glm::inverse (modelMatrix));
	glm::vec3 hitNormal = normalize (glm::vec3 (normalMatrix * glm::vec4 (normalize (unormalizedHitNormal), 1.0)));
	return std::pair(hitPosition, hitNormal);
}

float RayTracer::transmissionCoeff (const glm::vec3& wo, const glm::vec3& n, float refraction) {
	float n_r = refraction;
	float k = 1 - sqr(n_r)*(1 - sqr(dot(n, wo)));
	return k; // if k <= 0 -> no light transmission
}

Ray RayTracer::refractRay (const glm::vec3& wo, const glm::vec3& n, const glm::vec3& p, float refraction, float transmissionK) {
	float n_r = refraction;
	float k = transmissionK;
	glm::vec3 direction = -n_r * wo - (n_r - dot(n, wo) + sqrt(k)) * n;
	glm::vec3 origin = p;
	Ray refracted_ray(origin, direction);
	return refracted_ray; 
}

// Preparing for Monte Carlo Path Tracing...
glm::vec3 RayTracer::sample (const std::shared_ptr<Scene> scenePtr, const Ray & ray, size_t originMeshIndex, size_t originTriangleIndex, bool insideHit = false) {
	Hit hit;
	bool intersectionFound = rayTrace(ray, scenePtr, originMeshIndex, originTriangleIndex, hit, insideHit);
	bool single_refraction = false;
	if (intersectionFound && hit.m_distance > 0.f) {
		std::pair<glm::vec3,glm::vec3> PN = hitPN(scenePtr, hit);
		glm::vec3 hitPosition = PN.first;
		glm::vec3 hitNormal = PN.second;
		glm::vec3 wo = normalize(-ray.direction ());
		const std::shared_ptr<Material> materialPtr = scenePtr->material(scenePtr->mesh2material(hit.m_meshIndex));
		float transparency = materialPtr->transparency();

		if (transparency != 0) { // transparent
			float refraction = materialPtr->refraction();
			float transmissionK = transmissionCoeff(wo, hitNormal, refraction);
			bool inside = hit.m_meshIndex == originMeshIndex ? true : false;
			if (transmissionK > 0) { // refraction
				Ray refracted_ray = refractRay(wo, hitNormal, hitPosition, refraction, transmissionK);
				if (inside) { // hit inside an object
					return sample(scenePtr, refracted_ray, hit.m_meshIndex, hit.m_triangleIndex, true);
				}
				else{ // mix reflected and transmitied light
					return glm::mix(shade(scenePtr, ray, hit), sample(scenePtr, refracted_ray, hit.m_meshIndex, hit.m_triangleIndex, single_refraction), transparency);
				}
			}
			else { // no refraction
				if (inside) { // hit inside an object
					return glm::vec3(0.0); // approximation: suppose that ray lost all its energy inside the object
				}
				else { // reflect environment
					float reflectance = materialPtr->reflectance();
					return glm::mix(shade(scenePtr, ray, hit), scenePtr->backgroundColor (), reflectance);	
				}
			}
		}
		else { // opaque 
			return shade(scenePtr, ray, hit);
		}
	} 
	else {
		return scenePtr->backgroundColor ();
	}
}