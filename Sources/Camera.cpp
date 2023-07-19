// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#include "Camera.h"

Ray Camera::rayAt (float u, float v) const {
	glm::mat4 viewMat = inverse (computeViewMatrix());
	glm::vec3 viewRight = normalize (glm::vec3 (viewMat[0]));
	glm::vec3 viewUp = normalize (glm::vec3 (viewMat[1]));
	glm::vec3 viewDir = -normalize (glm::vec3 (viewMat[2]));
	glm::vec3 eye = glm::vec3 (viewMat[3]);
	float w = 2.0*float (tan (glm::radians (m_fov/2.0)));
	glm::vec3 rayDir = glm::normalize (viewDir + ((u - 0.5f) * m_aspectRatio * w) * viewRight + ((1.f-v) - 0.5f) * w * viewUp);  
	return Ray (eye, rayDir);
}