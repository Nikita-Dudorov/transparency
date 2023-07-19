// ----------------------------------------------
// Polytechnique - INF584 "Image Synthesis"
//
// Base code for practical assignments.
//
// Copyright (C) 2022 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------
#include "Ray.h"

#include <algorithm>
#include <limits>

using namespace std;

static float EPSILON = 0.00000001f;

bool Ray::triangleIntersect(
	const glm::vec3& p0,
	const glm::vec3& p1,
	const glm::vec3& p2,
	float& u,
	float& v,
	float& t) const {
	glm::vec3 e1 = p1 - p0, e2 = p2 - p0;
	glm::vec3 dxe2 = cross(m_direction, e2);
	float det = dot(e1, dxe2);
	if (fabs (det) < EPSILON)
		return false;
	float invDet = 1.f / det;
	glm::vec3 op0 = m_origin - p0;
	u = dot(op0, dxe2) * invDet;
	glm::vec3 op0xe1 = cross(op0, e1);
	v = dot(m_direction, op0xe1) * invDet;
	t = dot(e2, op0xe1) * invDet;
	if (u < 0.f || u > 1.f)
		return false;
	if (v >= 0.f && u + v <= 1.f)
		return true;
	return false;
}

bool Ray::boxIntersect(
	const glm::vec3& boxMin,
	const glm::vec3& boxMax,
	float& nearT,
	float& farT) const {
	nearT = std::numeric_limits<float>::min();
	farT = std::numeric_limits<float>::max();
	glm::vec3 dRcp(1.f / m_direction[0], 1.f / m_direction[1], 1.f / m_direction[2]);
	for (int i = 0; i < 3; i++) {
		const float direction = m_direction[i];
		const float origin = m_origin[i];
		const float minVal = boxMin[i], maxVal = boxMax[i];
		if (direction == 0) {
			if (origin < minVal || origin > maxVal)
				return false;
		}
		else {
			float t1 = (minVal - origin) * dRcp[i];
			float t2 = (maxVal - origin) * dRcp[i];
			if (t1 > t2)
				std::swap(t1, t2);
			nearT = std::max(t1, nearT);
			farT = std::min(t2, farT);
			if (!(nearT <= farT))
				return false;
		}
	}
	return true;
}
