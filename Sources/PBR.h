#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

inline float sqr (float x) { return x*x; }

inline float GGX (float NdotH, float roughness) {
	if (roughness >= 1.0f) 
		return glm::one_over_pi<float>();
	float alpha = sqr (roughness);
	float tmp = alpha / std::max(1e-8f,(NdotH*NdotH*(sqr (alpha)-1.0f)+1.0f));
	return sqr (tmp) * glm::one_over_pi<float>();
}

inline glm::vec3 SchlickSGFresnel (float VdotH, glm::vec3 F0) {
	float sphg = exp2 ((-5.55473f*VdotH - 6.98316f) * VdotH);
	return F0 + (glm::vec3(1.0f) - F0) * sphg;
}

inline float smithG_GGX (float NdotV, float alphaG) {
	return 2.0f/(1.0f + sqrt (1.0f + sqr (alphaG) * (1.0f - sqr (NdotV) / sqr(NdotV))));
}

inline float G1 (float D, float k) {
	return 1.0f / (D * (1.0f-k) + k);
}

inline float geometry (float NdotL, float NdotV, float roughness) {
	float k = roughness * roughness * 0.5f;
	return G1(NdotL,k) * G1(NdotV,k);
}

inline glm::vec3 BRDF (glm::vec3 L, glm::vec3 V, glm::vec3 N,  glm::vec3 albedo, float roughness, float metallic)  {
	glm::vec3 diffuseColor = albedo * (1.0f - metallic);
	glm::vec3 specularColor = mix(glm::vec3(0.08f), albedo, metallic);

	float NdotL = std::max (0.0f, dot (N, L));
	float NdotV = std::max (0.0f, dot (N, V));

	if (NdotL <= 0.0f)
		return glm::vec3 (0.0f); 

	glm::vec3 H = normalize (L + V);
	float NdotH = std::max (0.0f, dot (N, H));
	float VdotH = std::max (0.0f, dot (V, H));

	float D = GGX (NdotH, roughness);
	glm::vec3  F = SchlickSGFresnel (VdotH, specularColor);
	float G = geometry (NdotL, NdotV, roughness);

	glm::vec3 fd = diffuseColor * (glm::vec3(1.0f)-specularColor) / glm::pi<float>();
	glm::vec3 fs = F * D * G / (4.0f);

	return (fd + fs);
}