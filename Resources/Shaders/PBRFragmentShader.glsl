#version 410 core // Minimal GL version support expected from the GPU

const float PI = 3.14159265358979323846;
const float INV_PI = 0.31830988618379067153776752674503;

const float IOR = 1.31;

struct LightSource {
	vec3 direction;
	vec3 color;
	float intensity;
};

struct Material {
	vec3 albedo;
	float roughness;
	float metallicness;
};

const int MAX_NUM_OF_LIGHT_SOURCES = 8;

// Shader uniform variables, set from the CPU side code
uniform mat4 transInvViewMat;
uniform int numOfLightSources;
uniform LightSource lightSources[MAX_NUM_OF_LIGHT_SOURCES]; 
uniform Material material;

// Shader inputs, linearly interpolated by default from the previous stage outputs (here the vertex shader)
in vec3 fPosition; 
in vec3 fNormal;

out vec4 colorResponse; // Shader output: the color response attached to this fragment

float sqr (float x) { return x*x; }

float GGX (float NdotH, float roughness) {
	if (roughness >= 1.0) 
		return INV_PI;
	float alpha = sqr (roughness);
	float tmp = alpha / max(1e-8,(NdotH*NdotH*(sqr (alpha)-1.0)+1.0));
	return tmp * tmp * INV_PI;
}

vec3 SchlickSGFresnel (float VdotH, vec3 F0) {
	float sphg = exp2 ((-5.55473*VdotH - 6.98316) * VdotH);
	return F0 + (vec3(1.0) - F0) * sphg;
}

float smithG_GGX (float NdotV, float alphaG) {
    return 2.0/(1.0 + sqrt (1.0 + sqr (alphaG) * (1.0 - sqr (NdotV) / sqr(NdotV))));
}

float G1 (float D, float k) {
	return 1.0 / (D * (1.0-k) + k);
}

float geometry (float NdotL, float NdotV,	float roughness) {
	float k = roughness * roughness * 0.5;
	return G1(NdotL,k) * G1(NdotV,k);
}

vec3 BRDF (vec3 L, vec3 V, vec3 N, vec3 albedo, float roughness, float metallic) {
    vec3 diffuseColor = albedo * (1.0 - metallic);
	vec3 specularColor = mix(vec3(0.08), albedo, metallic);

    float NdotL = max (0.0, dot (N, L));
    float NdotV = max (0.0, dot (N, V));

    if (NdotL <= 0.0)
   	 	return vec3 (0.0); 
    
    vec3 H = normalize (L + V);
    float NdotH = max (0.0, dot (N, H));
    float VdotH = max (0.0, dot (V, H));
    
    float D = GGX (NdotH, roughness);
    vec3  F = SchlickSGFresnel (VdotH, specularColor);
    float G = geometry (NdotL, NdotV, roughness);

    vec3 fd = diffuseColor * (vec3(1.0)-specularColor) / PI;
    vec3 fs = F * D * G / (4.0);
  
    return (fd + fs);
}

vec3 computeReflectedRadiance (LightSource l, vec3 wo, vec3 n, vec3 albedo, float roughness, float metallic) {
	vec3 Li = l.color * l.intensity * PI;
	// Transform the world space light direction to view space
	vec3 wi = normalize ((transInvViewMat * vec4(-normalize (l.direction), 1.0)).xyz);
	vec3 fr = BRDF (wi, wo, n, albedo, roughness, metallic);
	return Li * fr * max (0.0, dot (wi, n));
}

vec3 pbrShading (vec3 wo, vec3 n, vec3 albedo, float roughness, float metallic) {
	vec3 radiance = vec3 (0.0);
	for (int i = 0; i < numOfLightSources; i++)
		radiance += computeReflectedRadiance (lightSources[i], wo, n, albedo, roughness, metallic);
   	return radiance; 
}

void main () {
	vec3 n = normalize (fNormal);// Linear barycentric interpolation does not preserve unit vectors, normal texture filtering
	vec3 wo = normalize (-fPosition);
	colorResponse = vec4 (pbrShading (wo, n, material.albedo, material.roughness, material.metallicness), 1.0);
}