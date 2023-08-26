#version 410 core // Minimal GL version support expected from the GPU

const float PI = 3.14159265358979323846;
const float INV_PI = 0.31830988618379067153776752674503;

const float DZ = 1.0; // approximate distance to opaque background

struct LightSource {
	vec3 direction;
	vec3 color;
	float intensity;
};

struct Material {
	vec3 albedo;
	float roughness;
	float metallicness;
	float transparency;
	float refraction;
};

const int MAX_NUM_OF_LIGHT_SOURCES = 8;

// Shader uniform variables, set from the CPU side code
uniform sampler2D imageTex;
uniform int width;
uniform int height;

uniform mat4 transInvViewMat, projectionMat;
uniform int numOfLightSources;
uniform LightSource lightSources[MAX_NUM_OF_LIGHT_SOURCES]; 
uniform Material material;

// Shader inputs, linearly interpolated by default from the previous stage outputs (here the vertex shader)
in vec3 fPosition; 
in vec3 fNormal;
in vec2 fTextCoord;

out vec4 colorResponse; // Shader output: the color response attached to this fragment

float sqr (float x) { return x*x; }

vec2 WorldToScreen (vec3 p) {
	vec4 p_projected = projectionMat * vec4 (p, 1.0);
	vec3 p_screen = p_projected.xyz / p_projected.w; // in [-1; 1] x [-1; 1]
	vec2 screenCoords = vec2 (0.5*(p_screen.x + 1), 0.5*(-p_screen.y + 1)); // in [0;1] x [0;1]
	return screenCoords;
}

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

	// optional grid texture
	if (fTextCoord != vec2(0.0)){
		float alpha = mod(floor (10 * fTextCoord.x) + floor (10 * fTextCoord.y), 2);
		fd = alpha * fd;
	}
  
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

float transmissionCoeff (vec3 wo, vec3 n, float refraction) {
	float n_r = refraction;
	float k = 1 - sqr(n_r)*(1 - sqr(dot(n, wo)));
	return k; // if k <= 0 -> no light transmission
}

vec3 SnellRefractedRay (vec3 wo, vec3 n, float refraction, float kTransmission) {
	float n_r = refraction;
	float k = kTransmission;
	vec3 refractedRay = -n_r * wo - (n_r - dot(n, wo) + sqrt(k)) * n;
	return refractedRay; 
}

vec2 RefractedRayScreen (vec3 wo, vec3 n, float refraction, float kTransmission) {
	vec3 ray_refracted = SnellRefractedRay(wo, n, refraction, kTransmission);
	vec3 ray_camera = vec3 (0.0, 0.0, -1.0);
	vec3 dest = fPosition + ray_refracted * DZ / dot(ray_refracted, ray_camera);
	vec2 dest_screen = WorldToScreen(dest);
	return dest_screen;
}

void main () {
	if (material.transparency == 0.0) { // Opaque background
		//colorResponse = vec4 (texture(imageTex, vec2 (gl_FragCoord.x/width, gl_FragCoord.y/height)).rgb, 1.0);
		colorResponse = vec4 (texture(imageTex, fTextCoord).rgb, 1.0);
		return;	
	}

	vec3 n = normalize (fNormal); // Linear barycentric interpolation does not preserve unit vectors, normal texture filtering
	vec3 wo = normalize (-fPosition);
	vec3 pbrColorResponse = pbrShading (wo, n, material.albedo, material.roughness, material.metallicness);

	float kTransmission = transmissionCoeff (wo, n, material.refraction);
	if (kTransmission > 0) { // there is a refracted ray
		vec2 offset = RefractedRayScreen(wo, n, material.refraction, kTransmission);
		vec3 transmissionColorResponse = texture(imageTex, offset).rgb; 
		colorResponse = vec4 (mix (pbrColorResponse, transmissionColorResponse, material.transparency), 1.0);
	}
	else { // no refraction
		colorResponse = vec4 (pbrColorResponse, 1.0);
	}
}