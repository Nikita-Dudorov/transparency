#version 410 core // Minimal GL version support expected from the GPU

layout(location=0) in vec3 vPosition; // The 1st input attribute is the position (CPU side: glVertexAttrib 0)
layout(location=1) in vec3 vNormal;
layout(location=2) in vec2 vTextCoord;

uniform mat4 modelMat, viewMat, projectionMat, normalMat;

out vec3 fPosition;
out vec3 fNormal;
out vec2 fTextCoord;

void main() {
    if (vNormal == vec3(0.0)) { // screen texture
        gl_Position =  vec4 (vPosition, 1.0);
        return;
    }

    vec4 p = viewMat * modelMat * vec4 (vPosition, 1.0);
    gl_Position =  projectionMat * p; 
    vec4 n = normalMat * vec4 (vNormal, 1.0);
    fPosition = p.xyz;
    fNormal = normalize (n.xyz);
    fTextCoord = vTextCoord;
}