#pragma once

#include <limits>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

/// A 3D ray with basic intersection primitives for ray tracing.
class Ray {
public:
    inline Ray (const glm::vec3 & origin, const glm::vec3 & direction) : m_origin (origin), m_direction (direction) {}
    
    inline ~Ray () {}

    inline const glm::vec3 & origin () const { return m_origin; }
    
    inline const glm::vec3 & direction () const { return m_direction; }
    
    bool triangleIntersect (const glm::vec3 &p0,
                            const glm::vec3 &p1,
                            const glm::vec3 &p2,
                            float & u,
                            float & v,
                            float & t) const;
    
    bool boxIntersect (const glm::vec3 & boxMin,
                       const glm::vec3 & boxMax,
                       float & nearT,
                       float & farT) const;

private:
    glm::vec3 m_origin;
    glm::vec3 m_direction;
};

