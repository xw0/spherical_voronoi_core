//
//  svMath.hpp
//  SphericalVoronoi
//
//  Created by Home on 2014-06-27.
//  Copyright (c) 2014 whenitsdone.org. All rights reserved.
//

#ifndef SphericalVoronoi_svMath_hpp
#define SphericalVoronoi_svMath_hpp

// cereal serialization

template <class Archive>
inline void save(Archive& archive, glm::vec3 const& v)
{
    archive(v.x, v.y, v.z);
}

template <class Archive>
inline void load(Archive& archive, glm::vec3& v)
{
    archive(v.x, v.y, v.z);
}

template <class Archive>
inline void save(Archive& archive, glm::vec2 const& v)
{
    archive(v.x, v.y);
}

template <class Archive>
inline void load(Archive& archive, glm::vec2& v)
{
    archive(v.x, v.y);
}

namespace sv
{

    inline AABB::AABB()
    {
        reset();
    }

    inline AABB::AABB(const Real3& p)
    : m_min(p), m_max(p)
    {
    }

    inline AABB::AABB(const Real3& min, const Real3& max)
    : m_min(min), m_max(max)
    {
        assert(m_min.x < m_max.x);
        assert(m_min.y < m_max.y);
        assert(m_min.z < m_max.z);
    }

    inline void AABB::reset()
    {
        m_min = Real3(FLT_MAX);
        m_max = Real3(-FLT_MAX);
    }

    inline bool AABB::isValid() const
    {
        bool result = m_min.x <= m_max.x;
        
        if (result)
        {
            assert(m_min.y <= m_max.y);
            assert(m_min.z <= m_max.z);
        }
        
        return result;
    }

    inline bool AABB::isEmpty() const
    {
        return m_min == m_max;
    }

    inline void AABB::unionWith(const Real3& p)
    {
        if (isValid())
        {
            m_min = glm::min(m_min, p);
            m_max = glm::max(m_max, p);
        }
        else
        {
            m_min = m_max = p;
        }
    }

    inline void AABB::unionWith(const AABB& aabb)
    {
        if (isValid())
        {
            m_min = glm::min(m_min, aabb.m_min);
            m_max = glm::max(m_max, aabb.m_max);
        }
        else
        {
            *this = aabb;
        }
    }

    inline bool AABB::contains(const Real3& p) const
    {
        if (!isValid()) return false;
        return
            (m_min.x <= p.x && p.x <= m_max.x) &&
            (m_min.y <= p.y && p.y <= m_max.y) &&
            (m_min.z <= p.z && p.z <= m_max.z);
    }
    
    inline bool AABB::contains(const AABB& aabb) const
    {
        if (!isValid() || !aabb.isValid()) return false;
        return
        (m_min.x <= aabb.m_min.x && aabb.m_max.x <= m_max.x) &&
        (m_min.y <= aabb.m_min.y && aabb.m_max.y <= m_max.y) &&
        (m_min.z <= aabb.m_min.z && aabb.m_max.z <= m_max.z);
    }
    
    inline void AABB::getMajorVertices(const Real3& direction, Real3& P, Real3& N) const
    {
        if (direction.x >= 0)
        {
            P.x = m_max.x;
            N.x = m_min.x;
        }
        else
        {
            P.x = m_min.x;
            N.x = m_max.x;
        }
        
        if (direction.y >= 0)
        {
            P.y = m_max.y;
            N.y = m_min.y;
        }
        else
        {
            P.y = m_min.y;
            N.y = m_max.y;
        }
        
        if (direction.z >= 0)
        {
            P.z = m_max.z;
            N.z = m_min.z;
        }
        else
        {
            P.z = m_min.z;
            N.z = m_max.z;
        }
    }
    
    inline bool AABB::operator == (const AABB& aabb) const
    {
        if (!isValid() || !aabb.isValid()) return false;
        return (m_min == aabb.m_min && m_max == aabb.m_max);
    }
    
    inline Ray::Ray()
    : Ray(Real3(0.0), Real3(1.0, 0, 0))
    {
    }
    
    inline Ray::Ray(const Real3& origin, const Real3& direction)
    {
        setOrigin(origin);
        setDirection(direction);
    }

    inline Plane::Plane() {}
    
    inline Plane::Plane(const Plane& other)
    : m_normal(other.m_normal),
    m_distance(other.m_distance)
    {}
    
    inline Plane::Plane(const Real3& normal, Real distance)
    : m_normal(normal),
    m_distance(distance)
    {}
    
    inline Plane::Plane(const Real4& vec)
    : m_normal(vec),
    m_distance(vec.w)
    {}
    
    inline Plane::Plane(const Real3& a, const Real3& b, const Real3& c)
    {
        m_normal = glm::normalize(glm::cross(b - a, c - a));
        m_distance = - glm::dot(a, m_normal);
    }

    
    inline Plane Plane::normalize() const {
        Plane res;
        
        Real denom = 1 / glm::length(m_normal);
        
        res.m_normal = m_normal * denom;
        res.m_distance = m_distance * denom;
        
        return res;
    }
    
    // The transform passed in should be the inverse transpose or be an orthogonal matrix.
    // http://stackoverflow.com/questions/7685495/transforming-a-3d-plane-by-4x4-matrix
    inline Plane Plane::transform(const glm::mat4 transform) const {
        Plane res;
        
        res.m_normal = Real3(Mat4(transform) * Real4(m_normal, 0));
        res.m_distance = m_distance - glm::dot(getTransformPosition(Mat4(transform)), res.m_normal);
        
        return res;
    }

    inline Real Plane::distance(const Real3& point) const {
        return glm::dot(m_normal, point) + m_distance;
    }
    
    inline bool Plane::pointOnSide(const Real3& point) const {
        return distance(point) >= 0;
    }
    
    // http://paulbourke.net/geometry/planeline/
    inline bool Plane::lineIntersection(const Real3& ptA, const Real3& ptB, Real3& resultDestination) const {
        Real3 vec = m_normal * (ptA - ptB);
        Real denom = vec.x + vec.y + vec.z;
        
        if(denom == 0) {
            //line is perpendicular to plane
            return false;
        }
        
        vec = m_normal * ptA;
        Real dist = (vec.x + vec.y + vec.z + m_distance) / denom;
        
        if(dist >= (Real)0 && dist <= (Real)1) {
            resultDestination = ptA + (ptB - ptA) * dist;
            return true;
        }
        else {
            return false;
        }
    }
    
    template <typename T>
    inline PositionT<T>::PositionT()
    : m_face(CF_INVALID), m_height(1.0f), m_surfacePoint(0), m_spacePosition(0)
    {
    }
    
    template <typename T>
    inline PositionT<T>::PositionT(ECubeFace face, T s, T t, T p)
    : m_face(face), m_height(p)
    {
        switch (face)
        {
            case CF_POSX:
                m_surfacePoint = Vec3(1.0, t, -s);
                break;
            case CF_NEGX:
                m_surfacePoint = Vec3(-1.0, t, s);
                break;
            case CF_POSY:
                m_surfacePoint = Vec3(t, 1.0, s);
                break;
            case CF_NEGY:
                m_surfacePoint = Vec3(-t, -1.0, s);
                break;
            case CF_POSZ:
                m_surfacePoint = Vec3(s, t, 1.0);
                break;
            case CF_NEGZ:
                m_surfacePoint = Vec3(-s, t, -1.0);
                break;
            default:
                assert(false);
        }
        m_spacePosition = glm::normalize(m_surfacePoint) * m_height;
    }
    
    template <typename T>
    inline PositionT<T>::PositionT(ECubeFace face, const typename PositionT<T>::Vec3& stp)
    : Position(face, stp.s, stp.t, stp.p)
    {
    }
    
    template <typename T>
    inline typename PositionT<T>::Vec3 PositionT<T>::stpCoords() const
    {
        switch (m_face)
        {
            case CF_POSX:
                return Vec3(-m_surfacePoint.z, m_surfacePoint.y, m_height);
            case CF_NEGX:
                return Vec3(m_surfacePoint.z, m_surfacePoint.y, m_height);
            case CF_POSY:
                return Vec3(m_surfacePoint.z, m_surfacePoint.x, m_height);
            case CF_NEGY:
                return Vec3(m_surfacePoint.z, -m_surfacePoint.x, m_height);
            case CF_POSZ:
                return Vec3(m_surfacePoint.x, m_surfacePoint.y, m_height);
            case CF_NEGZ:
                return Vec3(-m_surfacePoint.x, m_surfacePoint.y, m_height);
            default:
                assert(false);
                return Vec3(0, 0, 0);
        }
    }
}

#endif
