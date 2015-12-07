//
//  svpoint.h
//  SphericalVoronoi
//
//  Created by Home on 2014-04-17.
//  Copyright (c) 2014 whenitsdone.org. All rights reserved.
//

#ifndef SphericalVoronoi_svpoint_h
#define SphericalVoronoi_svpoint_h

#include <cassert>
#include <ostream>
#include <tuple>
#include <bitset>

#ifndef GLM_FORCE_RADIANS
#define GLM_FORCE_RADIANS
#endif

#include "glm/glm.hpp"
#include "xs_Float.h"  // http://stereopsis.com/sree/fpu2006.html

#include "svBase.h"
#include "svColor.h"

#define FLOOR_TO_INT(val)      xs_FloorToInt(val)
#define ROUND_TO_INT(val)      xs_RoundToInt(val)
#define CEIL_TO_INT(val)       xs_CeilToInt(val)

namespace sv
{
    typedef double Real;
    typedef glm::dvec2 Real2;
    typedef glm::dvec3 Real3;
    typedef glm::dvec4 Real4;
    typedef glm::dmat4 Mat4;
    
    typedef glm::vec2 F2;
    typedef glm::vec3 F3;
    typedef glm::vec4 F4;
    
    inline std::string stringFromReal3(const Real3& r3)
    {
        std::stringstream ss;
        ss << "<" << r3.x << "," << r3.y << "," << r3.z << ">";
        return ss.str();
    }
    
//    typedef float Real;
//    typedef glm::vec3 Real3;
//    typedef glm::vec4 Real4;
//    typedef glm::mat4 Mat4;
//
//    typedef glm::dvec3 Double3;
//    typedef glm::dvec4 Double4;
//    typedef glm::dmat4 DMat4;
    
    class Ray;
    class AABB;
    class Plane;
    class Point;
    
    inline Real3 getTransformPosition(const Mat4& transform) {
        return Real3(transform[3][0], transform[3][1], transform[3][2]);
    }
    
    enum ECubeFace
    {
        CF_POSX = 0,
        CF_FIRST = 0,
        CF_NEGX,
        CF_POSY,
        CF_NEGY,
        CF_POSZ,
        CF_NEGZ,
        
        CF_MAX,
        
        CF_POSX_BITMASK = (1 << CF_POSX),
        CF_NEGX_BITMASK = (1 << CF_NEGX),
        CF_POSY_BITMASK = (1 << CF_POSY),
        CF_NEGY_BITMASK = (1 << CF_NEGY),
        CF_POSZ_BITMASK = (1 << CF_POSZ),
        CF_NEGZ_BITMASK = (1 << CF_NEGZ),
        
        CF_INVALID = 0xFF
    };
    
    typedef std::bitset<CF_MAX> CubeFaceBitSet;
    
    class Point
    {
    public:
        Point() : theta(0), phi(0) { computePosition(); }
        
        Point(Real theta_, Real phi_)
        : theta(theta_), phi(phi_)
        {
            computePosition();
        }
        
        Point(const Real3& direction)
        {
            assignDirection(direction);
        }
        
        Point(Real x, Real y, Real z)
        {
            assignDirection(Real3(x, y, z));
        }
        
        Point(const F3& pos)
        : Point(pos.x, pos.y, pos.z)
        {
        }
        
        Real theta;
        Real phi;
        Real3 position;
        
        std::tuple<Real3, CubeFaceBitSet> cubeCoord() const;
        
        Real3 tangent() const;
        Real3 binormal() const;
        
        void assignDirection(const Real3& direction)
        {
            Real r = glm::length(direction);
            assert(r > 0);
            theta = glm::acos(glm::clamp<Real>(direction.z / r, -1.0, 1.0));
            phi = glm::atan(direction.y, direction.x);
            position = direction / r;
        }
        
        Real sphericalDistance(const Point& p2) const
        {
            Real dot = glm::dot(position, p2.position);
            Real result = glm::acos(glm::clamp<Real>(dot, -1.0, 1.0));
            return result;
        }
        
        bool equalWithEps(const Point& p2, float eps) const
        {
            using namespace glm;
            return abs(position.x - p2.position.x) < eps &&
                abs(position.y - p2.position.y) < eps &&
                abs(position.z - p2.position.z) < eps;
        }
        
        bool operator < (const Point& p2)
        {
            return (theta < p2.theta) || (theta == p2.theta && phi < p2.phi);
        }
        
        friend std::ostream& operator << (std::ostream& stream, const Point& p)
        {
            return stream << p.theta << "," << p.phi;
        }
        
        template <class Archive>
        void serialize(Archive& ar)
        {
            ar(theta, phi, position.x, position.y, position.z);
        }
        
    private:
        void computePosition()
        {
            using namespace glm;
            position = Real3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
        }
    };
    
    namespace Util
    {
        extern std::vector<Real3> splitSphericalLineSegment(const Point& start, const Point& end, Real deltaAngle = M_PI / 180.0);
        extern Real lagrangeInterpolate(Real x, const std::vector<Real>& xArray, const std::vector<Real>& yArray);
        extern Real interpolateSphericalSamples(const Point& p0, const std::vector<Point>& points, const std::vector<Real>& values);
        extern Real computeTriangleArea(const Real3& p0, const Real3& p1, const Real3& p2);
        
        extern void faceAxisDirection(ECubeFace face, Real3& s_dir, Real3& t_dir, Real3& p_dir);

        inline Real sqrDistance(Real2 a, Real2 b)
        {
            auto p = a - b;
            return glm::dot(p, p);
        }
    }
    
    class SphericalLine
    {
    public:
        SphericalLine() : direction(Real3(0, 0, 1)), xi(0) {}
        SphericalLine(const Real3& direction_, Real xi_) : direction(glm::normalize(direction_)), xi(xi_) {}
        
        Real3 direction;
        Real xi;
    };
    
    class AABB
    {
    public:
        AABB();
        AABB(const Real3& p);
        AABB(const Real3& min, const Real3& max);
        
        void reset();
        bool isValid() const;
        bool isEmpty() const;
        
        const Real3& min() const { assert(isValid()); return m_min; }
        const Real3& max() const { assert(isValid()); return m_max; }
        
        Real3 center() const { assert(isValid()); return (m_min + m_max) * (Real)0.5; }
        Real3 size() const { assert(isValid()); return m_max - m_min; }
        Real3 extent() const { return size() * (Real)0.5; }
        
        void getMajorVertices(const Real3& direction, Real3& P, Real3& N) const;
        
        void unionWith(const Real3& p);
        void unionWith(const AABB& aabb);
        bool contains(const Real3& p) const;
        bool contains(const AABB& aabb) const;
        
        bool operator == (const AABB& aabb) const;
    private:
        Real3 m_min, m_max;
    };
    
    class Ray
    {
    public:
        Ray();
        Ray(const Real3& origin, const Real3& direction);
        
        inline const Real3& origin() const { return m_origin; }
        inline void setOrigin(const Real3& origin) { m_origin = origin; }
        inline const Real3& direction() const { return m_direction; }
        inline void setDirection(const Real3& direction) { m_direction = glm::normalize(direction); }
        inline void setNormalizedDirection(const Real3& direction) { m_direction = direction; }
        
    private:
        Real3 m_origin;
        Real3 m_direction;
    };
    
    class Plane
    {
    public:
        Plane();
        Plane(const Plane& other);
        Plane(const Real3& normal, Real distance);
        Plane(const Real4& vec);
        Plane(const Real3& a, const Real3& b, const Real3& c);
        
        const Real3& normal() const { return m_normal; }
        void setNormal(const Real3& normal) { m_normal = normal; }
        Real distance() const { return m_distance; }
        void setDistance(Real distance) { m_distance = distance; }
        
        Plane normalize() const;
        Plane transform(const glm::mat4 transform) const;
        Real distance(const Real3& point) const;
        bool pointOnSide(const Real3& point) const;
        bool lineIntersection(const Real3& ptA, const Real3& ptB, Real3& resultDestination) const;

    private:
        Real3 m_normal;
        Real m_distance;
    };
    
    bool threePlanesIntersection(const Plane& planeA, const Plane& planeB, const Plane& planeC, Real3& result);
    
    bool rayAabbIntersection(const Ray& ray, const AABB& aabb);
    
    template <typename T>
    class PositionT
    {
    public:
        using Vec3 = glm::detail::tvec3<T>;
        PositionT();
        PositionT(ECubeFace face, T s, T t, T p);
        PositionT(ECubeFace face, const Vec3& stp);
        
        ECubeFace face() const { return m_face; }
        const Vec3& surfacePoint() const { return m_surfacePoint; }
        Vec3 stpCoords() const;
        const Vec3& spacePosition() const { return m_spacePosition; }
        
    private:
        ECubeFace m_face;
        T m_height;
        Vec3 m_surfacePoint;
        Vec3 m_spacePosition;
    };
    
    typedef PositionT<Real> Position;
    typedef PositionT<float> PositionF;
    typedef PositionT<double> PositionD;
}

#include "svMath.hpp"

#endif
