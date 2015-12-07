//
//  svMath.cpp
//  SphericalVoronoi
//
//  Created by Home on 2014-06-28.
//  Copyright (c) 2014 whenitsdone.org. All rights reserved.
//

#include "svPrefix.h"

#include "svMath.h"

namespace sv
{
    std::tuple<Real3, CubeFaceBitSet> Point::cubeCoord() const
    {
        Real3 absPos = glm::abs(position);
        Real3 result;
        CubeFaceBitSet faceSet;
        Real ratio;
        if (absPos.z >= glm::max(absPos.x, absPos.y))
        {
            ratio = 1.0 / absPos.z;
            faceSet |= (position.z > 0 ? CF_POSZ_BITMASK : CF_NEGZ_BITMASK);
            if (absPos.z == absPos.x)
            {
                faceSet |= (position.x > 0 ? CF_POSX_BITMASK : CF_NEGX_BITMASK);
            }
            if (absPos.z == absPos.y)
            {
                faceSet |= (position.y > 0 ? CF_POSY_BITMASK : CF_NEGY_BITMASK);
            }
        }
        else if (absPos.y >= glm::max(absPos.x, absPos.z))
        {
            ratio = 1.0 / absPos.y;
            faceSet |= (position.y > 0 ? CF_POSY_BITMASK : CF_NEGY_BITMASK);
            if (absPos.y == absPos.x)
            {
                faceSet |= (position.x > 0 ? CF_POSX_BITMASK : CF_NEGX_BITMASK);
            }
            if (absPos.y == absPos.z)
            {
                faceSet |= (position.z > 0 ? CF_POSZ_BITMASK : CF_NEGZ_BITMASK);
            }
        }
        else if (absPos.x >= glm::max(absPos.y, absPos.z))
        {
            ratio = 1.0 / absPos.x;
            faceSet |= (position.x > 0 ? CF_POSX_BITMASK : CF_NEGX_BITMASK);
            if (absPos.x == absPos.y)
            {
                faceSet |= (position.y > 0 ? CF_POSY_BITMASK : CF_NEGY_BITMASK);
            }
            if (absPos.x == absPos.z)
            {
                faceSet |= (position.z > 0 ? CF_POSZ_BITMASK : CF_NEGZ_BITMASK);
            }
        }
        result = position * ratio;
        return std::tie(result, faceSet);
    }
    
    Real3 Point::tangent() const
    {
        CubeFaceBitSet bitSet;
        
        Real3 binormal(1, 0, 0);
        
        if (bitSet.test(CF_POSX))
        {
            binormal = Real3(0, 1, 0);
        }
        else if (bitSet.test(CF_NEGX))
        {
            binormal = Real3(0, 1, 0);
        }
        else if (bitSet.test(CF_POSY))
        {
            binormal = Real3(1, 0, 0);
        }
        else if (bitSet.test(CF_NEGY))
        {
            binormal = Real3(-1, 0, 0);
        }
        else if (bitSet.test(CF_POSZ))
        {
            binormal = Real3(0, 1, 0);
        }
        else if (bitSet.test(CF_NEGZ))
        {
            binormal = Real3(0, 1, 0);
        }
        else
        {
            assert(false);
        }
        
        Real3 normal = position;
        Real3 result = glm::normalize(glm::cross(binormal, normal));
        return result;
    }
    
    Real3 Point::binormal() const
    {
        Real3 t = tangent();
        Real3 normal = position;
        Real3 result = glm::normalize(glm::cross(normal, t));
        return result;
    }
    
    // http://www.cgafaq.info/wiki/Intersection_of_three_planes
    bool threePlanesIntersection(const Plane& planeA, const Plane& planeB, const Plane& planeC, Real3& result)
    {
        Real3 bcCross = glm::cross(planeB.normal(), planeC.normal());
        Real denom = glm::dot(planeA.normal(), bcCross);
        
        if (denom == 0) {
            result = Real3(0);
            return false;
        }
        else {
            result = (-planeA.distance() * bcCross
                    - planeB.distance() * glm::cross(planeC.normal(), planeA.normal())
                    - planeC.distance() * glm::cross(planeA.normal(), planeB.normal())) / denom;
            return true;
        }
    }
    
    // http://tavianator.com/2011/05/fast-branchless-raybounding-box-intersections/
    bool rayAabbIntersection(const Ray& ray, const AABB& aabb)
    {
        Real3 n_inv = Real3(1.0) / ray.direction();
        
        double tx1 = (aabb.min().x - ray.origin().x)*n_inv.x;
        double tx2 = (aabb.max().x - ray.origin().x)*n_inv.x;
        
        double tmin = glm::min(tx1, tx2);
        double tmax = glm::max(tx1, tx2);
        
        double ty1 = (aabb.min().y - ray.origin().y)*n_inv.y;
        double ty2 = (aabb.max().y - ray.origin().y)*n_inv.y;
        
        tmin = glm::max(tmin, glm::min(ty1, ty2));
        tmax = glm::min(tmax, glm::max(ty1, ty2));
        
        double tz1 = (aabb.min().z - ray.origin().z)*n_inv.z;
        double tz2 = (aabb.max().z - ray.origin().z)*n_inv.z;
        
        tmin = glm::max(tmin, glm::min(tz1, tz2));
        tmax = glm::min(tmax, glm::max(tz1, tz2));
        
        return tmax >= glm::max(tmin, 0.0);
    }
    
    namespace Util
    {
    
        std::vector<Real3> splitSphericalLineSegment(const Point& start, const Point& end, Real deltaAngle)
        {
            std::vector<Real3> result;
            
            assert(start.position != -end.position);
            
            auto direction = glm::normalize(glm::cross(start.position, end.position));
            float distance = glm::acos(glm::dot(start.position, end.position));
            
            result.push_back(start.position);
            
            for (auto angle=deltaAngle; angle<distance; angle+=deltaAngle)
            {
                Mat4 rotation = glm::rotate(Mat4(1.0), angle, direction);
                Real3 pos = glm::normalize(Real3(rotation * Real4(start.position, 1.0)));
                
                result.push_back(pos);
            }
            
            result.push_back(end.position);
            
            return result;
        }
        
        Real lagrangeInterpolate(Real x, const std::vector<Real>& xArray, const std::vector<Real>& yArray)
        {
            assert(xArray.size() == yArray.size());
            
            Real sum = 0.0;
            for (unsigned int i = 0; i < xArray.size(); ++i)
            {
                Real Xi, Yi;
                Xi = xArray[i];
                Yi = yArray[i];
                Real factor = 1.0;
                for (unsigned int j = 0; j < xArray.size(); ++j)
                {
                    if (i != j)
                    {
                        Real Xj = xArray[j];
                        factor *= (x - Xj) / (Xi - Xj);
                    }
                }
                sum += factor * Yi;
            }
            return sum;
        }
        
        Real interpolateSphericalSamples(const Point& p0, const std::vector<Point>& points, const std::vector<Real>& values)
        {
            Real totalSqrDistance = std::accumulate(points.begin(), points.end(), 0.0, [p0](Real res, const Point& p) {
                Real d = p.sphericalDistance(p0);
                return res + d * d;
            });
            
            Real sum = 0.0;
            Real weight = 0.0;
            
            for (size_t i = 0; i < points.size(); ++i)
            {
                const Point& p = points[i];
                Real d = p.sphericalDistance(p0);
                Real w = (totalSqrDistance - d*d) / totalSqrDistance;
                sum += w * values[i];
                weight += w;
            }
            return sum / weight;
        }
        
        Real computeTriangleArea(const Real3& p0, const Real3& p1, const Real3& p2)
        {
            Real3 v12 = p2 - p1;
            Real3 v02 = p2 - p0;
            Real3 v12n = glm::normalize(v12);
            Real t = glm::dot(v02, v12n);
            Real3 c = p2 - v12n * t;
            Real d = glm::distance(p0, c);
            Real l12 = glm::length(v12);
            return l12 * d * 0.5;
        }
        
        void faceAxisDirection(ECubeFace face, Real3& s_dir, Real3& t_dir, Real3& p_dir)
        {
            switch (face)
            {
                case CF_POSX:
                    p_dir = Real3(1, 0, 0);
                    s_dir = Real3(0, 0, -1);
                    t_dir = Real3(0, 1, 0);
                    break;
                case CF_NEGX:
                    p_dir = Real3(-1, 0, 0);
                    s_dir = Real3(0, 0, 1);
                    t_dir = Real3(0, 1, 0);
                    break;
                case CF_POSY:
                    p_dir = Real3(0, 1, 0);
                    s_dir = Real3(0, 0, 1);
                    t_dir = Real3(1, 0, 0);
                    break;
                case CF_NEGY:
                    p_dir = Real3(0, -1, 0);
                    s_dir = Real3(0, 0, 1);
                    t_dir = Real3(-1, 0, 0);
                    break;
                case CF_POSZ:
                    p_dir = Real3(0, 0, 1);
                    s_dir = Real3(1, 0, 0);
                    t_dir = Real3(0, 1, 0);
                    break;
                case CF_NEGZ:
                    p_dir = Real3(0, 0, -1);
                    s_dir = Real3(-1, 0, 0);
                    t_dir = Real3(0, 1, 0);
                    break;
                default:
                    assert(0);
                    p_dir = Real3(1, 0, 0);
                    s_dir = Real3(0, 0, -1);
                    t_dir = Real3(0, 1, 0);
            }
        }
    }
    
}
