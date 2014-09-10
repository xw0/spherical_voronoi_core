//
//  svData.h
//  SphericalVoronoi
//
//  Created by Home on 2014-04-18.
//  Copyright (c) 2014 whenitsdone.org. All rights reserved.
//

#ifndef SphericalVoronoi_svdata_h
#define SphericalVoronoi_svdata_h

#include <memory>
#include <vector>
#include <cassert>
#include <set>
#include <ostream>
#include <iterator>
#include <algorithm>
#include "svMath.h"

namespace sv
{
    class half_edge;
    class cell;
    class vertex;
    class beach_arc;
    class site_event;
    class circle_event;
    
    typedef std::shared_ptr<vertex> vertex_ptr;
    typedef std::shared_ptr<half_edge> half_edge_ptr;
    typedef std::shared_ptr<cell> cell_ptr;
    
    class cell
    {
    public:
        cell (uint32_t i, const Point& p) : index(i), point(p)
        {
        }
        uint32_t index;
        Point point;
        std::vector<half_edge_ptr> halfEdges;
        
        void reset()
        {
            halfEdges.clear();
        }
        
        friend std::ostream& operator<< (std::ostream& stream, const cell& c)
        {
            return stream << "<" << c.index << "> " << "(" << c.point << ")";
        }
    };
    
    class vertex
    {
    public:
        vertex (const Point& p, std::shared_ptr<cell> c0, std::shared_ptr<cell> c1) : point(p)
        {
            cells.insert(c0);
            cells.insert(c1);
        }
        vertex (const Point& p, std::shared_ptr<cell> c0, std::shared_ptr<cell> c1, std::shared_ptr<cell> c2) : point(p)
        {
            cells.insert(c0);
            cells.insert(c1);
            cells.insert(c2);
        }
        uint32_t index;
        Point point;
        std::vector<half_edge_ptr> halfEdges;
        std::set<cell_ptr> cells;
        
        void reset()
        {
            halfEdges.clear();
            cells.clear();
        }
        
        friend std::ostream& operator<< (std::ostream& stream, const vertex& v)
        {
            stream << "point(" << v.point << ")";
            stream << " cells<";
            for (auto c : v.cells)
            {
                stream << c->index << ",";
            }
            stream << ">";
            return stream;
        }
    };
    
    class half_edge
    {
    public:
        half_edge(std::shared_ptr<vertex> s, std::shared_ptr<vertex> e)
        : start(s), end(e)
        {
        }
        uint32_t index;
        cell_ptr cell;
        cell_ptr otherCell;
        vertex_ptr start;
        vertex_ptr end;
        half_edge_ptr prev;
        half_edge_ptr next;
        half_edge_ptr twin;
        
        void reset()
        {
            cell.reset();
            start.reset();
            end.reset();
            prev.reset();
            next.reset();
            twin.reset();
        }
        
        friend std::ostream& operator<< (std::ostream& stream, const half_edge& e)
        {
            stream << "s: " << *e.start << "e: " << *e.end;
            return stream;
        }
    };
    
    class beach_arc
    {
    public:
        beach_arc(std::shared_ptr<cell> cell_)
        : cell(cell_)
        {
        }
        
        std::shared_ptr<cell> cell;
        
        std::shared_ptr<circle_event> circleEvent;      // the related circle event
        
        std::shared_ptr<vertex> startVertex;
        
        bool operator< (const beach_arc& ba) const
        {
            return cell->point.phi < ba.cell->point.phi;
        }
        
        friend std::ostream& operator<< (std::ostream& stream, const beach_arc& arc)
        {
            stream << "cell " << *arc.cell;
            if (arc.startVertex)
            {
                stream << "startVertex " << *arc.startVertex;
            }
            else
            {
                stream << "startVertex NONE";
            }
            return stream;
        }
    };
    
    typedef std::shared_ptr<beach_arc> beach_arc_ptr;
    typedef std::vector<beach_arc_ptr> beach_type;
    
    class site_event
    {
    public:
        site_event(std::shared_ptr<cell> cell_)
        : cell(cell_)
        {
            theta = cell->point.theta;
            phi = cell->point.phi;
        }
        
        std::shared_ptr<cell> cell;
        Real theta;
        Real phi;
        
        bool operator< (const site_event& right) const
        {
            return (theta < right.theta) || (theta == right.theta && phi < right.phi);
        }
        
        friend std::ostream& operator<< (std::ostream& stream, const site_event& e)
        {
            return stream << *e.cell;
        }
    };
    
    class circle_event
    {
    public:
        circle_event(const beach_arc_ptr& arc_i_, const beach_arc_ptr& arc_j_, const beach_arc_ptr& arc_k_)
        : arc_i(arc_i_), arc_j(arc_j_), arc_k(arc_k_)
        {
            using namespace glm;
            
            auto pij = cell_i()->point.position - cell_j()->point.position;
            auto pkj = cell_k()->point.position - cell_j()->point.position;
            auto direction = cross(pij, pkj);
            circle_center = Point(direction);
            circle_radius = acos(dot(circle_center.position, cell_i()->point.position));
            theta = acos(circle_center.position.z) + circle_radius;
        }

        beach_arc_ptr arc_i;
        beach_arc_ptr arc_j;
        beach_arc_ptr arc_k;
        
        cell_ptr cell_i() const { return arc_i->cell; }
        cell_ptr cell_j() const { return arc_j->cell; }
        cell_ptr cell_k() const { return arc_k->cell; }
        
        Point circle_center;
        Real circle_radius;
        
        Real theta;        // the lowest point on circle
        
        bool operator< (const circle_event& ce) const
        {
            return theta < ce.theta;
        }
        
        friend std::ostream& operator<< (std::ostream& stream, const circle_event& e)
        {
            stream << "[" << e.cell_i()->index << "," << e.cell_j()->index << "," << e.cell_k()->index << "] " << "theta " << e.theta;
            return stream;
        }
    };
    
    typedef std::shared_ptr<circle_event> circle_event_ptr;
    
    struct compare_circle_event_priority
    {
        bool operator()(const std::shared_ptr<circle_event>& left, const std::shared_ptr<circle_event>& right) const
        {
            return *left < *right;
        }
    };
    
}

#endif
