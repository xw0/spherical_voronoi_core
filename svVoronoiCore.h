//
//  svVoronoiCore.h
//  SphericalVoronoi
//
//  Created by Xiang Wei on 2014-05-03.
//  Copyright (c) 2014 whenitsdone.org. All rights reserved.
//

#ifndef SphericalVoronoi_svvoronoi_h
#define SphericalVoronoi_svvoronoi_h

#include "svMath.h"
#include "svData.h"

namespace sv
{
    class SphericalVoronoiCore
    {
    public:
        SphericalVoronoiCore(const std::vector<Real3>& directions);

        void setDebugMode(bool debugMode) { this->debugMode = debugMode; }

        bool isFinished() const;
        void step(Real maxDeltaXi);
        void solve(std::function<void(int)> cb = nullptr);       // step until finished

        const std::vector<half_edge_ptr>& getHalfEdges() const { return halfEdges; }
        const std::vector<vertex_ptr>& getVertices() const { return vertices; }
        const std::vector<cell_ptr>& getCells() const { return cells; }

    protected:
        bool debugMode;

        void dumpBeachState(std::ostream& stream);

        void finializeGraph();
        void cleanupMiddleVertices();
        void duplicateHalfEdges();
        void bindHalfEdgesToCells();

        beach_type::const_iterator getPrevArcOnBeach(beach_type::const_iterator it) const
        {
            if (it != beach.begin())
            {
                return std::prev(it);
            }
            else
            {
                return std::prev(beach.end());
            }
        }

        beach_type::const_iterator getNextArcOnBeach(beach_type::const_iterator it) const
        {
            auto next = std::next(it);
            if (next == beach.end())
            {
                next = beach.begin();
            }
            return next;
        }


        bool intersectWithNextArc(beach_type::const_iterator itArc, Real xi, Point& oPoint) const;
        bool intersectWithPrevArc(beach_type::const_iterator itArc, Real xi, Point& oPoint) const;
        void handleSiteEvent(site_event& event);
        void handleCircleEvent(const circle_event_ptr& event);

        static Point thetaToPoint(Real theta, bool positive, Real xi, Real theta1, Real phi1);
        static Point phiToPoint(Real phi, Real xi, Real theta1, Real phi1);
        static bool arcsIntersection(const beach_arc& arc1, const beach_arc& arc2, Real xi, Point& oPoint);

        int nbSteps;
        SphericalLine scanLine;

        constexpr static Real eps = 1e-5;

        std::vector<half_edge_ptr> halfEdges;
        std::vector<vertex_ptr> vertices;
        std::vector<cell_ptr> cells;

        beach_type beach;

        bool isArcOnBeach(const beach_arc_ptr& arc) const
        {
            return find(beach.begin(), beach.end(), arc) != beach.end();
        }

        std::vector<site_event> siteEventQueue;
        std::vector<circle_event_ptr> circleEventQueue;

        void addNewSiteEvent(const site_event& event)
        {
            using namespace std;
            auto it = lower_bound(siteEventQueue.begin(), siteEventQueue.end(), event);
            siteEventQueue.insert(it, event);
        }

        void addNewCircleEvent(const std::shared_ptr<circle_event>& event)
        {
            using namespace std;
            auto it = lower_bound(circleEventQueue.begin(), circleEventQueue.end(), event, compare_circle_event_priority());
            circleEventQueue.insert(it, event);
        }

        void removeCircleEvent(const std::shared_ptr<circle_event>& event)
        {
            using namespace std;
            auto it = find(circleEventQueue.begin(), circleEventQueue.end(), event);
            assert(it != circleEventQueue.end());
            circleEventQueue.erase(it);
        }
    };
}

#endif
