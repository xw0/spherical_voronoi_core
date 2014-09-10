//
//  svVoronoiCore.cpp
//  SphericalVoronoi
//
//  Created by Xiang Wei on 2014-05-03.
//  Copyright (c) 2014 whenitsdone.org. All rights reserved.
//
#include "svPrefix.h"

#include "svVoronoiCore.h"

#define SV_DEBUG(...)   do { if (debugMode) { __VA_ARGS__; } } while(0)

namespace sv
{
    namespace
    {
        struct IndexedDirection
        {
            IndexedDirection(const Real3& dir, int id)
            : direction(dir)
            , index(id)
            {
            }

            Real3 direction;
            int index;
        };
    }

    SphericalVoronoiCore::SphericalVoronoiCore(const std::vector<Real3>& directions)
    : scanLine(), nbSteps(0), debugMode(false)
    {
        using namespace glm;
        using namespace std;

        std::vector<std::pair<int, Real3>> sortedDirections;
        for (auto& dir : directions)
        {
            sortedDirections.push_back(std::pair<int, Real3>(sortedDirections.size(), dir));
        }

        sort(sortedDirections.begin(), sortedDirections.end(), [](const std::pair<int, Real3>& a, const std::pair<int, Real3>& b) { return a.second.z > b.second.z; });
        for (size_t i=0; i<sortedDirections.size(); ++i)
        {
            auto& d = sortedDirections[i].second;
            if (length(d) < eps) continue;
            Point p(d);
            for (auto& cell : cells)
            {
                if (cell->point.equalWithEps(p, eps))
                {
                    continue;
                }
            }
            auto c = shared_ptr<cell>(new cell(sortedDirections[i].first, p));
            cells.push_back(c);

            addNewSiteEvent(site_event(c));
        }
    }

    bool SphericalVoronoiCore::isFinished() const
    {
        return siteEventQueue.size() == 0 && circleEventQueue.size() == 0;
    }

    void SphericalVoronoiCore::dumpBeachState(std::ostream& stream)
    {
        using namespace std;

        stream << "  beach [";
        for (auto itArc=beach.begin(); itArc!=beach.end(); ++itArc)
        {
            auto arc = *itArc;
            stream << arc->cell->index;
            if (arc->circleEvent)
            {
                auto itPrevArc = getPrevArcOnBeach(itArc);
                auto itNextArc = getNextArcOnBeach(itArc);
                assert(arc->circleEvent->arc_j == arc);
                assert(arc->circleEvent->arc_i == *itPrevArc);
                assert(arc->circleEvent->arc_k == *itNextArc);
                stream << "*";
            }
            stream << ",";
        }
        stream << "]" << endl;
    }

    void SphericalVoronoiCore::step(Real maxDeltaXi)
    {
        using namespace std;

        if (!isFinished())
        {
            Real nextXi = scanLine.xi + maxDeltaXi;

            SV_DEBUG(cout << "step " << nbSteps << " " << scanLine.xi);

            if (siteEventQueue.size() > 0 && (circleEventQueue.size() == 0 || siteEventQueue[0].theta < circleEventQueue[0]->theta))
            {
                auto se = siteEventQueue[0];
                if (se.theta <= nextXi)
                {
                    scanLine.xi = se.theta;
                    SV_DEBUG(cout << " -> " << scanLine.xi << endl);
                    SV_DEBUG(dumpBeachState(cout));
                    handleSiteEvent(se);
                    siteEventQueue.erase(siteEventQueue.begin());
                    SV_DEBUG(dumpBeachState(cout));
                }
                else
                {
                    scanLine.xi = nextXi;
                    SV_DEBUG(cout << " -> " << scanLine.xi << endl);
                }
            }
            else if (circleEventQueue.size() > 0)
            {
                auto ce = circleEventQueue[0];
                if (ce->theta <= nextXi)
                {
                    scanLine.xi = ce->theta;
                    SV_DEBUG(cout << " -> " << scanLine.xi << endl);
                    SV_DEBUG(dumpBeachState(cout));
                    handleCircleEvent(ce);
                    circleEventQueue.erase(circleEventQueue.begin());
                    SV_DEBUG(dumpBeachState(cout));
                }
                else
                {
                    scanLine.xi = nextXi;
                    SV_DEBUG(cout << " -> " << scanLine.xi << endl);
                }
            }
            else
            {
                scanLine.xi = nextXi;
                SV_DEBUG(cout << " -> " << scanLine.xi << endl);
            }

            if (isFinished())
            {
                finializeGraph();
            }

            if (debugMode)
            {
                SV_DEBUG(cout << "=================" << endl);
            }

            ++ nbSteps;
        }
    }

    void SphericalVoronoiCore::solve(std::function<void(int)> cb)
    {
        int nbSteps = 0;
        while (!isFinished())
        {
            step(M_PI * 4);
            ++nbSteps;
            if (cb)
            {
                cb(nbSteps);
            }
        }
    }

    void SphericalVoronoiCore::finializeGraph()
    {
        using namespace std;
        using namespace glm;

        SV_DEBUG(cout << "Finalize graph" << endl);

        for (auto& edge: halfEdges)
        {
            auto v0 = edge->start;
            auto v1 = edge->end;
            v0->halfEdges.push_back(edge);
        }

        cleanupMiddleVertices();
        duplicateHalfEdges();
        bindHalfEdgesToCells();

        for (size_t i=0; i<halfEdges.size(); ++i)
        {
            halfEdges[i]->index = (uint32_t)i;
        }

        for (size_t i=0; i<vertices.size(); ++i)
        {
            vertices[i]->index = (uint32_t)i;
        }
    }

    void SphericalVoronoiCore::cleanupMiddleVertices()
    {
        using namespace std;

        SV_DEBUG(cout << "cleanupMiddleVertices" << endl);

        std::vector<vertex_ptr> deleteVertices;
        std::vector<half_edge_ptr> deleteHalfEdges;

        for (auto& v0 : vertices)
        {
            if (v0->cells.size() == 2)
            {
                assert(v0->halfEdges.size() == 2);
                auto v1 = v0->halfEdges[0]->end;
                auto v2 = v0->halfEdges[1]->end;

                auto newEdge = make_shared<half_edge>(v1, v2);
                v1->halfEdges.push_back(newEdge);
                halfEdges.push_back(newEdge);

                deleteVertices.push_back(v0);
                deleteHalfEdges.push_back(v0->halfEdges[0]);
                deleteHalfEdges.push_back(v0->halfEdges[1]);
            }
        }

        SV_DEBUG(cout << "Vertices to delete : " << deleteVertices.size() << endl);
        SV_DEBUG(cout << "Half edges to delete : " << deleteHalfEdges.size() << endl);

        for (auto& v : deleteVertices)
        {
            vertices.erase(find(vertices.begin(), vertices.end(), v));
        }

        for (auto& edge : deleteHalfEdges)
        {
            halfEdges.erase(find(halfEdges.begin(), halfEdges.end(), edge));
        }
    }

    void SphericalVoronoiCore::duplicateHalfEdges()
    {
        using namespace std;

        SV_DEBUG(cout << "duplicateHalfEdges" << endl);

        vector<half_edge_ptr> newHalfEdges;

        for (auto& edge: halfEdges)
        {
            auto newEdge = make_shared<half_edge>(edge->end, edge->start);
            edge->twin = newEdge;
            newEdge->twin = edge;
            edge->end->halfEdges.push_back(newEdge);
            newHalfEdges.push_back(newEdge);
        }

        copy(newHalfEdges.begin(), newHalfEdges.end(), back_inserter(halfEdges));
    }

    void SphericalVoronoiCore::bindHalfEdgesToCells()
    {
        using namespace std;

        SV_DEBUG(cout << "bindHalfEdgesToCells" << endl);

        for (auto& e : halfEdges)
        {
            e->cell.reset();
            e->next.reset();
            e->prev.reset();
        }

        for (auto& e : halfEdges)
        {
            vector<cell_ptr> common;
            set_intersection(e->start->cells.begin(), e->start->cells.end(), e->end->cells.begin(), e->end->cells.end(), back_inserter(common));
            assert(common.size() == 2);
            cell_ptr c0 = common[0];
            Real3 d0 = e->start->point.position - c0->point.position;
            Real3 d1 = e->end->point.position - c0->point.position;
            Real3 n = glm::cross(d0, d1);
            cell_ptr c = nullptr;
            if (glm::dot(n, c0->point.position) > 0)
            {
                c = c0;
            }
            else
            {
                c = common[1];
            }

            e->cell = c;
            c->halfEdges.push_back(e);
        }

        for (auto& c : cells)
        {
            vector<half_edge_ptr> potentialEdges(c->halfEdges);
            c->halfEdges.clear();

            for (auto& e : potentialEdges)
            {
                for (auto& e1 : potentialEdges)
                {
                    if (e1 == e) continue;

                    if (e->end == e1->start)
                    {
                        e->next = e1;
                        e1->prev = e;
                    }
                }
            }

            auto e = potentialEdges[0];
            while (potentialEdges.size() > 0)
            {
                c->halfEdges.push_back(e);
                potentialEdges.erase(find(potentialEdges.begin(), potentialEdges.end(), e));
                e = e->next;
            }

            assert(e == c->halfEdges[0]);
        }

        for (auto& e : halfEdges)
        {
            e->otherCell = e->twin->cell;
        }
    }


    bool SphericalVoronoiCore::arcsIntersection(const beach_arc& arc1, const beach_arc& arc2, Real xi, Point& oPoint)
    {
        Real theta1 = arc1.cell->point.theta;
        Real phi1 = arc1.cell->point.phi;

        Real theta2 = arc2.cell->point.theta;
        Real phi2 = arc2.cell->point.phi;

        if (theta1 >= xi)
        {
            if (theta2 >= xi)
            {
                return false;
            }
            else
            {
                Point pt = phiToPoint(phi1, xi, theta2, phi2);
                oPoint = pt;
                return true;
            }
        }
        if (theta2 >= xi)
        {
            if (theta1 >= xi)
            {
                return false;
            }
            else
            {
                Point pt = phiToPoint(phi2, xi, theta1, phi1);
                oPoint = pt;
                return true;
            }
        }

        Real cos_xi = glm::cos(xi);
        Real sin_xi = glm::sin(xi);
        Real cos_theta1 = glm::cos(theta1);
        Real sin_theta1 = glm::sin(theta1);
        Real cos_theta2 = glm::cos(theta2);
        Real sin_theta2 = glm::sin(theta2);
        Real cos_phi1 = glm::cos(phi1);
        Real sin_phi1 = glm::sin(phi1);
        Real cos_phi2 = glm::cos(phi2);
        Real sin_phi2 = glm::sin(phi2);
        Real a1 = (cos_xi - cos_theta2) * sin_theta1 * cos_phi1;
        Real a2 = (cos_xi - cos_theta1) * sin_theta2 * cos_phi2;
        Real a = a1 - a2;
        Real b1 = (cos_xi - cos_theta2) * sin_theta1 * sin_phi1;
        Real b2 = (cos_xi - cos_theta1) * sin_theta2 * sin_phi2;
        Real b = b1 - b2;
        Real c = (cos_theta1 - cos_theta2) * sin_xi;
        Real l = glm::sqrt(a*a + b*b);
        if (abs(a) > l || abs(c) > l)
        {
            return false;
        }
        else
        {
            auto gamma = glm::atan(a, b);
            auto sin_phi_int_plus_gamma_1 = c / l;
            //auto sin_phi_int_plus_gamma_2 = - c / l;
            auto phi_int_plus_gamma_1 = glm::asin(sin_phi_int_plus_gamma_1);
            //auto phi_int_plus_gamma_2 = asin(sin_phi_int_plus_gamma_2);
            auto pA = phi_int_plus_gamma_1 - gamma;
            //auto pB = phi_int_plus_gamma_2 - gamma;
            Point ptA_1 = phiToPoint(pA, xi, theta1, phi1);
//            point ptA_2 = phiToPoint(pA, xi, theta2, phi2);
//            assert(glm::distance(ptA_1.position, ptA_2.position) < eps);
            oPoint = ptA_1;
            if (oPoint.phi > M_PI)
            {
                oPoint.phi -= M_PI * 2;
            }
            if (oPoint.phi < -M_PI)
            {
                oPoint.phi += M_PI * 2;
            }
            return true;
        }
    }

    bool SphericalVoronoiCore::intersectWithNextArc(beach_type::const_iterator itArc, Real xi, Point& oPoint) const
    {
        auto itNextArc = getNextArcOnBeach(itArc);
        if (itNextArc == itArc)
        {
            return false;
        }
        bool result = arcsIntersection(**itArc, **itNextArc, xi, oPoint);
        return result;
    }

    bool SphericalVoronoiCore::intersectWithPrevArc(beach_type::const_iterator itArc, Real xi, Point& oPoint) const
    {
        auto itPrevArc = getPrevArcOnBeach(itArc);
        if (itPrevArc == itArc)
        {
            return false;
        }
        bool result = arcsIntersection(**itPrevArc, **itArc, xi, oPoint);
        return result;
    }

    void SphericalVoronoiCore::handleSiteEvent(site_event& event)
    {
        using namespace std;

        SV_DEBUG(cout << "HandleSiteEvent " << event << endl);

        if (beach.size() == 0)
        {
            beach.emplace_back(make_shared<beach_arc>(event.cell));
        }
        else if (beach.size() == 1)
        {
            beach.emplace_back(make_shared<beach_arc>(event.cell));
            auto arc = beach[0];
            auto newArc = beach[1];
            Point p = phiToPoint(event.phi, scanLine.xi, arc->cell->point.theta, arc->cell->point.phi);
            std::shared_ptr<vertex> newVertex = std::shared_ptr<vertex>(new vertex(p, event.cell, arc->cell));
            vertices.push_back(newVertex);
            arc->startVertex = newVertex;
            newArc->startVertex = newVertex;
        }
        else
        {
            bool intersectFound = false;
            for (beach_type::const_iterator itArc=beach.begin(); itArc!=beach.end(); ++itArc)
            {
                auto arc = *itArc;
                beach_type::const_iterator itPrevArc = getPrevArcOnBeach(itArc);
                auto prevArc = *itPrevArc;
                beach_type::const_iterator itNextArc = getNextArcOnBeach(itArc);
                auto nextArc = *itNextArc;
                /*auto& nextArc = */*itNextArc;

                Point pointPrev, pointNext;
                bool intPrev, intNext;
                intPrev = intersectWithPrevArc(itArc, scanLine.xi, pointPrev);
                intNext = intersectWithNextArc(itArc, scanLine.xi, pointNext);

                Real phi_start = arc->cell->point.phi - M_PI;
                if (intPrev)
                {
                    phi_start = pointPrev.phi;
                }

                Real phi_end = arc->cell->point.phi + M_PI;
                if (intNext)
                {
                    phi_end = pointNext.phi;
                }

                if (phi_start <= phi_end)
                {
                    intersectFound = phi_start <= event.phi && event.phi <= phi_end;
                }
                else
                {
                    intersectFound = event.phi < phi_end || event.phi > phi_start;
                }
                if (intersectFound)
                {
                    //TODO: handle event.phi == phi_start or phi_end
                    auto vertex1 = prevArc->startVertex;
                    auto vertex2 = arc->startVertex;

                    SV_DEBUG(cout << "Intersect with arc: " << *arc << endl);

                    if (arc->circleEvent)
                    {
                        SV_DEBUG(cout << "  clear circleEvent " << *arc->circleEvent << endl);
                        removeCircleEvent(arc->circleEvent);
                        arc->circleEvent.reset();
                    }
                    beach_type::const_iterator itArc2 = beach.insert(itArc, make_shared<beach_arc>(arc->cell));
                    itArc = std::next(itArc2);
                    auto arc2 = *itArc2;
                    beach_type::const_iterator itNewArc = beach.insert(itArc, make_shared<beach_arc>(event.cell));
                    itArc = std::next(itNewArc);
                    auto newArc = *itNewArc;

                    Point p = phiToPoint(event.phi, scanLine.xi, arc->cell->point.theta, arc->cell->point.phi);
                    std::shared_ptr<vertex> newVertex = std::shared_ptr<vertex>(new vertex(p, event.cell, arc->cell));
                    vertices.push_back(newVertex);

                    arc2->startVertex = newArc->startVertex = newVertex;

                    SV_DEBUG(cout << "  after insert arc: {" << *arc2 << "}, {" << *newArc << "}, {" << *arc << "}, {" << *nextArc << "}" << endl);

                    // refresh all iterators
                    itPrevArc = find(beach.begin(), beach.end(), prevArc);
                    itArc2 = getNextArcOnBeach(itPrevArc);
                    itNewArc = getNextArcOnBeach(itArc2);
                    itArc = getNextArcOnBeach(itNewArc);
                    itNextArc = getNextArcOnBeach(itNextArc);

                    auto ce1 = make_shared<circle_event>(prevArc, arc2, newArc);
                    if (ce1->theta >= event.theta)
                    {
                        arc2->circleEvent = ce1;
                        //itArc2->startVertex = newVertex;
                        SV_DEBUG(cout << "  create new circleEvent " << *ce1 << endl);
                        addNewCircleEvent(ce1);
                    }

                    auto ce2 = make_shared<circle_event>(newArc, arc, nextArc);
                    if (ce2->theta >= event.theta)
                    {
                        arc->circleEvent = ce2;
                        //itArc->startVertex = newVertex;
                        SV_DEBUG(cout << "  create new circleEvent " << *ce2 << endl);
                        addNewCircleEvent(ce2);
                    }

                    if (prevArc->circleEvent)
                    {
                        SV_DEBUG(cout << "  refresh circleEvent " << *prevArc->circleEvent << endl);
                        removeCircleEvent(prevArc->circleEvent);
                        prevArc->circleEvent.reset();
                        auto itPrevPrevArc = getPrevArcOnBeach(itPrevArc);
                        auto ce = make_shared<circle_event>(*itPrevPrevArc, *itPrevArc, *itArc2);
                        addNewCircleEvent(ce);
                        prevArc->circleEvent = ce;
                    }

                    break;
                }
            }
            assert(intersectFound);
        }
    }

    void SphericalVoronoiCore::handleCircleEvent(const circle_event_ptr& event)
    {
        using namespace std;

        SV_DEBUG(cout << "HandleCircleEvent " << *event << endl);

        auto arc_j = event->arc_j;
        auto arc_i = event->arc_i;
        auto arc_k = event->arc_k;

//        assert(isArcOnBeach(arc_j));
//        assert(isArcOnBeach(arc_i));
//        assert(isArcOnBeach(arc_k));

        assert(arc_j->circleEvent == event);
        arc_j->circleEvent.reset();
        if (arc_i->circleEvent)
        {
            SV_DEBUG(cout << "remove circleEvent " << *arc_i->circleEvent << " from arcI " << *arc_i << endl);
            removeCircleEvent(arc_i->circleEvent);
            arc_i->circleEvent.reset();
        }
        if (arc_k->circleEvent)
        {
            SV_DEBUG(cout << "remove circleEvent " << *arc_k->circleEvent << " from arcK " << *arc_k << endl);
            removeCircleEvent(arc_k->circleEvent);
            arc_k->circleEvent.reset();
        }

        auto newVertex = make_shared<vertex>(event->circle_center, arc_i->cell, arc_j->cell, arc_k->cell);
        vertices.push_back(newVertex);

        if (arc_i->startVertex)
        {
            auto edge = make_shared<half_edge>(arc_i->startVertex, newVertex);
            SV_DEBUG(cout << "  create half_edge for arcI " << *arc_i << "[" << *edge << "]" << endl);
            halfEdges.push_back(edge);
        }

        if (arc_j->startVertex)
        {
            auto edge = make_shared<half_edge>(arc_j->startVertex, newVertex);
            SV_DEBUG(cout << "  create half_edge for arcJ " << *arc_j << "[" << *edge << "]" << endl);
            halfEdges.push_back(edge);
        }
        //itArcK->startVertex = newVertex;

        {
            auto it = find(beach.begin(), beach.end(), arc_j);
            assert(it != beach.end());
            SV_DEBUG(cout << "  arc " << *arc_j << " removed from beach");
            beach.erase(it);
        }

        auto itArcI = find(beach.begin(), beach.end(), arc_i);
//        assert(itArcI != beach.end());
        auto itArcK = find(beach.begin(), beach.end(), arc_k);
//        assert(itArcK != beach.end());

        if (getPrevArcOnBeach(itArcI) == itArcK)
        {
            if (arc_k->startVertex)
            {
                auto edge = make_shared<half_edge>(arc_k->startVertex, newVertex);
                SV_DEBUG(cout << "  create half_edge for arcK " << *arc_k << "[" << *edge << "]" << endl);
                halfEdges.push_back(edge);
            }
            SV_DEBUG(cout << "  arc " << *arc_i << " removed from beach");
            beach.erase(itArcI); itArcI = beach.end();

            itArcK = find(beach.begin(), beach.end(), arc_k);
//            assert(itArcK != beach.end());
            SV_DEBUG(cout << "  arc " << *arc_k << " removed from beach");
            beach.erase(itArcK); itArcK = beach.end();
        }
        else
        {
            auto itArc1 = getPrevArcOnBeach(itArcI);
            auto itArc2 = getNextArcOnBeach(itArcK);
            auto& arc_1 = *itArc1;
            auto& arc_2 = *itArc2;

            if (arc_1->cell->index != arc_i->cell->index && arc_i->cell->index != arc_k->cell->index && arc_1->cell->index != arc_k->cell->index)
            {
//                assert(isArcOnBeach(arc_1));
//                assert(isArcOnBeach(arc_i));
//                assert(isArcOnBeach(arc_k));

                auto ceI = make_shared<circle_event>(arc_1, arc_i, arc_k);
                if (ceI->theta >= scanLine.xi)
                {
                    SV_DEBUG(cout << "  create new circleEvent " << *ceI << endl);
                    arc_i->circleEvent = ceI;
                    arc_i->startVertex = newVertex;
                    addNewCircleEvent(ceI);
                }
            }

            if (arc_i->cell->index != arc_k->cell->index && arc_k->cell->index != arc_2->cell->index && arc_i->cell->index != arc_2->cell->index)
            {
//                assert(isArcOnBeach(arc_i));
//                assert(isArcOnBeach(arc_k));
//                assert(isArcOnBeach(arc_2));

                auto ceK = make_shared<circle_event>(arc_i, arc_k, arc_2);
                if (ceK->theta >= scanLine.xi)
                {
                    SV_DEBUG(cout << "  create new circleEvent " << *ceK << endl);
                    arc_k->circleEvent = ceK;
                    //itArcK->startVertex = newVertex;
                    addNewCircleEvent(ceK);
                }
            }
        }
    }

    Point SphericalVoronoiCore::thetaToPoint(Real theta, bool positive, Real xi, Real theta1, Real phi1)
    {
        Real delta_phi = 0;

        if (theta != 0)
        {
            auto theta_p = theta1;
            auto cos_delta_phi = ((glm::cos(xi) - glm::cos(theta_p)) / glm::tan(theta) + glm::sin(xi)) / glm::sin(theta_p);
            cos_delta_phi = glm::clamp<Real>(cos_delta_phi, -1.0, 1.0);
            delta_phi = glm::acos(cos_delta_phi);
        }
        Real s = positive ? 1 : -1;
        Point p(theta, phi1 + delta_phi * s);
        return p;
    }

    Point SphericalVoronoiCore::phiToPoint(Real phi, Real xi, Real theta1, Real phi1)
    {
        if (theta1 >= xi)
        {
            assert(0);
            return Point(xi, phi);      // could be any point on the line segment
        }
        else
        {
            auto a = - (glm::sin(theta1) * glm::cos(phi - phi1) - glm::sin(xi));
            auto b = - (glm::cos(xi) - glm::cos(theta1));
            auto theta = glm::atan(b, a);
            return Point(theta, phi);
        }
    }


}
