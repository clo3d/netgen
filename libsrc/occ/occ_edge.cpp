#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"

#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomLProp_CLProps.hxx>

#pragma clang diagnostic pop

#include "occ_edge.hpp"
#include "occgeom.hpp"

namespace netgen
{
    OCCEdge::OCCEdge(TopoDS_Shape edge_, GeometryVertex & start_, GeometryVertex & end_)
        : GeometryEdge(start_, end_),
          edge(TopoDS::Edge(edge_))
    {
        curve = BRep_Tool::Curve(edge, s0, s1);
        // A degenerate edge (seam/pole) has no 3D curve. BRepGProp::LinearProperties builds a
        // BRepAdaptor_Curve internally, which dereferences the missing 3D curve and crashes on
        // some such edges. Skip it: a degenerate edge has ~zero length, so leaving props at its
        // zero default is correct.
        if(!curve.IsNull())
            BRepGProp::LinearProperties(edge, props);

        auto verts = GetVertices(edge);
        if(verts.size() != 2)
            throw Exception("OCC edge does not have 2 vertices");

        // A degenerate edge (seam/pole, common in imported surface models) has no 3D curve:
        // BRep_Tool::Curve returns null. GetPoint() would then dereference the null curve and
        // crash. Guard the start/end swap; GetPoint is also guarded so meshing can proceed.
        if(!curve.IsNull() && start != end)
        {
            // swap start/end if necessary
            double d00 = Dist(GetPoint(0), start->GetPoint());
            double d01 = Dist(GetPoint(0), end->GetPoint());
            if(d01 < d00)
                swap(start, end);
        }
    }

    double OCCEdge::GetLength() const
    {
        return props.Mass();
    }

    Point<3> OCCEdge::GetCenter() const
    {
        return occ2ng( props.CentreOfMass() );
    }

    Point<3> OCCEdge::GetPoint(double t) const
    {
        // Degenerate edge (no 3D curve): it collapses to its vertex.
        if(curve.IsNull())
            return start->GetPoint();
        return occ2ng( curve->Value(s0+t*(s1-s0)) );
    }

    double OCCEdge::CalcStep(double t, double sag) const
    {
        // Advance the (normalized) parameter by a curvature-adaptive step: small where the
        // curve bends sharply, larger where it is flat. Used by NetgenGeometry::Analyse to
        // sample the edge for local mesh-size restriction. The stock code left this a stub
        // that threw "not implemented", which aborts Analyse for any OCC model whose edges
        // reach this loop (e.g. large imported surface shells) — so meshing yields nothing.
        if(curve.IsNull())
            return 1.0; // degenerate edge collapses to a point: jump straight to the end
        const double range = s1 - s0;
        if(fabs(range) < 1e-30)
            return 1.0;
        const double u = s0 + t*range;
        double dt = 0.1; // fallback: ~10 uniform samples
        try
        {
            // Curvature()/D1() throw (LProp_NotDefined) at non-C2 points of B-spline curves,
            // which are common in imported CAD; fall back to uniform sampling there.
            GeomLProp_CLProps props(curve, u, 2, 1e-9);
            const double speed = props.D1().Magnitude();      // |dC/du|
            if(speed > 1e-12)
            {
                const double curvature = props.Curvature();
                if(curvature > 1e-12)
                {
                    const double radius = 1.0/curvature;
                    // arc length to keep the chord sag within sag*radius: ~2*sqrt(2*sag*radius)
                    const double darc = 2.0*sqrt(max2(0.0, 2.0*sag*radius));
                    dt = (darc / speed) / fabs(range);
                }
                else
                    dt = 0.2; // (near) straight segment: coarse steps are fine
            }
        }
        catch(...)
        {
            dt = 0.1;
        }
        // Never return <=0 (would hang the caller's while-loop) or an oversized step.
        return min2(0.34, max2(0.02, dt));
    }

    void OCCEdge::ProjectPoint(Point<3>& p, EdgePointGeomInfo* gi) const
    {
        auto pnt = ng2occ(p);
        // extend the projection parameter range, else projection might fail
        // for an endpoint
        // see discussion here: https://forum.ngsolve.org/t/how-to-apply-occidentification-correctly/2555
        // I do not see a better way using occ tolerances?
        double eps = 1e-7 * (s1-s0);
        GeomAPI_ProjectPointOnCurve proj(pnt, curve, s0-eps, s1+eps);
        pnt = proj.NearestPoint();
        if(gi)
            gi->dist = (proj.LowerDistanceParameter() - s0)/(s1-s0);
        p = occ2ng(pnt);
    }

    Vec<3> OCCEdge::GetTangent(double t) const
    {
        t = s0 + t*(s1-s0);
        gp_Pnt p;
        gp_Vec v;
        curve->D1(t, p, v);
        return occ2ng(v) * (s1-s0);
    }

}
