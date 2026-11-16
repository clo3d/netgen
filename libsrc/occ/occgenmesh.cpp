#ifdef OCCGEOMETRY

#include <mystdlib.h>
#include <meshing.hpp>
#include <map>

#include "occgeom.hpp"
#include "occ_face.hpp"
#include "occmeshsurf.hpp"

#include <BRepAdaptor_Curve.hxx>
#include <BRepGProp.hxx>
#include <BRepLProp_CLProps.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepTools.hxx>
#include <GProp_GProps.hxx>
#include <Quantity_Color.hxx>
#include <ShapeAnalysis.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopoDS_Edge.hxx>

namespace netgen
{


#define TCL_OK 0
#define TCL_ERROR 1

#define DIVIDEEDGESECTIONS 10000   // better solution to come soon
#define IGNORECURVELENGTH 0
#define VSMALL 1e-10


  DLL_HEADER bool merge_solids = false;


  // can you please explain what you intend to compute here (JS) !!!
  double Line :: Dist (Line l)
  {
    Vec<3> n = p1-p0;
    Vec<3> q = l.p1-l.p0;
    double nq = n*q;

    Point<3> p = p0 + 0.5*n;
    double lambda = (p-l.p0)*n / (nq + VSMALL);

    if (lambda >= 0 && lambda <= 1)
      {
        double d = (p-l.p0-lambda*q).Length();
        //        if (d < 1e-3) d = 1e99;
        return d;
      }
    else
      return 1e99;
  }


  double ComputeH (double kappa, const MeshingParameters & mparam)
  {
    kappa *= mparam.curvaturesafety;
    /*
    double hret;

    if (mparam.maxh * kappa < 1)
      hret = mparam.maxh;
    else
      hret = 1 / (kappa + VSMALL);

    if (mparam.maxh < hret)
      hret = mparam.maxh;

    return hret;
    */
    // return min(mparam.maxh, 1/kappa);
    return max2(mparam.minh, (mparam.maxh*kappa < 1) ? mparam.maxh : 1/kappa);
  }



  void RestrictHTriangle (gp_Pnt2d & par0, gp_Pnt2d & par1, gp_Pnt2d & par2,
                          BRepLProp_SLProps * prop, BRepLProp_SLProps * prop2, Mesh & mesh, int depth, double h, int layer, const MeshingParameters & mparam)
  {
    int ls = -1;

    gp_Pnt pnt0,pnt1,pnt2;

    prop->SetParameters (par0.X(), par0.Y());
    pnt0 = prop->Value();

    prop->SetParameters (par1.X(), par1.Y());
    pnt1 = prop->Value();

    prop->SetParameters (par2.X(), par2.Y());
    pnt2 = prop->Value();

    double aux;
    double maxside = pnt0.Distance(pnt1);
    ls = 2;
    aux = pnt1.Distance(pnt2);
    if(aux > maxside)
      {
        maxside = aux;
        ls = 0;
      }
    aux = pnt2.Distance(pnt0);
    if(aux > maxside)
      {
        maxside = aux;
        ls = 1;
      }



    gp_Pnt2d parmid;

    parmid.SetX( (par0.X()+par1.X()+par2.X()) / 3 );
    parmid.SetY( (par0.Y()+par1.Y()+par2.Y()) / 3 );

    if (depth%3 == 0)
      {
        double curvature = 0;

        prop2->SetParameters (parmid.X(), parmid.Y());
        if (!prop2->IsCurvatureDefined())
          {
            (*testout) << "curvature not defined!" << endl;
            return;
          }
        curvature = max(fabs(prop2->MinCurvature()),
                        fabs(prop2->MaxCurvature()));

        prop2->SetParameters (par0.X(), par0.Y());
        if (!prop2->IsCurvatureDefined())
          {
            (*testout) << "curvature not defined!" << endl;
            return;
          }
        curvature = max(curvature,max(fabs(prop2->MinCurvature()),
                                      fabs(prop2->MaxCurvature())));

        prop2->SetParameters (par1.X(), par1.Y());
        if (!prop2->IsCurvatureDefined())
          {
            (*testout) << "curvature not defined!" << endl;
            return;
          }
        curvature = max(curvature,max(fabs(prop2->MinCurvature()),
                                      fabs(prop2->MaxCurvature())));

        prop2->SetParameters (par2.X(), par2.Y());
        if (!prop2->IsCurvatureDefined())
          {
            (*testout) << "curvature not defined!" << endl;
            return;
          }
        curvature = max(curvature,max(fabs(prop2->MinCurvature()),
                                      fabs(prop2->MaxCurvature())));

        //(*testout) << "curvature " << curvature << endl;

        if (curvature < 1e-3)
          {
            //(*testout) << "curvature too small (" << curvature << ")!" << endl;
            return;
            // return war bis 10.2.05 auskommentiert
          }



        h = ComputeH (curvature+1e-10, mparam);

        if(h < 1e-4*maxside)
          return;


        // if (h > 30) return;
      }

    if (h < maxside && depth < 10)
      {
        //cout << "\r h " << h << flush;
        gp_Pnt2d pm;

        //cout << "h " << h << " maxside " << maxside << " depth " << depth << endl;
        //cout << "par0 " << par0.X() << " " << par0.Y()
        //<< " par1 " << par1.X() << " " << par1.Y()
        //   << " par2 " << par2.X() << " " << par2.Y()<< endl;

        if(ls == 0)
          {
            pm.SetX(0.5*(par1.X()+par2.X())); pm.SetY(0.5*(par1.Y()+par2.Y()));
            RestrictHTriangle(pm, par2, par0, prop, prop2, mesh, depth+1, h, layer, mparam);
            RestrictHTriangle(pm, par0, par1, prop, prop2, mesh, depth+1, h, layer, mparam);
          }
        else if(ls == 1)
          {
            pm.SetX(0.5*(par0.X()+par2.X())); pm.SetY(0.5*(par0.Y()+par2.Y()));
            RestrictHTriangle(pm, par1, par2, prop, prop2, mesh, depth+1, h, layer, mparam);
            RestrictHTriangle(pm, par0, par1, prop, prop2, mesh, depth+1, h, layer, mparam);
          }
        else if(ls == 2)
          {
            pm.SetX(0.5*(par0.X()+par1.X())); pm.SetY(0.5*(par0.Y()+par1.Y()));
            RestrictHTriangle(pm, par1, par2, prop, prop2, mesh, depth+1, h, layer, mparam);
            RestrictHTriangle(pm, par2, par0, prop, prop2, mesh, depth+1, h, layer, mparam);
          }

      }
    else
      {
        gp_Pnt pnt;
        Point3d p3d;

        prop->SetParameters (parmid.X(), parmid.Y());
        pnt = prop->Value();
        p3d = Point3d(pnt.X(), pnt.Y(), pnt.Z());
        mesh.RestrictLocalH (p3d, h, layer);

        p3d = Point3d(pnt0.X(), pnt0.Y(), pnt0.Z());
        mesh.RestrictLocalH (p3d, h, layer);

        p3d = Point3d(pnt1.X(), pnt1.Y(), pnt1.Z());
        mesh.RestrictLocalH (p3d, h, layer);

        p3d = Point3d(pnt2.X(), pnt2.Y(), pnt2.Z());
        mesh.RestrictLocalH (p3d, h, layer);

        //(*testout) << "p = " << p3d << ", h = " << h << ", maxside = " << maxside << endl;

      }
  }

  bool OCCMeshFace (const OCCGeometry & geom, Mesh & mesh, FlatArray<int, PointIndex> glob2loc,
                       const MeshingParameters & mparam, int nr, int projecttype, bool delete_on_failure,
                       double maxh_scale)
  {
    auto k = nr+1;
    if(1==0 && !geom.fvispar[k-1].IsDrawable())
      {
        (*testout) << "ignoring face " << k << endl;
        cout << "ignoring face " << k << endl;
        return true;
      }

    // if(master_faces[k]!=k)
    //     continue;

    (*testout) << "mesh face " << k << endl;
    multithread.percent = 100 * k / (mesh.GetNFD() + VSMALL);
    geom.facemeshstatus[k-1] = -1;

    // FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
    auto face = TopoDS::Face(geom.fmap(k));
    const auto& occface = dynamic_cast<const OCCFace&>(geom.GetFace(k-1));

    int oldnf = mesh.GetNSE();

    Box<3> bb = geom.GetBoundingBox();

    // int projecttype = PLANESPACE;
    // int projecttype = PARAMETERSPACE;
    
    static Timer tinit("init");
    tinit.Start();
    Meshing2OCCSurfaces meshing(geom, face, bb, projecttype, mparam);
    tinit.Stop();


    static Timer tprint("print");
    tprint.Start();
    if (meshing.GetProjectionType() == PLANESPACE)
      PrintMessage (2, "Face ", k, " / ", geom.GetNFaces(), " (plane space projection)");
    else
      PrintMessage (2, "Face ", k, " / ", geom.GetNFaces(), " (parameter space projection)");
    tprint.Stop();

    //      Meshing2OCCSurfaces meshing(f2, bb);
    // meshing.SetStartTime (starttime);
    //(*testout) << "Face " << k << endl << endl;


    auto segments = geom.GetFace(k-1).GetBoundary(mesh);

    if (meshing.GetProjectionType() == PLANESPACE)
      {
        static Timer t("MeshSurface: Find edges and points - Physical"); RegionTimer r(t);
        int cntp = 0;
        glob2loc = 0;

        for (Segment & seg : segments)
          // if (seg.si == k)
            for (int j = 0; j < 2; j++)
              {
                PointIndex pi = seg[j];
                if (glob2loc[pi] == 0)
                  {
                    meshing.AddPoint (mesh.Point(pi), pi);
                    cntp++;
                    glob2loc[pi] = cntp;
                  }
              }
        for(const auto& vert : geom.GetFaceVertices(geom.GetFace(k-1)))
          {
            PointIndex pi = vert->nr + IndexBASE<PointIndex>();
            if(glob2loc[pi] == 0)
              {
                auto gi = occface.Project(mesh[pi]);
                MultiPointGeomInfo mgi;
                mgi.AddPointGeomInfo(gi);
                meshing.AddPoint(mesh[pi], pi, &mgi);
                cntp++;
                glob2loc[pi] = cntp;
              }
          }


        /*
        for (int i = 1; i <= mesh.GetNSeg(); i++)
          {
            Segment & seg = mesh.LineSegment(i);
        */
        // for (Segment & seg : mesh.LineSegments())                
        for (Segment & seg : segments)
          //if (seg.si == k)
            {
              PointGeomInfo gi0, gi1;
              gi0.trignum = gi1.trignum = k;
              gi0.u = seg.epgeominfo[0].u;
              gi0.v = seg.epgeominfo[0].v;
              gi1.u = seg.epgeominfo[1].u;
              gi1.v = seg.epgeominfo[1].v;
              
              //if(orientation & 1)
              meshing.AddBoundaryElement (glob2loc[seg[0]], glob2loc[seg[1]], gi0, gi1);

            }
      }
    else
      {
        static Timer t("MeshSurface: Find edges and points - Parameter"); RegionTimer r(t);

        Array<PointGeomInfo> gis(2*segments.Size());
        gis.SetSize (0);
        glob2loc = 0;
        // int cntpt = 0;

        Box<2> uv_box(Box<2>::EMPTY_BOX);
        for(auto & seg : segments)
            for(auto i : Range(2))
                uv_box.Add( {seg.epgeominfo[i].u, seg.epgeominfo[i].v } );

        BoxTree<2> uv_tree(uv_box);
        double tol = 1e99;
        for(auto& seg : segments)
          {
            Point<2> p1 = { seg.epgeominfo[0].u, seg.epgeominfo[0].v };
            Point<2> p2 = { seg.epgeominfo[1].u, seg.epgeominfo[1].v };
            tol = min2(tol, Dist(p1, p2));
          }
        uv_tree.SetTolerance(0.9 * tol);
        Array<int> found_points;

        for(auto & seg : segments)
        {
            PointGeomInfo gi[2];
            gi[0].trignum = gi[1].trignum = k;
            gi[0].u = seg.epgeominfo[0].u;
            gi[0].v = seg.epgeominfo[0].v;
            gi[1].u = seg.epgeominfo[1].u;
            gi[1].v = seg.epgeominfo[1].v;

            int locpnum[2] = {0, 0};

            for (int j = 0; j < 2; j++)
            {
                Point<2> uv = {gi[j].u, gi[j].v};
                uv_tree.GetIntersecting(uv, uv, found_points);

                bool found = false;
                for(auto& fp : found_points)
                  {
                    if(meshing.GetGlobalIndex(fp - 1) == seg[j])
                      {
                        locpnum[j] = fp;
                        found = true;
                      }
                  }
                if(!found)
                {
                    PointIndex pi = seg[j];
                    locpnum[j] = meshing.AddPoint (mesh.Point(pi), pi) + 1;
                    glob2loc[pi] = locpnum[j];
                    gis.Append (gi[j]);
                    uv_tree.Insert(uv, locpnum[j]);
                }
            }

            meshing.AddBoundaryElement (locpnum[0], locpnum[1], gi[0], gi[1]);
        }
        for(const auto& vert : geom.GetFaceVertices(geom.GetFace(k-1)))
          {
            PointIndex pi = vert->nr + IndexBASE<PointIndex>();
            if(glob2loc[pi] == 0)
              {
                auto gi = occface.Project(mesh[pi]);
                MultiPointGeomInfo mgi;
                mgi.AddPointGeomInfo(gi);
                glob2loc[pi] = meshing.AddPoint(mesh[pi], pi, &mgi) + 1;
                gis.Append(gi);
                Point<2> uv = { gi.u, gi.v };
                uv_tree.Insert(uv, glob2loc[pi]);
              }
          }
      }


    // Philippose - 15/01/2009
    auto& props = occface.properties;
    double maxh = min2(geom.face_maxh[k-1], props.maxh);
    if (maxh_scale != 1.0 && mparam.maxh < 1e98)
      maxh = min2(maxh, mparam.maxh * maxh_scale);
    //double maxh = mparam.maxh;
    //      int noldpoints = mesh->GetNP();
    int noldsurfel = mesh.GetNSE();
    int layer = props.layer;

    static Timer tsurfprop("surfprop");
    tsurfprop.Start();
    GProp_GProps sprops;
    BRepGProp::SurfaceProperties(TopoDS::Face(geom.fmap(k)),sprops);
    tsurfprop.Stop();
    const double facearea = sprops.Mass();
    meshing.SetMaxArea(2.*facearea);

    // A degenerate face (sliver / fully trimmed away) has ~zero surface area, so
    // maxarea = 2*area is non-positive. Meshing2's overshoot guard is `if(maxarea > 0 ...)`,
    // where maxarea <= 0 is the "no area limit" sentinel — so a zero-area face silently
    // disables the guard and the advancing front runs away (millions of elements, front
    // never empties: observed as an endless hang on imported surface shells). Such a face
    // has nothing to mesh; skip it cleanly (0 elements, treated as handled, not failed).
    if (facearea <= 0.0)
      {
        geom.facemeshstatus[k-1] = 1;
        return false;
      }

    MESHING2_RESULT res;

    // TODO: check overlap not correctly working here
    MeshingParameters mparam_without_overlap = mparam;
    mparam_without_overlap.checkoverlap = false;
    
    try {
      static Timer t("GenerateMesh"); RegionTimer reg(t);
      res = meshing.GenerateMesh (mesh, mparam_without_overlap, maxh, k, layer);
    }

    catch (SingularMatrixException)
      {
        // (*myerr) << "Singular Matrix" << endl;
        res = MESHING2_GIVEUP;
      }

    catch (UVBoundsException)
      {
        // (*myerr) << "UV bounds exceeded" << endl;
        res = MESHING2_GIVEUP;
      }

    static Timer t1("rest of loop"); RegionTimer reg1(t1);
      
    bool meshing_failed = res != MESHING2_OK;
    if(meshing_failed && delete_on_failure)
    {
        for (SurfaceElementIndex sei = noldsurfel; sei < mesh.GetNSE(); sei++)
            mesh.Delete(sei);

        mesh.Compress();
    }

    for (SurfaceElementIndex sei = oldnf; sei < mesh.GetNSE(); sei++)
      mesh[sei].SetIndex (k);

    auto n_illegal_trigs = mesh.FindIllegalTrigs();
    PrintMessage (3, n_illegal_trigs, " illegal triangles");
    return meshing_failed;
  }


  // Mesh a face the advancing front gave up on by ear-clipping its existing
  // boundary loop. Only boundary vertices already in the mesh are used, so the
  // result stays conformal with the neighbouring faces. Meant for small
  // degenerate (near-planar) faces; multi-loop boundaries are not handled.
  bool OCCMeshFaceFromBoundary (const OCCGeometry & geom, Mesh & mesh, int nr)
  {
    int k = nr+1;
    auto segments = geom.GetFace(k-1).GetBoundary(mesh);
    size_t nseg = segments.Size();
    if (nseg < 3 || nseg > 200) return false;

    // chain segments into a single closed loop (multi-loop faces are skipped)
    std::map<int,int> nxt;
    for (size_t i = 0; i < nseg; i++)
      {
        int s = int(segments[i][0]);
        if (nxt.count(s)) return false;          // non-manifold boundary
        nxt[s] = int(i);
      }
    NgArray<PointIndex> loop;
    PointIndex start = segments[0][0], cur = start;
    for (size_t s = 0; s < nseg; s++)
      {
        loop.Append (cur);
        auto it = nxt.find (int(cur));
        if (it == nxt.end()) return false;
        cur = segments[it->second][1];
      }
    if (cur != start) return false;              // open chain
    int n = loop.Size();
    if (n != (int)nseg) return false;            // more than one loop

    // Newell normal + 2D basis
    Vec<3> nrm(0,0,0);
    for (int i = 0; i < n; i++)
      {
        Point<3> a = mesh[loop[i]], b = mesh[loop[(i+1)%n]];
        nrm(0) += (a(1)-b(1))*(a(2)+b(2));
        nrm(1) += (a(2)-b(2))*(a(0)+b(0));
        nrm(2) += (a(0)-b(0))*(a(1)+b(1));
      }
    double nl = nrm.Length();
    if (nl < 1e-30) return false;
    nrm *= 1.0/nl;
    Vec<3> u = (fabs(nrm(0)) < 0.9) ? Vec<3>(1,0,0) : Vec<3>(0,1,0);
    u -= (u * nrm) * nrm;  u *= 1.0/u.Length();
    Vec<3> v = Cross (nrm, u);

    NgArray<double> px(n), py(n);
    Point<3> orig = mesh[loop[0]];
    for (int i = 0; i < n; i++)
      {
        Vec<3> d = mesh[loop[i]] - orig;
        px[i] = d * u;  py[i] = d * v;
      }

    auto addTri = [&](PointIndex a, PointIndex b, PointIndex c)
      {
        Element2d el(3);
        el.PNum(1) = a;  el.PNum(2) = b;  el.PNum(3) = c;
        el.SetIndex (k);
        mesh.AddSurfaceElement (el);
      };

    // ear clipping in the projected plane (loop is CCW w.r.t. Newell normal by construction)
    NgArray<int> idx(n);
    for (int i = 0; i < n; i++) idx[i] = i;
    int remaining = n, guard = 0;
    while (remaining > 3 && guard < 2*n*n)
      {
        bool clipped = false;
        for (int i = 0; i < remaining && !clipped; i++, guard++)
          {
            int i0 = idx[(i+remaining-1)%remaining], i1 = idx[i], i2 = idx[(i+1)%remaining];
            double ax=px[i0], ay=py[i0], bx=px[i1], by=py[i1], cx=px[i2], cy=py[i2];
            double cr = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
            if (cr <= 0) continue;               // reflex corner
            bool inside = false;
            for (int j = 0; j < remaining && !inside; j++)
              {
                int p = idx[j];
                if (p==i0 || p==i1 || p==i2) continue;
                double d1 = (bx-ax)*(py[p]-ay) - (by-ay)*(px[p]-ax);
                double d2 = (cx-bx)*(py[p]-by) - (cy-by)*(px[p]-bx);
                double d3 = (ax-cx)*(py[p]-cy) - (ay-cy)*(px[p]-cx);
                if (d1 >= 0 && d2 >= 0 && d3 >= 0) inside = true;
              }
            if (inside) continue;
            addTri (loop[i0], loop[i1], loop[i2]);
            for (int j = i; j < remaining-1; j++) idx[j] = idx[j+1];
            remaining--;
            clipped = true;
          }
        if (!clipped) break;                     // no ear found (degenerate) -> fan the rest
      }
    if (remaining == 3)
      addTri (loop[idx[0]], loop[idx[1]], loop[idx[2]]);
    else
      for (int i = 1; i < remaining-1; i++)      // fallback fan (near-degenerate slivers)
        addTri (loop[idx[0]], loop[idx[i]], loop[idx[i+1]]);

    PrintMessage (3, "face ", k, ": meshed from boundary loop (", n, " points)");
    return true;
  }


  void OCCSetLocalMeshSize(const OCCGeometry & geom, Mesh & mesh,
                           const MeshingParameters & mparam, const OCCParameters& occparam)
  {
    static Timer t1("OCCSetLocalMeshSize");
    RegionTimer regt(t1);
    mesh.SetGlobalH (mparam.maxh);
    mesh.SetMinimalH (mparam.minh);

    NgArray<double> maxhdom;
    maxhdom.SetSize (geom.NrSolids());
    maxhdom = mparam.maxh;
    int maxlayer = 1;

    for(auto dom : Range(geom.GetNSolids()))
    {
      auto & props = geom.GetSolid(dom).properties;
      maxhdom[dom] = min2(maxhdom[dom], props.maxh);
      maxlayer = max2(maxlayer, props.layer);
    }

    for(auto & f : geom.Faces())
      maxlayer = max2(maxlayer, f->properties.layer);

    for(auto & e : geom.Edges())
      maxlayer = max2(maxlayer, e->properties.layer);

    mesh.SetMaxHDomain (maxhdom);

    Box<3> bb = geom.GetBoundingBox();
    bb.Increase (bb.Diam()/10);

    if (mparam.uselocalh)
      {
        const char * savetask = multithread.task;
        multithread.percent = 0;

        for(auto layer : Range(1, maxlayer+1))
            mesh.SetLocalH (bb.PMin(), bb.PMax(), mparam.grading, layer);

        for(auto& v : geom.Vertices())
          {
            auto& props = v->properties;
            if(props.maxh < 1e99)
                mesh.GetLocalH(props.layer)->SetH(v->GetPoint(), props.maxh);
          }

        int nedges = geom.emap.Extent();

        double mincurvelength = IGNORECURVELENGTH;
        double maxedgelen = 0;
        double minedgelen = 1e99;

        if(occparam.resthminedgelenenable) 
          {
            mincurvelength = occparam.resthminedgelen;
            if(mincurvelength < IGNORECURVELENGTH) mincurvelength = IGNORECURVELENGTH;
          }

        multithread.task = "Setting local mesh size (elements per edge)";

        // setting elements per edge
        for (int i = 1; i <= nedges && !multithread.terminate; i++)
          {
            TopoDS_Edge e = TopoDS::Edge (geom.emap(i));
            multithread.percent = 100 * (i-1)/double(nedges);
            if (BRep_Tool::Degenerated(e)) continue;

            double len = Mass(e);

            if (len < mincurvelength)
              {
                (*testout) << "ignored" << endl;
                continue;
              }

            bool is_identified_edge = false;
            // TODO: change to use hash value
            const auto& gedge = geom.GetEdge(e);
            auto& v0 = gedge.GetStartVertex();
            auto& v1 = gedge.GetEndVertex();
            for(auto & ident : v0.identifications)
            {
                auto other = ident.from == &v0 ? ident.to : ident.from;
                if(other->nr == v1.nr && ident.type == Identifications::CLOSESURFACES)
                {
                    is_identified_edge = true;
                    break;
                }
            }

            if(is_identified_edge)
              continue;

            double localh = len/mparam.segmentsperedge;
            double s0, s1;

            Handle(Geom_Curve) c = BRep_Tool::Curve(e, s0, s1);
            if (c.IsNull()) continue; // no 3D curve (not flagged degenerate): c->Value would crash

            const auto & props = gedge.properties;
            localh = min2(localh, props.maxh);
            localh = max2(localh, mparam.minh);
            maxedgelen = max (maxedgelen, len);
            minedgelen = min (minedgelen, len);
            int maxj = max((int) ceil(len/localh)*2, 2);

            for (int j = 0; j <= maxj; j++)
              {
                gp_Pnt pnt = c->Value (s0+double(j)/maxj*(s1-s0));
                mesh.RestrictLocalH (Point3d(pnt.X(), pnt.Y(), pnt.Z()), localh, props.layer);
              }
          }

        multithread.task = "Setting local mesh size (edge curvature)";

        // setting edge curvature

        int nsections = 20;

        for (int i = 1; i <= nedges && !multithread.terminate; i++)
          {
            double maxcur = 0;
            multithread.percent = 100 * (i-1)/double(nedges);
            TopoDS_Edge edge = TopoDS::Edge (geom.emap(i));
            if (BRep_Tool::Degenerated(edge)) continue;
            double s0, s1;
            Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);
            if (c.IsNull()) continue; // no 3D curve (not flagged degenerate): c->Value would crash
            BRepAdaptor_Curve brepc(edge);
            BRepLProp_CLProps prop(brepc, 2, 1e-5);
            auto layer = geom.GetEdge(edge).properties.layer;

            for (int j = 1; j <= nsections; j++)
              {
                double s = s0 + j/(double) nsections * (s1-s0);
                prop.SetParameter (s);
                double curvature = 0;
                if(prop.IsTangentDefined())
                  curvature = prop.Curvature();
                if(curvature> maxcur) maxcur = curvature;

                if (curvature >= 1e99)
                  continue;

                gp_Pnt pnt = c->Value (s);

                mesh.RestrictLocalH (Point3d(pnt.X(), pnt.Y(), pnt.Z()), ComputeH (fabs(curvature), mparam), layer);
              }
          }

        multithread.task = "Setting local mesh size (face curvature)";

        // setting face curvature

        int nfaces = geom.fmap.Extent();

        BuildTriangulation(geom.shape);
        for (int i = 1; i <= nfaces && !multithread.terminate; i++)
          {
            multithread.percent = 100 * (i-1)/double(nfaces);
            TopoDS_Face face = TopoDS::Face(geom.fmap(i));
            TopLoc_Location loc;
            Handle(Geom_Surface) surf = BRep_Tool::Surface (face);
            Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation (face, loc);

            if(triangulation.IsNull())
              {
                if (geom.shape.Infinite())
                  throw Exception("Cannot generate mesh for an infinite geometry");
                // BRepMesh failed to triangulate this one face. Skip its curvature-based
                // local mesh-size contribution instead of aborting the whole model — the
                // face is still meshed later by MeshSurface. Throwing here made a single
                // un-triangulatable face (common in large imported surface shells) discard
                // the entire body's mesh.
                continue;
              }
            
            BRepAdaptor_Surface sf(face, Standard_True);
            // one prop for evaluating and one for derivatives
            BRepLProp_SLProps prop(sf, 0, 1e-5);
            BRepLProp_SLProps prop2(sf, 2, 1e-5);
            auto layer = geom.GetFace(face).properties.layer;

            int ntriangles = triangulation -> NbTriangles();
            for (int j = 1; j <= ntriangles; j++)
              {
                gp_Pnt p[3];
                gp_Pnt2d par[3];

                for (int k = 1; k <=3; k++)
                  {
                    // int n = triangulation->Triangles()(j)(k);
                    // p[k-1] = triangulation->Nodes()(n).Transformed(loc);
                    // par[k-1] = triangulation->UVNodes()(n);
                    // fix for OCC7.6.0-dev
                    int n = triangulation->Triangle(j)(k);
                    p[k-1] = triangulation->Node(n).Transformed(loc);
                    par[k-1] = triangulation->UVNode(n);
                  }

                //double maxside = 0;
                //maxside = max (maxside, p[0].Distance(p[1]));
                //maxside = max (maxside, p[0].Distance(p[2]));
                //maxside = max (maxside, p[1].Distance(p[2]));
                //cout << "\rFace " << i << " pos11 ntriangles " << ntriangles << " maxside " << maxside << flush;

                RestrictHTriangle (par[0], par[1], par[2], &prop, &prop2, mesh, 0, 0, layer, mparam);
                //cout << "\rFace " << i << " pos12 ntriangles " << ntriangles << flush;
              }
          }

        // setting close edges

        if (mparam.closeedgefac.has_value())
          {
            multithread.task = "Setting local mesh size (close edges)";

            int sections = 100;

            NgArray<Line> lines(sections*nedges);

            /*
            BoxTree<3> * searchtree =
              new BoxTree<3> (bb.PMin(), bb.PMax());
            */
            BoxTree<3> searchtree(bb.PMin(), bb.PMax());
            
            int nlines = 0;
            Array<int> edgenumber;
            for (int i = 1; i <= nedges && !multithread.terminate; i++)
              {
                TopoDS_Edge edge = TopoDS::Edge (geom.emap(i));
                int layer = geom.GetEdge(edge).properties.layer;
                if (BRep_Tool::Degenerated(edge)) continue;

                double s0, s1;
                Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);
                BRepAdaptor_Curve brepc(edge);
                BRepLProp_CLProps prop(brepc, 1, 1e-5);
                prop.SetParameter (s0);

                gp_Vec d0 = prop.D1();
                if(d0.Magnitude() > gp::Resolution())
                    d0 = d0.Normalized();
                double s_start = s0;
                // int count = 0;
                for (int j = 1; j <= sections; j++)
                  {
                    double s = s0 + (s1-s0)*(double)j/(double)sections;
                    prop.SetParameter (s);
                    gp_Vec d1 = prop.D1();
                    if(d1.Magnitude() > gp::Resolution())
                        d1 = d1.Normalized();
                    double cosalpha = fabs(d0*d1);
                    if ((j == sections) || (cosalpha < cos(10.0/180.0*M_PI)))
                      {
                        // count++;
                        gp_Pnt p0 = c->Value (s_start);
                        gp_Pnt p1 = c->Value (s);
                        lines[nlines].p0 = Point<3> (p0.X(), p0.Y(), p0.Z());
                        lines[nlines].p1 = Point<3> (p1.X(), p1.Y(), p1.Z());
                        lines[nlines].layer = layer;

                        Box3d box;
                        box.SetPoint (Point3d(lines[nlines].p0));
                        box.AddPoint (Point3d(lines[nlines].p1));

                        searchtree.Insert (box.PMin(), box.PMax(), nlines+1);
                        nlines++;
                        edgenumber.Append(i);

                        s_start = s;
                        d0 = d1;
                      }
                  }
              }

            NgArray<int> linenums;
            auto is_identified_edge = [&](int e0, int e1) {
                const auto& edge0 = geom.GetEdge(e0-1);
                const auto& edge1 = geom.GetEdge(e1-1);

                if(edge0.primary == edge1.primary)
                    return true;

                Array<const GeometryVertex *> v0 = { &edge0.GetStartVertex(), &edge0.GetEndVertex() };
                Array<const GeometryVertex *> v1 = { &edge1.GetStartVertex(), &edge1.GetEndVertex() };
                for(auto i : Range(2))
                    for(auto j : Range(2))
                        if(v0[i]->primary == v1[j]->primary)
                            return true;

                return false;
            };

            for (int i = 0; i < nlines; i++)
              {
                multithread.percent = (100*i)/double(nlines);
                Line & line = lines[i];

                Box3d box;
                box.SetPoint (Point3d(line.p0));
                box.AddPoint (Point3d(line.p1));
                double maxhline = max (mesh.GetH(box.PMin(), line.layer),
                                       mesh.GetH(box.PMax(), line.layer));
                box.Increase(maxhline);

                double mindist = 1e99;
                linenums.SetSize(0);
                searchtree.GetIntersecting(box.PMin(),box.PMax(),linenums);

                for (int j = 0; j < linenums.Size(); j++)
                  {
                    int num = linenums[j]-1;
                    if (i == num) continue;
                    if (line.layer != lines[num].layer) continue;
                    if( is_identified_edge(edgenumber[i], edgenumber[num]) ) continue;
                    if ((line.p0-lines[num].p0).Length2() < 1e-15) continue;
                    if ((line.p0-lines[num].p1).Length2() < 1e-15) continue;
                    if ((line.p1-lines[num].p0).Length2() < 1e-15) continue;
                    if ((line.p1-lines[num].p1).Length2() < 1e-15) continue;
                    mindist = min (mindist, line.Dist(lines[num]));
                  }

                mindist /= (*mparam.closeedgefac + VSMALL);

                if (mindist < 1e-3 * bb.Diam())
                  {
                    (*testout) << "extremely small local h: " << mindist
                               << " --> setting to " << 1e-3 * bb.Diam() << endl;
                    (*testout) << "somewhere near " << line.p0 << " - " << line.p1 << endl;
                    mindist = 1e-3 * bb.Diam();
                  }

                mesh.RestrictLocalHLine(line.p0, line.p1, mindist, line.layer);
              }
          }

        for (auto mspnt : mparam.meshsize_points)
          mesh.RestrictLocalH(mspnt.pnt, mspnt.h, mspnt.layer);

        multithread.task = savetask;

      }

    mesh.LoadLocalMeshSize (mparam.meshsizefilename);
  }
}

#endif
