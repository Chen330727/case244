#if AXI 

/**
On trees we need refinement functions. */
scalar cm_css_test2[];
face vector fm_fss_test2[];

#if TREE
static void refine_cm_css_test2_axi (Point point, scalar cm_css_test2)
{
#if !EMBED
  fine(cm_css_test2,0,0) = fine(cm_css_test2,1,0) = y - Delta/4.;
  fine(cm_css_test2,0,1) = fine(cm_css_test2,1,1) = y + Delta/4.;
#else // EMBED
  if (css_test2[] > 0. && css_test2[] < 1.) {
    coord n = interface_normal (point, css_test2);
    // Better? involve fss_test2 (troubles w prolongation)
    // coord n = facet_normal (point, css_test2, fss_test2);
    foreach_child() {
      if (css_test2[] > 0. && css_test2[] < 1.) {
	coord p;
    	double alpha = plane_alpha (css_test2[], n);
	plane_center (n, alpha, css_test2[], &p);
	cm_css_test2[] = (y + Delta*p.y)*css_test2[];
      }
      else
	cm_css_test2[] = y*css_test2[];
    }
  }
  else
    foreach_child()
      cm_css_test2[] = y*css_test2[];
#endif // EMBED
}

static void refine_face_fss_test2_x_axi (Point point, scalar fm_fss_test2)
{
#if !EMBED
  if (!is_refined(neighbor(-1))) {
    fine(fm_fss_test2,0,0) = y - Delta/4.;
    fine(fm_fss_test2,0,1) = y + Delta/4.;
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors) {
    fine(fm_fss_test2,2,0) = y - Delta/4.;
    fine(fm_fss_test2,2,1) = y + Delta/4.;
  }
  fine(fm_fss_test2,1,0) = y - Delta/4.;
  fine(fm_fss_test2,1,1) = y + Delta/4.;
#else // EMBED
  double sig = 0., ff = 0.;
  if (css_test2[] > 0. && css_test2[] < 1.) {
    coord n = facet_normal (point, css_test2, fss_test2);
    sig = sign(n.y)*Delta/4.;
  }
  if (!is_refined(neighbor(-1))) {
    ff = fine(fss_test2.x,0,0);
    fine(fm_fss_test2,0,0) = (y - Delta/4. - sig*(1. - ff))*ff;
    ff = fine(fss_test2.x,0,1);
    fine(fm_fss_test2,0,1) = (y + Delta/4. - sig*(1. - ff))*ff;
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors) {
    ff = fine(fss_test2.x,2,0);
    fine(fm_fss_test2,2,0) = (y - Delta/4. - sig*(1. - ff))*ff;
    ff = fine(fss_test2.x,2,1);
    fine(fm_fss_test2,2,1) = (y + Delta/4. - sig*(1. - ff))*ff;
  }
  ff = fine(fss_test2.x,1,0);
  fine(fm_fss_test2,1,0) = (y - Delta/4. - sig*(1. - ff))*ff;
  ff = fine(fss_test2.x,1,1);
  fine(fm_fss_test2,1,1) = (y + Delta/4. - sig*(1. - ff))*ff;
#endif // EMBED
}

static void refine_face_fss_test2_y_axi (Point point, scalar fm_fss_test2)
{
#if !EMBED
  if (!is_refined(neighbor(0,-1)))
    fine(fm_fss_test2,0,0) = fine(fm_fss_test2,1,0) = max(y - Delta/2., 1e-20);
  if (!is_refined(neighbor(0,1)) && neighbor(0,1).neighbors)
    fine(fm_fss_test2,0,2) = fine(fm_fss_test2,1,2) = y + Delta/2.;
  fine(fm_fss_test2,0,1) = fine(fm_fss_test2,1,1) = y;
#else // EMBED
  if (!is_refined(neighbor(0,-1))) {
    fine(fm_fss_test2,0,0) = (max(y - Delta/2., 1e-20))*fine(fss_test2.y,0,0) ;
    fine(fm_fss_test2,1,0) = (max(y - Delta/2., 1e-20))*fine(fss_test2.y,1,0);
  }
  if (!is_refined(neighbor(0,1)) && neighbor(0,1).neighbors) {
    fine(fm_fss_test2,0,2) = (y + Delta/2.)*fine(fss_test2.y,0,2);
    fine(fm_fss_test2,1,2) = (y + Delta/2.)*fine(fss_test2.y,1,2);
  }
  fine(fm_fss_test2,0,1) = y*fine(fss_test2.y,0,1);
  fine(fm_fss_test2,1,1) = y*fine(fss_test2.y,1,1);
#endif // EMBED
}
#endif

/**
If embedded solids are presents, *cm_css_test2*, *fm_fss_test2* and the fluxes need to be
updated consistently with the axisymmetric cylindrical coordinates and
the solid fractions. */

#if EMBED
// double axi_factor (Point point, coord p) {
//   return (y + p.y*Delta);
// }

void cm_css_test2_update (scalar cm_css_test2, scalar css_test2, face vector fss_test2)
{
  foreach() {
    if (css_test2[] > 0. && css_test2[] < 1.) {
      coord p, n = facet_normal (point, css_test2, fss_test2);
      double alpha = plane_alpha (css_test2[], n);
      plane_center (n, alpha, css_test2[], &p);
      cm_css_test2[] = (y + Delta*p.y)*css_test2[];
    }
    else
      cm_css_test2[] = y*css_test2[];
  }
  cm_css_test2[top] = dirichlet(y*css_test2[]);
  cm_css_test2[bottom] = dirichlet(y*css_test2[]);
}

void fm_fss_test2_update (face vector fm_fss_test2, scalar css_test2, face vector fss_test2)
{
  foreach_face(x) {
    double sig = 0.;
    if (css_test2[] > 0. && css_test2[] < 1.) {
      coord n = facet_normal (point, css_test2, fss_test2);
      sig = sign(n.y)*Delta/2.;
    }
    fm_fss_test2.x[] = (y - sig*(1. - fss_test2.x[]))*fss_test2.x[];
  }
  foreach_face(y)
    fm_fss_test2.y[] = max(y, 1e-20)*fss_test2.y[];
  fm_fss_test2.t[top] = dirichlet(y*fss_test2.t[]);
  fm_fss_test2.t[bottom] = dirichlet(y*fss_test2.t[]);
}
#endif // EMBED

// event metric (i = 0) {
void metric_css_test2(){
  /**
  The volume/area of a cell is proportional to $r$ (i.e. $y$). We need
  to set boundary conditions at the top and bottom so that *cm_css_test2* is
  interpolated properly when refining/coarsening the mesh. */

  scalar cm_css_test2v = cm_css_test2;
  foreach()
    cm_css_test2v[] = y;
  cm_css_test2[top] = dirichlet(y);
  cm_css_test2[bottom] = dirichlet(y);

  /**
  We do the same for the length scale factors. The "length" of faces
  on the axis of revolution is zero ($y=r=0$ on the axis). To avoid
  division by zero we set it to epsilon (note that mathematically the
  limit is well posed). */

  face vector fm_fss_test2v = fm_fss_test2;
  foreach_face()
    fm_fss_test2v.x[] = max(y, 1e-20);
  fm_fss_test2.t[top] = dirichlet(y);
  fm_fss_test2.t[bottom] = dirichlet(y);
  
  /**
  We set our refinement/prolongation functions on trees. */

#if TREE
  cm_css_test2.refine = cm_css_test2.prolongation = refine_cm_css_test2_axi;
  fm_fss_test2.x.prolongation = refine_face_fss_test2_x_axi;
  fm_fss_test2.y.prolongation = refine_face_y_axi;
#endif
}

#endif //#if AXI

