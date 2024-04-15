/**
# Volume-Of-Fluid advection

We want to approximate the solution of the advection equations
$$
\partial_tc_i + \mathbf{u}_f\cdot\nabla c_i = 0
$$
where $c_i$ are volume fraction fields describing sharp interfaces.

This can be done using a conservative, non-diffusive geometric VOF
scheme.

We also add the option to transport diffusive tracers (aka "VOF
concentrations") confined to one side of the interface i.e. solve the
equations
$$
\partial_tt_{i,j} + \nabla\cdot(\mathbf{u}_ft_{i,j}) = 0
$$
with $t_{i,j} = c_if_j$ (or $t_{i,j} = (1 - c_i)f_j$) and $f_j$ is a
volumetric tracer concentration (see [Lopez-Herrera et al, 2015](#lopez2015)).

The list of tracers associated with the volume fraction is stored in
the *tracers* attribute. For each tracer, the "side" of the interface
(i.e. either $c$ or $1 - c$) is controlled by the *inverse*
attribute). */

attribute {
  scalar * tracers, c;
  bool inverse;
  bool khaki;
}

/**
We will need basic functions for volume fraction computations. */

#include "fractions.h"

/**
The list of volume fraction fields `interfaces`, will be provided by
the user.

The face velocity field `ulf` will be defined by a solver as well
as the timestep. */

extern scalar * interfaces;
extern face vector ulf,ugf; 
extern double dt;
//add for Temperature by cxb, when vof_phase change from 0 to 1 or 1 to 0 
extern double Tsat0;
extern scalar Tlff,Tgff;
extern scalar ff;
extern scalar flux_x,flux_y;
extern vector flux_show;

// extern scalar aiml;
int vof_flag;
int globalii;
// extern scalar flux_out;

//extern scalar css_test3_n,cs;


extern face vector uf;  //fir advect if tracer gas side.
extern double tracex,tracey;
extern double delta_min;

#define EPS3 0.0000000001
//#define CCMIN CCCCMIN 
#define CCMIN CCCCMIN 
// in mycenter.h 0.4

/**
The gradient of a VOF-concentration `t` is computed using a standard
three-point scheme if we are far enough from the interface (as
controlled by *cmin*), otherwise a two-point scheme biased away from
the interface is used. */

foreach_dimension()
static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
{
  //static const double cmin = 0.5;
  static const double cmin = 0.5;
  double cl = c[-1], cc = c[], cr = c[1];
  if (t.inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && t.gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin)
	return t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
      else
	return (t[1]/cr - t[]/cc)/Delta;
    }
    else if (cl >= cmin)
      return (t[]/cc - t[-1]/cl)/Delta;
  }
  return 0.;
}

foreach_dimension()
static double vof_concentration_gradient_f_x (Point point, scalar c, scalar t)
{
  //static const double cmin = 0.5;
  static const double cmin = 0.5;
  double cl = c[-1], cc = c[], cr = c[1];
  double result=0.0;
  int flag=0;
  if (t.inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin) {
    // if(t.gradient != zero){
          // printf("t.gradient !=zero\n");
          if (cr >= cmin) {
            if (cl >= cmin){
        // return t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
              // return (t[1]/cr-t[-1]/cl)/(2.0*Delta);
              result = (t[1]/cr-t[-1]/cl)/(2.0*Delta);
              flag = 1;
            }else{
              // return (t[1]/cr - t[]/cc)/Delta;
               result = (t[1]/cr - t[]/cc)/Delta;
               flag = 2;
            }
          }
          else if (cl >= cmin){
            // return (t[]/cc - t[-1]/cl)/Delta;
             result = (t[]/cc - t[-1]/cl)/Delta;
             flag = 3;
          }
    // }else{
    //     // printf("t.gradient ==zero\n");
    // }
  }
  if(fabs(result)>1e+29){
    printf("result very big: result=%g flag=%d t[]=%g cc=%g t[-1]=%g cl=%g t[1]=%g cr=%g\n",result,flag,t[],cc,t[-1],cl,t[1],cr);
  }
  return result;
}

/**
On trees, VOF concentrations need to be refined properly i.e. using
volume-fraction-weighted linear interpolation of the concentration. */

#if TREE
static void vof_concentration_refine (Point point, scalar s)
{
  scalar f = s.c;
  if (cm[] == 0. || (!s.inverse && f[] <= 0.) || (s.inverse && f[] >= 1.))
    foreach_child()
      s[] = 0.;
  else {
    coord g;
    foreach_dimension()
      g.x = Delta*vof_concentration_gradient_x (point, f, s);
    double sc = s.inverse ? s[]/(1. - f[]) : s[]/f[], cmc = 4.*cm[];
    foreach_child() {
      s[] = sc;
      foreach_dimension()
	s[] += child.x*g.x*cm[-child.x]/cmc;
      s[] *= s.inverse ? 1. - f[] : f[];
    }
  }
}

/**
On trees, we need to setup the appropriate prolongation and
refinement functions for the volume fraction fields. */

event defaults (i = 0)
{
  for (scalar c in interfaces) {
    c.refine = c.prolongation = fraction_refine;
    c.dirty = true;
    scalar * tracers = c.tracers;
    for (scalar t in tracers) {
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = vof_concentration_refine;
      t.dirty = true;
      t.c = c;
    }
  }
}
#endif // TREE

/**
Boundary conditions for VOF-advected tracers usually depend on
boundary conditions for the VOF field. */

event defaults (i = 0)
{
  for (scalar c in interfaces) {
    scalar * tracers = c.tracers;
    for (scalar t in tracers)
      t.depends = list_add (t.depends, c);
  }
}

/**
We need to make sure that the CFL is smaller than 0.5 to ensure
stability of the VOF scheme. */

event stability (i++) {
  if (CFL > 0.5)
    CFL = 0.5;
}

/**
## One-dimensional advection

The simplest way to implement a multi-dimensional VOF advection scheme
is to use dimension-splitting i.e. advect the field along each
dimension successively using a one-dimensional scheme.

We implement the one-dimensional scheme along the x-dimension and use
the [foreach_dimension()](/Basilisk C#foreach_dimension) operator to
automatically derive the corresponding functions along the other
dimensions. */


//add for Temperature by cxb, when vof_phase change from 0 to 1 or 1 to 0 


foreach_dimension()
static void sweep_x (scalar c, scalar cc, scalar * tcl)
{
  vector n[];
  scalar alpha[], flux[];
  double cfl = 0.;
   
//add for Temperature by cxb, when vof_phase change from 0 to 1 or 1 to 0 

/*
  scalar cold[];
  foreach(){
    cold[] = c[];
  }
*/


  /**
  If we are also transporting tracers associated with $c$, we need to
  compute their gradient i.e. $\partial_xf_j = \partial_x(t_j/c)$ or
  $\partial_xf_j = \partial_x(t_j/(1 - c))$ (for higher-order
  upwinding) and we need to store the computed fluxes. We first
  allocate the corresponding lists. */

  scalar * tracers = c.tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    for (scalar t in tracers) {
      scalar gf = new scalar, flux = new scalar;
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }

    /**
    The gradient is computed using the "interface-biased" scheme above. */

    foreach() {
      scalar t, gf;
      for (t,gf in tracers,gfl)
	gf[] = vof_concentration_gradient_x (point, c, t);
    }
  }
  
  /**
  We reconstruct the interface normal $\mathbf{n}$ and the intercept
  $\alpha$ for each cell. Then we go through each (vertical) face of
  the grid. */

  reconstruction (c, n, alpha);
  foreach_face(x, reduction (max:cfl)) {

    /**
    To compute the volume fraction flux, we check the sign of the velocity
    component normal to the face and compute the index `i` of the
    corresponding *upwind* cell (either 0 or -1). */

    double un = ulf.x[]*dt/(Delta*fm.x[] + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;

    /**
    We also check that we are not violating the CFL condition. */

#if EMBED
//#if 1==1
    printf("vof EMBED=true\n");
     //if (css_test3_n[] >= 1.)
    if (cs[] >= 1.)
#endif
    if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
      cfl = un*fm.x[]*s/(cm[] + SEPS);

    /**
    If we assume that `un` is negative i.e. `s` is -1 and `i` is 0, the
    volume fraction flux through the face of the cell is given by the dark
    area in the figure below. The corresponding volume fraction can be
    computed using the `rectangle_fraction()` function.
    
    ![Volume fraction flux](figures/flux.svg)
    
    When the upwind cell is entirely full or empty we can avoid this
    computation. */

    double cf; // fixme: ternary operator not properly detected by qcc stencil
    if (c[i] <= 0. || c[i] >= 1.)
      cf = c[i];
    else
      cf = rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
			       (coord){-0.5, -0.5, -0.5},
			       (coord){s*un - 0.5, 0.5, 0.5});
    
    /**
    Once we have the upwind volume fraction *cf*, the volume fraction
    flux through the face is simply: */

    flux[] = cf*ulf.x[];

    /**
    If we are transporting tracers, we compute their flux using the
    upwind volume fraction *cf* and a tracer value upwinded using the
    Bell--Collela--Glaz scheme and the gradient computed above. */
    
    scalar t, gf, tflux;
    for (t,gf,tflux in tracers,gfl,tfluxl) {
      double cf1 = cf, ci = c[i];
      if (t.inverse)
	cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
	double fff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
	tflux[] = fff*cf1*ulf.x[];
      }
      else
	tflux[] = 0.;
    }
  }
  delete (gfl); free (gfl);
  
  /**
  We warn the user if the CFL condition has been violated. */

  if (cfl > 0.5 + 1e-6)
    fprintf (ferr, 
	     "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n", 
	     cfl - 0.5), fflush (ferr);

  /**
  Once we have computed the fluxes on all faces, we can update the
  volume fraction field according to the one-dimensional advection
  equation
  $$
  \partial_tc = -\nabla_x\cdot(\mathbf{u}_f c) + c\nabla_x\cdot\mathbf{u}_f
  $$
  The first term is computed using the fluxes. The second term -- which is
  non-zero for the one-dimensional velocity field -- is approximated using
  a centered volume fraction field `cc` which will be defined below. 

  For tracers, the one-dimensional update is simply
  $$
  \partial_tt_j = -\nabla_x\cdot(\mathbf{u}_f t_j)
  $$
  */

#if !EMBED
//#if 1==0
  foreach() {
    c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta);
    scalar t, tc, tflux;
    for (t, tc, tflux in tracers, tcl, tfluxl)
      t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta);
//add for Temperature by cxb, when vof_phase change from 0 to 1 or 1 to 0 

    /*
    for (t in {Tlff,Tgff}){
      if((c[]-0.5)*(cold[]-0.5)<=0){
           if(t.inverse){
                t[] = Tsat0*(1.0-c[]);
           }else{            
                t[] = Tsat0*c[];
           }
       }
    }
   */

  }
#else // EMBED
  /**
  When dealing with embedded boundaries, we simply ignore the fraction
  occupied by the solid. This is a simple approximation which has the
  advantage of ensuring boundedness of the volume fraction and
  conservation of the total tracer mass (if it is computed also
  ignoring the volume occupied by the solid in partial cells). */
  
  foreach()
    if (cs[] > 0.) {
    //if (css_test3_n[] > 0.) {
      // c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/Delta; // ulf only used in liquid side, other wise, it is not divergence free
      c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/Delta;
      scalar t, tc, tflux;
      for (t, tc, tflux in tracers, tcl, tfluxl)
	t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/Delta;
    }
#endif // EMBED

  delete (tfluxl); free (tfluxl);
}

extern int level_interface;
void get_topomask2(scalar topo_mask_vof){

     
     foreach(){
        int phase_sign = (ff[]>=0.5-EPS3) ? 1 : -1;
	topo_mask_vof[] = 3*phase_sign;
        if(ff[]<1.0-EPS3 && ff[]>EPS3){
 		topo_mask_vof[] = 0;
        }
     }
     foreach(){
        if(level==level_interface){
          bool is1= false;
          foreach_dimension(){
                  int temp1=topo_mask_vof[1];
                  int temp2=topo_mask_vof[-1];
                  if(temp1==0 || temp2==0){
                      is1 = true;
                  }
          }
          int temp3=topo_mask_vof[];
          if (is1 && temp3!=0){
             topo_mask_vof[] = (ff[]>=0.5-EPS3) ? 1 : -1 ;
          }
        }
     } 
  for(int phase=0;phase<=1;phase++){
     foreach(){
        if(level==level_interface){
            int phase_sign = ff[]>=0.5-EPS3 ? 1 : -1;
            bool is1= false;
            foreach_dimension(){
              int temp1=topo_mask_vof[1];
              int temp2=topo_mask_vof[-1];
              if( temp1==(2*phase-1) || temp2==(2*phase-1)){
                is1 = true;
              }
            }
            int temp3 = topo_mask_vof[];
            if (is1 && temp3==3*phase_sign){
              topo_mask_vof[] = 2*(2*phase-1) ;
            }
        }
     }
  }


       
}



///////////////////////////////////////////////////////////////////////////////
//// for advection of Tgff
//// Xiangbin 20220812
////////////////////////////////////////////////////////////////////////////////


#include "./reconstruction4.h"


// foreach_dimension()
// static void sweep2_x (scalar c, scalar cc, scalar flux_temp, scalar * tcl)
foreach_dimension()
static void sweep2_x (scalar c, scalar cc, scalar * tcl)
{
  vector n[];
  scalar alpha[], flux[];
  double cfl = 0.;
   
//add for Temperature by cxb, when vof_phase change from 0 to 1 or 1 to 0 

/*
  scalar cold[];
  foreach(){
    cold[] = c[];
  }
*/


fprintf(stderr,"sweep2 vof 434\n");
  /**
  If we are also transporting tracers associated with $c$, we need to
  compute their gradient i.e. $\partial_xf_j = \partial_x(t_j/c)$ or
  $\partial_xf_j = \partial_x(t_j/(1 - c))$ (for higher-order
  upwinding) and we need to store the computed fluxes. We first
  allocate the corresponding lists. */

  scalar flux2[];//tflux2[];

  scalar fluxg[];//tfluxg[];

  // scalar flux_temp[];

  scalar * tracers = c.tracers, * gfl = NULL, * tfluxl = NULL;
  scalar * gfl2 = NULL, * gflg = NULL, * tfluxl2 = NULL, * tfluxlg = NULL;

  fprintf(stderr,"sweep2 tracers 449\n");
  if (tracers) {
    for (scalar t in tracers) {
      fprintf(stderr,"sweep2 tracers 450\n");
      scalar gf = new scalar, flux = new scalar;
      scalar gf2 = new scalar, flux2 = new scalar;
      scalar gfg = new scalar, fluxg = new scalar;
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
      gfl2 = list_append (gfl2, gf2);
      tfluxl2 = list_append (tfluxl2, flux2);
      gflg = list_append (gflg, gfg);
      tfluxlg = list_append (tfluxlg, fluxg);
       fprintf(stderr,"sweep2 tracers 462\n");
    }

    /**
    The gradient is computed using the "interface-biased" scheme above. */
fprintf(stderr,"sweep2 tracers 467\n");
    foreach() {
      scalar t, gf, gf2, gfg;
      for (t, gf, gf2, gfg in tracers, gfl, gfl2, gflg){
	      gf[] = vof_concentration_gradient_f_x (point, c, t);
        // printf("  %g",gf[]);
        if(fabs(gf[])>1e+29){
            printf("gf[]=%g, is true\n",gf[]);
        }
        gf2[] = gf[];
        gfg[] = gf[];
      }
    }
    fprintf(stderr,"sweep2 tracers 474\n");
  }
  boundary(gfl);
  boundary(gfl2);
  boundary(gflg);
  /**
  We reconstruct the interface normal $\mathbf{n}$ and the intercept
  $\alpha$ for each cell. Then we go through each (vertical) face of
  the grid. */

fprintf(stderr,"sweep2 vof 479\n");

if(1==1){
  reconstruction (c, n, alpha);
}else{
  reconstruction4 (c, n, alpha);
}
   scalar topo_mask2[];
   get_topomask2(topo_mask2);

   scalar topo_mask2_g[];
   foreach(){
      topo_mask2_g[] = -topo_mask2[];
   }
fprintf(stderr,"sweep2 vof 493\n");


  foreach_face(x, reduction (max:cfl)) {

    /**
    To compute the volume fraction flux, we check the sign of the velocity
    component normal to the face and compute the index `i` of the
    corresponding *upwind* cell (either 0 or -1). */

   
   
   double khaki_velocity=0.0;
//   if(((topo_mask2[-1]==-1) && (topo_mask2[]==0))  || (((topo_mask2[-1]==0) && (topo_mask2[]==1)) )){
//       if()
//    }


    double un = ulf.x[]*dt/(Delta*fm.x[] + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;

    double ung = ugf.x[]*dt/(Delta*fm.x[] + SEPS), sg = sign(ung);
    int ig = -(sg + 1.)/2.;

    double un2 = uf.x[]*dt/(Delta*fm.x[] + SEPS), s2 = sign(un2);
    int i2 = -(s2 + 1.)/2.;

    /**
    We also check that we are not violating the CFL condition. */

#if EMBED
//#if 1==1
    if (cs[] >= 1.){
   // if (css_test3_n[] >= 1.){  
        if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
             cfl = un*fm.x[]*s/(cm[] + SEPS);

        if (ung*fm.x[]*sg/(cm[] + SEPS) > cfl)
             cfl = ung*fm.x[]*sg/(cm[] + SEPS);

        if (un2*fm.x[]*s2/(cm[] + SEPS) > cfl)
             cfl = un2*fm.x[]*s2/(cm[] + SEPS);
    }
#else
     if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
      cfl = un*fm.x[]*s/(cm[] + SEPS);

    if (ung*fm.x[]*sg/(cm[] + SEPS) > cfl)
      cfl = ung*fm.x[]*sg/(cm[] + SEPS);

    if (un2*fm.x[]*s2/(cm[] + SEPS) > cfl)
      cfl = un2*fm.x[]*s2/(cm[] + SEPS);
#endif
   

    /**
    If we assume that `un` is negative i.e. `s` is -1 and `i` is 0, the
    volume fraction flux through the face of the cell is given by the dark
    area in the figure below. The corresponding volume fraction can be
    computed using the `rectangle_fraction()` function.
    
    ![Volume fraction flux](figures/flux.svg)
    
    When the upwind cell is entirely full or empty we can avoid this
    computation. */

    double cf; // fixme: ternary operator not properly detected by qcc stencil
    double cf2;
     double cfg;
    if (c[i] <= 0. || c[i] >= 1.){
      cf = c[i];
    }else{
      cf = rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
			       (coord){-0.5, -0.5, -0.5},
			       (coord){s*un - 0.5, 0.5, 0.5});
    }
    /**
    Once we have the upwind volume fraction *cf*, the volume fraction
    flux through the face is simply: */

    if (c[i2] <= 0. || c[i2] >= 1.){
      cf2 = c[i2];
    }else{  
      cf2 = rectangle_fraction ((coord){-s2*n.x[i2], n.y[i2], n.z[i2]}, alpha[i2],
			       (coord){-0.5, -0.5, -0.5},
			       (coord){s2*un2 - 0.5, 0.5, 0.5});
    }
    
   if (c[ig] <= 0. || c[ig] >= 1.){
      cfg = c[ig];
    }else{
      cfg = rectangle_fraction ((coord){-sg*n.x[ig], n.y[ig], n.z[ig]}, alpha[ig],
			       (coord){-0.5, -0.5, -0.5},
			       (coord){sg*ung - 0.5, 0.5, 0.5});
    }

    flux[] = cf*ulf.x[];
    flux2[] = cf2*uf.x[];
    fluxg[] = cfg*ugf.x[];
    /**
    If we are transporting tracers, we compute their flux using the
    upwind volume fraction *cf* and a tracer value upwinded using the
    Bell--Collela--Glaz scheme and the gradient computed above. */
    
    scalar t, gf, tflux, gf2, tflux2, gfg, tfluxg;
    //scalar tflux2[];
    // double ddd;
    for (t,gf,tflux,gf2,tflux2,gfg,tfluxg in tracers,gfl,tfluxl,gfl2,tfluxl2,gflg,tfluxlg) {
      double cf1 = cf, ci = c[i];
      double cf21 = cf2, ci2 = c[i2];
      double cf1g = cfg, cig = c[ig];
      
      if (t.inverse){
            cf1 = 1. - cf1, ci = 1. - ci;
       }
      if (ci > 1e-10) {
          double fff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
          tflux[] = fff*cf1*ulf.x[];
          if(!t.inverse){
            flux_show.x[]=tflux[];
              // flux_out[] = tflux[];
              // flux_temp[] = tflux[];
              // if(fabs(x-tracex+delta_min/2.0*2)<1e-6 && fabs(y-tracey)<1e-6){
              if(fabs(x-tracex+delta_min/2.0*2)<1e-6 && (fabs(y-tracey-delta_min/2.0)<1e-6 || fabs(y-tracey+delta_min/2.0)<1e-6)){
                  printf("x face: tflux=%g fff=%g cf1=%g ulf.x[]=%g gf[i]=%g intx=%g inty=%g\n",tflux[],fff,cf1,ulf.x[],gf[i],x/delta_min,y/delta_min);
              }
              // if(fabs(x-tracex)<1e-6 && fabs(y-tracey+delta_min/2.0*2)<1e-6){
                if((fabs(x-tracex-delta_min/2.0)<1e-6 || fabs(x-tracex+delta_min/2.0)<1e-6) && fabs(y-tracey+delta_min/2.0*2)<1e-6){
                  printf("y face: tflux=%g fff=%g cf1=%g ulf.x[]=%g gf[i]=%g  intx=%g inty=%g\n",tflux[],fff,cf1,ulf.x[],gf[i],x/delta_min,y/delta_min);
              }

              if(fabs(tflux[]>1e16)){
                printf("tflux too big, x=%g y=%g  tflux=%g fff=%g cf1=%g ulf.x[]=%g t[i]=%g gf[i]=%g\n",x,y,tflux[],fff,cf1,ulf.x[],t[i],gf[i]);
              }
          }
      }
      else{
	        tflux[] = 0.;
          if(!t.inverse){
              // flux_temp[] =-1;
              //  flux_out[] =-1;
              // ddd = -1;
          }
      }

      if (t.inverse){
        cf21 = 1. - cf21, ci2 = 1. - ci2;
       }
      if (ci2 > 1e-10) {
        double fff2 = t[i2]/ci2 + s2*min(1., 1. - s2*un2)*gf[i2]*Delta/2.;
        tflux2[] = fff2*cf21*uf.x[];
      }
      else{
        tflux2[]=0.;
      }

       if (t.inverse){
            cf1g = 1. - cf1g, cig = 1. - cig;
       }
      if (cig > 1e-10) {
          double fffg = t[ig]/cig + sg*min(1., 1. - sg*ung)*gf[ig]*Delta/2.;
          tfluxg[] = fffg*cf1g*ugf.x[];
      }
      else{
	        tfluxg[] = 0.;
      }
    }
    // flux_temp[] = ddd;
  }

  #if TREE
  scalar * fluxl = list_concat (NULL, tfluxl);
  fluxl = list_concat (fluxl, tfluxl2);
  fluxl = list_concat (fluxl, tfluxlg);
  fluxl = list_append (fluxl, flux);
  for (int l = depth() - 1; l >= 0; l--)
    foreach_halo (prolongation, l) {
#if dimension == 1
      if (is_refined (neighbor(-1)))
	for (scalar fl in fluxl)
	  fl[] = fine(fl);
      if (is_refined (neighbor(1)))
	for (scalar fl in fluxl)
	  fl[1] = fine(fl,2);
#elif dimension == 2
      if (is_refined (neighbor(-1)))
	for (scalar fl in fluxl)
	  fl[] = (fine(fl,0,0) + fine(fl,0,1))/2.;
      if (is_refined (neighbor(1)))
	for (scalar fl in fluxl)
	  fl[1] = (fine(fl,2,0) + fine(fl,2,1))/2.;
#else // dimension == 3
      if (is_refined (neighbor(-1)))
	for (scalar fl in fluxl)
	  fl[] = (fine(fl,0,0,0) + fine(fl,0,1,0) +
		  fine(fl,0,0,1) + fine(fl,0,1,1))/4.;
      if (is_refined (neighbor(1)))
	for (scalar fl in fluxl)
	  fl[1] = (fine(fl,2,0,0) + fine(fl,2,1,0) +
		   fine(fl,2,0,1) + fine(fl,2,1,1))/4.;
#endif
    }
  free (fluxl);
#endif

  // if(vof_flag==1){
    // // scalar ttt,tt;
    // // for(ttt, tt in tracers, tfluxl){
    // //             if(!ttt.inverse){
    // //         char name93[80]; 
    // //         sprintf(name93,"vofout-dimension%d-%d.dat",vof_flag,pid());
    // //         FILE * fp93 = fopen(name93,"w");
    // //         int num2=0;
    // //         foreach(){
    // //           //fprintf (fp93,"%g %g %g %g %g %g\n",  x, y, z, T[], Tgff[],Tlff[]);
    // //               fprintf (fp93,"%g %g %g\n",  x, y, tt[]);
    // //         }
    // //         fclose(fp93);
    // //         MPI_Barrier (MPI_COMM_WORLD);
    // //         if(pid()==0){
    // //                   char command1[150];
    // //                   sprintf(command1, "LC_ALL=C cat vofout-dimension%d*.dat > outfacets/vofout-dimension%d-%d",vof_flag,vof_flag,globalii);
    // //                   system(command1);

    // //                   char command7[150];
    // //                   sprintf(command7, "LC_ALL=C rm -rf vofout-dimension*.dat");
    // //                   system(command7);
    // //               }
    // //        }
    // //    }
  // }
// foreach(){
//   aiml[] = flux_temp[];
// }
  // boundary(tfluxl);
  // restriction(tfluxl);



  // int l = maxl-1;
  // foreach_halo (prolongation, l){
  //   scalar t1;
  //   for(t1 in tfluxl){
  //     printf("halo_l=%d: x = %g y = %g a = %g\n",l, x, y, t1[]);
  //   }
  // }


  fprintf(stderr,"sweep2 vof 635\n");
  delete (gfl); free (gfl);
  delete (gfl2); free (gfl2);
  delete (gflg); free (gflg);
  
  /**
  We warn the user if the CFL condition has been violated. */

  if (cfl > 0.5 + 1e-6)
    fprintf (ferr, 
	     "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n", 
	     cfl - 0.5), fflush (ferr);

  /**
  Once we have computed the fluxes on all faces, we can update the
  volume fraction field according to the one-dimensional advection
  equation
  $$
  \partial_tc = -\nabla_x\cdot(\mathbf{u}_f c) + c\nabla_x\cdot\mathbf{u}_f
  $$
  The first term is computed using the fluxes. The second term -- which is
  non-zero for the one-dimensional velocity field -- is approximated using
  a centered volume fraction field `cc` which will be defined below. 

  For tracers, the one-dimensional update is simply
  $$
  \partial_tt_j = -\nabla_x\cdot(\mathbf{u}_f t_j)
  $$
  */

#if !EMBED
fprintf(stderr,"665 vof\n");
    // printf("my_vof exit\n");
    // exit(1);
    foreach() {
     // if(cs[]>0.0){
     // if(css_test3_n[]>0.0){
            //if(topo_mask2[]>=-1) // add 20230103
              c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta);

            scalar t, tc, tflux;
            scalar tflux2, tfluxg;
            for (t, tc, tflux, tflux2, tfluxg in tracers, tcl, tfluxl, tfluxl2, tfluxlg){
              int khaki_flag=0;
              if(t.inverse){ 
                  if(topo_mask2_g[]>=-1){ //becase vof method in paris is different from here.
                       // t[] += dt*(tfluxg[] - tfluxg[1] + tc[]*(ugf.x[1] - ugf.x[]))/(cm[]*Delta);  

                       // to be more accurate Tgff shoud be advected by ugf, here using ulf   
                        t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta);            
                  }
              } else{
                if(topo_mask2[]>=-1){ 
                        t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta);
                       // t[] += dt*(tflux2[] - tflux2[1] + tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
                }
              }

            }
        //add for Temperature by cxb, when vof_phase change from 0 to 1 or 1 to 0 

            /*
            for (t in {Tlff,Tgff}){
              if((c[]-0.5)*(cold[]-0.5)<=0){
                  if(t.inverse){
                        t[] = Tsat0*(1.0-c[]);
                  }else{            
                        t[] = Tsat0*c[];
                  }
              }
            }
          */
     // }
  }
//#if 1==0
// // // // //   foreach() {
// // // // //     if(topo_mask2[]>=-1) // add 20230103
// // // // //       c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta);

// // // // //     scalar t, tc, tflux,t2;
// // // // //     for (t, tc, tflux in tracers, tcl, tfluxl){
// // // // //        int khaki_flag=0;
// // // // //        if(t.khaki){
// // // // //            if(topo_mask2[]<=-1){ //becase vof method in paris is different from here.
// // // // //                 t[] += dt*(tflux2[] - tflux2[1] + tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta); 
// // // // //                //t[] += dt*(tflux2[] - tflux2[1]+ tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
// // // // //                // t[] += dt*(tflux2[] - tflux2[1])/(cm[]*Delta); 
// // // // //                 khaki_flag = 1;              
// // // // //            }
// // // // //        }
// // // // //        if(!khaki_flag){
// // // // //           t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta);
// // // // //           //t[] += dt*(tflux2[] - tflux2[1]+ tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
// // // // //           //t[] += dt*(tflux[] - tflux[1])/(cm[]*Delta);
// // // // //        }

// // // // //     }
// // // // // //add for Temperature by cxb, when vof_phase change from 0 to 1 or 1 to 0 

// // // // //     /*
// // // // //     for (t in {Tlff,Tgff}){
// // // // //       if((c[]-0.5)*(cold[]-0.5)<=0){
// // // // //            if(t.inverse){
// // // // //                 t[] = Tsat0*(1.0-c[]);
// // // // //            }else{            
// // // // //                 t[] = Tsat0*c[];
// // // // //            }
// // // // //        }
// // // // //     }
// // // // //    */

// // // // //   }
#else // EMBED
  /**
  When dealing with embedded boundaries, we simply ignore the fraction
  occupied by the solid. This is a simple approximation which has the
  advantage of ensuring boundedness of the volume fraction and
  conservation of the total tracer mass (if it is computed also
  ignoring the volume occupied by the solid in partial cells). */
  
  // foreach()
  //   if (cs[] > 0.) {
  //     // c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/Delta; // ulf only used in liquid side, other wise, it is not divergence free
  //     c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/Delta;
  //     scalar t, tc, tflux;
  //     for (t, tc, tflux in tracers, tcl, tfluxl)
	// t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/Delta;
  //   }


  // foreach() {
  //     if(cs[]>0.0){
  //    // if(css_test3_n[]>0.0){
  //           if(topo_mask2[]>=-1) // add 20230103
  //             c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta);

  //           scalar t, tc, tflux,t2;
  //           for (t, tc, tflux in tracers, tcl, tfluxl){
  //             int khaki_flag=0;
  //             if(t.khaki){
  //                 if(topo_mask2[]<=-1){ //becase vof method in paris is different from here.
  //                       t[] += dt*(tflux2[] - tflux2[1] + tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta); 
  //                     //t[] += dt*(tflux2[] - tflux2[1]+ tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
  //                     // t[] += dt*(tflux2[] - tflux2[1])/(cm[]*Delta); 
  //                       khaki_flag = 1;              
  //                 }
  //             }
  //             if(!khaki_flag){
  //                 t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta);
  //               //  t[] += dt*(tflux2[] - tflux2[1]+ tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
  //                 //t[] += dt*(tflux[] - tflux[1])/(cm[]*Delta);
  //             }

  //           }
  //       //add for Temperature by cxb, when vof_phase change from 0 to 1 or 1 to 0 

  //           /*
  //           for (t in {Tlff,Tgff}){
  //             if((c[]-0.5)*(cold[]-0.5)<=0){
  //                 if(t.inverse){
  //                       t[] = Tsat0*(1.0-c[]);
  //                 }else{            
  //                       t[] = Tsat0*c[];
  //                 }
  //             }
  //           }
  //         */
  //     }
  // }
fprintf(stderr,"802 vof\n");
  foreach() {
      if(cs[]>0.0){
     //  if(cs[]>1e-6){
     // if(css_test3_n[]>0.0){
            //if(topo_mask2[]>=-1) // add 20230103
          //     c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/(cm[]/cs[]*Delta);
                // c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta); //before 20230614
         //   if(topo_mask2[]>=-1)
             c[] += dt*(flux[] - flux[1] + cc[]*(ulf.x[1] - ulf.x[]))/(max(y,1e-20)*Delta); // 20230614

            scalar t, tc, tflux;
            scalar tflux2, tfluxg;
            // scalar gf;
            // for (t, tc, tflux, tflux2, tfluxg, gf in tracers, tcl, tfluxl, tfluxl2, tfluxlg, gfl){
              for (t, tc, tflux, tflux2, tfluxg in tracers, tcl, tfluxl, tfluxl2, tfluxlg){
              int khaki_flag=0;
              if(t.inverse){
                  if(topo_mask2_g[]>=-1){ //becase vof method in paris is different from here.
                       // t[] += dt*(tfluxg[] - tfluxg[1] + tc[]*(ugf.x[1] - ugf.x[]))/(cm[]*Delta);  

                       // to be more accurate Tgff shoud be advected by ugf, here using ulf   
                     //    t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]/cs[]*Delta);   
                          // t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta);   //before 20230614

                      //origin
                      // t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/Delta);   
                      //modify
                       t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(max(y,1e-20)*Delta);    //20230614 
                        //  t[] += dt*(tflux2[] - tflux2[1])/(max(y,1e-20)*Delta);     
                  }
              } else{
                if(topo_mask2[]>=-1){ 

        if(fabs(x-tracex)<1e-6 && fabs(y-tracey)<1e-6 && !t.inverse){
           printf("1myvof.h t[]=%g dt=%g tflux[]=%g tflux[1]=%g tc[]=%g ulf.x[]=%g ulf.x[1]=%g c[]=%g\n",t[],dt,tflux[],tflux[1],tc[],ulf.x[],ulf.x[1],c[]);
           printf("1flux[]=%g flux[1]=%g\n",flux[],flux[1]);
          //  printf("gf[-1]=%g,gf[]=%g,gf[1]=%g",gf[-1],gf[],gf[1]);
          //  printf("fine(tflux,0,0)=%g,fine(tflux,0,1)=%g\n",fine(tflux,0,0),fine(tflux,0,1));
        }
        if(fabs(x-1.70996e-5)<1e-6 && fabs(y-2.18262e-5)<1e-6 && !t.inverse){
           printf("2myvof.h t[]=%g dt=%g tflux[]=%g tflux[1]=%g tc[]=%g ulf.x[]=%g ulf.x[1]=%g c[]=%g\n",t[],dt,tflux[],tflux[1],tc[],ulf.x[],ulf.x[1],c[]);
           printf("2flux[]=%g flux[1]=%g\n",flux[],flux[1]);
          //  printf("gf[-1]=%g,gf[]=%g,gf[1]=%g",gf[-1],gf[],gf[1]);
          //  printf("fine(tflux,0,0)=%g,fine(tflux,0,1)=%g\n",fine(tflux,0,0),fine(tflux,0,1));
        }
        if(fabs(x-1.70996e-5)<1e-6 && fabs(y-1.8916e-5)<1e-6 && !t.inverse){
           printf("3myvof.h t[]=%g dt=%g tflux[]=%g tflux[1]=%g tc[]=%g ulf.x[]=%g ulf.x[1]=%g c[]=%g\n",t[],dt,tflux[],tflux[1],tc[],ulf.x[],ulf.x[1],c[]);
           printf("3flux[]=%g flux[1]=%g\n",flux[],flux[1]);
          //  printf("gf[-1]=%g,gf[]=%g,gf[1]=%g",gf[-1],gf[],gf[1]);
          //  printf("fine(tflux,0,0)=%g,fine(tflux,0,1)=%g\n",fine(tflux,0,0),fine(tflux,0,1));
        }
         if(fabs(x-2.00098e-5)<1e-6 && fabs(y-1.60059e-5)<1e-6 && !t.inverse){
           printf("4myvof.h t[]=%g dt=%g tflux[]=%g tflux[1]=%g tc[]=%g ulf.x[]=%g ulf.x[1]=%g c[]=%g\n",t[],dt,tflux[],tflux[1],tc[],ulf.x[],ulf.x[1],c[]);
           printf("4flux[]=%g flux[1]=%g\n",flux[],flux[1]);
          //  printf("gf[-1]=%g,gf[]=%g,gf[1]=%g",gf[-1],gf[],gf[1]);
          //  printf("fine(tflux,0,0)=%g,fine(tflux,0,1)=%g\n",fine(tflux,0,0),fine(tflux,0,1));
        }
         if(fabs(x-2.29199e-5)<1e-6 && fabs(y-1.60059-5)<1e-6 && !t.inverse){
           printf("5myvof.h t[]=%g dt=%g tflux[]=%g tflux[1]=%g tc[]=%g ulf.x[]=%g ulf.x[1]=%g c[]=%g\n",t[],dt,tflux[],tflux[1],tc[],ulf.x[],ulf.x[1],c[]);
           printf("5flux[]=%g flux[1]=%g\n",flux[],flux[1]);
          //  printf("gf[-1]=%g,gf[]=%g,gf[1]=%g",gf[-1],gf[],gf[1]);
          //  printf("fine(tflux,0,0)=%g,fine(tflux,0,1)=%g\n",fine(tflux,0,0),fine(tflux,0,1));
        }

// myvof.h t[]=1.54964e+09 dt=1.14604e-08 tflux[]=0 tflux[1]=3588.36 tc[]=1.54964e+09 ulf.x[]=1.96766e-06 ulf.x[1]=2.3156e-06
// myvof.h t[]=1.54935e+09 dt=1.14604e-08 tflux[]=0 tflux[1]=3284.94 tc[]=1.54964e+09 ulf.x[]=2.46815e-06 ulf.x[1]=2.12021e-06
                       //  t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]/cs[]*Delta);
                        //  t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta); //before 20230614
                       t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(max(y,1e-20)*Delta); //20230614
                      //  t[] += dt*(tflux2[] - tflux2[1] )/(max(y,1e-20)*Delta);
                      //  t[] += dt*(tflux2[] - tflux2[1] + tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
                }
              }

            }
        //add for Temperature by cxb, when vof_phase change from 0 to 1 or 1 to 0 

            /*
            for (t in {Tlff,Tgff}){
              if((c[]-0.5)*(cold[]-0.5)<=0){
                  if(t.inverse){
                        t[] = Tsat0*(1.0-c[]);
                  }else{            
                        t[] = Tsat0*c[];
                  }
              }
            }
          */
      }
  }
#endif // EMBED
//  delete (gfl); free (gfl);

  delete (tfluxl); free (tfluxl);
   delete (tfluxl2); free (tfluxl2);
   delete (tfluxlg); free (tfluxlg);
   

}






















/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */

void vof_advection (scalar * interfaces, int i)
{
  // vof_flag = 1;
  globalii = i;
  for (scalar c in interfaces) {

    /**
    We first define the volume fraction field used to compute the
    divergent term in the one-dimensional advection equation above. We
    follow [Weymouth & Yue, 2010](/src/references.bib#weymouth2010) and use a
    step function which guarantees exact mass conservation for the
    multi-dimensional advection scheme (provided the advection velocity
    field is exactly non-divergent). */

    scalar cc[], * tcl = NULL, * tracers = c.tracers; 
    // scalar cc[], flux_temp[], * tcl = NULL, * tracers = c.tracers; 
       
    for (scalar t in tracers) {
      scalar tc = new scalar;
      tcl = list_append (tcl, tc);
#if TREE
      if (t.refine != vof_concentration_refine) {
	t.refine = t.prolongation = vof_concentration_refine;
	t.restriction = restriction_volume_average;
	t.dirty = true;
	t.c = c;
      }
#endif // TREE
    }
    foreach() {
      cc[] = (c[] > 0.5);
      scalar t, tc;
      for (t, tc in tracers, tcl) {
	if (t.inverse)
	  tc[] = c[] < 0.5 ? t[]/(1. - c[]) : 0.;
	else
	  tc[] = c[] > 0.5 ? t[]/c[] : 0.;
      }
    }
 fprintf(stderr,"even vof 926\n");
    /**
    We then apply the one-dimensional advection scheme along each
    dimension. To try to minimise phase errors, we alternate dimensions
    according to the parity of the iteration index `i`. */

    void (* sweep[dimension]) (scalar, scalar, scalar *);
    // void (* sweep[dimension]) (scalar, scalar, scalar, scalar *);
    int d = 0;
    foreach_dimension(){
      //sweep[d++] = sweep_x;
      sweep[d++] = sweep2_x;
    }
    for (d = 0; d < dimension; d++){
      vof_flag = (i + d) % dimension;
      fprintf(stderr,"sweep2 vof 939\n");
      // printf("dimension=%d\n",dimension);
      sweep[(i + d) % dimension] (c, cc, tcl);
      // sweep[(i + d) % dimension] (c, cc, flux_temp, tcl);
      fprintf(stderr,"sweep2 vof 940\n");
      
    }
    // foreach(){
    //   aiml[] = flux_temp[];
    // }
    delete (tcl), free (tcl);
  }
}


extern face vector usf;
extern double thick_init,L0;

event vof (i++){ 
/*
#ifdef trace_source
  foreach(){
	 ff_old[] = ff[];
  }
#endif
*/
  scalar temp_cold[];
  foreach(){
	temp_cold[] = ff[];
  }
 fprintf(stderr,"even vof 960\n");
  vof_advection (interfaces, i);
  fprintf(stderr,"even vof 962\n");
  boundary({ff});
    //  20220812-output
/*
  char name69[80];
  sprintf(name69,"ulf-check-%g.dat", t);
  FILE * fp69 = fopen(name69,"a");
  fprintf(fp69,"#i j x y ulf uf usf\n");
  foreach_face(x){
          if(fabs(x-(L0-thick_init))<4*Delta){
               fprintf(fp69,"%d %d %g %g %g %g %g %g\n",point.i,point.j,x,y,ulf.x[],uf.x[],usf.x[]);         
          }
  }
  fclose(fp69);
*/
  scalar div2[];
  foreach(){
       div2[] =0.0;
       double ff_temp = ff[];
       foreach_dimension(){
         if(ff_temp>0.0)
          div2[] = div2[] + ulf.x[1] - ulf.x[] ;
       }
  }

/*
  char name70[80];
  sprintf(name70,"ff-check-%g.dat", t);
  FILE * fp70 = fopen(name70,"a");
  fprintf(fp70,"#i j x y ffold ff_after_vof cha\n");
  foreach(){
     //if(fabs(x-(L0-thick_init))<4*Delta)
     if(ff[]<1.0 && ff[]>0)
     fprintf(fp70,"%d %d %g %g %g %g %g %g\n", point.i,point.j,x,y,temp_cold[],ff[],temp_cold[]-ff[],div2[]);
  } 
  fclose(fp70);
*/


}

/**
## References

~~~bib
@Article{lopez2015,
  title = {A VOF numerical study on the electrokinetic effects in the 
           breakup of electrified jets},
  author = {J. M. Lopez-Herrera and A. M. Ganan-Calvo and S. Popinet and
            M. A. Herrada},
  journal = {International Journal of Multiphase Flows},
  pages = {14-22},
  volume = {71},
  year = {2015},
  doi = {doi.org/10.1016/j.ijmultiphaseflow.2014.12.005},
  url = {http://gfs.sf.net/papers/lopez2015.pdf}
}
~~~
*/
