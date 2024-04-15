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
#include "./myc2d2.h"
extern scalar css_test;
extern scalar modify_near_region;//, cc_css_test;
extern bool energy_advecting_flag;

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
// extern scalar flux_x,flux_y;
extern vector flux_show;

// extern scalar aiml;
int vof_flag;
int globalii;
// extern scalar flux_out;

//extern scalar css_test3_n,cs;


extern face vector uf;  //fir advect if tracer gas side.
extern double tracex,tracey;
extern double delta_min;
extern scalar topo_mask_s;

#define EPS3 0.0000000001
//#define CCMIN CCCCMIN 
#define CCMIN CCCCMIN 
// in mycenter.h 0.4

scalar cc_css_test[];
// scalar modify_near_region[];

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
    // printf("vof EMBED=true\n");
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
  scalar alpha[], flux[],flux_css_test[];

  vector solid_n[];
  scalar solid_neg_alpha[];

  vector  css_test_n[];
  scalar css_test_alpha[];
  double cfl = 0.;
   
  scalar cs_neg[];
  foreach(){
    cs_neg[] = 1.0-cs[];
  }
//add for Temperature by cxb, when vof_phase change from 0 to 1 or 1 to 0 

/*
  scalar cold[];
  foreach(){
    cold[] = c[];
  }
*/


// fprintf(stderr,"sweep2 vof 434\n");
  /**
  If we are also transporting tracers associated with $c$, we need to
  compute their gradient i.e. $\partial_xf_j = \partial_x(t_j/c)$ or
  $\partial_xf_j = \partial_x(t_j/(1 - c))$ (for higher-order
  upwinding) and we need to store the computed fluxes. We first
  allocate the corresponding lists. */

  scalar flux2[];//tflux2[];

  scalar fluxg[];//tfluxg[];

  // scalar flux_temp[];

  scalar * tracers = c.tracers, * gfl = NULL, * tfluxl = NULL, * tfluxl_css_test = NULL;
  scalar * gfl2 = NULL, * gflg = NULL, * tfluxl2 = NULL, * tfluxlg = NULL;

  // fprintf(stderr,"sweep2 tracers 449\n");
  if (tracers) {
    for (scalar t in tracers) {
      // fprintf(stderr,"sweep2 tracers 450\n");
      scalar gf = new scalar, flux = new scalar;
      scalar gf2 = new scalar, flux2 = new scalar;
      // fprintf(stderr,"sweep2 tracers 451\n");
      scalar gfg = new scalar, fluxg = new scalar;
      // fprintf(stderr,"sweep2 tracers 452\n");
      scalar flux_css_test = new scalar;
      //  fprintf(stderr,"sweep2 tracers 453\n");
      gfl = list_append (gfl, gf);
      // fprintf(stderr,"sweep2 tracers 454\n");
      tfluxl = list_append (tfluxl, flux);
      // fprintf(stderr,"sweep2 tracers 455\n");
      gfl2 = list_append (gfl2, gf2);
      tfluxl2 = list_append (tfluxl2, flux2);
      // fprintf(stderr,"sweep2 tracers 456\n");
      gflg = list_append (gflg, gfg);
      tfluxlg = list_append (tfluxlg, fluxg);
      //  fprintf(stderr,"sweep2 tracers 462\n");
        tfluxl_css_test = list_append (tfluxl_css_test, flux_css_test);
    }

    /**
    The gradient is computed using the "interface-biased" scheme above. */
// fprintf(stderr,"sweep2 tracers 467\n");
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
    // fprintf(stderr,"sweep2 tracers 474\n");
  }
  boundary(gfl);
  boundary(gfl2);
  boundary(gflg);
  /**
  We reconstruct the interface normal $\mathbf{n}$ and the intercept
  $\alpha$ for each cell. Then we go through each (vertical) face of
  the grid. */

// fprintf(stderr,"sweep2 vof 479\n");

if(1==1){
  reconstruction (c, n, alpha);
}else{
  reconstruction4 (c, n, alpha);
}
  reconstruction (cs_neg, solid_n, solid_neg_alpha);

   scalar topo_mask2[];
   get_topomask2(topo_mask2);

   scalar topo_mask2_g[];
   foreach(){
      topo_mask2_g[] = -topo_mask2[];
   }
// fprintf(stderr,"sweep2 vof 493\n");


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


double cf_css_test=0.0;
if(modify_near_region[]>0 && cs[]>0.0){
      if(cs[i]>=1.0){
          cf_css_test = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
            rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
              (coord){-0.5, -0.5, -0.5},
              (coord){s*un - 0.5, 0.5, 0.5});
          // set flux limit between big cell to small cell; then there will be volume inconservative.
      }else if((cs[i]>0.0 && cs[i]<1.0) && c[i]>=1.0){
        // }else if((cs[i]>0.0 && cs[i]<1.0) && css_test[i]>=cs[i]){
          double moving_gas=0.0;
          double moving_liquid=0.0;
          coord nsf = (coord){solid_n.x[i],solid_n.y[i]}; //css_test3
          double alphasf = solid_neg_alpha[i];

            coord face_n = (coord){sign(un),0};
            if(fabs(un)>1e-30){
                // nlg.x = n.x[i];
                // nlg.y = n.y[i];
                // origin: left bottom
                // un is + :1-un  i=-1 -->  fabs(i)-un  x=1-un x-(1-un)=0 x-1+un=0  x-fabs(i)+fabs(un)=0 a=sign(un) c=-fabs(i)+fabs(un)
                // un is - : -un   i=0  -->  fabs(i)-un x=-un x+un -x-un=0 -x-un=0  -x-fabs(i)+fabs(un)=0
                // double alphalg = alpha[i];
                // double alphalg = fabs(i)-un;
                double face_alpha = -fabs(i)+fabs(un); // ax+b+c=0 if a<0,then c=<0
                double a1,b1,c1,a2,b2,c2;
                a1 = nsf.x,b1=nsf.y,c1=-0.5*(nsf.x+nsf.y)-alphasf;
                // a2 = nlg.x,b2=nlg.y,c2=-0.5*(nlg.x+nlg.y)-alphalg; 
                a2 = face_n.x,b2=face_n.y;//c2=-0.5*(face_n.x+face_n.y)-face_alpha;
                c2=face_alpha; 
                        //a1*x+b1*y+c1=0,a2*x+b2*y+c2=0;
                // printf("a1=%g,b1=%g,c1=%g,a2=%g,b2=%g,c2=%g\n",a1,b1,c1,a2,b2,c2);
                double data2[13];
                // cut_line_test_in_basilisk(a1,b1,c1,a2,b2,c2,data2);
                cut_line_test_in_basilisk((coord){a1,b1},c1,(coord){a2,b2},c2,data2);
                moving_liquid = data2[0];
                // cf = gas_volume/cs[i];
                // cf_css_test = moving_liquid/fabs(un);
                
                if(fs.x[]>0.0){
                    cf_css_test = moving_liquid/(fabs(un)*fs.x[]);
                    if(cf_css_test>1){
                        cf_css_test=1.0;
                    }
                }else{
                    cf_css_test=0.0;
                }
                // printf("1111 cf=%g moving_liquid=%g un=%g f[]=%g css_test[]=%g\n",cf_css_test,moving_liquid,un,c[],css_test[]);
          // cf=1.0;
        }else{
            cf_css_test = 0.0;
        }
        // cf = 1.0;
        }else if((cs[i]<1.0 && cs[i]>0.0) && (c[i]<1.0 && c[i]>0.0)){
          // printf("22222222\n");
      // }else if((cs[i]<1.0 && cs[i]>0.0) && (css_test[i]<cs[i] && css_test[i]>0.0)){
          double gas_volume=0.0;
          double liquid_volume=0.0;
          coord nsf = (coord){solid_n.x[i],solid_n.y[i]};
          double alphasf = solid_neg_alpha[i];

          coord nlg = (coord){n.x[i],n.y[i]};
          
          // coord nlg = (coord){css_test_n.x[i],css_test_n.y[i]};
          double alphalg = alpha[i];

          // double alphalg = css_test_alpha[i];

          coord face_n = (coord){sign(un),0};
          if(fabs(un)>1e-30){
              double face_alpha = -fabs(i)+fabs(un);
              double a1,b1,c1,a2,b2,c2,a3,b3,c3;
              a1 = nsf.x,b1=nsf.y,c1=-0.5*(nsf.x+nsf.y)-alphasf;
              a2 = nlg.x,b2=nlg.y,c2=-0.5*(nlg.x+nlg.y)-alphalg; 
              a3 = face_n.x,b3=face_n.y; //c3=-0.5*(face_n.x+face_n.y)-face_alpha; 
              c3=face_alpha;
                      //a1*x+b1*y+c1=0,a2*x+b2*y+c2=0;
              // printf("a1=%g,b1=%g,c1=%g,a2=%g,b2=%g,c2=%g,a3=%g,b3=%g,c3=%g\n",a1,b1,c1,a2,b2,c2,a3,b3,c3);
              cut_line_test_in_basilisk_3_times((coord){a1,b1},c1,(coord){a2,b2},c2,(coord){a3,b3},c3,&liquid_volume,&gas_volume);
              // cf = liquid_volume/(liquid_volume+gas_volume);
              // cf_css_test = liquid_volume/fabs(un);
              if(fs.x[]>0.0){
                    cf_css_test = liquid_volume/(fabs(un)*fs.x[]);
                    if(cf_css_test>1){
                        cf_css_test=1.0;
                    }
                }else{
                    cf_css_test=0.0;
                }
              // printf("2222 cf=%g liquid_volume=%g un=%g f[]=%g css_test[]=%g\n",cf_css_test,liquid_volume,un,c[],css_test[]);
          }else{
              cf_css_test = 0;//1;
          }
      }else{
          cf_css_test = 0.0;
      }

      //limit
          if(((cs[]>=1 && cs[i]<1) || (cs[i]>=1 && cs[]<1)) && fabs(un)>0.0){
                double flux1 = cf_css_test*(fabs(un)*fs.x[]);
                // double limit = min(cs[],cs[i])*0.8;
                double limit = min(cs[],cs[i])*1.0;
                if(flux1>limit && fs.x[]>0.0){
                    cf_css_test = limit/(fabs(un)*fs.x[]);
                }
          }
  }


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

    flux_css_test[] = cf_css_test*ulf.x[];
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
              // // if(fabs(x-tracex+delta_min/2.0*2)<1e-6 && (fabs(y-tracey-delta_min/2.0)<1e-6 || fabs(y-tracey+delta_min/2.0)<1e-6)){
              // //     printf("x face: tflux=%g fff=%g cf1=%g ulf.x[]=%g gf[i]=%g intx=%g inty=%g\n",tflux[],fff,cf1,ulf.x[],gf[i],x/delta_min,y/delta_min);
              // // }
              // // // if(fabs(x-tracex)<1e-6 && fabs(y-tracey+delta_min/2.0*2)<1e-6){
              // //   if((fabs(x-tracex-delta_min/2.0)<1e-6 || fabs(x-tracex+delta_min/2.0)<1e-6) && fabs(y-tracey+delta_min/2.0*2)<1e-6){
              // //     printf("y face: tflux=%g fff=%g cf1=%g ulf.x[]=%g gf[i]=%g  intx=%g inty=%g\n",tflux[],fff,cf1,ulf.x[],gf[i],x/delta_min,y/delta_min);
              // // }

              // // if(fabs(tflux[]>1e16)){
              // //   printf("tflux too big, x=%g y=%g  tflux=%g fff=%g cf1=%g ulf.x[]=%g t[i]=%g gf[i]=%g\n",x,y,tflux[],fff,cf1,ulf.x[],t[i],gf[i]);
              // // }
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

  // fprintf(stderr,"sweep2 vof 865\n");
  #if TREE
  scalar * fluxl = list_concat (NULL, tfluxl);
  fluxl = list_concat (fluxl, tfluxl2);
  fluxl = list_concat (fluxl, tfluxlg);
  fluxl = list_append (fluxl, flux);
  fluxl = list_append (fluxl, flux_css_test);
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
// boundary_flux(tfluxl_css_test);
// boundary_flux(tfluxl);

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


  // fprintf(stderr,"sweep2 vof 635\n");
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
                       //  t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]/cs[]*Delta);
                        //  t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(cm[]*Delta); //before 20230614
                       t[] += dt*(tflux[] - tflux[1] + tc[]*(ulf.x[1] - ulf.x[]))/(max(y,1e-20)*Delta); //20230614
                      //  t[] += dt*(tflux2[] - tflux2[1] )/(max(y,1e-20)*Delta);
                      //  t[] += dt*(tflux2[] - tflux2[1] + tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
                }
              }

            }
      }
  }


  // scalar modify_near_region[];
  foreach(){
     if(modify_near_region[]>0 && cs[]>0.0){
        if (cs[] >= 1.0) { 
          // c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/Delta;
          css_test[] += dt*(flux_css_test[] - flux_css_test[1] + (cc_css_test[])*(ulf.x[1] - ulf.x[]))/(max(y,1e-20)*Delta);
        }else if(cs[]<1.0 && cs[]>0.0){
          css_test[] += dt*(flux_css_test[] - flux_css_test[1] + (cc_css_test[])*(ulf.x[1] - ulf.x[]))/(max(y,1e-20)*Delta);
        }
     }
      // }
  }


#endif // EMBED
//  delete (gfl); free (gfl);

  delete (tfluxl); free (tfluxl);
   delete (tfluxl2); free (tfluxl2);
   delete (tfluxlg); free (tfluxlg);
   delete (tfluxl_css_test); free (tfluxl_css_test);
   

}






















/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */

void vof_advection (scalar * interfaces, int i)
{
  // vof_flag = 1;
  globalii = i;
  for (scalar c in interfaces) {

foreach(){
    modify_near_region[] = 0.0;
    if(level==level_interface){
      bool not_near_triple=true;
      bool c_me = (c[]>=0.5);
      // double aa = css_test[];
        foreach_neighbor(1){
          bool c_nei = (c[]>=0.5);
          // double bb = css_test[];
          // if(((cs[]<1.0 && cs[]>0.0)) && ((c[]<1.0 && c[]>0.0) || (c_me ^ c_nei))){
          // bool flag =( (aa>=cs[] && bb<=0) || (aa<=0 && bb>=cs[]));
          // if( (cs[]<1.0 && cs[]>0.0) && ((c[]<1.0 && c[]>0.0) || flag )){
            if(((cs[]<1.0 && cs[]>0.0)) && (c[]<1.0 && c[]>0.0) && (level==level_interface)){
          // if(((cs[]<1.0 && cs[]>0.0)) && ((c[]<1.0 && c[]>0.0) || (css_test[]>0 && css_test[]<cs[]))){
          // if(cs[]<1.0 && cs[]>0.0 && ((c[]<1.0 && c[]>0.0) ||  (css_test[]>0 && css_test[]<cs[]))){
              not_near_triple = false;
          }
        }
        // if(((cs[]<1.0 && cs[]>0.0)) && (c[]<1.0 && c[]>0.0)){
        //       modify_near_region[] = 1;
        // }
        if(!not_near_triple){
            modify_near_region[] = 1;
        }
    }
  }

  foreach(){
    if((level==level_interface) && (modify_near_region[]==0)){
        bool find_flag=false;
        Point me = point;
        foreach_neighbor(2){
          if(level==level_interface){
            if(!((me.i==point.i) && (me.j==point.j)) && (!find_flag)){
                if(modify_near_region[]==1){
                    find_flag = true;
                }
            }
          }
        }
        if(find_flag){
            modify_near_region[] = 2;
        }
    }
  }

  foreach(){
    if((level==level_interface) && (modify_near_region[]==0)){
        bool find_flag=false;
        Point me = point;
        foreach_neighbor(2){
          if(level==level_interface){
            if(!((me.i==point.i) && (me.j==point.j)) && (!find_flag)){
                if(modify_near_region[]==2){
                    find_flag = true;
                }
            }
          }
        }
        if(find_flag){
            modify_near_region[] = 3;
        }
    }
  }
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
//  fprintf(stderr,"even vof 926\n");
    /**
    We then apply the one-dimensional advection scheme along each
    dimension. To try to minimise phase errors, we alternate dimensions
    according to the parity of the iteration index `i`. */

    // foreach(){
    //     if(css_test[]>cs[]){
    //         css_test[]=cs[];
    //     }else if(css_test[]<0.0){
    //         css_test[] = 0.0;
    //     }
    //   }

     foreach(){
      // cc_css_test[] = (css_test[] > (1-css_test[]));
      cc_css_test[] = (css_test[] > (cs[]-css_test[]));
      if(cs[]<1.0 && cs[]>0.0 && ff[]<1.0 && ff[]>0.0){
          cc_css_test[] = 0;
      }
      // cc_css_test[] = css_test[]>0.0;
    }

    void (* sweep[dimension]) (scalar, scalar, scalar *);
    // void (* sweep[dimension]) (scalar, scalar, scalar, scalar *);
    int d = 0;
    foreach_dimension(){
      //sweep[d++] = sweep_x;
      sweep[d++] = sweep2_x;
    }
    for (d = 0; d < dimension; d++){
      vof_flag = (i + d) % dimension;
      // fprintf(stderr,"sweep2 vof 939\n");
      // printf("dimension=%d\n",dimension);
      sweep[(i + d) % dimension] (c, cc, tcl);
      // sweep[(i + d) % dimension] (c, cc, flux_temp, tcl);
      // fprintf(stderr,"sweep2 vof 940\n");
      
    }
    // foreach(){
    //   aiml[] = flux_temp[];
    // }
    delete (tcl), free (tcl);


      // foreach(){
      //   if(css_test[]>cs[]){
      //       css_test[]=cs[];
      //   }else if(css_test[]<0.0){
      //       css_test[] = 0.0;
      //   }
      // }

    vector n[];
    scalar alpha[];
    vector solid_n[];
    scalar solid_neg_alpha[];
    scalar cs_neg[];
    foreach(){
      cs_neg[] = 1.0-cs[];
    }
    reconstruction(cs_neg,solid_n,solid_neg_alpha);
    // reconstruction_contact(c,n,alpha);
    reconstruction(c,n,alpha);
  foreach(){
      // bool not_near_triple=true;
      // foreach_neighbor(1){
      //   if(cs[]<1.0 && cs[]>0.0 && c[]<1.0 && c[]>0.0){
      //       not_near_triple = false;
      //   }
      // }
      if(modify_near_region[]==1){ //change f based on css_test
      double c_old2 = c[];
          if(cs[]>=1){
              c[] = css_test[];
          }else if(cs[]<=0){
              c[] = c[];
          }else{ //0<cs<1
              if(css_test[]>=cs[]){
                    c[] = 1.0;
              }else if(css_test[]<=0){
                    c[] = 0.0;
              }else{ //actrually move in one direction allow volume source
                    double a1,b1,c1,a2,b2,c2;
                    double gas_volume=0.0;
                    double liquid_volume=0.0;
                    coord nsf = (coord){solid_n.x[],solid_n.y[]};
                    double alphasf = solid_neg_alpha[];
                    a1 = nsf.x,b1=nsf.y,c1=-0.5*(nsf.x+nsf.y)-alphasf;
                    //  printf("nsf.x=%g nsf.y=%g alphasf=%g c1=%g\n",a1,b1,alphasf,c1);

                    coord nlg = (coord){n.x[],n.y[]};
                    // coord nlg = (coord){css_test_n.x[],css_test_n.y[]};
                    // double alphalg=alpha[];
                    double alphalg;
                    a2 = nlg.x,b2=nlg.y; 

                //  if((fabs(a2)<1e-20 && fabs(b2)<1e-20)){
                //           //       coord ns = mycs (point, cs);
                //           coord nf;
                //           double distance = HUGE;
                //           double aimx,aimy;
                //           aimx = x;
                //           aimy = y;
                //           coord nf_temp;
                //           foreach_neighbor(1){
                //               if((!(fabs(aimx-x)<1e-29 && fabs(aimy-y)<1e-29)) && (c[]<1.0 && c[]>0.0 && cs[]>=1.0)){
                //                     double distance_temp = sqrt(sq(aimx-x)+sq(aimy-y));
                //                     if(distance > distance_temp){
                //                         distance = distance_temp;
                //                         nf_temp.x = n.x[];
                //                         nf_temp.y = n.y[];
                //                     }
                //               }
                //           }
                //           if(distance<HUGE/2.0){
                //               nf.x = nf_temp.x;
                //               nf.y = nf_temp.y;
                //               coord nlg2 = (coord){nf.x,nf.y};
                //               a2= nlg2.x,b2=nlg2.y;
                //           }
                //           fprintf(stderr,"gougou a2=%g b2=%g\n",a2,b2);
          
                //  }
                if(!(fabs(a2)<1e-20 && fabs(b2)<1e-20)){
                    double c_min,c_max;
                    line_intersection_range(nlg, &c_min, &c_max);
                    // printf("c_min=%g,c_max=%g\n",c_min,c_max);
                    c_min = c_min+1e-10;
                    c_max = c_max-1e-10;
                    int iter = 0;
                    double c_mid, area_mid;
                    double max_iter =20;
                    double c_mid_alpha;
                    double target_area = css_test[];
                    bool find_target=false;
                    double tol=1e-8;
                    double data3[13];
                    // cut_line_test_in_basilisk((coord){a1,b1},c1,(coord){a2,b2},-0.5*(nlg.x+nlg.y)-c_min,data3);
                    cut_line_test_in_basilisk((coord){a1,b1},c1,(coord){a2,b2},c_min,data3);
                    double area_cmin = data3[1];
                    double data4[13];
                    // cut_line_test_in_basilisk((coord){a1,b1},c1,(coord){a2,b2},-0.5*(nlg.x+nlg.y)-c_max,data4);
                    cut_line_test_in_basilisk((coord){a1,b1},c1,(coord){a2,b2},c_max,data4);
                    double area_cmax = data4[1];
                    bool min_max_increase=true;
                    min_max_increase = area_cmax>area_cmin;
                    while(iter < max_iter && (!find_target)) {
                        // printf("iter=%d !!!!!!!!!!!!!!!!!!\n",iter);
                        c_mid = (c_min + c_max) / 2.0;
                        // c_mid=-0.5*(nlg.x+nlg.y)-c_mid_temp; //c2
                        // area_mid = area_func(a, b, c_mid);
                        double data2[13];
                        // printf("iter=%d 00000!!!!!!!!!!!!!!!!!!\n",iter);
                        cut_line_test_in_basilisk((coord){a1,b1},c1,(coord){a2,b2},c_mid,data2);
                        area_mid = data2[1];
                        if (fabs(area_mid - target_area) < tol) {
                            // return c_mid; // Found c satisfying the condition
                            find_target = true;
                        }
                        // printf("iter=%d 111111!!!!!!!!!!!!!!!!!!\n",iter);
                        if ((area_mid < target_area && min_max_increase) || (area_mid>target_area && (!min_max_increase))) {
                            c_min = c_mid;
                        } else {
                            c_max = c_mid;
                        }
                        
                        iter++;
                        // fprintf(stderr,"vof iter=%d\n",iter);
                        // if(i%10==0)
                        //   fprintf(stderr,"vof iter=%d target=%g get=%g\n",iter,target_area,area_mid);
                    }
                    
                    // if(find_target){
                      c_mid_alpha = -0.5*(a2+b2) - c_mid; //alpha after transfer
                      // c1=-0.5*(nsf.x+nsf.y)-alphasf;
                     
        // ax + by + c =0 
        // ax1 + by1 = alpha
        // x1 = x - 0.5
        //       a(x1+0.5) + b(y1+0.5) + c =0 ----> ax1+by1= -c-0.5(a+b) = alpha 

        //ax1 + by1 = alpha
        // x = x1 + 0.5 ----->  a(x-0.5)+b(y-0.5)=alpha ax+by-0.5(a+b)-alpha=0     
                        double c_old = c[];
                        c[] = plane_volume (nlg, c_mid_alpha);

                       

                        if(c[]==0){
                            printf("target=%g result=%g error=%g cs[]=%g a2=%g b2 =%g alpha[]=%g c_mid_alpha=%g area_cmin=%g area_cmax=%g c_old=%g c[]===%g c_mid=%g\n",
                            css_test[],area_mid,css_test[]-area_mid,cs[],a2,b2,alpha[],c_mid_alpha,area_cmin,area_cmax,c_old,c[],c_mid);
                            printf("nsf.x=%g nsf.y=%g alphasf=%g c1=%g\n",a1,b1,alphasf,c1);
                            printf("lelele\n");
                        }

                  }else{
                        // c[] = css_test[];
                        printf("xiaolaji\n");
                  }
                        
                    // }
              }
          }
         //energy Tlff and Tgff is based on ff[], so when ff changed, these value should also be changed
         if(cs[]>0.0){
                        for(scalar tt in tracers){
                           double T_temp = Tsat00;
                           double value_l,value_g;
                          //  bool energy_advecting_flag = true;
                              if(energy_advecting_flag){
                                  value_l = Trhol*Tcpl;
                                  value_g = Trhog*Tcpg;
                              }else{
                                  value_l = 1.0;
                                  value_g = 1.0;
                              }
                           if(!tt.inverse){ //liquid
                              if(c[]>0.0){
                                 if(c_old2>0.0){
                                    T_temp = tt[]/(c_old2*value_l);
                                 }else{
                                    T_temp = Tsat00;
                                 }
                                  tt[] = c[]*T_temp*value_l;
                              }
                           }else{ //gas
                              if(1.0-c[]>0.0){
                                 if(1.0-c_old2>0.0){
                                    T_temp = tt[]/((1.0-c_old2)*value_g);
                                 }else{
                                    T_temp = Tsat00;
                                 }
                                   tt[] = (1.0-c[])*T_temp*value_g;
                              }
                           }
                        }
         }

     }else{ //far away modified css_test based on ff
            if(cs[]>=1){
                css_test[]=c[];
            }else if(cs[]<=0){
                css_test[]=0.0;
            }else{
                if(c[]>=1){
                    css_test[] = cs[];
                }else if(c[]<=0){
                    css_test[] = 0.0;
                }else{
                    coord nsf = (coord){solid_n.x[],solid_n.y[]};
                    double alphasf = solid_neg_alpha[];
                    coord nlg = (coord){n.x[],n.y[]};
                    double alphalg=alpha[];
                    double a1,b1,c1,a2,b2,c2;
                    a1 = nsf.x,b1=nsf.y,c1=-0.5*(nsf.x+nsf.y)-alphasf;
                    a2 = nlg.x,b2=nlg.y,c2=-0.5*(nlg.x+nlg.y)-alphalg; 
                    // printf("a1=%g,b1=%g,c1=%g,a2=%g,b2=%g,c2=%g\n",a1,b1,c1,a2,b2,c2);
                    double data2[13];
                    cut_line_test_in_basilisk((coord){a1,b1},c1,(coord){a2,b2},c2,data2);
                    css_test[]  = data2[1]; //=css_test_value;
                }
            }
     }
  }// foreach

    foreach(){
      if(c[]>0.0 && c[]<1e-8){
            Tlff[]=0.0;
            Tgff[]=(1.0-0.0)/(1.0-c[])*Tgff[];
            c[]=0.0;
      }
      if(c[]<1.0 && c[]>1-1e-8){
            Tgff[]=0.0;
            Tlff[]=(1.0)/c[]*Tlff[];
            c[]=1.0;
      }
    }
// event ("contact");

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
//  fprintf(stderr,"even vof 960\n");
  vof_advection (interfaces, i);
  // fprintf(stderr,"even vof 962\n");
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
