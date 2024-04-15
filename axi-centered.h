/**
# Incompressible Navier--Stokes solver (centered formulation)

We wish to approximate numerically the incompressible,
variable-density Navier--Stokes equations
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} = 0
$$
with the deformation tensor 
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$.

The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2011](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver. If embedded boundaries are used, a different scheme is used
for viscosity. */

#include "run.h"
// #include "timestep.h"
#include "./my-timestep.h"
#include "bcg.h"
#if EMBED
//# include "./my-viscosity-embed.h"
// #include "./my-viscosity-embed-official-update.h"
#include "./newest-viscosity.h"
//#include "./my-viscosity.h"
#else
# include "./my-viscosity.h"
#endif

/**
The primary variables are the centered pressure field $p$ and the
centered velocity field $\mathbf{u}$. The centered vector field
$\mathbf{g}$ will contain pressure gradients and acceleration terms.

We will also need an auxilliary face velocity field $\mathbf{u}_f$ and
the associated centered pressure field $p_f$. */

// extern scalar ps;
// extern face vector usf;

scalar p[];
vector u[], g[];
scalar pf[];
face vector uf[];
extern scalar ff;
extern scalar ps;
#define is_phase(v) (v>=0.5)
extern face vector ulf, ugf, usf, usfg;
extern scalar T,Tlff,Tgff;
extern scalar source_pc,source_pc2;
extern double Tsat00;
extern int level_interface;
extern bool sigma_in_project;
#define EPS 0.0000000001
/**
In the case of variable density, the user will need to define both the
face and centered specific volume fields ($\alpha$ and $\alpha_c$
respectively) i.e. $1/\rho$. If not specified by the user, these
fields are set to one i.e. the density is unity.

Viscosity is set by defining the face dynamic viscosity $\mu$; default
is zero.

The face field $\mathbf{a}$ defines the acceleration term; default is
zero.

The statistics for the (multigrid) solution of the pressure Poisson
problems and implicit viscosity are stored in *mgp*, *mgpf*, *mgu*
respectively. 

If *stokes* is set to *true*, the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$ is omitted. This is a
reference to [Stokes flows](http://en.wikipedia.org/wiki/Stokes_flow)
for which inertia is negligible compared to viscosity. */

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
mgstats mgp, mgpf, mgu;
bool stokes = false;

/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{a}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{a}$, this can be written */

#if EMBED
      #if  AXI

    //  # define neumann_pressure(i) (fm.x[i] ? a.n[i]*fm.n[i]/alpha.n[i] :	\
    //             a.n[i]*rho[]/(cm[] + SEPS))
    # define neumann_pressure(i) (fm.n[i] ? a.n[i]*fm.n[i]/alpha.n[i] :	\
                a.n[i]*rho[]/(cm[] + SEPS))
      #else
    # define neumann_pressure(i) (alpha.n[i] ? a.n[i]*fm.n[i]/alpha.n[i] :	\
                a.n[i]*rho[]/(cm[] + SEPS))
      #endif
#else
# define neumann_pressure(i) (a.n[i]*fm.n[i]/alpha.n[i])
#endif

p[right] = neumann (neumann_pressure(ghost));
p[left]  = neumann (- neumann_pressure(0));

#if AXI
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
p[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
p[top]    = neumann (neumann_pressure(ghost));
p[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
p[front]  = neumann (neumann_pressure(ghost));
p[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

/**
For [embedded boundaries on trees](/src/embed-tree.h), we need to
define the pressure gradient for prolongation of pressure close to
embedded boundaries. */

#if TREE && EMBED
void pressure_embed_gradient (Point point, scalar p, coord * g)
{
  foreach_dimension()
    g->x = rho[]/(cm[] + SEPS)*(a.x[] + a.x[1])/2.;
}
#endif // TREE && EMBED

/**
## Initial conditions */

event defaults (i = 0)
{

  CFL = 0.8;

  /**
  The pressures are never dumped. */

  p.nodump = pf.nodump = false;
  
  /**
  The default density field is set to unity (times the metric). */

  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = fm.x[];
  }

  /**
  On trees, refinement of the face-centered velocity field needs to
  preserve the divergence-free condition. */

#if TREE
  uf.x.refine = refine_face_solenoidal;

  /**
  When using [embedded boundaries](/src/embed.h), the restriction and
  prolongation operators need to take the boundary into account. */

#if EMBED
  uf.x.refine = refine_face;
  ulf.x.refine = refine_face; ///

  foreach_dimension(){
    uf.x.prolongation = refine_embed_face_x;
    ulf.x.prolongation = refine_embed_face_x;
  }
  for (scalar s in {p, pf, u, g}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
  for (scalar s in {p, pf})
    s.embed_gradient = pressure_embed_gradient;


#endif // EMBED
#endif // TREE
}


/**
We had some objects to display by default. */

event default_display (i = 0)
  display ("squares (color = 'u.x', spread = -1);");

/**
After user initialisation, we initialise the face velocity and fluid
properties. */

double dtmax;
void delete_small_fraction(){
    //double small_fraction_limit_T = 1e-12; //EPS;
    double small_fraction_limit_T = 1e-12; //EPS;
    foreach(){
        if(ff[]<small_fraction_limit_T && ff[]>0.0){ //this is just for evaporate, if cooling add >1e-6 situation
           Tlff[]=Tsat00;
           Tgff[]=Tgff[];
           //Tgff[]=Tsat00*(1.0-0.0);
           ff[]=0.0;
        }else if(ff[]>1-small_fraction_limit_T && ff[]<1.0){
           Tlff[]=Tlff[];
          // Tlff[]=Tsat00*1.0;
           Tgff[]=Tsat00;
           ff[]=1.0;
        }
    }
    //double small_fraction_limit_T = EPS; //EPS;

}

event init (i = 0)
{
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0);

  /**
  We update fluid properties. */
 //delete_small_fraction();
  event ("properties");

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");
}

/**
## Time integration

The timestep for this iteration is controlled by the CFL condition,
applied to the face centered velocity field $\mathbf{u}_f$; and the
timing of upcoming events. */
extern double CFL_evap;
double timestep_evap (const face vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL_evap;
  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt = Delta/fabs(u.x[]);
#if EMBED
      assert (fm.x[]);
      dt *= fm.x[];
#else
      dt *= cm[];
#endif
      if (dt < dtmax) dtmax = dt;
    }
  dtmax *= CFL_evap;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}

extern int case_number;
extern bool restart_Tsat,surface_heat;


extern double Tks,Tkl,Tkg;
extern double Trhos,Trhol,Trhog;
extern double Tcps,Tcpl,Tcpg,delta_min;
double Tkslg_stability(double ks,double kl, double kg, double delta_min){
    double result;
    double maxk;
    maxk = max(ks,kl);
    maxk = max(maxk,kg);
    result = 4.0*sq(delta_min)/(4.0*maxk);
    return result;
}

event set_dtmax (i++,last) dtmax = DT;

// event stability (i++,last) {
//   dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
// }

extern int globali,outstep;
extern bool out_flag;
event stability (i++,last) {
  // printf("pid %d: event stability\n",pid());

  // if(t<0.5){
  //     dt_plot = 0.01;
  // }else{
  //   dt_plot=0.1;
  // }
  globali=i;
  out_flag =((globali%outstep)==0);// || (globali>=3050 && globali<=4000)) ; //globali>120;//((globali%outstep)==0);
  //fprintf(stderr,"before: timestep dtmax=%g\n",dtmax);

  if(1==1){
      
           // dtmax = dtnext (stokes ? dtmax : timestep (uf, dtmax));
            // fprintf(stderr,"after: itime=%d t=%g dtmax=%g dt=%g\n",i,t,dtmax,dt);
            // dt = dtnext (stokes ? dtmax : timestep (ulf, dtmax));
            // if(dt>DT){
            //     dt = DT;
            // }
            double dtmax1,dtmax2;
            // if(case_number==3 && (!restart_Tsat) && surface_heat){
            //       dtmax = DT;
            // }
            // if(case_number==2){
            //       dtmax = DT;
            // }
            dtmax1 = timestep (uf, dtmax);
            dtmax2 = timestep (ulf, dtmax);
            dt = dtmax1 < dtmax2 ? (dtmax1 < DT ? dtmax1 : DT) : (dtmax2 < DT ? dtmax2 : DT);
            double dt_T_stability = Tkslg_stability(Tks/(Trhos*Tcps),Tkl/(Trhol*Tcpl),Tkg/(Trhog*Tcpg),delta_min);
            if(dt>1e-20){
              dt = min(dt,dt_T_stability);
            }else{
              fprintf(stderr,"not normal\n");
              dt = dt_T_stability;
            }
            // fprintf(stderr,"dt=%g\n",dt);
            dt = dtnext(dt);
          //  min = a < b ? (a < c ? a : c) : (b < c ? b : c);
            fprintf(stderr,"after: itime=%d t=%g dtmax=%g dtmax1=%g dtmax2=%g DT=%g dt_T_stability=%g dt=%g\n",i,t,dtmax,dtmax1,dtmax2,DT,dt_T_stability,dt);
            
  }else{
       dtmax = DT;
       dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
       fprintf(stderr,"after: itime=%d t=%g dtmax=%g dt=%g\n",i,t,dtmax,dt);
  }

  
  //dt = dtnext (stokes ? dtmax : timestep_evap (ulf, dtmax));
  
  //dtmax = dt;
  //dt = dtnext (stokes ? dtmax : timestep (ulf, dtmax));
  //dt = dtnext(DT);
  //fprintf(stderr,"force timestep=0.0005,dt=%g\n",dt);
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */
event before_vof(i++,last){
  // printf("t=%g dt=%g\n",t,dt);
  // printf("pid %d: event before_vof\n",pid());
    event("contact");
};

event vof (i++,last){
  // printf("pid %d: event vof\n",pid());
};
// event tracer_advection (i++,last);
// event tracer_diffusion (i++,last);
event after_vof(i++,last){
  //  printf("pid %d: event after_vof\n",pid());
     event("contact");
};

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last);
event diffusionT_one (i++,last){
  MPI_Barrier(MPI_COMM_WORLD);
  // printf("pid %d: event diffusionT\n",pid());
}; 
// event extract1 (i++) {
//   printf("pid %d: event extract1\n",pid());
// };
// event pictures (i++) {
//   printf("pid %d: event picture1\n",pid());
// };

/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz [advection scheme](/src/bcg.h) and the
pressure gradient and acceleration terms at time $t$ (stored in vector
$\mathbf{g}$). */

void prediction()
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

   if (u.x.gradient){
   // if (1==0){
   // printf("woshi shabi \n");
    foreach(){
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
          du.x[] = 0.;
        else{
             du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
           
             //for embed_face_gradient_x
           // du.x = u.x.gradient(u.x,0);
           
           
           
            // du.x[] = face_gradient_x (u.x,0);
            // for(scalar a in {du.x}){
            // a[] = (a.third && fs.x[] < 1. && fs.x[] > 0. ?			
            //      embed_face_gradient_x (point, a, 0) :			
            //     (a[] - a[-1])/Delta);
            //  printf("i=%d: du.x[]=%g u[-1]=%g u[]=%g u[1]=%g\n",globali,du.x[],u.x[-1], u.x[], u.x[1]);
            // }
        }
#else
	  du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
#endif
      }
    }// foreach
  }else{ //u.x.graidnet
    foreach(){
        foreach_dimension() {
  #if EMBED
          if (!fs.x[] || !fs.x[1])
              du.x[] = 0.;
          else{
              du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
            //  printf("i=%d du.x[]=%g u[-1]=%g u[]=%g u[1]=%g\n",globali,du.x[],u.x[-1], u.x[], u.x[1]);
          }
  #else
      du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
  #endif
      }//foreach dimension
    }//foreach
  }//u.x.gradient 
  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + (g.x[] + g.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;
    #if dimension > 1
    if (fm.y[i,0] && fm.y[i,1]) {
     // if (fs.y[i,0] && fs.y[i,1]) {
      double fyy = u.y[i] < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
    //  printf("i=%d before-uf.x[]=%g u.y[i,0]=%g dt=%g fyy=%g g.x[]=%g,g.x[-1]=%g\n",globali,uf.x[],u.y[i],dt,fyy,g.x[],g.x[-1]);
      uf.x[] -= dt*u.y[i]*fyy/(2.*Delta);
    //  printf("i=%d uf.x[]=%g u.y[i,0]=%g dt=%g\n",globali,uf.x[],u.y[i],dt);
    }
    #endif
    #if dimension > 2
    if (fm.z[i,0,0] && fm.z[i,0,1]) {
      double fzz = u.z[i] < 0. ? u.x[i,0,1] - u.x[i] : u.x[i] - u.x[i,0,-1];
      uf.x[] -= dt*u.z[i]*fzz/(2.*Delta);
    }
    #endif
    uf.x[] *= fm.x[];
  }

  delete ((scalar *){du});
}

/**
### Advection term

We predict the face velocity field $\mathbf{u}_f$ at time $t+\Delta
t/2$ then project it to make it divergence-free. We can then use it to
compute the velocity advection term, using the standard
Bell-Collela-Glaz advection scheme for each component of the velocity
field. */

#include "./my-bcg2.h"
event advection_term (i++,last)
{

  //  #if EMBED
  //               foreach_dimension(){
  //                 //u.x.gradient = embed_face_gradient_x ;
  //                 u.x.gradient = minmod;// gradients ;
  //               }
  //   #endif
  if (!stokes) {
fprintf(stderr,"mypoisson1 begin\n");
    prediction(); 
    //mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax); 
     // mgpf = project_source (uf, pf, alpha, source_pc, dt/2.0, mgpf.nrelax);
     mgpf = project_source (uf, pf, alpha, source_pc2, dt/2.0, mgpf.nrelax);
      advection ((scalar *){u}, uf, dt, (scalar *){g});
    // advection2 ((scalar *){u}, uf, dt, (scalar *){g});
fprintf(stderr,"mypoisson1 end\n");
  }


 
}

/**
### Viscous term

We first define a function which adds the pressure gradient and
acceleration terms. */

static void correction (double dt)
{
  foreach()
    foreach_dimension()
      u.x[] += dt*g.x[];
}

/**
The viscous term is computed implicitly. We first add the pressure
gradient and acceleration terms, as computed at time $t$, then call
the implicit viscosity solver. We then remove the acceleration and
pressure gradient terms as they will be replaced by their values at
time $t+\Delta t$. */

event viscous_term (i++,last)
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax);
    correction (-dt);
  }

  /**
  We reset the acceleration field (if it is not a constant). */

  if (!is_constant(a.x)) {
    face vector af = a;
    trash ({af});
    foreach_face()
      af.x[] = 0.;
  }
}

/**
### Acceleration term

The acceleration term $\mathbf{a}$ needs careful treatment as many
equilibrium solutions depend on exact balance between the acceleration
term and the pressure gradient: for example Laplace's balance for
surface tension or hydrostatic pressure in the presence of gravity.

To ensure a consistent discretisation, the acceleration term is
defined on faces as are pressure gradients and the centered combined
acceleration and pressure gradient term $\mathbf{g}$ is obtained by
averaging. 

The (provisionary) face velocity field at time $t+\Delta t$ is
obtained by interpolation from the centered velocity field. The
acceleration term is added. */

event acceleration (i++,last)
{
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*a.x[]);
}

/**
## Approximate projection

This function constructs the centered pressure gradient and
acceleration field *g* using the face-centered acceleration field *a*
and the cell-centered pressure field *p*. */

void centered_gradient (scalar p, vector g)
{

  /**
  We first compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[]*a.x[] - alpha.x[]*(p[] - p[-1])/Delta;

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
}

/**
To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). Then compute the centered gradient field *g*. */

event projection (i++,last)
{
//   //mgp = project (uf, p, alpha, dt, mgp.nrelax);
// char name33[80];
//   sprintf(name33,"mass_record5-pid%d.dat",pid());
//   FILE * fp33 = fopen(name33,"w");
// #if EMBED
//   foreach(){
//     //fprintf(fp33,"%g %g %g %g\n",x,y,T[],cs[]);
//     fprintf(fp33,"%g %g %g %g %g\n",x,y,ff[],T[],source_pc[]);
//   }
//   // foreach_face(){
//   //   fprintf(fp33,"%g %g %g %g %g %g %g %d\n",x,y,ff[],T[],fs.x[],fm.x[],source_pc[],topomask[]);
//   // }
// #else
//   foreach(){
//     //fprintf(fp33,"%g %g %g\n",x,y,T[]);
//    // fprintf(fp33,"%g %g %g %g %g %g %d %g %g\n",x,y,ff[],T[],cm[],source_pc[],topo_mask[],phase0Tgrad[],phase1Tgrad[]);
//    fprintf(fp33,"%g %g %g %g\n",x,y,ff[],T[]);
//   }
//   //   foreach_face(){
//   //   fprintf(fp33,"%g %g %g %g\n",x,y,fm.x[],alpha.x[]);
//   // }
// #endif

//   fclose(fp33);
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(pid()==0){
//                   char command1[150];
//                   sprintf(command1, "LC_ALL=C cat mass_record5-pid*.dat > outfacets/mass_record5-%g",t);
//                   system(command1);

//                   char command7[150];
//                   sprintf(command7, "LC_ALL=C rm -rf mass_record5-pid*.dat");
//                   system(command7);
//     }
//   //mgp = project_source (uf, p, alpha, source_pc, dt, mgp.nrelax);
fprintf(stderr,"mypoisson2 begin\n");
  mgp = project_source (uf, p, alpha, source_pc2, dt, mgp.nrelax);

  
  centered_gradient (p, g);

  /**
  We add the gradient field *g* to the centered velocity field. */

  correction (dt);
fprintf(stderr,"mypoisson2 finish\n");
}

event poisson_ps(i++,last){
   // get ps and usf
    // printf("pid %d: event poisson_ps\n",pid());
}


event get_ulf(i++,last){
   
    // printf("pid %d: event get_ulf\n",pid());
   // get ps and usf


  //  char name33[80];
  // sprintf(name33,"mass_record4-pid%d.dat",pid());
  // FILE * fp33 = fopen(name33,"w");
  // foreach(){
  //   fprintf(fp33,"%g %g %g %g %g %g %g\n",x,y,ff[],T[],source_pc[],p[],ps[]);
  // }
  //        foreach_boundary(left)
  //   for (int i = -BGHOSTS; i < 0; i++)
  //     fprintf (fp33, "%g %g %g %g %g %g %g\n", x + (i + 0.5)*Delta, y, ff[i],T[i],source_pc[i],p[i],ps[i]);
  
  // foreach_boundary(right)
  //   for (int i = 1; i <= BGHOSTS; i++)
  //     fprintf (fp33, "%g %g %g %g %g %g %g\n", x + (i - 0.5)*Delta, y,ff[i], T[i],source_pc[i],p[i],ps[i]);
    
  // foreach_boundary(bottom)
  //   for (int i = -BGHOSTS; i < 0; i++)
  //     fprintf (fp33, "%g %g %g %g %g %g %g\n", x , y + (i + 0.5)*Delta,ff[0,i], T[0,i],source_pc[0,i],p[0,i],ps[0,i]);
    
  // foreach_boundary(top)
  //   for (int i = 1; i <= BGHOSTS; i++)
  //     fprintf (fp33, "%g %g %g %g %g %g %g\n", x, y + (i - 0.5)*Delta, ff[0,i],T[0,i],source_pc[0,i],p[0,i],ps[0,i]);


  // fclose(fp33);
  //  MPI_Barrier(MPI_COMM_WORLD);
  //  if(pid()==0){
  //                 char command1[150];
  //                 sprintf(command1, "LC_ALL=C cat mass_record4-pid*.dat > outfacets/mass_record4-%g",t);
  //                 system(command1);

  //                 char command7[150];
  //                 sprintf(command7, "LC_ALL=C rm -rf mass_record4-pid*.dat");
  //                 system(command7);
  //   }
}

/**
Some derived solvers need to hook themselves at the end of the
timestep. */

event end_timestep (i++, last);

/**
## Adaptivity

After mesh adaptation fluid properties need to be updated. When using
[embedded boundaries](/src/embed.h) the fluid fractions and face
fluxes need to be checked for inconsistencies. */

#if TREE

event update_Tl_Tg(i++,last){
   
  //  printf("pid %d: update_Tl_Tg\n",pid());
}
event adapt (i++,last) {
#if EMBED
  fractions_cleanup (cs, fs);

  //
   boundary({cs,fs});
               #if AXI     
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});
                  
              #endif
  foreach_face()
    if ((uf.x[] || ulf.x[] || usf.x[]) && !fs.x[]){ ////??????? velocity jump in triple point
      uf.x[] = 0.;
      ulf.x[]=0.0;
      usf.x[] = 0.0;
      ugf.x[] = 0.0;
      usfg.x[] = 0.0;
    }

  boundary({uf,ulf,usf});
#endif
  event ("properties");
}
#endif

/**
## See also

* [Double projection](double-projection.h)
* [Performance monitoring](perfs.h)
*/
