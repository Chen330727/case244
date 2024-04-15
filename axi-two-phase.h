/**
# Two-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h). 

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

// #include "./my-vof.h"
#include "./my-vof-css-test.h"
//include "fractions.h"

scalar ff[], * interfaces = {ff};
scalar ff_oppo[];
scalar f_height[];
double Trhol = 1., mu1 = 0., Trhog = 1., mu2 = 0.;
double Tkg,Tkl,Tcpg,Tcpl,hfg;
double k1,k0;
scalar rhocp[];
double Trhos;
double Tcps;
double Tks;

extern scalar T;
extern scalar css_test3, css_test3_n;
extern scalar fs_solid;
extern face vector fss_test3_n;
extern scalar topo_mask_s;
extern scalar intersect_true;
extern bool flag_rho;
extern bool sigma_in_project;
#define is_solid(vv) vv>0.0
#define is_solid2(vv) vv>0.5
/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

// attribute {
//   bool inverse;
//   bool khaki;
// }

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  
  if (mu1 || mu2)
    mu = new face vector;

  /**
  We add the interface to the default display. */

  display ("draw_vof (c = 'ff');");
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(v) (clamp(v,0.,1.)*(Trhol - Trhog) + Trhog)
#endif
//function rhocp_v is similar to rho() 
#define rhocp_v(v) (clamp(v,0.,1.)*(Trhol*Tcpl - Trhog*Tcpg) + Trhog*Tcpg)
#define energy_v(v,Tl,Tg) (clamp(v,0.,1.)*(Trhol*Tcpl*Tl - Trhog*Tcpg*Tg) + Trhog*Tcpg*Tg)/(clamp(v,0.,1.)*(Trhol*Tcpl - Trhog*Tcpg) + Trhog*Tcpg)

#ifndef mu
# define mu(ff)  (clamp(ff,0.,1.)*(mu1 - mu2) + mu2)
#endif

//function cond_v is similar to mu
#define cond_v(v) (clamp(v,0.,1.)*(Tkl - Tkg) + Tkg)
/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
#else
# define sf ff
#endif

event tracer_advection (i++)
{
  // printf("i=%d t=%g\n",i,t);
  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
   if(1==0){
      #if dimension <= 2
        foreach()
          sf[] = (4.*ff[] + 
            2.*(ff[0,1] + ff[0,-1] + ff[1,0] + ff[-1,0]) +
            ff[-1,-1] + ff[1,-1] + ff[1,1] + ff[-1,1])/16.;
      #else // dimension == 3
        foreach()
          sf[] = (8.*ff[] +
            4.*(ff[-1] + ff[1] + ff[0,1] + f[0,-1] + ff[0,0,1] + ff[0,0,-1]) +
            2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
          ff[0,1,1] + ff[0,1,-1] + ff[0,-1,1] + ff[0,-1,-1] +
          ff[1,1] + ff[1,0,1] + ff[1,-1] + ff[1,0,-1]) +
            ff[1,-1,1] + ff[-1,1,1] + f[-1,1,-1] + ff[1,1,1] +
            ff[1,1,-1] + ff[-1,-1,-1] + ff[1,-1,-1] + ff[-1,-1,1])/64.;
      #endif
   }else{
      foreach(){
        sf[] = ff[];
      }
   }
#endif // !sf

#if TREE
 // sf.prolongation = refine_bilinear;
 // sf.dirty = true; // boundary conditions need to be updated
 ff.refine = ff.prolongation = fraction_refine;
 ff.dirty = true;
#endif
}

scalar true_interface[];

event properties (i++)
{

  // scalar true_interface[];
  foreach(){
      true_interface[]=0;
      bool flag=false;
      // foreach_neighbor(1){
      //       if(intersect_true[]==1){
      //           flag = true;
      //       }
      // }
       if(intersect_true[]==1){
                flag = true;
      }
      if(topo_mask_s[]==0 && flag && (ff[]<1.0 && ff[]>0.0)){
          true_interface[]=1;
      }
  }

  // bool flag_rho=true;
#if EMBED
if(sigma_in_project){
          // viscosity_ghost_fluid1();
        //rho for momemtum equation, in solid rho=0
          foreach(){
              if(flag_rho){
                  if(topo_mask_s[]==0 && (ff[]>0.0 && ff[]<1.0) && (true_interface[]!=1)){
                      rhov[] = cm[]*rho(1.0);
                  }else{
                      rhov[] = cm[]*rho(sf[]); 
                  }
              }else{
                  rhov[] = cm[]*rho(sf[]); 
              }
          }
    }else{
      if(flag_rho){
        foreach_face() {    
                  //alphav and mu for momemtum equation, in solid value=0
                  // alphav and mu is used in momemtum euqation, mu and alpha in solid =0, using css_test3_n and fss_test3_n
                double sf0,sf1;
                if(topo_mask_s[]==0 && (ff[]<1.0 && ff[]>0.0) && (true_interface[]!=1)){
                    sf1 = 1.0;
                }else{
                    sf1 = sf[];
                }
                if(topo_mask_s[-1]==0 && (ff[-1]<1.0 && ff[-1]>0.0) && (true_interface[-1]!=1)){
                    sf0 = 1.0;
                }else{
                    sf0 = sf[-1];
                }

                double fff1 = (sf1 + sf0)/2.;
                alphav.x[] = fm.x[]/rho(fff1);
                if (mu1 || mu2) {
                    face vector muv = mu;
                    //muv.x[] = fss_test3_n.x[]*mu(fff1);
                    muv.x[] = fm.x[]*mu(fff1);
                    // mu_f.x[] = fm.x[]*mu2(fff1);
                  }
        }
      }else{
          foreach_face() {
            //alphav and mu for momemtum equation, in solid value=0
            // alphav and mu is used in momemtum euqation, mu and alpha in solid =0, using css_test3_n and fss_test3_n
          double fff1 = (sf[] + sf[-1])/2.;
            //double ffs = (css_test3[] + css_test3[-1])/2.;
          double ffs2 = (fs_solid[] + fs_solid[-1])/2.;
          //alphav.x[] = fss_test3_n.x[]/rho(fff1);
          alphav.x[] = fm.x[]/rho(fff1);
          if (mu1 || mu2) {
              face vector muv = mu;
              //muv.x[] = fss_test3_n.x[]*mu(fff1);
              muv.x[] = fm.x[]*mu(fff1);
            }
        }
      }
        
    //rho for momemtum equation, in solid rho=0
      foreach(){
          if(flag_rho){
              if(topo_mask_s[]==0 && (ff[]>0.0 && ff[]<1.0) && (true_interface[]!=1)){
                  rhov[] = cm[]*rho(1.0);
              }else{
                  rhov[] = cm[]*rho(sf[]); 
              }
          }else{
              rhov[] = cm[]*rho(sf[]); 
          }
      }
    }
        




foreach(){
    if(css_test3_n[]>=1.0){
    //  	rhov[] = cm[]*rho(sf[]);
      //rhocp[]= cm[]*rhocp_v(sf[]); 

    //	rhocp[]= rhocp_v(ff[]);
      // rhocp[]= rhocp_v(ff[])*max(y,1e-20)*css_test3_n[];  
       rhocp[]= rhocp_v(ff[])*cm[];  
      //T[] = (Tlff[]*Trhol*Tcpl*sf[] + Tgff[]*Trhog*Tcpg*(1.0-sf[]))/(Trhol*Tcpl*sf[] + Trhog*Tcpg*(1.0-sf[]));

      //T[] = energy_v(ff[],Tlff[],Tgff[]);
    }else if(css_test3_n[]<=0.0){
    //      rhov[] = cm[]*Trhos;
      //rhocp[] = cm[]*Trhos*Tcps;

      //rhocp[] = Trhos*Tcps;
    //  rhocp[] = Trhos*Tcps*max(y,1e-20)*(1-css_test3_n[]);
       rhocp[] = Trhos*Tcps*cm[];
      //printf("properties rhocp=%g\n",rhocp[]);
    }else{
      //rhov[] = cm[]*Trhos;
  //		rhov[] = css_test3[]*Trhos + (1.0-css_test3[])*rho(sf[]);
      //rhocp[] = cm[]*Trhos*Tcps;
      //rhocp[] = (1.0-css_test3_n[])*Trhos*Tcps + css_test3_n[]*rhocp_v(ff[]);
      rhocp[] = ((1.0-css_test3_n[])*Trhos*Tcps + css_test3_n[]*rhocp_v(ff[]))*cm[];
  // // fprintf(stderr,"i=%d,error in rhov and rhocp, contain solid and fluid in one cell;css_test3_n=%g,css_test3=%g\n,x=%g,y=%g,z=%g\n",i,css_test3_n[],css_test3[],x,y,z);
      //fflush(stderr);
      //exit(1);
    }
  }
  //rho for momemtum equation, in solid rho=0
  foreach(){
      if(flag_rho){
          if(topo_mask_s[]==0 && (ff[]>0.0 && ff[]<1.0) && (true_interface[]!=1)){
              rhov[] = cm[]*rho(1.0);
          }else{
              rhov[] = cm[]*rho(sf[]); 
          }
      }else{
          rhov[] = cm[]*rho(sf[]); 
      }
  }


#else
  foreach_face() {
    double fff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(fff);  //alphav.x = fm.x* is important for pressure equation
    if (mu1 || mu2) {
      face vector muv = mu;
      muv.x[] = fm.x[]*mu(fff);
    }
  }
  
  foreach()
    rhov[] = cm[]*rho(sf[]);



#endif
// foreach(){
//     if(css_test3[]<=0.0){
//   //  	rhov[] = cm[]*rho(sf[]);
// 		//rhocp[]= cm[]*rhocp_v(sf[]); 

// 		rhocp[]= rhocp_v(sf[]); 
// 		//T[] = (Tlff[]*Trhol*Tcpl*sf[] + Tgff[]*Trhog*Tcpg*(1.0-sf[]))/(Trhol*Tcpl*sf[] + Trhog*Tcpg*(1.0-sf[]));

// 		//T[] = energy_v(ff[],Tlff[],Tgff[]);
// 	}else if(css_test3[]>=1){
//   //      rhov[] = cm[]*Trhos;
// 		//rhocp[] = cm[]*Trhos*Tcps;

// 		rhocp[] = Trhos*Tcps;
// 	}else{
// 		//rhov[] = cm[]*Trhos;
// //		rhov[] = css_test3[]*Trhos + (1.0-css_test3[])*rho(sf[]);
// 		//rhocp[] = cm[]*Trhos*Tcps;
// 		rhocp[] = css_test3[]*Trhos*Tcps + (1.0-css_test3[])*rhocp_v(sf[]);
// // // fprintf(stderr,"i=%d,error in rhov and rhocp, contain solid and fluid in one cell;css_test3_n=%g,css_test3=%g\n,x=%g,y=%g,z=%g\n",i,css_test3_n[],css_test3[],x,y,z);
// 		//fflush(stderr);
// 		//exit(1);
// 	}
//   }
#if TREE  
 // sf.prolongation = fraction_refine;
 // sf.dirty = true; // boundary conditions need to be updated
  ff.refine = ff.prolongation = fraction_refine;
 ff.dirty = true;
#endif
}
