/**
# Multigrid Poisson--Helmholtz solvers

We want to solve Poisson--Helmholtz equations of the general form
$$
L(a) = \nabla\cdot (\alpha\nabla a) + \lambda a = b
$$
This can be done efficiently using a multigrid solver. 

An important aspect of Poisson--Helmholtz equations is that the
operator $L()$ is linear. This property can be used to build better
estimates of a solution by successive *corrections* to an initial
guess. If we define an approximate solution $\tilde{a}$ as
$$
\tilde{a} + da = a
$$
where $a$ is the exact (unknown) solution, using the linearity of the
operator we find that $da$ verifies
$$
L(da) = b - L(\tilde{a})
$$
where the right-hand-side is often called the *residual* of the
approximate solution $\tilde{a}$.

## Multigrid cycle

Here we implement the multigrid cycle proper. Given an initial guess
*a*, a residual *res*, a correction field *da* and a relaxation
function *relax*, we will provide an improved guess at the end of the
cycle. */
extern double k0,k1; //k0=k/rho_cp (gas);
extern double k01; //depend on phase_flag
extern face vector modphase0;
extern face vector modphase1;
extern int phase_flag3;
extern double Tsat00;
extern scalar css_test;
extern face vector fss_test;
scalar tempa[];
extern int level_interface;
extern int minl;
extern double total_area_diff,flux_total_diff;
extern double t;
extern int globali;
extern scalar ff,rhocp;
bool flag_v0_nodata;

bool embed3_flux_flag = false; //false ; // false; //true;

scalar res2[];
bool res2_flag = false; //true; //false; //false; //true;

bool relax_embed2 = true;

extern int choice_T_adap; // embed_flux2 different methods for choice_T_adap
extern double Tkl,Tkg;
extern double Trhog,Tcpg;


extern bool new_flag;
#define CCCCCMIN CCCCMIN

#include "gradient_poisson.h" //for face_gradient
#include "fractions.h"

//scalar res_T[];
double maxres_T;
//double lim_cut = 1e-12;
// double point_to_line(coord nn, double alpha_n, coord pp){
//     if(sqrt(nn.x*nn.x+nn.y*nn.y)<lim_cut){
//          return 1.0/sqrt(2.0);
//      }

//      double a = nn.x/sqrt(nn.x*nn.x+nn.y*nn.y), b=nn.y/sqrt(nn.x*nn.x+nn.y*nn.y);
//      double val1=fabs(a*pp.x+b*pp.y-alpha_n)/fabs(sqrt(a*a+b*b));
    
//     if(val1<=1.0/sqrt(2.0))
//       return val1;
//     else if(val1>1.0/sqrt(2.0)){
//       return 1.0/sqrt(2.0);
//     }
//      //return 1.0;
// }



#define face_gradient_no_cs_x(a,i) ((a[i] - a[i-1])/Delta)
#define face_gradient_no_cs_y(a,i) ((a[0,i] - a[0,i-1])/Delta)
#define face_gradient_no_cs_z(a,i) ((a[0,0,i] - a[0,0,i-1])/Delta)
#define face_value_no_cs(a,i)      ((a[i] + a[i-1])/2.)
#define center_gradient_no_cs(a)   ((a[1] - a[-1])/(2.*Delta))



static inline double bilinear_no_cs (Point point, scalar s)
{
  #if dimension == 1
    return (3.*coarse(s) + coarse(s,child.x))/4.;
  #elif dimension == 2
    return (9.*coarse(s) + 
	    3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
	    coarse(s,child.x,child.y))/16.;
  #else // dimension == 3
    return (27.*coarse(s) + 
	    9.*(coarse(s,child.x) + coarse(s,0,child.y) +
		coarse(s,0,0,child.z)) + 
	    3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
		coarse(s,0,child.y,child.z)) + 
	    coarse(s,child.x,child.y,child.z))/64.;
  #endif
}

//for multigrid-bilinear
#if MULTIGRID
//static inline double bilinear_embed2 (Point point, scalar s)
double bilinear_embed2 (Point point, scalar s)
{
  if (!coarse(css_test) || !coarse(css_test,child.x))
    return coarse(s);
  #if dimension >= 2
  if (!coarse(css_test,0,child.y) || !coarse(css_test,child.x,child.y))
    return coarse(s);
  #endif
  #if dimension >= 3
  if (!coarse(css_test,0,0,child.z) || !coarse(css_test,child.x,0,child.z) ||
      !coarse(css_test,0,child.y,child.z) ||
      !coarse(css_test,child.x,child.y,child.z))
    return coarse(s);  
  #endif
  return bilinear_no_cs (point, s);
  //return bilinear (point, s);
}

//#define bilinear(point, s) bilinear_embed2(point, s)
#endif // MULTIGRID

void mg_cycle2 (scalar * a, scalar * res, scalar * da,
	       void (* relax2) (scalar * da, scalar * res, 
			       int depth, void * data),
	       void * data,
	       int nrelax, int minlevel, int maxlevel)
{

// printf("pid=%d: begin mg_cycle 110\n",pid());
  /**
  We first define the residual on all levels. */
  for(scalar ss in res){
        bool flag7=true;
        ss.flag_Tsat0_0 = flag7; 
        restriction ({ss}); 
        //version 1: no limit restriction
        //version 2: with limit restriction 
  }
 //restriction (res);
//printf("pid=%d: mg_cycle 119\n",pid());
  /**
  We then proceed from the coarsest grid (*minlevel*) dcss_testown to the
  finest grid. */

  minlevel = min (minlevel, maxlevel);
  // for(scalar s in da){
  //     T_tree_function_T_poisson(s,css_test,Tsat00);  // I am here
  // }
 // printf("pid=%d: mg_cycle 128\n",pid());


  for (int l = minlevel; l <= maxlevel; l++) {

    /**
    On the coarsest grid, we take zero as initial guess. */

    if (l == minlevel)
      foreach_level_or_leaf (l)
	        for (scalar s in da)
	             foreach_blockf (s)
	                 s[] = 0.;

    /**
    On all other grids, we take as initial guess the approximate solution
    on the coarser grid bilinearly interpolated onto the current grid. */

    else
      foreach_level (l)
          for (scalar s in da){
            foreach_blockf (s){
              //s[] = bilinear (point, s);
            // s[] = bilinear_embed2 (point,s);
              //  int choice1 = 2;
              //  if(choice1==1){  //choice1=1 is for choise_T_adapt=2 get_T_from3
              //         bool flag7=true;
              //         s.flag_Tsat0_0 = flag7; 
              //         //s.ff3 = css_test;
              //         s[] = Tl_refine33_T(point,s);    //checked - 3D                    
              //  }else if(choice1==2){  //if T is at the center, then use these;
              //          s[] = bilinear_embed2 (point,s);                    
              //  }
              s[] = bilinear_no_cs (point, s);  //without any embed implimentation at the interface of solid and fluid
             //  s[] = bilinear (point, s);  //without any embed implimentation at the interface of solid and fluid
              //version 1: no limit refine
              //version 2: with limit refine 
            }
          } 
    
    /**
    We then apply homogeneous boundary conditions and do several
    iterations of the relaxation function to refine the initial guess. */

    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax2 (da, res, l, data);
      boundary_level (da, l);
    }
  //  printf("pid=%d: mg_cycle 171\n",pid());
  }

  /**
  And finally we apply the resulting correction to *a*. */

  foreach() {
    scalar s, ds;
    for (s, ds in a, da)
      foreach_blockf (s)
        	s[] += ds[];
  }
}

/**
## Multigrid solver

The multigrid solver itself uses successive calls to the multigrid
cycle to refine an initial guess until a specified tolerance is
reached. 

The maximum number of iterations is controlled by *NITERMAX2* and the
tolerance by *TOLERANCE* with the default values below. */

int NITERMAX2 = 800, NITERMIN2 = 1;
double TOLERANCE2 = 1000.0; //1e-3; //1e-3; //1e-6; tolemrace2=1 ---> 1um

/**
Information about the convergence of the solver is returned in a structure. */

typedef struct {
  int i;              // number of iterations
  double resb, resa;  // maximum residual before and after the iterations
  double sum;         // sum of r.h.s.
  int nrelax;         // number of relaxations
  int minlevel;       // minimum level of the multigrid hierarchy
} mgstats2;

/**
The user needs to provide a function which computes the residual field
(and returns its maximum) as well as the relaxation function. The
user-defined pointer *data* can be used to pass arguments to these
functions. The optional number of relaxations is *nrelax* (default is
one) and *res* is an optional list of fields used to store the
residuals. The minimum level of the hierarchy can be set (default is
zero i.e. the root cell). */

struct MGSolve2 {
  scalar * a, * b;
  double (* residual2) (scalar * a, scalar * b, scalar * res,
		       void * data);
  void (* relax2) (scalar * da, scalar * res, int depth, 
		  void * data);
  void * data;
  
  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats2 mg_solve2 (struct MGSolve2 p)
{
  //  printf("pid=%d: begin mg_solve 234\n",pid());
  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);

  /**
  The boundary conditions for the correction fields are the
  *homogeneous* equivalent of the boundary conditions applied to
  *a*. */

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];
  
  /**
  We initialise the structure storing convergence statistics. */

  mgstats2 s = {0};
  double sum = 0.;
  foreach (reduction(+:sum))
    for (scalar s in p.b)
      sum += s[];
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;
  
  /**
  Here we compute the initial residual field and its maximum. */
//  printf("pid=%d: begin mg_solve 265\n",pid());
  double resb;
  // for(scalar sss in p.a){
  //     sss.flag_Tsat0_0 = false;
  // }
  resb = s.resb = s.resa = p.residual2 (p.a, p.b, res, p.data);   // I am here 
  new_flag = false;
//  printf("pid=%d: begin mg_solve 271\n",pid());
  /**
  We then iterate until convergence or until *NITERMAX2* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE2;
  for (s.i = 0;
       s.i < NITERMAX2 && (s.i < NITERMIN2 || s.resa > p.tolerance);
       s.i++) {
    //  for (s.i = 0;
    //    s.i < NITERMAX2 && (s.i < NITERMIN2 || maxres_T > p.tolerance);
    //    s.i++) {
        // printf("pid=%d: begin mg_solve 282\n",pid());
    mg_cycle2 (p.a, res, da, p.relax2, p.data,
	      s.nrelax,
	      p.minlevel,
	      grid->maxdepth);  // I am here

    // for(scalar tempaa in p.a){
    //   foreach(){
    //     if(css_test[]<CCCCCMIN){
    //         tempaa[] = Tsat00;
    //     }
    //   }
    // }
// printf("pid=%d,i+%d, mg_solve2:291,cccccmin=%g\n",pid(),s.i,CCCCCMIN);
    // for(scalar tempaaa in p.a){
    //   foreach(){
    //    // if(css_test[]<=CCCCCMIN){
    //    // if(css_test[]<=CCCCCMIN){
    //     if(topo_mask[]==0 && css_test2[]<=1.0-1e-6 && css_test2[]>=1e-6){
    //         tempaaa[] = Tsat00;
    //     }
    //   }
    // }
//   for(scalar sss in p.a){
//       sss.flag_Tsat0_0 = false;
//   }
    s.resa = p.residual2 (p.a, p.b, res, p.data);
 //   printf("pid%d: resa=%g resb=%g\n",pid(),s.resa,s.resb);
    /**
    We tune the number of relaxations so that the residual is reduced
    by between 2 and 20 for each cycle. This is particularly useful
    for stiff systems which may require a larger number of relaxations
    on the finest grid. */

#if 1
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
	s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
	s.nrelax--;
    }
#else
    if (s.resa == resb) /* convergence has stopped!! */
      break;
    if (s.resa > resb/1.1 && p.minlevel < grid->maxdepth)
      p.minlevel++;
#endif

    resb = s.resa;
  }
  s.minlevel = p.minlevel;
  
  /**
  If we have not satisfied the tolerance, we warn the user. */

  if (s.resa > p.tolerance) {
  // if (maxres_T > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr, 
	     "WARNING: convergence for %s not reached after %d iterations\n"
	     "  res: %g sum: %g nrelax: %d\n", v.name,
	     s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }
    
  /**
  We deallocate the residual and correction fields and free the lists. */

  if (!p.res)
    delete (res), free (res);
  delete (da), free (da);

  return s;
}









// teperature transfer
 attribute {
    bool flag_Tsat1_0;
 }

void get_tempa_from_Tlff_Tgff_forTg(scalar aa, scalar cc){ //from temperature of fraction of liquid and gass to  temperature of center;
//cc should be original ff

// this should be get from each level??????????

///clamp(ff,0,1)  //////////////////////////////////////////////////////////////////////>??????????????????????????????????????
//?????????????????????????????????????????
//scalar tempa = lista[0];
bool flag_Tsat1_0 = aa.flag_Tsat1_0;

     foreach(){
double eps_new = CCCCMIN;
double eps_new2 = CCCCMIN;//1e-3; // CCCCMIN;
double Tsat1;
if(flag_Tsat1_0){
      Tsat1 = 0.0;
}else{
      Tsat1 = Tsat00;
}



          //if(cc[]>1.0 - eps_new && cc[]<1.0){
          //       tempa[] = aa[];
         //  }else if(cc[]<eps_new && cc[]>0.0){
           if(cc[]<eps_new && cc[]>0.0){
                  tempa[] = Tsat1;
           }else if(cc[]>eps_new && cc[]<(1.0-eps_new)){
          
              //coord na= interface_normal3 (point, ff,hhh),p1a;
              coord na= interface_normal (point, cc),p1a;
	           double alphaa = plane_alpha (cc[], na);
             
              double lengtha = line_length_center (na, alphaa, &p1a);//2-dimension,p is centriod of the interface
              //p1a.x = p1a.x*Delta, p1a.y= p1a.y*Delta;
              coord pla;
              plane_center(na,alphaa,cc[],&pla); //center of liquid
              //pla.x = pla.x*Delta, pla.y = pla.y*Delta;
              coord pga;
             // pga.x = -cc[]*pla.x/(1.0-cc[]);
             // pga.y = -cc[]*pla.y/(1.0-cc[]);
              coord nna;
              nna.x = na.x/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              nna.y = na.y/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              double d1a,d2a,d3a;
              d2a = (p1a.x - 0.0)*nna.x +  (p1a.y - 0.0)*nna.y; // directional distance, toward oppsite normal of interface
              d1a = (p1a.x - pla.x)*nna.x + (p1a.y - pla.y)*nna.y; // normal distance from p1 to pl
            //  d3a = fabs((p1a.x - pga.x)*nna.x + (p1a.y - pga.y)*nna.y); // normal distance from p1 to pg
              if(fabs(d2a)<eps_new2){  
                  tempa[] = Tsat1;

              }else if(fabs(d1a)<eps_new2){ //delet_small is needed, this is taken as some cc[] =0 case
                        // don't know how to hand this, for don't calculate this
                        tempa[] = Tsat1; // in poinson equation, we don't calculate this tempeartue.
              //}else if(!( fabs(nna.x)<eps_new2 || fabs(nna.y)<eps_new2)){
                 }else{

                        tempa[]= (d2a/d1a)*aa[] + (d1a - d2a)/d1a*Tsat1; 
              //      }

              // }else{
              //        T[] = aa[];
              }
         }else if(cc[]>=(1.0-eps_new)){ //pure liquid cell

                 tempa[] = aa[]; //ff[]=1.0

         }else if(cc[]<=0.0){

              tempa[] = Tsat1; //(1-ff[])=1.0

         }

   }

}

void get_tempa_from_Tlff_Tgff_forTg2(scalar aa, scalar cc){ //from temperature of fraction of liquid and gass to  temperature of center;
//cc should be original ff

// this should be get from each level??????????

///clamp(ff,0,1)  //////////////////////////////////////////////////////////////////////>??????????????????????????????????????
//?????????????????????????????????????????
//scalar tempa = lista[0];
bool flag_Tsat1_0 = aa.flag_Tsat1_0;

     foreach(){
double eps_new = CCCCMIN;
double eps_new2 = CCCCMIN;//1e-3; // CCCCMIN;
double Tsat1;
if(flag_Tsat1_0){
      Tsat1 = 0.0;
}else{
      Tsat1 = Tsat00;
}



          //if(cc[]>1.0 - eps_new && cc[]<1.0){
          //       tempa[] = aa[];
         //  }else if(cc[]<eps_new && cc[]>0.0){
           //if(cc[]<eps_new && cc[]>0.0){
            if(cc[]<0.5){
                  tempa[] = Tsat1;
           }else if(cc[]<(1.0-eps_new)){
          
              //coord na= interface_normal3 (point, ff,hhh),p1a;
              coord na= interface_normal (point, cc),p1a;
	           double alphaa = plane_alpha (cc[], na);
             
              double lengtha = line_length_center (na, alphaa, &p1a);//2-dimension,p is centriod of the interface
              //p1a.x = p1a.x*Delta, p1a.y= p1a.y*Delta;
              coord pla;
              plane_center(na,alphaa,cc[],&pla); //center of liquid
              //pla.x = pla.x*Delta, pla.y = pla.y*Delta;
              coord pga;
             // pga.x = -cc[]*pla.x/(1.0-cc[]);
             // pga.y = -cc[]*pla.y/(1.0-cc[]);
              coord nna;
              nna.x = na.x/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              nna.y = na.y/sqrt(pow(na.x,2.0)+pow(na.y,2.0));
              double d1a,d2a,d3a;
              d2a = (p1a.x - 0.0)*nna.x +  (p1a.y - 0.0)*nna.y; // directional distance, toward oppsite normal of interface
              d1a = (p1a.x - pla.x)*nna.x + (p1a.y - pla.y)*nna.y; // normal distance from p1 to pl
            //  d3a = fabs((p1a.x - pga.x)*nna.x + (p1a.y - pga.y)*nna.y); // normal distance from p1 to pg
              if(fabs(d2a)<eps_new2){  
                  tempa[] = Tsat1;

              }else if(fabs(d1a)<eps_new2){ //delet_small is needed, this is taken as some cc[] =0 case
                        // don't know how to hand this, for don't calculate this
                        tempa[] = Tsat1; // in poinson equation, we don't calculate this tempeartue.
              //}else if(!( fabs(nna.x)<eps_new2 || fabs(nna.y)<eps_new2)){
                 }else{

                        tempa[]= (d2a/d1a)*aa[] + (d1a - d2a)/d1a*Tsat1; 
              //      }

              // }else{
              //        T[] = aa[];
              }
         }else if(cc[]>=(1.0-eps_new)){ //pure liquid cell

                 tempa[] = aa[]; //ff[]=1.0

         }
        //  else if(cc[]<=0.0){

        //       tempa[] = Tsat1; //(1-ff[])=1.0

        //  }

   }

}

void get_tempa_from_Tlff_Tgff_forTg3(scalar aa, scalar cc){ //from temperature of fraction of liquid and gass to  temperature of center;
//cc should be original ff

// this should be get from each level??????????

///clamp(ff,0,1)  //////////////////////////////////////////////////////////////////////>??????????????????????????????????????
//?????????????????????????????????????????
//scalar tempa = lista[0];
bool flag_Tsat1_0 = aa.flag_Tsat1_0;

     foreach(){
double eps_new = CCCCMIN;
double eps_new2 = CCCCMIN;//1e-3; // CCCCMIN;
double Tsat1;
if(flag_Tsat1_0){
      Tsat1 = 0.0;
}else{
      Tsat1 = Tsat00;
}

      if(cc[]>eps_new){
         tempa[] = aa[];
      }else{
         tempa[] = Tsat1;
      }
   }

}

/**
## Application to the Poisson--Helmholtz equation

We now apply the generic multigrid solver to the Poisson--Helmholtz equation
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$
We first setup the data structure required to pass the extra
parameters $\alpha$ and $\lambda$. We define $\alpha$ as a face
vector field because we need values at the face locations
corresponding to the face gradients of field $a$. 

*alpha* and *lambda* are declared as *(const)* to indicate that the
function works also when *alpha* and *lambda* are constant vector
(resp. scalar) fields. If *tolerance* is set, it supersedes the
default *TOLERANCE* of the multigrid solver, *nrelax* controls the
initial number of relaxations (default is one), *minlevel* controls
the minimum level of the hierarchy (default is one) and *res* is an
optional list of fields used to store the final residual (which can be
useful to monitor convergence). */

struct Poisson2 {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};

/**
We can now write the relaxation function. We first recover the extra
parameters from the data pointer. */

// static void relax2 (scalar * al, scalar * bl, int l, void * data)
// {
//   scalar a = al[0], b = bl[0];
//   struct Poisson2 * p = (struct Poisson2 *) data;
//   (const) face vector alpha = p->alpha;
//   (const) scalar lambda = p->lambda;
// #if EMBED
//   bool embedded = (a.boundary[embed] != symmetry);
// #endif
  
//   /**
//   We use either Jacobi (under)relaxation or we directly reuse values
//   as soon as they are updated. For Jacobi, we need to allocate space
//   for the new field *c*. Jacobi is useful mostly as it gives results
//   which are independent of the order in which the cells are
//   traversed. This is not the case for the simple traversal, which
//   means for example that results will depend on whether a tree or
//   a multigrid is used (because cells will be traversed in a different
//   order). The same comment applies to OpenMP or MPI parallelism. In
//   practice however Jacobi convergence tends to be slower than simple
//   reuse. */
  
// #if JACOBI
//   scalar c[];
// #else
//   scalar c = a;
// #endif
  
//   if(res2_flag){
//       if(l==level_interface){
//             foreach(){
//                  if(level == level_interface && (css_test[]>CCCCCMIN) && (css_test[]<1.0-CCCCCMIN)){
//                            b[] = b[] + res2[];
//                  }

//             }

//       }

//   }


//   /**
//   We use the face values of $\alpha$ to weight the gradients of the
//   5-points Laplacian operator. We get the relaxation function. */

//   //int choice_T_adap =2;

//   if(choice_T_adap==2){
//     tempa.flag_Tsat1_0 = true;
//     //get_tempa_from_Tlff_Tgff_forTg(a,css_test); //choice_T_adap=2 in diffusion
//   }

//   foreach_level_or_leaf (l) {


// //if(css_test[]>0){
// //if(css_test[]>CCCCMIN){ //contain small in residual, but not in relax
// if((l<level_interface) || ((css_test[]>CCCCMIN) && (l==level_interface))){ //contain small in residual, but not in relax

//     double n = - sq(Delta)*b[], d = - lambda[]*sq(Delta);
//     foreach_dimension() {
//       n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
//       d += alpha.x[1] + alpha.x[];
//       if((phase_flag3==1) && (globali==37 || globali==36 ) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
//                     printf("globali:%d line 637: alpha.x[1]=%g alpha.x[]=%g a[1]=%g a[]=%g\n",globali,alpha.x[1], alpha.x[],a[1],a[]);
//               }
//     }
// #if EMBED
//     if (embedded) {
//       double c, e = embed_flux (point, a, alpha, &c);
//       n -= c*sq(Delta);
//       d += e*sq(Delta);
//     }
//     if (!d)
//       c[] = b[] = 0.;
//     else
// #endif // EMBED
//    // if (embedded) {
//       //homogeneous boundary for da
//    //if (css_test[]>1e-6 && (css_test[]<(1.0-1e-6)) &&  (is_leaf(cell))) {  //????????????????????? in foreach_level cell, should we contain source;
//                                                 // in residual,we don't have source in level cell
   
   
//    if(relax_embed2){
//     //if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN) &&  (is_leaf(cell))) {
//       if (css_test[]>1e-6 && (css_test[]<1.0-1e-6) &&  (is_leaf(cell))) {
//  //     if (css_test[]>1e-6 && (css_test[]<1.0-1e-6)) {
//     //if (1) {
      
//     //if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN)) {
//    // if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN) && (l==level_interface) && (is_leaf(cell))) {
//    //  if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN) &&  (is_leaf(cell)) && (l==level_interface)) {
//       double bc1;
//      // if(a.boundary[0] != a.boundary_homogeneous[0])
//      //          bc1 = Tsat00;
//       //else
//                bc1 = 0.0;
//        double c1,e;
      
//       // if(choice_T_adap==2){
          
//       //     e = embed_flux2 (point, tempa, alpha, &c1,bc1);
//       // }else{
//      //      e = embed_flux2 (point, a, alpha, &c1,bc1);
//       // }
//   //   if(embed3_flux_flag){
//   //          e = embed_flux3 (point, tempa, alpha, &c1,bc1);
//   //   }else{
//            e = embed_flux2 (point, a, alpha, &c1,bc1);
//   //   }

//       n -= c1*sq(Delta);
//       d += e*sq(Delta);

//        if((phase_flag3==0) && (globali==37 || globali==36 ) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
//                     printf("globali:%d embed relax2: c1=%g e=%g c[]=%g n=%g d=%g\n",globali,c1,e,c[],n,d);
//         }
//    }
//   }




//    //  }
//     if (!d)
//       c[] = b[] = 0.;
//     else{

//         // c[] = n/d;
//      // if(css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN) &&  (is_leaf(cell))){
//        //       c[] = n/d/css_test[]; 
//        //}else{
//                c[] = n/d; 
//       //}
//               //  if(css_test[] < CCCCCMIN){
//               //        c[] = 0.0;
//               //  }
//               // if(globali==37){
//               //     if(c[]>2.0){
//               //         printf("ff[]=%g n=%g, d=%g c[]=%g\n",css_test[],n,d,c[]);
//               //     }
//               // }
//               if((phase_flag3==0) && (globali==37 || globali==36 ) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
//                     printf("globali:%d line relax2 final: c[]=%g n=%g d=%g\n",globali,c[],n,d);
//               }

//     }

//     }else{//css_test[]
//         c[] =0.0;
//     }
//   }

//   /**
//   For weighted Jacobi we under-relax with a weight of 2/3. */
  
// #if JACOBI
//   foreach_level_or_leaf (l)
//     a[] = (a[] + 2.*c[])/3.;
// #endif
  
// #if TRASH
//   scalar a1[];
//   foreach_level_or_leaf (l)
//     a1[] = a[];
//   trash ({a});
//   foreach_level_or_leaf (l)
//     a[] = a1[];
// #endif
// }

static void relax2 (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson2 * p = (struct Poisson2 *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
#if EMBED3
  bool embedded = (a.boundary[embed] != symmetry);
#endif
  
  /**
  We use either Jacobi (under)relaxation or we directly reuse values
  as soon as they are updated. For Jacobi, we need to allocate space
  for the new field *c*. Jacobi is useful mostly as it gives results
  which are independent of the order in which the cells are
  traversed. This is not the case for the simple traversal, which
  means for example that results will depend on whether a tree or
  a multigrid is used (because cells will be traversed in a different
  order). The same comment applies to OpenMP or MPI parallelism. In
  practice however Jacobi convergence tends to be slower than simple
  reuse. */
  
#if JACOBI
  scalar c[];
#else
  scalar c = a;
#endif
  
  if(res2_flag){
      if(l==level_interface){
            foreach(){
                 if(level == level_interface && (css_test[]>CCCCCMIN) && (css_test[]<1.0-CCCCCMIN)){
                           b[] = b[] + res2[];
                 }

            }

      }

  }


  /**
  We use the face values of $\alpha$ to weight the gradients of the
  5-points Laplacian operator. We get the relaxation function. */

  //int choice_T_adap =2;

  if(choice_T_adap==2){
    tempa.flag_Tsat1_0 = true;
    //get_tempa_from_Tlff_Tgff_forTg(a,css_test); //choice_T_adap=2 in diffusion
  }

  foreach_level_or_leaf (l) {


//if(css_test[]>0){
//if(css_test[]>CCCCMIN){ //contain small in residual, but not in relax
//if((l<level_interface) || ((css_test[]>CCCCMIN) && (l==level_interface))){ //contain small in residual, but not in relax
if(1==1){ 
    double n = - sq(Delta)*b[], d = - lambda[]*sq(Delta);
    foreach_dimension() {
      n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
      d += alpha.x[1] + alpha.x[];
      // if((phase_flag3==1) && (globali==37 || globali==36 ) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
      //               printf("globali:%d line 637: alpha.x[1]=%g alpha.x[]=%g a[1]=%g a[]=%g\n",globali,alpha.x[1], alpha.x[],a[1],a[]);
      //         }
    }
#if EMBED3
    if (embedded) {
      double c, e = embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      c[] = b[] = 0.;
    else
#endif // EMBED
   // if (embedded) {
      //homogeneous boundary for da
   //if (css_test[]>1e-6 && (css_test[]<(1.0-1e-6)) &&  (is_leaf(cell))) {  //????????????????????? in foreach_level cell, should we contain source;
                                                // in residual,we don't have source in level cell
   
   
   if(relax_embed2){
    //if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN) &&  (is_leaf(cell))) {
     // if (css_test[]>1e-6 && (css_test[]<1.0-1e-6) &&  (is_leaf(cell))) {
    //  if (css_test[]>1e-6 && (css_test[]<1.0-1e-6) &&  (is_leaf(cell))) {
    // if (css_test3_n[]>0.0 && (css_test[]>1e-6) && (css_test[]<1.0-1e-6) &&  (is_leaf(cell))) {
      if (1==0) {
   // if (css_test3_n[]>0.0 && (ff[]>1e-6 && ff[]<1.0-1e-6) && (is_leaf(cell))) {// && (css_test[]>1e-6 && css_test[]<1.0-1e-6)){ //  &&  (is_leaf(cell))) {
 //     if (css_test[]>1e-6 && (css_test[]<1.0-1e-6)) {
    //if (1) {
      
    //if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN)) {
   // if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN) && (l==level_interface) && (is_leaf(cell))) {
   //  if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN) &&  (is_leaf(cell)) && (l==level_interface)) {
      double bc1;
     // if(a.boundary[0] != a.boundary_homogeneous[0])
     //          bc1 = Tsat00;
      //else
               bc1 = 0.0;
       double c1,e;
      
      // if(choice_T_adap==2){
          
      //     e = embed_flux2 (point, tempa, alpha, &c1,bc1);
      // }else{
     //      e = embed_flux2 (point, a, alpha, &c1,bc1);
      // }
  //   if(embed3_flux_flag){
  //          e = embed_flux3 (point, tempa, alpha, &c1,bc1);
  //   }else{
           e = embed_flux2 (point, a, alpha, &c1,bc1);
  //   }

      n -= c1*sq(Delta);
      d += e*sq(Delta);
      double c1_2,e_2;
      e_2 = embed_flux2_css_test2 (point, a, alpha, &c1_2,bc1);
      n -= c1_2*sq(Delta);
      d += e_2*sq(Delta);

      //  if((phase_flag3==0) && (globali==37 || globali==36 ) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
      //               printf("globali:%d embed relax2: c1=%g e=%g c[]=%g n=%g d=%g\n",globali,c1,e,c[],n,d);
      //   }
   }
  }




   //  }
    if (!d)
      c[] = b[] = 0.;
    else{

        // c[] = n/d;
     // if(css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN) &&  (is_leaf(cell))){
       //       c[] = n/d/css_test[]; 
       //}else{
               c[] = n/d; 
      //}
              //  if(css_test[] < CCCCCMIN){
              //        c[] = 0.0;
              //  }
              // if(globali==37){
              //     if(c[]>2.0){
              //         printf("ff[]=%g n=%g, d=%g c[]=%g\n",css_test[],n,d,c[]);
              //     }
              // }
              // if((phase_flag3==0) && (globali==37 || globali==36 ) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
              //       printf("globali:%d line relax2 final: c[]=%g n=%g d=%g\n",globali,c[],n,d);
              // }

    }

    }else{//css_test[]
        c[] =0.0;
    }
  }

  /**
  For weighted Jacobi we under-relax with a weight of 2/3. */
  
#if JACOBI
  foreach_level_or_leaf (l)
    a[] = (a[] + 2.*c[])/3.;
#endif
  
#if TRASH
  scalar a1[];
  foreach_level_or_leaf (l)
    a1[] = a[];
  trash ({a});
  foreach_level_or_leaf (l)
    a[] = a1[];
#endif
}


/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree meshrelax2 
requires more careful consideration... */
//#if TGRAD_LEON
#if 1
   #include "Tgrad-leon.h"
#endif

// static double residual2 (scalar * al, scalar * bl, scalar * resl, void * data)
// {
//   scalar a = al[0], b = bl[0], res = resl[0];
//   struct Poisson2 * p = (struct Poisson2 *) data;
//   (const) face vector alpha = p->alpha;
//   (const) scalar lambda = p->lambda;
// #if EMBED
//   bool embedded = (a.boundary[embed] != symmetry);
// #endif
//   double maxres = 0.;
// #if TREE

// printf("pid=%d: residual2 770\n",pid());

// //int choice_T_adap =2;
// //scalar temp_a[];
// if( choice_T_adap == 2){
//     tempa.flag_Tsat1_0 = false;
//     //get_tempa_from_Tlff_Tgff_forTg(a,css_test); //choice_T_adap=2 in diffusion

//    // get_tempa_from_Tlff_Tgff_forTg(a,css_test); //choice_T_adap=2 in diffusion //comment for Tlff

//           // char name93[80]; 
//           // sprintf(name93,"tempa-poisson-check%g.dat", t);
//           // FILE * fp93 = fopen(name93,"w");
//           // int num2=0;
//           // foreach(){
//           //   //if(ff[]>EPS && (1.0-ff[])>EPS)
//           //   fprintf (fp93,"%g %g %g %g\n",  x, y,a[], tempa[]);
//           // }
//           // fclose(fp93);
//     get_tempa_from_Tlff_Tgff_forTg3(a,css_test);

//     // useless
//     //foreach(){
//     //  a[] = tempa[];
//     //}
//     boundary({a});
//     //
// }
// printf("pid=%d: residual2 798\n",pid());
//   /* conservative coarse/fine discretisation (2nd order) */
//   face vector g[];
// //  foreach_face()
// //    g.x[] = alpha.x[]*face_gradient_x (a, 0);
//   foreach_face(){
//       //g.x[] = alpha.x[]*embed_face_gradient2_x(point,a,0);
//       g.x[] = alpha.x[]*face_gradient2_x(a,0);  // I am here 
// //       g.x[] = alpha.x[]*face_gradient2_x(tempa,0);
//       // g.x = alpha.x * face_gradient2_x();
//   }
// printf("pid=%d: residual2 809\n",pid());
// //#if TGRAD_LEON 
// #if 1
//     scalar phase0Tg[],phase1Tg[];
//     Tgrad_leon(css_test,phase0Tg,phase1Tg);
// #endif
// printf("pid=%d: residual2 815\n",pid());

//   total_area_diff = 0.0;
//   flux_total_diff = 0.0;

//   foreach (reduction(max:maxres)) {

//  // if(css_test[]>0.0){
//    if(css_test[]>CCCCCMIN){
//          if((phase_flag3==0) && (globali==37 || globali==36) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
//                     printf("globali:%d b lambda1 residual2: res[]=%g\n",globali,res[]);
//         }
//     res[] = b[] - lambda[]*a[];
//          if((phase_flag3==0) && (globali==37 || globali==36) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
//                     printf("globali:%d b lambda2 residual2: b[]=%g lambda[]=%g a[]=%g res[]=%g\n",globali,b[1],lambda[],a[],res[]);
//         }
//     foreach_dimension(){
//       res[] -= (g.x[1] - g.x[])/Delta;
//           if((phase_flag3==0) && (globali==37 || globali==36) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
//                         printf("globali:%d gx residual2: g.x[1]=%g g.x[]=%g res[]=%g\n",globali,g.x[1],g.x[],res[]);
//             }
//     }
// #if EMBED
//     if (embedded) {
//       double c, e = embed_flux (point, a, alpha, &c);
//       res[] += c - e*a[];
//     }
// #endif // EMBED 


//     //res without embed_flux;
    

//       if(!res2_flag){
//         // if (embedded) {
//             //dirichlet boundary for a
//         if (css_test[]>1e-6 && (css_test[]<(1.0-1e-6))) {
//         //if (css_test[]>1e-6) {

//         // if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN)) {

//       // even ff<ccccmin, I should calculate the diffusion flux across the interface

//       //  if (css_test[]>1e-6 && (css_test[]<1.0-1e-6)) {
//             double bc1;
//           // if(a.boundary[0] != a.boundary_homogeneous[0]){
//                     bc1 = Tsat00;
//           // }
//             //else{
//             //         bc1 = 0.0;
//           // }


//           /*
//             char name152[80];
//             sprintf(name152,"check-boundary-homogeneous.dat");
//             FILE * fp152 = fopen(name152,"a");
//             fprintf(fp152,"residual,bc1= %g\n",bc1);
//             fclose(fp152);
//             */
            
//             double c1, e;
//       #if 0
//           // c1 = 0.0;
//             //here ff is chnanged by phase 
//             //Tgrad_leon(ff,phase0Tgrad,phase1Tgrad);
//             double mua = 0., fa = 0.;
//             foreach_dimension() {
//                 mua += alpha.x[] + alpha.x[1];
//                 fa  += fss_test.x[] + fss_test.x[1];
//               }
//             // total_area_diff += area;
//             // flux_total_diff += area * grad;
//             // *val = - mua/(fa + SEPS)*grad*area/Delta;
//             // return - mua/(fa + SEPS)*coef*area/Delta;
//             // because css_test is changed by phase, here we 
//             //    adopt grad=phase1Tg[],another is grad=-phase0Tg[];
//             coord n = interface_normal7(point,css_test,hhh); //interface_normal7
//             coord p;
//             double alpha = plane_alpha (css_test[], n);
//             double area = plane_area_center (n, alpha, &p);
//             double grad;
//             grad = phase1Tg[];
//             c1 =  - mua/(fa + SEPS)*grad*area/Delta;
//             e = -mua/(fa+SEPS)*0.0*area/Delta;
//             res[] += c1 - e*a[];
//       #else
//             if(!embed3_flux_flag){
//               if(choice_T_adap==2){
//               //    e = embed_flux2 (point, tempa, alpha, &c1,bc1);
//                   flag_v0_nodata = false;
//                   e = embed_flux2 (point, a, alpha, &c1,bc1);
//                   if(flag_v0_nodata){
//                       double mua = 0., fa = 0.;
//                       foreach_dimension() {
//                           mua += alpha.x[] + alpha.x[1];
//                           fa  += fss_test.x[] + fss_test.x[1];
//                         }
//                        coord n = interface_normal7(point,css_test,hhh); // interface_normal7(point,css_test,hhh); //interface_normal7
//                         coord p;
//                         double alpha = plane_alpha (css_test[], n);
//                         double area = plane_area_center (n, alpha, &p);
//                       double grad = phase1Tg[];
//                       c1 =  - mua/(fa + SEPS)*grad*area/Delta;
//                       e = -mua/(fa+SEPS)*0.0*area/Delta;
//                   }
//                 // e = embed_flux2 (point, tempa, alpha, &c1,bc1);
//               }else{
//                   e = embed_flux2 (point, a, alpha, &c1,bc1);
//               }
//             }else{
//                   e = embed_flux3 (point, tempa, alpha, &c1,bc1);
//             }
//             res[] += c1 - e*a[];
//             if((phase_flag3==0) && (globali==37 || globali==36) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
//                     printf("globali:%d line embed residual2: c1=%g e=%g a[]=%g res[]=%g\n",globali,c1,e,a[],res[]);
//         }
//       #endif


            

//             /////////////////////////////////////
//             //add this 2022115; T*c = left + right---T=(left+right)/c
//             // require small dt ?????????
//             /////////////////////////////////////
//             //res[] = res[]/css_test[];
//         }

//   }else{  //res[] without embed_flux; res2[] with embed_flux;
//        res2[] = 0.0;
//          if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN)) {
//       double bc1;
//                bc1 = Tsat00;      
//       double c1, e;
//       if(!embed3_flux_flag){
//         if(choice_T_adap==2){
//             e = embed_flux2 (point, tempa, alpha, &c1,bc1);
//            // e = embed_flux2 (point, a, alpha, &c1,bc1);
//         }else{
//             e = embed_flux2 (point, a, alpha, &c1,bc1);
//         }
//       }else{
//             e = embed_flux3 (point, tempa, alpha, &c1,bc1);
//       }
//      // res[] += c1 - e*a[];
//       res2[] = (c1 - e*a[]);

//    }
    
//   }

//  }else{

//       res[] = 0.0;
//   }
//    // }
//     //  if(css_test[]<CCCCCMIN){
//     //         res[] = 0.0;
//     //   }
   
//     if (fabs (res[]) > maxres){
//     // if ((fabs (res[]) > maxres) && (css_test[]>CCCCCMIN)){
//       maxres = fabs (res[]);
//     }
//   }  //foreach

//    //printf("t=%g diff: total_area_diff=%g flux_total_diff=%g\n",t,total_area_diff,flux_total_diff);
// #else // !TREE
//   /* "naive" discretisation (only 1st order on trees) */
//   foreach (reduction(max:maxres)) {
//     res[] = b[] - lambda[]*a[];
//     foreach_dimension(){
//       // res[] += (alpha.x[0]*face_gradient_x (a, 0) -
//       //		alpha.x[1]*face_gradient_x (a, 1))/Delta;  
//      // res[] += (alpha.x[0]*embed_face_gradient2_x (point,a, 0) -
// 	//	alpha.x[1]*embed_face_gradient2_x (point,a, 1))/Delta;  
//        res[] += (alpha.x[0]*face_gradient2_x (a, 0) -
// 		alpha.x[1]*face_gradient2_x (a, 1))/Delta;  

//      }
// #if EMBED
//     if (embedded) {
//       double c, e = embed_flux (point, a, alpha, &c);
//       res[] += c - e*a[];
//     }
// #endif // EMBED

//      //if (embedded) {
//        //dirichlet boundary for a
//        double bc1 = Tsat00;
//        double c1, e;
//       // if(choice_T_adap==2){
//       //    // tempa.flag_Tsat1_0 = false;
//       //      e = embed_flux2 (point, tempa, alpha, &c1,bc1);
//       // }else{
//       //     e = embed_flux2 (point, a, alpha, &c1,bc1);
//       // }
//       if(!embed3_flux_flag){
//         if(choice_T_adap==2){
//             e = embed_flux2 (point, tempa, alpha, &c1,bc1);
//         }else{
//             e = embed_flux2 (point, a, alpha, &c1,bc1);
//         }
//       }else{
//             e = embed_flux3 (point, tempa, alpha, &c1,bc1);
//       }
//       res[] += c1 - e*a[];
//      //}
     
//     if (fabs (res[]) > maxres)
//       maxres = fabs (res[]);
//   }
  
// #endif // !TREE
//   return maxres;




// }

//val1-0.5deltax-val-0.5deltax-val2-deltax-val3

// double thirdorder(double val1,double val2,double val3){
//       return ()
// }

foreach_dimension()
static double vof_concentration_gradient_forT_x (Point point, scalar c, scalar t, bool inver,double cmin)
{
  //static const double cmin = 0.5;//0.5;
  double cl = c[-1], cc = c[], cr = c[1];
  if (inver)
     cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  //if (cc >= cmin && t.gradient != zero) {
  // // // if (cc >= cmin) {  
  // // //   if (cr >= cmin) {
  // // //     if (cl >= cmin) {
	// // // // if (t.gradient)
	// // // //   return t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
  // // //       if (1==1)
  // // //         return thirdorder(t[-1], t[], t[1])/Delta;
  // // //       else
  // // //         return (t[1] - t[-1])/(2.*Delta);
  // // //     }
  // // //     else{  // tl - val - t - tr //extrapolate
	// // //         return (t[1] - t[])/Delta;
  // // //     }
  // // //   }
  // // //   else if (cl >= cmin)
  // // //     return (t[] - t[-1])/Delta;
  // // // }
  if (cc >= cmin) {  
    if (cl >= cmin) {
      if (cr >= cmin) {
	// if (t.gradient)
	//   return t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
        if (1==1){
          //return thirdorder(t[-1], t[], t[1])/Delta;
        }else{
          //return (t[1] - t[-1])/(2.*Delta);
          return (t[] - t[-1])/Delta;
        }
      }
      else{  // tl - val - t - tr //extrapolate
	        return (t[] - t[-1])/Delta;
      }
    }
    else if (cr >= cmin)
      return (t[1] - t[])/Delta;
  }
  return 0.;
}


extern scalar phase0Tgrad,phase1Tgrad;

static double residual2 (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson2 * p = (struct Poisson2 *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
#if EMBED3
  bool embedded = (a.boundary[embed] != symmetry);
#endif
  double maxres = 0.;
  maxres_T = 0.0;
#if TREE

// printf("pid=%d: residual2 770\n",pid());

// //int choice_T_adap =2;
// //scalar temp_a[];
// if( choice_T_adap == 2){
//     tempa.flag_Tsat1_0 = false;
//     //get_tempa_from_Tlff_Tgff_forTg(a,css_test); //choice_T_adap=2 in diffusion

//    // get_tempa_from_Tlff_Tgff_forTg(a,css_test); //choice_T_adap=2 in diffusion //comment for Tlff

//           // char name93[80]; 
//           // sprintf(name93,"tempa-poisson-check%g.dat", t);
//           // FILE * fp93 = fopen(name93,"w");
//           // int num2=0;
//           // foreach(){
//           //   //if(ff[]>EPS && (1.0-ff[])>EPS)
//           //   fprintf (fp93,"%g %g %g %g\n",  x, y,a[], tempa[]);
//           // }
//           // fclose(fp93);
//     get_tempa_from_Tlff_Tgff_forTg3(a,css_test);

//     // useless
//     //foreach(){
//     //  a[] = tempa[];
//     //}
//     boundary({a});
//     //
// }
// printf("pid=%d: residual2 798\n",pid());

  /* conservative coarse/fine discretisation (2nd order) */
  // face vector g[];

  // foreach_face(){
  //    // g.x[] = alpha.x[]*face_gradient2_x(a,0);  // I am here 
  //     g.x[] = alpha.x[]*face_gradient_s_x(a,0);  // temperature for gas and fluid is based on css_test3_n; temperature for solid based on css_test3
  // }











//20230119 I am here
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////












//caculate gradient at the gas and liquid interface for both side;
#if 1
    // scalar phase0Tg[],phase1Tg[];
    Tgrad_leon(ff,phase0Tg,phase1Tg);  //can not use temperature outside of gas and liquid, can not use temperature in solid(css_test3_n<0)
    //phase0Tg for gas and phase1Tg for liquid, direction of both are from gas to liquid
    // following res[] , inside is positive; so, phase1Tg for liquid is out ,then res[]-=phase1Tg; phase0Tg for gas is out, then res[]-=phase0Tg; 
    //still we need to calculate gradient at the multiphase interface, but on both side, attedtion
    //Tgrad_leon_s(css_test_n,phase0Tg,phase1Tg); // I write code to calculate gradient at the solid interface
#endif

face vector g[];
// foreach_face()
//    g.x[] = alpha.x[]*face_gradient_x (a, 0);

//based on css_test

if(1==1){
    foreach_face(){
      // both is ok, because cs=1.0
      g.x[] = alpha.x[]*face_gradient_x (a, 0);
      //g.x[] = alpha.x[]*face_gradient_no_cs_x (a, 0);
    }
}else{
    foreach_face(){
        double val;
        double cmin = 0.3;//1e-12;
        if(css_test3_n[]<=0){ // pure solid cell
              bool flag=true;
              val = vof_concentration_gradient_forT_x(point,css_test3_n,a,flag,cmin);
        }else if(css_test3_n[]>=1){
              if(css_test[]>=1){ //pure liquid cell
                    bool flag=false;
                    val = vof_concentration_gradient_forT_x(point,css_test,a,flag,cmin);
              }else if(css_test[]<=0){ // pure gas cell
                    bool flag=false;
                    val = vof_concentration_gradient_forT_x(point,css_test2,a,flag,cmin);
              }else{
                  double val1=0.0,val2=0.0,weight1=0.0,weight2=0.0;
                  if(css_test[]>=cmin){
                    bool flag = false;
                    val1 = vof_concentration_gradient_forT_x(point,css_test,a,flag,cmin);
                    weight1 = css_test[];
                  } 
                  if(css_test2[]>=cmin){
                    bool flag = false;
                    val2 = vof_concentration_gradient_forT_x(point,css_test2,a,flag,cmin); 
                    weight2 = css_test2[];
                  }
                  if((weight1+weight2>0.01)){
                      val = (val1*weight1 + val2*weight2)/(weight1+weight2);
                  }else{
                      val = 0.0;
                  }
              }
        }else{ // solid-fluid(gas and liquid) mixed-cell
              double val1=0.0,val2=0.0,weight1=0.0,weight2=0.0;
              double val_fluid;
                  if(css_test[]>=cmin){
                    bool flag = false;
                    val1 = vof_concentration_gradient_forT_x(point,css_test,a,flag,cmin);
                    weight1 = css_test[];
                  } 
                  if(css_test2[]>=cmin){
                    bool flag = false;
                    val2 = vof_concentration_gradient_forT_x(point,css_test2,a,flag,cmin); 
                    weight2 = css_test2[];
                  }
                  if((weight1+weight2>0.001)){
                      val_fluid = (val1*weight1 + val2*weight2)/(weight1+weight2);
                  }else{
                      weight1=0.0;
                      weight2=0.0;
                      val_fluid = 0.0;
                  }
              double val_solid,weight_solid=0.0;
              if((1.0-css_test3_n[])>cmin){
                  bool flag = true;
                  val_solid = vof_concentration_gradient_forT_x(point,css_test3_n,a,flag,cmin);
                  weight_solid = 1.0-css_test3_n[];
              }
              if(weight1+weight2>0.0 && weight_solid>0.0){
                    val = val_solid*weight_solid   + val_fluid*(1-weight_solid);
              }else if(weight1+weight2>0.0){
                    val = val_fluid;
              }else{
                    val = val_solid;
              }
        }
        g.x[] = alpha.x[]*val;
      }

}



//n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
//      d += alpha.x[1] + alpha.x[];

// I am here

  total_area_diff = 0.0;
  flux_total_diff = 0.0;

  foreach (reduction(max:maxres) reduction(max:maxres_T)) {
bool flag2=false;
 // if(css_test[]>0.0){
 //  if(css_test[]>CCCCCMIN){
  if(1==1){ // without any specification        
    res[] = b[] - lambda[]*a[];        
    foreach_dimension(){
      res[] -= (g.x[1] - g.x[])/Delta;
    }
#if EMBED3
    if (embedded) {
      double c, e = embed_flux (point, a, alpha, &c);
      res[] += c - e*a[];
    }
#endif // EMBED 


    //res without embed_flux;
    bool one_temp=false;

      if(!res2_flag){
        // if (embedded) {
            //dirichlet boundary for a
        //if (css_test[]>1e-6 && (css_test[]<(1.0-1e-6))) {
           if(css_test3_n[]>0.0 && (!one_temp)){



  //     if(css_test3_n[]>0.0 && ((ff[]>1e-6) && (ff[]<1.0-1e-6))){  //caculate phase flux only in two-phase interface
         // if ((ff[]>1e-6) && (ff[]<(1.0-1e-6))) {
            if (1==1) {
          //if (css_test[]>1e-6) {
          // if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN)) {
        // even ff<ccccmin, I should calculate the diffusion flux across the interface
        //  if (css_test[]>1e-6 && (css_test[]<1.0-1e-6)) {
              double bc1;
            // if(a.boundary[0] != a.boundary_homogeneous[0]){
                      bc1 = Tsat00;
            // }
              //else{
              //         bc1 = 0.0;
            // }
            /*
              char name152[80];
              sprintf(name152,"check-boundary-homogeneous.dat");
              FILE * fp152 = fopen(name152,"a");
              fprintf(fp152,"residual,bc1= %g\n",bc1);
              fclose(fp152);
              */
              
              double c1, e=0.0;
              double c1_2,e_2=0.0;
        #if 0
            // c1 = 0.0;
              //here ff is chnanged by phase 
              //Tgrad_leon(ff,phase0Tgrad,phase1Tgrad);
              double mua = 0., fa = 0.;
              foreach_dimension() {
                  mua += alpha.x[] + alpha.x[1];
                  fa  += fss_test.x[] + fss_test.x[1];
                }
              // total_area_diff += area;
              // flux_total_diff += area * grad;
              // *val = - mua/(fa + SEPS)*grad*area/Delta;
              // return - mua/(fa + SEPS)*coef*area/Delta;
              // because css_test is changed by phase, here we 
              //    adopt grad=phase1Tg[],another is grad=-phase0Tg[];
              coord n = interface_normal7(point,css_test,hhh); //interface_normal7
              coord p;
              double alpha = plane_alpha (css_test[], n);
              double area = plane_area_center (n, alpha, &p);
              double grad;
              grad = phase1Tg[];
              c1 =  - mua/(fa + SEPS)*grad*area/Delta;
              e = -mua/(fa+SEPS)*0.0*area/Delta;
              res[] += c1 - e*a[];
        #else 
              if(!embed3_flux_flag){
                if(choice_T_adap==2){
                //    e = embed_flux2 (point, tempa, alpha, &c1,bc1);
                    flag_v0_nodata = false;//true;//false;
         //          e = embed_flux2 (point, a, alpha, &c1,bc1);
                      //  if((!flag_v0_nodata) && (fabs(c1)>1000.0 || fabs(e)>1000.0) ){ //&& new_flag){
                      //        printf("classic:c1=%g,e=%g,x=%g,y=%g,z=%g,ff=%g,css_test3=%g,css_test3_n=%g\n",c1,e,x,y,z,ff[],css_test3[],css_test3_n[]);
                      //   }
                    double area=0.0;
                  //  if(flag_v0_nodata){
                     if(1==1){
                        double mua = 0., fa = 0.;
                        foreach_dimension() {
                           mua += alpha.x[] + alpha.x[1];
                           fa  += fss_test.x[] + fss_test.x[1];
                         }
                         //mua = k1;//Tkl;
                        coord n = interface_normal(point, ff);//interface_normal7(point,ff,hhh); //
                          double alpha = plane_alpha (ff[], n);
                          coord pp;
                          area = plane_area_center (n, alpha, &pp);
                        double grad = phase1Tg[];// for fluid this is value for from water to out direction
                      //  c1 =  - mua/(fa + SEPS)*grad*area/Delta;

                        mua = Tkl;
                        fa =1.0;
                        c1 =  - mua/(fa + SEPS)*grad*area/Delta;
                        //c1 =  - mua*grad*area/Delta;
                        // if(fabs(c1)>1000.0 ){ //&& new_flag ){
                        //      printf("pid=%d,leon:c1=%g,x=%g,y=%g,z=%g,grad=%g,ff=%g,css_test3=%g,css_test3_n=%g\n",pid(),c1,x,y,z,grad,ff[],css_test3[],css_test3_n[]);
                        // }
                        //e = -mua/(fa+SEPS)*0.0*area/Delta;
                        e = 0.0;
                    }
                    flag_v0_nodata = false;//false;
           //         e_2 = embed_flux2_css_test2 (point, a, alpha, &c1_2, bc1);
                    // if((!flag_v0_nodata) && (fabs(c1_2)>1000.0 || fabs(e_2)>1000.0) ){ //&& new_flag){
                    //          printf("classic:c1_2=%g,e_2=%g,x=%g,y=%g,z=%g,ff=%g,css_test3=%g,css_test3_n=%g\n",c1_2,e_2,x,y,z,ff[],css_test3[],css_test3_n[]);
                    //     }
                    double area2=0.0;
             //       if(flag_v0_nodata){
                      if(1==1){
                        double mua = 0., fa = 0.;
                        foreach_dimension() {
                            mua += alpha.x[] + alpha.x[1];
                            fa  += fss_test2.x[] + fss_test2.x[1];
                          }
                        //mua= k0;//Tkg;
                        coord n = interface_normal(point, ff_oppo);// interface_normal7(point,ff_oppo,hhh_ff_oppo); // interface_normal7(point,css_test,hhh); //interface_normal7
                          coord pp;
                          double alpha = plane_alpha (ff_oppo[], n);
                          area2 = plane_area_center (n, alpha, &pp);
   //direction of phase1Tg and phase0Tg is from fluid to out; for liquid it is out, for gas it is in; so add one -, finally it is +
                        double grad = -phase0Tg[];  
                        // if(fabs(area-area2)<1e-4){
                        //     printf("area1-area2 is too big:%g\n",area-area2);
                        // }
                   //     c1_2 =  - mua/(fa + SEPS)*grad*area2/Delta;
                        
                        mua = Tkg;
                        fa= 1.0;
                        c1_2 =  - mua/(fa + SEPS)*grad*area2/Delta;
                        //c1_2 =  - mua*grad*area2/Delta;
                        // if(fabs(c1_2)>1000.0 ){ //&& new_flag){
                        //      printf("leon:c1_2=%g,x=%g,y=%g,z=%g,grad=%g,ff=%g,css_test3=%g,css_test3_n=%g\n",c1_2,x,y,z,grad,ff[],css_test3[],css_test3_n[]);
                        // }
                       // e_2 = -mua/(fa+SEPS)*0.0*area/Delta;
                       e_2 = 0.0;
                    }
                    // if(fabs(area-area2)>1e-4){
                    //         printf("area1-area2 is too big:%g; area=%g,area2=%g,%g,%g,%g\n",area-area2,area,area2,x,y,z);
                    //     }
                  // e = embed_flux2 (point, tempa, alpha, &c1,bc1);

                  if(fabs(c1)>1e-1 || fabs(c1_2)>1e-1){
                    flag2=true;
                  }
                }else{
                    e = embed_flux2 (point, a, alpha, &c1,bc1);
                }
              }else{
                    e = embed_flux3 (point, tempa, alpha, &c1,bc1);
              }
              res[] = res[] + ((c1_2 + c1) - (e + e_2)*a[]);
              // if(!((fabs(c1)<1e-3) && (fabs(c1_2)<1e-3))){
              //   if(fabs(res[])>1000.0){
              //         //printf("t=%g:c1=%g,c1_2=%g,-phase0=%g,phase1=%g,%g,%g,%g\n",t,c1,c1_2,-phase0Tg[],phase1Tg[],x,y,z);
              //         printf("     x=%g,y=%g,z=%g,ff=%g,c1_2=%g,c1=%g,css_test3=%g,css_test3_n=%g,css_test[]=%g\n",x,y,z,ff[],c1_2,c1,css_test3[],css_test3_n[],css_test[]);
              //         foreach_neighbor(){
              //              printf("     x=%g,y=%g,z=%g,ff=%g,c1_2=%g,c1=%g,css_test3=%g,css_test3_n=%g,css_test[]=%g\n",x,y,z,ff[],c1_2,c1,css_test3[],css_test3_n[],css_test[]);
              //         }
              //  }
              // if((phase_flag3==0) && (globali==37 || globali==36) && fabs(x-0.636719)<0.001 &&  fabs(y-0.496094)<0.001 && fabs(z-0.496094)<0.001){
              //         printf("globali:%d line embed residual2: c1=%g e=%g a[]=%g res[]=%g\n",globali,c1,e,a[],res[]);
              
         // }
        #endif


              

              /////////////////////////////////////
              //add this 2022115; T*c = left + right---T=(left+right)/c
              // require small dt ?????????
              /////////////////////////////////////
              //res[] = res[]/css_test[];
          }
          // if((1.0-css_test3_n[]>1e-6) && (css_test3_n[]>1e-6)){

          // }
      }
  }else{  //res[] without embed_flux; res2[] with embed_flux;
       printf("no mpi no mpi no mpi\n");
       res2[] = 0.0;
         if (css_test[]>CCCCCMIN && (css_test[]<1.0-CCCCCMIN)) {
      double bc1;
               bc1 = Tsat00;      
      double c1, e;
      if(!embed3_flux_flag){
        if(choice_T_adap==2){
            e = embed_flux2 (point, tempa, alpha, &c1,bc1);
           // e = embed_flux2 (point, a, alpha, &c1,bc1);
        }else{
            e = embed_flux2 (point, a, alpha, &c1,bc1);
        }
      }else{
            e = embed_flux3 (point, tempa, alpha, &c1,bc1);
      }
     // res[] += c1 - e*a[];
      res2[] = (c1 - e*a[]);

   }
    
  }

 }else{

      res[] = 0.0;
  }
   // }
    //  if(css_test[]<CCCCCMIN){
    //         res[] = 0.0;
    //   }
   
    if (fabs (res[]) > maxres){
    // if ((fabs (res[]) > maxres) && (css_test[]>CCCCCMIN)){
      maxres = fabs (res[]);
    }
    if(flag2){
        if(fabs(res[]/(Trhog*Tcpg)) > maxres_T){
          maxres_T = fabs(res[]/(Trhog*Tcpg));
        }
    }else{
        if(fabs(res[]/rhocp[]) > maxres_T){
          maxres_T = fabs(res[]/rhocp[]);
        }
    }
  }  //foreach

   //printf("t=%g diff: total_area_diff=%g flux_total_diff=%g\n",t,total_area_diff,flux_total_diff);
#else // !TREE  // For three phase problem, I did not modified anything

  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension(){
      // res[] += (alpha.x[0]*face_gradient_x (a, 0) -
      //		alpha.x[1]*face_gradient_x (a, 1))/Delta;  
     // res[] += (alpha.x[0]*embed_face_gradient2_x (point,a, 0) -
	//	alpha.x[1]*embed_face_gradient2_x (point,a, 1))/Delta;  
       res[] += (alpha.x[0]*face_gradient2_x (a, 0) -
		alpha.x[1]*face_gradient2_x (a, 1))/Delta;  

     }
#if EMBED3
    if (embedded) {
      double c, e = embed_flux (point, a, alpha, &c);
      res[] += c - e*a[];
    }
#endif // EMBED

     //if (embedded) {
       //dirichlet boundary for a
       double bc1 = Tsat00;
       double c1, e;
      // if(choice_T_adap==2){
      //    // tempa.flag_Tsat1_0 = false;
      //      e = embed_flux2 (point, tempa, alpha, &c1,bc1);
      // }else{
      //     e = embed_flux2 (point, a, alpha, &c1,bc1);
      // }
      if(!embed3_flux_flag){
        if(choice_T_adap==2){
            e = embed_flux2 (point, tempa, alpha, &c1,bc1);
        }else{
            e = embed_flux2 (point, a, alpha, &c1,bc1);
        }
      }else{
            e = embed_flux3 (point, tempa, alpha, &c1,bc1);
      }
      res[] += c1 - e*a[];
     //}
     
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  
#endif // !TREE
  return maxres;




}






















/**
## User interface

Finally we provide a generic user interface for a Poisson--Helmholtz
equation of the form
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$ */

mgstats2 poisson4 (struct Poisson2 p)
{

  printf("poisson4 begin\n");
 // printf("pid=%d: begin poisson4\n",pid());
  p.minlevel = minl;
  /**
  If $\alpha$ or $\lambda$ are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda.i)
    p.lambda = zeroc;

  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  face vector alpha = p.alpha;
  scalar lambda = p.lambda;
  //printf("pid=%d:  poisson4 1060\n",pid());
  restriction ({alpha,lambda});
  //printf("pid=%d:  poisson4 1062\n",pid());
  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  printf("poisson4 1814 %d\n",pid());


  double defaultol = TOLERANCE2;
  if (p.tolerance)
    TOLERANCE2 = p.tolerance;

//  printf("pid=%d:  poisson4 1071\n",pid());
  scalar a = p.a, b = p.b;

//  printf("pid=%d, poisson4:before mg_solves\n",pid());
  mgstats2 s = mg_solve2 ({a}, {b}, residual2, relax2,
			&p, p.nrelax, p.res, minlevel = max(1, p.minlevel));  // I am here
//  printf("pid=%d, poisson4:after mg_solves\n",pid());

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE2 = defaultol;

  return s;
}

/**
## Projection of a velocity field

The function below "projects" the velocity field *u* onto the space of
divergence-free velocity fields i.e.
$$
\mathbf{u}_f^{n+1} \leftarrow \mathbf{u}_f - \Delta t\alpha\nabla p
$$
so that
$$
\nabla\cdot\mathbf{u}_f^{n+1} = 0
$$
This gives the Poisson equation for the pressure
$$
\nabla\cdot(\alpha\nabla p) = \frac{\nabla\cdot\mathbf{u}_f}{\Delta t}
$$ */

struct Project2 {
  face vector uf;
  scalar p;
  face vector alpha; // optional: default unityf
  double dt;         // optional: default one
  int nrelax;        // optional: default four
};

trace
mgstats2 project2 (struct Project2 q)
{
  face vector uf = q.uf;
  scalar p = q.p;
  (const) face vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;
  
  /**
  We allocate a local scalar field and compute the divergence of
  $\mathbf{u}_f$. The divergence is scaled by *dt* so that the
  pressure has the correct dimension. */

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += uf.x[1] - uf.x[];
    div[] /= dt*Delta;
  }

  /**
  We solve the Poisson problem. The tolerance (set with *TOLERANCE*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats2 mgp = poisson4 (p, div, alpha,
			 tolerance = TOLERANCE2/sq(dt), nrelax = nrelax);

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face()
    uf.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}


// foreach_dimension()
// static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
// {
//   static const double cmin = 0.5;
//   double cl = c[-1], cc = c[], cr = c[1];
//   if (t.inverse)
//     cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
//   if (cc >= cmin && t.gradient != zero) {
//     if (cr >= cmin) {
//       if (cl >= cmin) {
// 	if (t.gradient)
// 	  return t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
// 	else
// 	  return (t[1]/cr - t[-1]/cl)/(2.*Delta);
//       }
//       else
// 	return (t[1]/cr - t[]/cc)/Delta;
//     }
//     else if (cl >= cmin)
//       return (t[]/cc - t[-1]/cl)/Delta;
//   }
//   return 0.;
// }

