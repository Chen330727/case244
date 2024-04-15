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
extern double Tkl,Tkg;
extern face vector modphase0;
extern face vector modphase1;
extern int phase_flag3;
extern scalar ff;
extern int minl;
extern int level_interface;
extern double Tsat00;
extern scalar cs,T;
extern face vector fs;
//extern scalar rhocp;

double k_value;

extern int maxl;

#define EPS_FS 0

extern  vector smallmodl,bigmodl,smallmodg,bigmodg;
extern scalar resl,resg;

void mg_cycle1 (scalar * a, scalar * res, scalar * da,
	       void (* relax1) (scalar * da, scalar * res, 
			       int depth, void * data),
	       void * data,
	       int nrelax, int minlevel, int maxlevel)
{

  /**
  We first define the residual on all levels. */

  restriction (res);

  /**
  We then proceed from the coarsest grid (*minlevel*) down to the
  finest grid. */

  minlevel = min (minlevel, maxlevel);
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
	    s[] = bilinear (point, s);
    }
  }
    /**
    We then apply homogeneous boundary conditions and do several
    iterations of the relaxation function to refine the initial guess. */

    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax1 (da, res, l, data);
      boundary_level (da, l);
    }
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

The maximum number of iterations is controlled by *NITERMAX* and the
tolerance by *TOLERANCE* with the default values below. */

int NITERMAX1 = 10, NITERMIN1 = 1;
double TOLERANCE1 = 1e-2; 

/**
Information about the convergence of the solver is returned in a structure. */

typedef struct {
  int i;              // number of iterations
  double resb, resa;  // maximum residual before and after the iterations
  double sum;         // sum of r.h.s.
  int nrelax;         // number of relaxations
  int minlevel;       // minimum level of the multigrid hierarchy
} mgstats1;

/**
The user needs to provide a function which computes the residual field
(and returns its maximum) as well as the relaxation function. The
user-defined pointer *data* can be used to pass arguments to these
functions. The optional number of relaxations is *nrelax* (default is
one) and *res* is an optional list of fields used to store the
residuals. The minimum level of the hierarchy can be set (default is
zero i.e. the root cell). */

struct MGSolve1 {
  scalar * a, * b;
  double (* residual1) (scalar * a, scalar * b, scalar * res,
		       void * data);
  void (* relax1) (scalar * da, scalar * res, int depth, 
		  void * data);
  void * data;
  
  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats1 mg_solve1 (struct MGSolve1 p)
{

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

  mgstats1 s = {0};
  double sum = 0.;
  foreach (reduction(+:sum))
    for (scalar s in p.b)
      sum += s[];
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;
  
  /**
  Here we compute the initial residual field and its maximum. */

  double resb;
  resb = s.resb = s.resa = p.residual1 (p.a, p.b, res, p.data);

  /**
  We then iterate until convergence or until *NITERMAX* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE1;
  for (s.i = 0;
       s.i < NITERMAX1 && (s.i < NITERMIN1 || s.resa > p.tolerance);
       s.i++) {
    mg_cycle1 (p.a, res, da, p.relax1, p.data,
	      s.nrelax,
	      p.minlevel,
	      grid->maxdepth);
    s.resa = p.residual1 (p.a, p.b, res, p.data);

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
  
//trace the final res[]
for(scalar s in res){
    if(phase_flag3==0){
      foreach(){
          resg[] = s[];
      }
    }else{
        foreach(){
          resl[] = s[];
      }
    }
}
  /**
  If we have not satisfied the tolerance, we warn the user. */

  if (s.resa > p.tolerance) {


//  char name34[80];
//   sprintf(name34,"mass_record17-pid%d.dat",pid());
//   FILE * fp34 = fopen(name34,"w");
// #if EMBED
//   foreach(noauto){
//     for(scalar a in res)
//     fprintf(fp34,"%g %g %g %g %g %g %g\n",x,y,ff[],T[],rhocp[],a[],cm[]);    
//      //        fprintf(fp34,"%g %g %g\n",x,y,u.y[]);
//       }

//  for(scalar a in res){
//          foreach_boundary(left)
//     for (int i = -BGHOSTS; i < 0; i++)
//       fprintf (fp34, "%g %g %g %g %g %g %g\n", x + (i + 0.5)*Delta, y, ff[i],T[i],rhocp[i],a[i],cm[i]);
  
//   foreach_boundary(right)
//     for (int i = 1; i <= BGHOSTS; i++)
//       fprintf (fp34, "%g %g %g %g %g %g %g\n", x + (i - 0.5)*Delta, y,ff[i], T[i],rhocp[i],a[i],cm[i]);
    
//   foreach_boundary(bottom)
//     for (int i = -BGHOSTS; i < 0; i++)
//       fprintf (fp34, "%g %g %g %g %g %g %g\n", x , y + (i + 0.5)*Delta,ff[0,i], T[0,i],rhocp[0,i],a[0,i],cm[0,i]);
    
//   foreach_boundary(top)
//     for (int i = 1; i <= BGHOSTS; i++)
//       fprintf (fp34, "%g %g %g %g %g %g %g\n", x, y + (i - 0.5)*Delta, ff[0,i],T[0,i],rhocp[0,i],a[0,i],cm[0,i]);
// }
// //  foreach_face(){
// //     fprintf(fp34,"%g %g %g %g %g %g %g %g %g %g\n",x,y,ff[],T[],fs.x[],fm.x[],source_pc[],ulf.x[],uf.x[],usf.x[]);
// //   }
// #else
//   // foreach(){ 
//   //   //fprintf(fp33,"%g %g %g\n",x,y,T[]);
//   //  // fprintf(fp33,"%g %g %g %g %g %g %d %g %g\n",x,y,ff[],T[],cm[],source_pc[],topo_mask[],phase0Tgrad[],phase1Tgrad[]);
//   //  fprintf(fp33,"%g %g %g %g %g\n",x,y,ff[],T[],source_pc[]);
//   // }
//   //   foreach_face(){
//   //   fprintf(fp33,"%g %g %g %g\n",x,y,fm.x[],alpha.x[]);
//   // }
// #endif

//   fclose(fp34);
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(pid()==0){
//                   char command1[150];
//                   sprintf(command1, "LC_ALL=C cat mass_record17-pid*.dat > outfacets/mass_record17-%g-%d",t,phase_flag3);
//                   system(command1);

//                   char command7[150];
//                   sprintf(command7, "LC_ALL=C rm -rf mass_record17-pid*.dat");
//                   system(command7);
//     }




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

struct Poisson1 {
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

static void relax1 (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson1 * p = (struct Poisson1 *) data;
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
  
  /**
  We use the face values of $\alpha$ to weight the gradients of the
  5-points Laplacian operator. We get the relaxation function. */

  foreach_level_or_leaf (l) {
    //  if(css_test3_n[]>=0.0){ //save4
     if(css_test3_n[]>=1.0){ //save5
         // bool solid_flag= false;
          //solid_flag = (css_test3_n[]<1.0);
          bool ff_flag = ((((ff[]>=0.5 ) && (phase_flag3==1)) || ((ff[]<0.5 ) && (phase_flag3==0))) && (fabs(ff[]-0.5)>1e-6)); //ff not changed
          bool solid_flag = (css_test3_n[]<1.0);
          if(ff_flag || solid_flag){
              double n = - sq(Delta)*b[], d = - lambda[]*sq(Delta);
              foreach_dimension() {
                if(1==0){
                      n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
                      d += alpha.x[1] + alpha.x[];
                }else{
                    //    foreach_dimension() {
                        double rightd;
                        double leftd;
                        double value_left;
                        double value_right;
                        bool left_s = false;
                        bool right_s =false;

                        bool center_mix = false;
                        if(solid_flag){
                            center_mix = true;
                        }
                        if(css_test3_n[1]<1.0){
                            right_s = true;
                        }
                        if(css_test3_n[-1]<1.0){
                            left_s = true;
                        }

                        //if(1==0 && l==level_interface){
                          if( l==level_interface){
                              if(phase_flag3==1){
                                    rightd = modphase1.x[1];
                                    leftd  = modphase1.x[];
                              }else{
                                    rightd = modphase0.x[1];
                                    leftd  = modphase0.x[];
                              }
                            if(fabs(rightd)<1e-2){
                                  rightd = 1.0;
                              }else if(fabs(rightd)>1.0){
                                  rightd = 1.0;
                              }else if(right_s | center_mix ){ // distance between onefluid and cell and fluid cell keeps 1
                                  rightd = 1.0;
                              }
                              if(fabs(leftd)<1e-2){
                                  leftd = 1.0;
                              }else if(fabs(leftd)>1.0){
                                  leftd =1.0;
                              }else if(left_s | center_mix){
                                  leftd = 1.0;
                              }
                                    if(((ff[]>0.5)^(ff[-1]>0.5)) && fabs(leftd)<1.0){ //1 0or 0 1 
                                          value_left = 0.0; //Tsat00
                                    }else{
                                          value_left = a[-1];
                                    }
                                    
                                    // if(left_s){
                                    //       value_left = 0.0;
                                    // }

                                    if(((ff[]>0.5)^(ff[1]>0.5)) && fabs(rightd)<1.0 ){ //1 0or 0 1
                                          value_right = 0.0; //Tsat00
                                    }else{
                                          value_right = a[1];
                                    }

                                    // if(right_s){
                                    //        value_right = 0.0;
                                    // }
                          }else{
                              rightd =  1.0;
                              leftd = 1.0;
                              value_right = a[1];
                                value_left = a[-1];
                                //  if(left_s){ // neibor is solid one-fluid cell
                                //         value_left = 0.0;
                                //   }
                                //   if(right_s){ // neibor is solid one-fluid cell
                                //         value_right = 0.0;
                                //   }
                          }
                      double left_k,right_k;// = k_value*fm.x[];
                          //  #if EMBED
                          //     double temp1,temp2,temp3,temp4;
                          //     temp1 = fs.x[];
                          //       if(temp1<EPS_FS){
                          //         temp1 = EPS_FS;
                          //       }
                          //       temp2 = max(y-Delta/2.0,1e-20);
                          //     left_k = k_value*temp1*temp2;//Tkg*fss_test.x[];

                          //     temp3 = fs.x[1];
                          //       if(temp3<EPS_FS){
                          //         temp3 = EPS_FS;
                          //       }
                          //       temp4 = max(y+Delta/2.0,1e-20);
                          //     right_k = k_value*temp3*temp4;//Tkg*fss_test.x[];
                          // #else
                                left_k = k_value*fm.x[];//max(y,1e-20);//fm.x[];
                                right_k = k_value*fm.x[1];//max(y,1e-20); //;fm.x[1];
                          // #endif
                      
                      if(left_s | center_mix){ // neibor is solid one-fluid cell
                              left_k = alpha.x[];
                              value_left = a[-1]; //0.0;
                        }
                        if(right_s | center_mix){ // neibor is solid one-fluid cell
                              right_k = alpha.x[1];
                              value_right = a[1];// 0.0;
                        }
                        // if(left_k<k_value*Delta*0.4){
                        //       left_k = k_value*Delta*0.4;
                        // }
                        // if(right_k<k_value*Delta*0.4){
                        //       right_k = k_value*Delta*0.4;
                        // }
                      n += 2.0*(right_k)/((rightd+leftd)*(rightd))*value_right + 2.0*(left_k)/((rightd+leftd)*(leftd))*value_left;
                      d += 2.0*(right_k)/((rightd+leftd)*(rightd)) + 2.0*(left_k)/((rightd+leftd)*(leftd));

                      // n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
                      // d += alpha.x[1] + alpha.x[];
                  //  } //foreach_dimension
                }
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
                c[] = n/d;
      }else{ //!ff_flag
          c[] =0.0;
      }
    }else{
       c[] = 0.0;
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
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

extern bool poisson_check;
static double residual1 (scalar * al, scalar * bl, scalar * resl, void * data)
{ 
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson1 * p = (struct Poisson1 *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
#if EMBED3
  bool embedded = (a.boundary[embed] != symmetry); 
#endif
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*face_gradient_x (a, 0);
  foreach (reduction(max:maxres)) {
    //bool ff_flag = ((((ff[]>0.5 ) && (phase_flag3==1)) || ((ff[]<0.5 ) && (phase_flag3==0))) && (fabs(ff[]-0.5)>1e-6)); //ff not changed
    //bool ff_flag = ((((ff[]>0.5 ) && (phase_flag3==1)) || ((ff[]<0.5 ) && (phase_flag3==0))) && (fabs(ff[]-0.5)>1e-6)); //ff not changed
    bool solid_flag = false;
    if(css_test3_n[]>=1){
           bool ff_flag = ((((ff[]>0.5 ) && (phase_flag3==1)) || ((ff[]<0.5 ) && (phase_flag3==0))) && (fabs(ff[]-0.5)>1e-6)); //ff not changed
           if(ff_flag){
                res[] = b[] - lambda[]*a[];
                foreach_dimension(){
                  if(1==0){
                    res[] -= (g.x[1] - g.x[])/Delta;
                  }else{
            //////////////////////////////////////////////
            double rightd;
                    double leftd;
                    double value_left;
                    double value_right;
                    bool left_s = false;
                    bool right_s = false;
                     if(css_test3_n[1]<1.0){
                        right_s = true;
                     }
                     if(css_test3_n[-1]<1.0){
                        left_s = true;
                     }
                  // if(1==0 && level==level_interface){
                    if(level==level_interface){
                          if(phase_flag3==1){
                                rightd = modphase1.x[1];
                                leftd  = modphase1.x[];
                          }else{
                                rightd = modphase0.x[1];
                                leftd  = modphase0.x[];
                          }
                          if(fabs(rightd)<1e-2){
                              rightd = 1.0;
                          }else if(fabs(rightd)>1.0){
                              rightd = 1.0;
                          }else if(right_s || solid_flag){ // distance between onefluid and cell and fluid cell keeps 1
                              rightd = 1.0;
                          }
                          if(fabs(leftd)<1e-2){
                              leftd = 1.0;
                          }else if(fabs(leftd)>1.0){
                              leftd =1.0;
                          }else if(left_s || solid_flag){
                              leftd = 1.0;
                          }
                              if(((ff[]>0.5)^(ff[-1]>0.5)) && fabs(leftd)<1.0 ){ //1 0or 0 1
                                    value_left = Tsat00;
                                  // printf("leftd=%g,value_left=%g\n",leftd,value_left);
                              }else{
                                    value_left = a[-1];
                              }

                                if(((ff[]>0.5)^(ff[1]>0.5)) && fabs(rightd)<1.0 ){ //1 0or 0 1
                                    value_right =  Tsat00;
                                  //  printf("rightd=%g,value_right=%g\n",rightd,value_right);
                                }else{
                                    value_right = a[1];
                                }
                      }else{
                          rightd =  1.0;
                          leftd = 1.0;
                            value_left = a[-1];
                            value_right = a[1];
                      }

                  
                  // n += 2.0*(k_value)/((rightd+leftd)*(rightd))*value_right + 2.0*(k_value)/((rightd+leftd)*(leftd))*value_left;
                  // d += 2.0*(k_value)/((rightd+leftd)*(rightd)) + 2.0*(k_value)/((rightd+leftd)*(leftd));

                  double value_center = a[];
                  //double left_k = k_value*fm.x[], right_k = k_value*fm.x[1];
                  double left_k,right_k;
                          // #if EMBED
                          //     double temp1,temp2,temp3,temp4;
                          //     temp1 = fs.x[];
                          //       if(temp1<EPS_FS){
                          //         temp1 = EPS_FS;
                          //       }
                          //       temp2 = max(y-Delta/2.0,1e-20);
                          //     left_k = k_value*temp1*temp2;//Tkg*fss_test.x[];

                          //     temp3 = fs.x[1];
                          //       if(temp3<EPS_FS){
                          //         temp3 = EPS_FS;
                          //       }
                          //       temp4 = max(y+Delta/2.0,1e-20);
                          //     right_k = k_value*temp3*temp4;//Tkg*fss_test.x[];
                          // #else
                                left_k = k_value*fm.x[];//*max(y,1e-20);//fm.x[];
                                right_k = k_value*fm.x[1];//max(y,1e-20); //;fm.x[1];
                          // #endif
                  if(left_s || solid_flag){ // neibor is solid one-fluid cell
                        left_k = alpha.x[];
                        value_left = a[-1]; //0.0;
                  }
                  if(right_s || solid_flag){ // neibor is solid one-fluid cell
                        right_k = alpha.x[1];
                        value_right = a[1];
                        
                  }
                  // if(left_k<k_value*Delta*0.4){
                  //       left_k = k_value*Delta*0.4;
                  // }
                  // if(right_k<k_value*Delta*0.4){
                  //       right_k = k_value*Delta*0.4;
                  // }
                if(phase_flag3==0){
                    smallmodg.x[]=leftd;
                    bigmodg.x[]=rightd;
                }else{
                    smallmodl.x[]=leftd;
                    bigmodl.x[]=rightd;
                }
                  res[] -= (2.0*(right_k)*(value_right-value_center)/(Delta*rightd) - 2.0*(left_k)*(value_center-value_left)/(Delta*leftd))/((rightd+leftd)*Delta);

            /////////////////////////////////////////////
                  }
                }
            }else if(fabs(ff[]-0.5)<=1e-6){
                  a[] = Tsat00;
                  res[] = 0.0;
             }else{
                  a[] = Tsat00;
                    res[] = 0.0;
                }
          
                if (fabs (res[]) > maxres){
                  maxres = fabs (res[]);
                  // if(fabs(maxres-376528)<2.0){
                  //       printf("maxres=%g x=%g y=%g f=%g cs=%g\n",maxres,x,y,ff[],cs[]);
                  // }
                }
    }else{
       //a[] = Tsat00;
       res[] = 0.0;
    }
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] += (alpha.x[0]*face_gradient_x (a, 0) -
		alpha.x[1]*face_gradient_x (a, 1))/Delta;  
#if EMBED3
    if (embedded) {
      double c, e = embed_flux (point, a, alpha, &c);
      res[] += c - e*a[];
    }
#endif // EMBED
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif // !TREE
  

  if(poisson_check){
      if(pid()==0){
            char name93[80]; 
            sprintf(name93,"poisson_check.dat");
            FILE * fp93 = fopen(name93,"a");
            int num2=0;
              fprintf (fp93," %g\n",  maxres);
            fclose(fp93);
      }
  }
  return maxres;
}

/**
## User interface

Finally we provide a generic user interface for a Poisson--Helmholtz
equation of the form
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$ */

mgstats1 poisson1 (struct Poisson1 p)
{
  if(phase_flag3==1){
      k_value = Tkl;
  }else{
      k_value = Tkg;
  }
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
  restriction ({alpha,lambda});

  // for(int l=minl;l<maxl;l++){
  // char name93[150];
  // sprintf(name93,"outfacets/alpha_%d_level%d_pid%d.dat",phase_flag3,l,pid());
  // FILE * fp93 = fopen(name93,"a");
  //   foreach_level(l){
  //       fprintf(fp93,"%g %g %d %g %g %g %g %g\n",x,y,l,ff[],alpha.x[],alpha.x[1],alpha.y[],alpha.y[1]);
  //  }
  //  fclose(fp93);
  //   MPI_Barrier(MPI_COMM_WORLD);
  //  if(pid()==0){
  //                 char command1[150];
  //                 sprintf(command1, "LC_ALL=C cat outfacets/alpha_0_level%d_pid*.dat > outfacets/alpha_0_level%d_%d",l,l,globali);
  //                 system(command1);

  //                 char command7[150];
  //                 sprintf(command7, "LC_ALL=C rm -rf outfacets/alpha_0_level%d_pid*.dat",l);
  //                 system(command7);

  //                  char command2[150];
  //                sprintf(command2, "LC_ALL=C cat outfacets/alpha_1_level%d_pid*.dat > outfacets/alpha_1_level%d_%d",l,l,globali);
  //                 system(command2);

  //                 char command3[150];
  //                 sprintf(command3, "LC_ALL=C rm -rf outfacets/alpha_1_level%d_pid*.dat",l);
  //                 system(command3);
  //   }
  // }
 


  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE1;
  if (p.tolerance)
    TOLERANCE1 = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats1 s = mg_solve1 ({a}, {b}, residual1, relax1,
			&p, p.nrelax, p.res, minlevel = max(1, p.minlevel));

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE1 = defaultol;

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

struct Project1 {
  face vector uf;
  scalar p;
  face vector alpha; // optional: default unityf
  double dt;         // optional: default one
  int nrelax;        // optional: default four
};

trace
mgstats1 project1 (struct Project1 q)
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

  //for lubomir,maxl = 12

  // if(maxl==12){
  //   if(phase_flag3==0){
  //       // TOLERANCE1 = 1e-5; //for the two temperature
  //       TOLERANCE1 = 1e-7; //for the two temperature
  //   }else{
  //       TOLERANCE1 = 5e-2;
  //   }
  // // }

  // mgstats1 mgp = poisson1 (p, div, alpha,
	// 		 tolerance = TOLERANCE1/sq(dt), nrelax = nrelax);


  TOLERANCE1 = 1e-6;
  mgstats1 mgp = poisson1 (p, div, alpha,
			 tolerance = TOLERANCE1, nrelax = nrelax);

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face()
    uf.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}
