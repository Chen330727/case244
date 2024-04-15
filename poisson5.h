/**
# Multigrid Poisson5--Helmholtz solvers

We want to solve Poisson5--Helmholtz equations of the general form
$$
L(a) = \nabla\cdot (\alpha\nabla a) + \lambda a = b
$$
This can be done efficiently using a multigrid solver. 

An important aspect of Poisson5--Helmholtz equations is that the
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
where the right-hand-side is often called the *residual5* of the
approximate solution $\tilde{a}$. 

## Multigrid cycle

Here we implement the multigrid cycle proper. Given an initial guess
*a*, a residual5 *res*, a correction field *da* and a relaxation
function *relax*, we will provide an improved guess at the end of the
cycle. */

extern scalar poisson_source2;
extern scalar ff,css_test3_n;
extern scalar T;
extern int maxl;
void mg_cycle5 (scalar * a, scalar * res, scalar * da,
	       void (* relax5) (scalar * da, scalar * res, 
			       int depth, void * data),
	       void * data,
	       int nrelax, int minlevel, int maxlevel)
{

  /**
  We first define the residual5 on all levels. */

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
	for (scalar s in da)
	  foreach_blockf (s)
	    s[] = bilinear (point, s);
    
    /**
    We then apply homogeneous boundary conditions and do several
    iterations of the relaxation function to refine the initial guess. */

    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax5 (da, res, l, data); 
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

The maximum number of iterations is controlled by *NITERMAX5* and the
tolerance by *TOLERANCE5* with the default values below. */

int NITERMAX5 = 100, NITERMIN5 = 1;
double TOLERANCE5 = 1e-3; //1e-3;

/**
Information about the convergence of the solver is returned in a structure. */

typedef struct {
  int i;              // number of iterations
  double resb, resa;  // maximum residual5 before and after the iterations
  double sum;         // sum of r.h.s.
  int nrelax;         // number of relaxations
  int minlevel;       // minimum level of the multigrid hierarchy
} mgstats5;

/**
The user needs to provide a function which computes the residual5 field
(and returns its maximum) as well as the relaxation function. The
user-defined pointer *data* can be used to pass arguments to these
functions. The optional number of relaxations is *nrelax* (default is
one) and *res* is an optional list of fields used to store the
residual5s. The minimum level of the hierarchy can be set (default is
zero i.e. the root cell). */

struct MGSolve5 {
  scalar * a, * b;
  double (* residual5) (scalar * a, scalar * b, scalar * res,
		       void * data);
  void (* relax5) (scalar * da, scalar * res, int depth, 
		  void * data);
  void * data;
  
  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats5 mg_solve5 (struct MGSolve5 p)
{

  /**
  We allocate a new correction and residual5 field for each of the scalars
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

  mgstats5 s = {0};
  double sum = 0.;
  foreach (reduction(+:sum))
    for (scalar s in p.b)
      sum += s[];
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;
  
  /**
  Here we compute the initial residual5 field and its maximum. */

  double resb;
  resb = s.resb = s.resa = p.residual5 (p.a, p.b, res, p.data);

  /**
  We then iterate until convergence or until *NITERMAX5* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual5 is lower than *TOLERANCE5*. */

  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE5;
  for (s.i = 0;
       s.i < NITERMAX5 && (s.i < NITERMIN5 || s.resa > p.tolerance);
       s.i++) {
    mg_cycle5 (p.a, res, da, p.relax5, p.data,
	      s.nrelax,
	      p.minlevel,
	      grid->maxdepth);
    s.resa = p.residual5 (p.a, p.b, res, p.data);

    /**
    We tune the number of relaxations so that the residual5 is reduced
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



//  char name34[80];
//   sprintf(name34,"mass_record18-pid%d.dat",pid());
//   FILE * fp34 = fopen(name34,"w");
// #if EMBED
//   foreach(noauto){
//     for(scalar a in res)
//     fprintf(fp34,"%g %g %g %g %g\n",x,y,ff[],T[],a[]);    
//      //        fprintf(fp34,"%g %g %g\n",x,y,u.y[]);
//       }

//  for(scalar a in res){
//          foreach_boundary(left)
//     for (int i = -BGHOSTS; i < 0; i++)
//       fprintf (fp34, "%g %g %g %g %g\n", x + (i + 0.5)*Delta, y, ff[i],T[i],a[i]);
  
//   foreach_boundary(right)
//     for (int i = 1; i <= BGHOSTS; i++)
//       fprintf (fp34, "%g %g %g %g %g\n", x + (i - 0.5)*Delta, y,ff[i], T[i],a[i]);
    
//   foreach_boundary(bottom)
//     for (int i = -BGHOSTS; i < 0; i++)
//       fprintf (fp34, "%g %g %g %g %g\n", x , y + (i + 0.5)*Delta,ff[0,i], T[0,i],a[0,i]);
    
//   foreach_boundary(top)
//     for (int i = 1; i <= BGHOSTS; i++)
//       fprintf (fp34, "%g %g %g %g %g\n", x, y + (i - 0.5)*Delta, ff[0,i],T[0,i],a[0,i]);
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
//                   sprintf(command1, "LC_ALL=C cat mass_record18-pid*.dat > outfacets/mass_record18-poisson5-%g",t);
//                   system(command1);

//                   char command7[150];
//                   sprintf(command7, "LC_ALL=C rm -rf mass_record18-pid*.dat");
//                   system(command7);
//     }


    scalar v = p.a[0];
    fprintf (ferr, 
	     "poisson5--WARNING: convergence for %s not reached after %d iterations\n"
	     "  res: %g sum: %g nrelax: %d\n", v.name,
	     s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }
    
  /**
  We deallocate the residual5 and correction fields and free the lists. */

  if (!p.res)
    delete (res), free (res);
  delete (da), free (da);

  return s;
}

/**
## Application to the Poisson5--Helmholtz equation

We now apply the generic multigrid solver to the Poisson5--Helmholtz equation
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
default *TOLERANCE5* of the multigrid solver, *nrelax* controls the
initial number of relaxations (default is one), *minlevel* controls
the minimum level of the hierarchy (default is one) and *res* is an
optional list of fields used to store the final residual5 (which can be
useful to monitor convergence). */

struct Poisson5 {
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

static void relax5 (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson5 * p = (struct Poisson5 *) data;
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
      // if((css_test3_n[]<=0.0 && l<level_interface) || (css_test3_n[]<1 && l==level_interface)){
      bool no_solid = (css_test3_n[]>=1.0);
        if( (css_test3_n[]<1.0) || no_solid){
      //  if( (css_test3_n[]<1.0)){
                  double n = - sq(Delta)*b[], d = - lambda[]*sq(Delta);
                  foreach_dimension() {
                    n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
                    d += alpha.x[1] + alpha.x[];
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
                      c[] = n/d;
              #else 
                  if (!d)
                    c[] = b[] = 0.;
                  else
                      c[] = n/d;
              #endif
       }else{
             //c[] = b[]= 0.;
             c[] = 0.;
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
The equivalent residual5 function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

static double residual5 (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson5 * p = (struct Poisson5 *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
#if EMBED3
  bool embedded = (a.boundary[embed] != symmetry);
#endif
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  // foreach_face()
  //   g.x[] = alpha.x[]*face_gradient_x (a, 0);
  foreach (reduction(max:maxres)) {
      if(css_test3_n[]<1.0){
              res[] = b[] - lambda[]*a[];
              foreach_dimension(){ 



        res[] -= (2.0*(alpha.x[1])*(a[1]-a[])/(Delta) - 2.0*(alpha.x[])*(a[]-a[-1])/(Delta))/(2.0*Delta);
              
             //   res[] -= (g.x[1] - g.x[])/Delta;
              }
                      #if EMBED3
                          if (embedded) {
                            double c, e = embed_flux (point, a, alpha, &c);
                            res[] += c - e*a[];
                          }
          
          #endif // EMBED 
      }else{
        res[] = 0.0;
      }  
              if (fabs (res[]) > maxres)
                maxres = fabs (res[]);
  }
#else // !TREE
          printf("!TREE in poisson5 \n");
          exit(1);
#endif // !TREE
  return maxres;
}

/**
## User interface

Finally we provide a generic user interface for a Poisson5--Helmholtz
equation of the form
$$
\nabla\cdot (\alpha\nabla a) + \lambda a = b
$$ */

mgstats5 poisson5 (struct Poisson5 p)
{

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

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE5;
  if (p.tolerance)
    TOLERANCE5 = p.tolerance;

  scalar a = p.a, b = p.b;

  // foreach(){
  //   if((css_test3_n[]<1.0-1e-6 && css_test3_n[]>1e-6) && (ff[]>1e-6 && ff[]<1.0-1e-6))
  //     b[] = b[] - poisson_source2[];
  // }


  //for lubomir,maxl = 12
  //  if(maxl==12){
      TOLERANCE5 = 2e-2;//2e-3;//1e-3;
  //  }
  mgstats5 s = mg_solve5 ({a}, {b}, residual5, relax5,
			&p, p.nrelax, p.res, minlevel = max(1, p.minlevel));

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE5 = defaultol;

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
This gives the Poisson5 equation for the pressure
$$
\nabla\cdot(\alpha\nabla p) = \frac{\nabla\cdot\mathbf{u}_f}{\Delta t}
$$ */

struct Project5 {
  face vector uf;
  scalar p;
  face vector alpha; // optional: default unityf
  double dt;         // optional: default one
  int nrelax;        // optional: default four
};

trace
mgstats5 project5 (struct Project5 q)
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
  We solve the Poisson5 problem. The tolerance (set with *TOLERANCE5*) is
  the maximum relative change in volume of a cell (due to the divergence
  of the flow) during one timestep i.e. the non-dimensional quantity 
  $$
  |\nabla\cdot\mathbf{u}_f|\Delta t 
  $$ 
  Given the scaling of the divergence above, this gives */

  mgstats5 mgp = poisson5 (p, div, alpha,
			 tolerance = TOLERANCE5/sq(dt), nrelax = nrelax);

  /**
  And compute $\mathbf{u}_f^{n+1}$ using $\mathbf{u}_f$ and $p$. */

  foreach_face()
    uf.x[] -= dt*alpha.x[]*face_gradient_x (p, 0);

  return mgp;
}
