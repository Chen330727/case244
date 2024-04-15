/**
# Time-implicit discretisation of reaction--diffusion equations

We want to discretise implicitly the reaction--diffusion equation
$$
\theta\partial_tf = \nabla\cdot(D\nabla f) + \beta f + r
$$ 
where $\beta f + r$ is a reactive term,  $D$ is the diffusion
coefficient and $\theta$ can be a density term.

Using a time-implicit backward Euler discretisation, this can be
written
$$
\theta\frac{f^{n+1} - f^{n}}{dt} = \nabla\cdot(D\nabla f^{n+1}) + \beta
f^{n+1} + r
$$
Rearranging the terms we get
$$
\nabla\cdot(D\nabla f^{n+1}) + (\beta - \frac{\theta}{dt}) f^{n+1} =
- \frac{\theta}{dt}f^{n} - r
$$
This is a Poisson--Helmholtz problem which can be solved with a
multigrid solver. */

#include "./axi-poisson1-2.h"

/**
The parameters of the `diffusion()` function are a scalar field `f`,
scalar fields `r` and $\beta$ defining the reactive term, the timestep
`dt` and a face vector field containing the diffusion coefficient
`D`. If `D` or $\theta$ are omitted they are set to one. If $\beta$ is
omitted it is set to zero. Both `D` and $\beta$ may be constant
fields.

Note that the `r`, $\beta$ and $\theta$ fields will be modified by the solver.

The function returns the statistics of the Poisson solver. */

extern int phase_flag3;

struct Diffusion1 {
  // mandatory
  scalar f;
  double dt;
  // optional
  face vector D;  // default 1
  scalar r, beta; // default 0
  scalar theta;   // default 1
};

trace
mgstats1 diffusion1 (struct Diffusion1 p) 
{

  /**
  If *dt* is zero we don't do anything. */

  if (p.dt == 0.) {
    mgstats1 s = {0};
    return s;
  }

  /**
  We define $f$ and $r$ for convenience. */

  scalar f = p.f, r = automatic (p.r);

  /**
  We define a (possibly constant) field equal to $\theta/dt$. */

  const scalar idt[] = - 1./p.dt;
  (const) scalar theta_idt = p.theta.i ? p.theta : idt;
  
  if (p.theta.i) {
    scalar theta_idt = p.theta;
    foreach()
      theta_idt[] *= idt[];
  }

  /**
  We use `r` to store the r.h.s. of the Poisson--Helmholtz solver. */

  if (p.r.i)
    foreach()
      r[] = theta_idt[]*f[] - r[];
  else {// r was not passed by the user
    foreach()
      r[] = theta_idt[]*f[];
  }


   
  /**
  If $\beta$ is provided, we use it to store the diagonal term $\lambda$. */

  scalar lambda = theta_idt;
  if (p.beta.i) {
    scalar beta = p.beta;
    foreach()
      beta[] += theta_idt[];
    lambda = beta;
  }

  /**
  Finally we solve the system. */
double toleranceT;
if(phase_flag3==0){
      toleranceT=1e-5;//  5e-6;mesu3-14     //toleranceT=1e-6;
}else{
      toleranceT=1e-8;//2e-2;//5e-3;//1e-3;
}
  return poisson1 (f, r, p.D, lambda, toleranceT);
}
