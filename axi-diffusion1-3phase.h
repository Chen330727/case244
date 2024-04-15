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

#include "./axi-poisson1-3phase.h"

/**
The parameters of the `diffusion()` function are a scalar field `f`,
scalar fields `r` and $\beta$ defining the reactive term, the timestep
`dt` and a face vector field containing the diffusion coefficient
`D`. If `D` or $\theta$ are omitted they are set to one. If $\beta$ is
omitted it is set to zero. Both `D` and $\beta$ may be constant
fields.

Note that the `r`, $\beta$ and $\theta$ fields will be modified by the solver.

The function returns the statistics of the Poisson solver. */


struct Diffusion1_3phase {
  // mandatory
  scalar f1,f2,f3;
  double dt;
  // optional
  face vector D;  // default 1
  scalar theta;   // default 1
  scalar r1, beta1; // default 0
  scalar r2, beta2; // default 0
  scalar r3, beta3; // default 0
};

void diffusion1_3phase (struct Diffusion1 p) 
{

  /**
  We define $f$ and $r$ for convenience. */

  scalar f1 = p.f1, r1 = automatic (p.r1);
  scalar f2 = p.f2, r1 = automatic (p.r2);
  scalar f1 = p.f3, r1 = automatic (p.r3);

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

  if (p.r1.i)
    foreach()
      r1[] = theta_idt[]*f1[] - r1[];
  else {// r was not passed by the user
    foreach()
      r1[] = theta_idt[]*f1[];
  }

  if (p.r2.i)
    foreach()
      r2[] = theta_idt[]*f2[] - r2[];
  else {// r was not passed by the user
    foreach()
      r2[] = theta_idt[]*f2[];
  }

    if (p.r3.i)
    foreach()
      r3[] = theta_idt[]*f3[] - r3[];
  else {// r was not passed by the user
    foreach()
      r3[] = theta_idt[]*f3[];
  }

   
  /**
  If $\beta$ is provided, we use it to store the diagonal term $\lambda$. */

  scalar lambda1 = theta_idt;
  if (p.beta1.i) {
    scalar beta1 = p.beta1;
    foreach()
      beta1[] += theta_idt[];
    lambda1 = beta1;
  }
    scalar lambda2 = theta_idt;
  if (p.beta2.i) {
    scalar beta2 = p.beta2;
    foreach()
      beta2[] += theta_idt[];
    lambda2 = beta2;
  }
    scalar lambda3 = theta_idt;
  if (p.beta3.i) {
    scalar beta3 = p.beta3;
    foreach()
      beta3[] += theta_idt[];
    lambda3 = beta3;
  }

  /**
  Finally we solve the system. */
double toleranceT;
toleranceT=1e-8;
// if(phase_flag3==0){
//       toleranceT=1e-5;//  5e-6;mesu3-14     //toleranceT=1e-6;
// }else{
//       toleranceT=1e-8;//2e-2;//5e-3;//1e-3;
// }
  return poisson1_3phase (f1, r1, f2,r2,f3,r3, p.D, lambda1,lambda2,lambda3, toleranceT);
    // return poisson1 (f, r, p.D, lambda, toleranceT);
}
