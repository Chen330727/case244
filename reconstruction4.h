#include "fractions.h"

#include "curvature.h"
coord interface_normal8 (Point point, scalar c, vector h)
{
  coord n;
  if (!h.x.i || (n = height_normal (point, c, h)).x == nodata)
    n = mycs (point, c);
  return n;
}

void reconstruction4 (const scalar c, vector n, scalar alpha)
{
  vector h_temp[];
  heights(c,h_temp);
  foreach() {

    /**
    If the cell is empty or full, we set $\mathbf{n}$ and $\alpha$ only to
    avoid using uninitialised values in `alpha_refine()`. */

    if (c[] <= 0. || c[] >= 1.) {
      alpha[] = 0.;
      foreach_dimension()
	n.x[] = 0.;
    }
    else {

      /**
      Otherwise, we compute the interface normal using the
      Mixed-Youngs-Centered scheme, copy the result into the normal field
      and compute the intercept $\alpha$ using our predefined function. */

      coord m = interface_normal8 (point, c, h_temp); //mycs
      foreach_dimension()
	n.x[] = m.x;
      alpha[] = plane_alpha (c[], m);
    }
  }

#if TREE

  /**
  On a tree grid, for the normal to the interface, we don't use
  any interpolation from coarse to fine i.e. we use straight
  "injection". */

  foreach_dimension()
    n.x.refine = n.x.prolongation = refine_injection;

  /**
  We set our refinement function for *alpha*. */

  alpha.n = n;
  alpha.refine = alpha.prolongation = alpha_refine;
#endif
}

