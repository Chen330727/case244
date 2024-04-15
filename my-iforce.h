/**
# Interfacial forces

We assume that the interfacial acceleration can be expressed as
$$
\phi\mathbf{n}\delta_s/\rho
$$
with $\mathbf{n}$ the interface normal, $\delta_s$ the interface Dirac
function, $\rho$ the density and $\phi$ a generic scalar field. Using
a CSF/Peskin-like approximation, this can be expressed as
$$
\phi\nabla f/\rho
$$
with $f$ the volume fraction field describing the interface.

The interfacial force potential $\phi$ is associated to each VOF
tracer. This is done easily by adding the following [field
attributes](/Basilisk C#field-attributes). */
extern scalar true_interface,topo_mask_s,ff;
extern bool sigma_in_project;
extern bool flag_rho;
extern scalar arealg;
extern vector flux_show;

#define is_phase2(v) (v>=0.5)

attribute {
  scalar phi;
}

/**
Interfacial forces are a source term in the right-hand-side of the
evolution equation for the velocity of the [centered Navier--Stokes
solver](navier-stokes/centered.h) i.e. it is an acceleration. If
necessary, we allocate a new vector field to store it. */

event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
  }
}

/**
The calculation of the acceleration is done by this event, overloaded
from [its definition](navier-stokes/centered.h#acceleration-term) in
the centered Navier--Stokes solver. */

event acceleration (i++)
{
  
  /**
  We check for all VOF interfaces for which $\phi$ is allocated. The
  corresponding volume fraction fields will be stored in *list*. */

  scalar * list = NULL;
  for (scalar f in interfaces)
    if (f.phi.i) {
      list = list_add (list, f);

      /**
      To avoid undeterminations due to round-off errors, we remove
      values of the volume fraction larger than one or smaller than
      zero. */

      foreach()
	f[] = clamp (f[], 0., 1.);
    }

  /**
  On trees we need to make sure that the volume fraction gradient
  is computed exactly like the pressure gradient. This is necessary to
  ensure well-balancing of the pressure gradient and interfacial force
  term. To do so, we apply the same prolongation to the volume
  fraction field as applied to the pressure field. */
  
#if TREE
  for (scalar f in list) {
    f.prolongation = p.prolongation;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif

  /**
  Finally, for each interface for which $\phi$ is allocated, we
  compute the interfacial force acceleration
  $$
  \phi\mathbf{n}\delta_s/\rho \approx \alpha\phi\nabla f
  $$ 
  */

  face vector ia = a;
  foreach_face(){
    for (scalar f in list)
      if (f[] != f[-1] && fm.x[] > 0.) {

	/**
	We need to compute the potential *phif* on the face, using its
	values at the center of the cell. If both potentials are
	defined, we take the average, otherwise we take a single
	value. If all fails we set the potential to zero: this should
	happen only because of very pathological cases e.g. weird
	boundary conditions for the volume fraction. */
	
	scalar phi = f.phi;
	double phif =
	  (phi[] < nodata && phi[-1] < nodata) ?
	  (phi[] + phi[-1])/2. :
	  phi[] < nodata ? phi[] :
	  phi[-1] < nodata ? phi[-1] :
	  0.;

double sf0,sf1;
  if(flag_rho){
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
  }else{
        sf1 = sf[];
        sf0 = sf[-1];
  }
 
//	ia.x[] += alpha.x[]/fm.x[]*phif*(f[] - f[-1])/Delta;
	// ia.x[] += alpha.x[]/(fm.x[] + 
double area_total = arealg[-1] + arealg[];
double left_w=arealg[-1]/area_total;
double right_w=arealg[]/area_total;
coord nn_f_l,nn_f_r;
nn_f_l.x=flux_show.x[-1];
nn_f_l.y=flux_show.y[-1];
nn_f_r.x=flux_show.x[];
nn_f_r.y=flux_show.y[];
normalize(&nn_f_l);
normalize(&nn_f_r);
coord nn_f;
nn_f.x = (nn_f_l.x*left_w+nn_f_r.x*right_w);
nn_f.y = (nn_f_l.y*left_w+nn_f_r.y*right_w);
normalize(&nn_f);

        if(sigma_in_project){
            if(topo_mask_s[]!=0){
                if(is_phase(sf1)){
                    sf1=1;
                }else{
                    sf1=0;
                }
            }else{ //topo_mask_s[]==1
                if(ff[]>0 && ff[]<1.0){
                    sf1=1;
                }else{ //ff[]<=0 or ff[]>=1
                    if(is_phase(sf1)){
                        sf1=1;
                    }else{
                        sf1=0;
                    }
                }
            }
            if(topo_mask_s[-1]!=0){
                if(is_phase(sf0)){
                    sf0=1;
                }else{
                    sf0=0;
                }
            }else{ //topo_mask_s[-1]==1
                if(ff[-1]>0 && ff[-1]<1.0){
                    sf0=1;
                }else{ //ff[]<=0
                    if(is_phase(sf0)){
                        sf0=1;
                    }else{
                        sf0=0;
                    }
                }
            }
            //=sigma*curvature*(H_i+1 - H_i)/Delta  pressure jump, phif=sigma*curvature
            
            // jump_p_f.x[] += phif*(sf1 - sf0)*fabs(nn_f.x)/Delta;
            // jump_p_f.x[] += phif*(sf1 - sf0)*fabs(nn_f.x);
            // jump_p_f.x[] += phif*(sf1 - sf0);
        }else{
            if(1==0){
              if(is_phase(sf1)){
                  sf1=1;
              }else{
                  sf1=0;
              }
              if(is_phase(sf0)){
                  sf0=1;
              }else{
                  sf0=0;
              }
            }else{
            //      if(topo_mask_s[]!=0){
            //           if(is_phase(sf1)){
            //               sf1=1;
            //           }else{
            //               sf1=0;
            //           }
            //       }else{ //topo_mask_s[]==1
            //           if(ff[]>0 && ff[]<1.0){
            //               sf1=1;
            //           }else{ //ff[]<=0 or ff[]>=1
            //               if(is_phase(sf1)){
            //                   sf1=1;
            //               }else{
            //                   sf1=0;
            //               }
            //           }
            //       }
            //       if(topo_mask_s[-1]!=0){
            //           if(is_phase(sf0)){
            //               sf0=1;
            //           }else{
            //               sf0=0;
            //           }
            //       }else{ //topo_mask_s[-1]==1
            //           if(ff[-1]>0 && ff[-1]<1.0){
            //               sf0=1;
            //           }else{ //ff[]<=0
            //               if(is_phase(sf0)){
            //                   sf0=1;
            //               }else{
            //                   sf0=0;
            //               }
            //           }
            //       }
            } 
            ia.x[] += alpha.x[]/(fm.x[] + SEPS)*phif*(sf1 - sf0)/Delta;
        }

      }
      }

  /**
  On trees, we need to restore the prolongation values for the
  volume fraction field. */
  
#if TREE
  for (scalar f in list) {
    f.prolongation = fraction_refine;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif
  
  /**
  Finally we free the potential fields and the list of volume
  fractions. */

  for (scalar f in list) {
    scalar phi = f.phi;
    delete ({phi});
    f.phi.i = 0;
  }
  free (list);
}

/**
## References

See Section 3, pages 8-9 of:

~~~bib
@hal{popinet2018, hal-01528255}
~~~
*/
