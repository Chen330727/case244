// note: u is weighted by fm
double timestep (const face vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt;
#if EMBED
    // #if AXI
    //   // if(fs.x[]>1e-3){
    //   //     assert (fm.x[]);
    //   //     dt *= fm.x[];
    //   // }else{
    //   //     dt *= 1e-3*y;
    //   // }
    //   if(fs.x[]>1e-6){
    //       dt *= fs.x[];
    //   }else{
    //       dt *= 1e-3;
    //   }
    // #else
    if(fs.x[]>1e-1){
           dt = Delta/fabs(u.x[]);
           dt *=  fm.x[];
    }
    // #endif
#else
      dt = Delta/fabs(u.x[]);
      dt *= cm[];
#endif
      if (dt < dtmax) dtmax = dt;
    }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
