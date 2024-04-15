scalar solidflux[];
scalar solidsurfaceT[];
int NITERMAX1_3phase = 10, NITERMIN1_3phase = 1;
double TOLERANCE1_3phase = 1e-2; 

typedef struct {
  int i;              // number of iterations
  double resb1, resa1;  // maximum residual before and after the iterations
  double sum1;         // sum of r.h.s.
  int nrelax1;         // number of relaxations
  

 
  double resb2, resa2;  // maximum residual before and after the iterations
  double sum2;         // sum of r.h.s.
  int nrelax2;         // number of relaxations
 

  double resb3, resa3;  // maximum residual before and after the iterations
  double sum3;         // sum of r.h.s.
  int nrelax3;         // number of relaxations

  int minlevel;       // minimum level of the multigrid hierarchy
} mgstats1_3phase;


struct Poisson1_3phase {
  scalar a1, b1;
  scalar a2, b2;
  scalar a3, b3;
  (const) face vector alpha;
  (const) scalar lambda1;
  (const) scalar lambda2;
  (const) scalar lambda3;
  double tolerance;
  int nrelax1,nrelax2,nrelax3, minlevel;
  scalar * res1;
  scalar * res2;
  scalar * res3;
};

struct MGSolve1_phase3 {
  scalar * a1, * b1;
  double (* residual1) (scalar * a1, scalar * b1, scalar * res1,
		       void * data1);
  void (* relax1) (scalar * da1, scalar * res1, int depth1, 
		  void * data1);
  scalar * a2, * b2;
  double (* residual2) (scalar * a2, scalar * b2, scalar * res2,
		       void * data2);
  void (* relax2) (scalar * da2, scalar * res2, int depth2, 
		  void * data2);
  scalar * a3, * b3;
  double (* residual3) (scalar * a3, scalar * b3, scalar * res3,
		       void * data3);
  void (* relax3) (scalar * da3, scalar * res3, int depth3, 
		  void * data3);
  void * data;
  
  int nrelax1,nrelax2,nrelax3;
  scalar * res1;
  scalar * res2;
  scalar * res3;
  int minlevel;
  double tolerance;
};


void mg_cycle1_3phase (scalar * a1, scalar * res1, scalar * da1,
           scalar * a2, scalar * res2, scalar * da2,
           scalar * a3, scalar * res3, scalar * da3,
	       void (* relax1) (scalar * da1, scalar * res1, 
			       int depth1, void * data1),
           void (* relax2) (scalar * da2, scalar * res2, 
			       int depth2, void * data2),
           void (* relax3) (scalar * da3, scalar * res3, 
			       int depth3, void * data3),
	       void * data,int nrelax1,int nrelax2,
	       int nrelax3, int minlevel, int maxlevel)
{

  /**
  We first define the residual on all levels. */

  restriction (res1);
  restriction (res2);
  restriction (res3);

  /**
  We then proceed from the coarsest grid (*minlevel*) down to the
  finest grid. */

  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {

    /**
    On the coarsest grid, we take zero as initial guess. */

    if (l == minlevel){
      foreach_level_or_leaf (l){
        for (scalar s in da1){
            foreach_blockf (s){
                s[] = 0.;
            }
        }
      }

       foreach_level_or_leaf (l){
        for (scalar s in da2){
            foreach_blockf (s){
                s[] = 0.;
            }
        }
      }

       foreach_level_or_leaf (l){
        for (scalar s in da3){
            foreach_blockf (s){
                s[] = 0.;
            }
        }
      }
    /**
    On all other grids, we take as initial guess the approximate solution
    on the coarser grid bilinearly interpolated onto the current grid. */

    }else{
      foreach_level (l){
            for (scalar s in da1){
                foreach_blockf (s){
                    s[] = bilinear (point, s);
                }
            }
      }

      foreach_level (l){
            for (scalar s in da2){
                foreach_blockf (s){
                    s[] = bilinear (point, s);
                }
            }
      }

      foreach_level (l){
            for (scalar s in da3){
                foreach_blockf (s){
                    s[] = bilinear (point, s);
                }
            }
      }
    }
    /**
    We then apply homogeneous boundary conditions and do several
    iterations of the relaxation function to refine the initial guess. */

    boundary_level (da1, l);
    boundary_level (da2, l);
    boundary_level (da3, l);



    for (int i = 0; i < nrelax1; i++) {
      relax1 (da1, res1, l, data);
      boundary_level (da1, l);
    }
    for (int i = 0; i < nrelax2; i++) {
      relax1 (da2, res2, l, data);
      boundary_level (da2, l);
    }
    for (int i = 0; i < nrelax3; i++) {
      relax1 (da3, res3, l, data);
      boundary_level (da3, l);
    }
  }

  /**
  And finally we apply the resulting correction to *a*. */

  foreach() {
    scalar s1, ds1;
    for (s1, ds1 in a1, da1)
      foreach_blockf (s1)
	s1[] += ds1[];

    scalar s2, ds2;
    for (s2, ds2 in a2, da2)
      foreach_blockf (s2)
	s2[] += ds2[];

    scalar s3, ds3;
    for (s3, ds3 in a3, da3)
      foreach_blockf (s3)
	s3[] += ds3[];
  }
}


mgstats1_3phase mg_solve1_3phase (struct MGSolve1_3phase p)
{

  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  scalar * da1 = list_clone (p.a1), * res1 = p.res1;
  if (!res1)
  res1 = list_clone (p.b1);

    scalar * da2 = list_clone (p.a2), * res2 = p.res2;
  if (!res2)
  res2 = list_clone (p.b2);

    scalar * da3 = list_clone (p.a3), * res3 = p.res3;
  if (!res3)
  res3 = list_clone (p.b3);

  /**
  The boundary conditions for the correction fields are the
  *homogeneous* equivalent of the boundary conditions applied to
  *a*. */

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da1)
      s.boundary[b] = s.boundary_homogeneous[b];

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da2)
      s.boundary[b] = s.boundary_homogeneous[b];

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da3)
      s.boundary[b] = s.boundary_homogeneous[b];
  
  /**
  We initialise the structure storing convergence statistics. */

  mgstats1_3phase s = {0};

  double sum1 = 0.;
  foreach (reduction(+:sum1))
    for (scalar s in p.b1)
      sum1 += s[];
  s.sum1 = sum1;
  s.nrelax1 = p.nrelax1 > 0 ? p.nrelax1 : 4;

  double sum2 = 0.;
  foreach (reduction(+:sum2))
    for (scalar s in p.b2)
      sum2 += s[];
  s.sum2 = sum2;
  s.nrelax2 = p.nrelax2 > 0 ? p.nrelax2 : 4;

  double sum3 = 0.;
  foreach (reduction(+:sum3))
    for (scalar s in p.b3)
      sum3 += s[];
  s.sum3 = sum3;
  s.nrelax3 = p.nrelax3 > 0 ? p.nrelax3 : 4;
  
  /**
  Here we compute the initial residual field and its maximum. */

  double resb1;
  resb1 = s.resb1 = s.resa1 = p.residual1 (p.a1, p.b1, res1, p.data);

  double resb2;
  resb2 = s.resb2 = s.resa2 = p.residual2 (p.a2, p.b2, res2, p.data);

  double resb3;
  resb3 = s.resb3 = s.resa3 = p.residual3 (p.a3, p.b3, res3, p.data);

  /**
  We then iterate until convergence or until *NITERMAX* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE1_3phase;
  for (s.i = 0;
       s.i < NITERMAX1_3phase && (s.i < NITERMIN1_3phase || 
       ((s.resa1 > p.tolerance) && (s.resa2 > p.tolerance) && (s.resa3 > p.tolerance)) );
       s.i++) {
    mg_cycle1_3phase (p.a1, res1, da1, 
    p.a2, res2, da2,
    p.a3, res3, da3, p.relax1, p.relax2,p.relax3, 
    p.data,
	s.nrelax1,s.nrelax2,s.nrelax3,
	p.minlevel,
	grid->maxdepth);




    s.resa1 = p.residual1 (p.a1, p.b1, res1, p.data);
    s.resa2 = p.residual2 (p.a2, p.b2, res2, p.data);
    s.resa3 = p.residual3 (p.a3, p.b3, res3, p.data);

    /**
    We tune the number of relaxations so that the residual is reduced
    by between 2 and 20 for each cycle. This is particularly useful
    for stiff systems which may require a larger number of relaxations
    on the finest grid. */

#if 1
    if (s.resa1 > p.tolerance) {
      if (resb1/s.resa1 < 1.2 && s.nrelax1 < 100)
	s.nrelax1++;
      else if (resb1/s.resa1 > 10 && s.nrelax1 > 2)
	s.nrelax1--;
    }

     if (s.resa2 > p.tolerance) {
      if (resb2/s.resa2 < 1.2 && s.nrelax2 < 100)
	s.nrelax2++;
      else if (resb2/s.resa2 > 10 && s.nrelax2 > 2)
	s.nrelax2--;
    }

     if (s.resa3 > p.tolerance) {
      if (resb3/s.resa3 < 1.2 && s.nrelax3 < 100)
	s.nrelax3++;
      else if (resb3/s.resa3 > 10 && s.nrelax3 > 2)
	s.nrelax3--;
    }
#else
    if (s.resa == resb) /* convergence has stopped!! */
      break;
    if (s.resa > resb/1.1 && p.minlevel < grid->maxdepth)
      p.minlevel++;
#endif

    resb1 = s.resa1;
    resb2 = s.resa2;
    resb3 = s.resa3;
  }
  s.minlevel = p.minlevel;
  
//trace the final res[]
// for(scalar s in res){
//     if(phase_flag3==0){
//       foreach(){
//           resg[] = s[];
//       }
//     }else{
//         foreach(){
//           resl[] = s[];
//       }
//     }
// }
  /**
  If we have not satisfied the tolerance, we warn the user. */

  if ((s.resa1 > p.tolerance) || (s.resa2 > p.tolerance) || (s.resa3 > p.tolerance)) {
    fprintf (ferr, 
	     "WARNING: convergence for T not reached after %d iterations\n", 
	     s.i), fflush (ferr);

    if((s.resa1 > p.tolerance)) {
            scalar v = p.a1[0];
            fprintf (ferr, 
                "WARNING: convergence for %s not reached after %d iterations\n"
                "  res: %g sum: %g nrelax: %d\n", v.name,
                s.i, s.resa1, s.sum1, s.nrelax1), fflush (ferr);
    }
    if((s.resa2 > p.tolerance)) {
            scalar v = p.a2[0];
            fprintf (ferr, 
                "WARNING: convergence for %s not reached after %d iterations\n"
                "  res: %g sum: %g nrelax: %d\n", v.name,
                s.i, s.resa2, s.sum2, s.nrelax2), fflush (ferr);
    }
    if((s.resa3 > p.tolerance)) {
            scalar v = p.a3[0];
            fprintf (ferr, 
                "WARNING: convergence for %s not reached after %d iterations\n"
                "  res: %g sum: %g nrelax: %d\n", v.name,
                s.i, s.resa3, s.sum3, s.nrelax3), fflush (ferr);
    }
  }
    
  /**
  We deallocate the residual and correction fields and free the lists. */

  if ((!p.res1) && (!p.res2) && (!p.res3)){
    delete (res1), free (res1);
     delete (res2), free (res2);
      delete (res3), free (res3);
  }
  delete (da), free (da);

  return s;
}



/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

extern bool poisson_check;
static double residual1 (scalar * al, scalar * bl, scalar * resl, void * data) //for gas
{ 
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson1_3phase * p = (struct Poisson1_3phase *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda1;
#if EMBED3
  bool embedded = (a.boundary[embed] != symmetry); 
#endif
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */

  get_solidsurfaceT();
  get_solidflux():


  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*face_gradient_x (a, 0);
  foreach (reduction(max:maxres)) {
    //bool ff_flag = ((((ff[]>0.5 ) && (phase_flag3==1)) || ((ff[]<0.5 ) && (phase_flag3==0))) && (fabs(ff[]-0.5)>1e-6)); //ff not changed
    //bool ff_flag = ((((ff[]>0.5 ) && (phase_flag3==1)) || ((ff[]<0.5 ) && (phase_flag3==0))) && (fabs(ff[]-0.5)>1e-6)); //ff not changed
    bool solid_flag = false;
    if(css_test3_n[]>=1 && ff[]<0.5){
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


mgstats1_3phase poisson1_3phase (struct Poisson1_3phase p)
{

   p.minlevel = minl;
  /**
  If $\alpha$ or $\lambda$ are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda1.i)
    p.lambda1 = zeroc;
   if (!p.lambda2.i)
    p.lambda2 = zeroc;
 if (!p.lambda3.i)
    p.lambda3 = zeroc;
  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  face vector alpha = p.alpha;
  scalar lambda1 = p.lambda1;
  scalar lambda2 = p.lambda2;
  scalar lambda3 = p.lambda3;
  restriction ({alpha,lambda1,lambda2,lambda3});

  double defaultol = TOLERANCE1_3phase;
  if (p.tolerance)
    TOLERANCE1_3phase = p.tolerance;

  scalar a1 = p.a1, b = p.b1;
  scalar a2 = p.a2, b = p.b2;
  scalar a3 = p.a3, b = p.b3;
  mgstats1_3phase s = mg_solve1_3phase ({a1}, {b1},residual1, relax1,
            {a2}, {b2},residual2, relax2,
            {a3}, {b3},residual3, relax3,
			&p, p.nrelax11, p.nrelax12, p.nrelax13, 
            p.res1,p.res2,p.res3, minlevel = max(1, p.minlevel));

//   mgstats1 s = mg_solve1 ({a}, {b}, residual1, relax1,
// 			&p, p.nrelax, p.res, minlevel = max(1, p.minlevel));

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE1_3phase = defaultol;

  return s;
}
