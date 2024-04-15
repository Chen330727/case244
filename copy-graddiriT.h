

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

//extern scalar css_test; //scalar cs[];
//extern face vector fss_test; //face vector fs[];


// foreach_dimension()
// double dirichlet_gradient2_x (Point point, scalar s, scalar css_test,
// 					   coord n, coord p, double bc,
// 					   double * coef , face vector fss_test)
foreach_dimension()
static inline double dirichlet_gradient2_x (Point point, scalar s, scalar css_test,
					   coord n, coord p, double bc,
					   double * coef , face vector fss_test)
{
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fss_test.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fss_test.x[i + (i < 0),j] && fss_test.y[i,j] && fss_test.y[i,j+1] &&
	  css_test[i,j-1] && css_test[i,j] && css_test[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fss_test.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fss_test.y[i,j,k+m] || !fss_test.y[i,j+1,k+m] ||
	    !fss_test.z[i,j+m,k] || !fss_test.z[i,j+m,k+1] ||
	    !css_test[i,j+m,k-1] || !css_test[i,j+m,k] || !css_test[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
  if (v[0] == nodata) {

    /**
    This is a degenerate case, we use the boundary value and the
    cell-center value to define the gradient. */
	
    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);

    if(fabs(bc/(d[0]*Delta))>60.0){
       // char name151[80];
       // sprintf(name151,"gradient-check-test.dat");
       // FILE * fp151 = fopen(name151,"a");
       // fprintf(fp151,"here!!!!!\n");
       // fclose(fp151);

    }
    return bc/(d[0]*Delta);
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
  *coef = 0.;
  if (v[1] != nodata){ // third-order gradient
    if(fabs((d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta))>60.0){
        double temp9 = fabs((d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta));
       // char name151[80];
       // sprintf(name151,"gradient-check-test.dat");
       // FILE * fp151 = fopen(name151,"a");
       // fprintf(fp151,"there1!!!!! d0=%g d1=%g v0=%g v1=%g grad=%g\n",d[0],d[1],v[0],v[1],temp9);
       // fclose(fp151);
    }else{
        double temp9 = fabs((d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta));
      //  char name151[80];
      //  sprintf(name151,"gradient-check-test.dat");
      //  FILE * fp151 = fopen(name151,"a");
      //  fprintf(fp151,"nice1!!!!! d0=%g d1=%g v0=%g v1=%g grad=%g\n",d[0],d[1],v[0],v[1],temp9);
      //  fclose(fp151);

    }
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  
  }
  if(fabs(bc - v[0])/(d[0]*Delta)>60.0){     
      //  char name151[80];
      //  sprintf(name151,"gradient-check-test.dat");
      //  FILE * fp151 = fopen(name151,"a");
      //  fprintf(fp151,"there2!!!!!\n");
      //  fclose(fp151);
  }
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double dirichlet_gradient2 (Point point, scalar s, scalar css_test,
			   coord n, coord p, double bc, double * coef, face vector fss_test)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient2_x (point, s, css_test, n, p, bc, coef, fss_test);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient2_x (point, s, css_test, n, p, bc, coef, fss_test);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient2_y (point, s, css_test, n, p, bc, coef, fss_test);
  return dirichlet_gradient2_z (point, s, css_test, n, p, bc, coef, fss_test);
#endif // dimension == 3
  return nodata;
}


