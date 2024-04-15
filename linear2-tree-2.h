#include "curvature.h"
coord interface_normal3 (Point point, scalar c, vector h)
{
  coord n;
  if (!h.x.i || (n = height_normal (point, c, h)).x == nodata)
    n = mycs (point, c);
  return n;
}


attribute {
  scalar ff6;
}

void restriction_zero (Point point, scalar s)
{
  scalar ff6 = s.ff6;
  double sum = 0.;
  double weight= 0.;
  foreach_child(){
    if(ff6[]>0.0){
        sum += s[]*ff6[];
        weight += ff6[];
    }    
  }
  // s[] = sum/(1 << dimension)/(cm[] + 1e-30);
  if(weight>0.0){
    s[] = sum/weight;
  }
}


void restriction_flux_sum (Point point, scalar s)
{
  double sum = 0.;
  foreach_child(){
        sum += s[];
  }
  s[] = sum;
}

double bilinear_no_cs3 (Point point, scalar s)
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
    return 0.0;

}

static inline double bilinear_embed_css_test3 (Point point, scalar s)
{
  if (!coarse(css_test3) || !coarse(css_test3,child.x))
    return coarse(s);
  #if dimension >= 2
  if (!coarse(css_test3,0,child.y) || !coarse(css_test3,child.x,child.y))
    return coarse(s);
  #endif
  #if dimension >= 3
  if (!coarse(css_test3,0,0,child.z) || !coarse(css_test3,child.x,0,child.z) ||
      !coarse(css_test3,0,child.y,child.z) ||
      !coarse(css_test3,child.x,child.y,child.z))
    return coarse(s);  
  #endif
  return bilinear_no_cs3 (point, s);
}

static inline double bilinear_embed_css_test2 (Point point, scalar s)
{
  if (!coarse(css_test2) || !coarse(css_test2,child.x))
    return coarse(s);
  #if dimension >= 2
  if (!coarse(css_test2,0,child.y) || !coarse(css_test2,child.x,child.y))
    return coarse(s);
  #endif
  #if dimension >= 3
  if (!coarse(css_test2,0,0,child.z) || !coarse(css_test2,child.x,0,child.z) ||
      !coarse(css_test2,0,child.y,child.z) ||
      !coarse(css_test2,child.x,child.y,child.z))
    return coarse(s);  
  #endif
  return bilinear_no_cs3 (point, s);
}

static inline double bilinear_embed_css_test (Point point, scalar s)
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
  return bilinear_no_cs3 (point, s);
}

  
double bilinear_embed_ff6 (Point point, scalar s)
{
  scalar ff6 = s.ff6;
  if (!coarse(ff6) || !coarse(ff6,child.x)){
     if(ff6[]>=0.5){ 
        return coarse(s);
     }else{ 
        return 0.0;
     }
  }
  #if dimension >= 2
  if (!coarse(ff6,0,child.y) || !coarse(ff6,child.x,child.y)){
    if(ff6[]>=0.5){ 
        return coarse(s);
     }else{ 
        return 0.0;
     }
  }
  #endif
  #if dimension >= 3
  if (!coarse(ff6,0,0,child.z) || !coarse(ff6,child.x,0,child.z) ||
      !coarse(ff6,0,child.y,child.z) ||
      !coarse(ff6,child.x,child.y,child.z)){
          if(ff6[]>=0.5){ 
              return coarse(s);
          }else{ 
              return 0.0;
          } 
      }
  #endif
  return bilinear_no_cs3 (point, s);
}