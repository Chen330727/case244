// void get_solidsurfaceT(){
//     1. search interface_cell

//     " we do not knoe the T and gradientT and direction at the solid surface, so maybe we could first approximate that the directoin of gradient T of the surface is
//       approximate to the direction of gradient T of the nerest cell", then using the direction to get a line and T2 and T2 and distance.
//      when solid surface is far away from the gas liquid interface, the difference is small
//                 3. get Tg2 Tg1 Tboundary Ts1 Ts2, get Tl2 Tl1 Tboundary Ts1 Ts2
//                 3.1 first use embed gradient to get T2 T1, for g and l and s
//                 3.2 if no T2 and T1, then 
      
//     4. Get solidsurfaceT[] from kf*gradient(Tf) 
// }

// #include "fractions.h"

extern scalar merge_to_me_s_c,merge_to_me_s_energy;
extern vector merge_to_me_s_position;

extern scalar merge_to_me_l_c,merge_to_me_l_energy;
extern vector merge_to_me_l_position;

extern scalar merge_to_me_g_c,merge_to_me_g_energy;
extern vector merge_to_me_g_position;
extern scalar intersect_ture;
extern scalar T;

extern bool Rcc_flag_l, Rcc_flag_g; // for heat resistance;
extern scalar aiml_s,aimg_s;
extern double Rcc;
extern double delta_min;
// extern vector number_s_f;

#define EPS32 0.0
extern face vector modphase0,modphase1,modphase_s_0;
bool structure_normal_flag = true;
struct Topo_m3 {
    scalar topo_mask;
    scalar ff;
    int level_interface;
};

void get_topomask_soilid(struct Topo_m3 q) {
     scalar topo_mask = q.topo_mask;
     scalar ff=q.ff;
     int level_interface=q.level_interface;

     foreach(){
        int phase_sign = (ff[]>=0.5-EPS32) ? 1 : -1;
	topo_mask[] = 3*phase_sign;
        //if(ff[]<1.0-EPS32 && ff[]>EPS32){
         if(ff[]<1.0-EPS32 && ff[]>EPS32){
 		topo_mask[] = 0;
        }
     }
     boundary({topo_mask});
     foreach(){
       // if(level==level_interface){
         if(level==level_interface  ){
            bool is1= false;
            foreach_dimension(){
               int temp1=topo_mask[1];
               int temp2=topo_mask[-1];
            if(temp1==0 || temp2==0){
               is1 = true;
               }
            }
            int temp3=topo_mask[];
            if (is1 && temp3!=0){
              // topo_mask[] = (ff[]>=0.5-EPS32) ? 1 : -1 ;
	      if(temp3==3){
	          topo_mask[] = 1;
	      }else if(topo_mask[]==-3){
	          topo_mask[] = -1;
	      }
            }
        }
        
     } 
  for(int phase=0;phase<=1;phase++){
     foreach(){
        // if(level==level_interface){
          if(level==level_interface ){
               int phase_sign = ff[]>=0.5-EPS32 ? 1 : -1;
            bool is1= false;
               foreach_dimension(){
                  int temp1=topo_mask[1];
                  int temp2=topo_mask[-1];
               if( temp1==(2*phase-1) || temp2==(2*phase-1)){
                  is1 = true;
                  }
               }
               int temp3 = topo_mask[];
               if (is1 && temp3==3*phase_sign){
                  topo_mask[] = 2*(2*phase-1) ;
               }
         }
     }
  }


       
}

void gradient_s (scalar * f, vector * g, scalar c, scalar topo_mask, bool inverse)
{
  //c is css_test3, inverse is for css_test3(c)
  //topo_mask is for c;
  assert (list_len(f) == vectors_len(g));
  double cmin = 0.5;
  foreach() {
    scalar s; vector v;
    for (s,v in f,g) {
        foreach_dimension() {
            v.x[] = 0.0;
            double cl = c[-1], cc = c[], cr = c[1];
            if (inverse)
                        cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
            if(fabs(topo_mask[])<3 && cc>=0.5){
                if (cr >= cmin) {
                    if (cl >= cmin){
                        v.x[] =(s[1] - s[-1])/(2.0*Delta);// (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
                    }else{
                        
                        v.x[] = (s[1] - s[])/Delta;
                    }
                }else if (cl >= cmin)
                    v.x[] =  (s[] - s[-1])/Delta;
             }
        }
    }
  }
}
extern double Tsat00;

static double gradient_distacew(double T1, double T2, double T3, double a, double b) {
    double F = (T2 - T1) / a;
    double B = (T3 - T2) / b;

    double dT_dx = F + ((B - F) * 2 / (a + b) * a / 2);

    return dT_dx;
}


static double gradient_quadratic(double T1, double T2, double T3, double a, double b) {
    double A = ((b*T1) - (a*T3) + ((a+b)*T2)) / (a*b*(a+b));
    double B = ((T3-T2)/b - (T2-T1)/a) / (b-a);
    
    double dT_dx = 2*A*0 + B; // At x = 0

    return dT_dx;
}

void gradient_ff (scalar * f, vector * g, scalar c, scalar c2, scalar topo_mask, scalar topo_mask_s, bool inverse)
{
  //c is ff, inverse is for css_test3(c)
  //c2 is css_test3
  //topo_mask is for c;
  assert (list_len(f) == vectors_len(g));
  double cmin = 0.5;
  foreach() {
    scalar s; vector v;
    for (s,v in f,g) {
        foreach_dimension() {
            v.x[] = 0.0;
            double cl = c[-1], cc = c[], cr = c[1];
            bool is_solid_l=c2[-1]>0.5;
            bool is_solid_c=c2[]>0.5;
            bool is_solid_r=c2[1]>0.5; 
            if (inverse)
                        cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
            if(!is_solid_c){ //not solid
                 if((!is_solid_l) && (!is_solid_r)){ //all not solid
                    // if(fabs(topo_mask[])<3 && cc>=0.5){
                      if(((fabs(topo_mask[])<3) || (topo_mask_s[]<=0 && level==level_interface)) && cc>=0.5){
                        if (cr >= cmin) {
                            if (cl >= cmin){
                                        v.x[] =(s[1] - s[-1])/(2.0*Delta);// (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
                            }else{
                              double leftd=1.0;
                              if(phase_flag3==1){
                                    leftd  = modphase1.x[];
                              }else{
                                    leftd  = modphase0.x[];
                              }
                              if(leftd<1e-2){
                                  leftd=1e-2;
                              }
                                v.x[] = gradient_distacew (Tsat00,s[],s[1],leftd,1)/Delta;//(s[1] - s[])/Delta;
                            }
                        }else if (cl >= cmin){
                           double rightd=1.0;
                              if(phase_flag3==1){
                                    rightd  = modphase1.x[1];
                              }else{
                                    rightd  = modphase0.x[1];
                              }
                              if(rightd<1e-2){
                                  rightd=1e-2;
                              }
                            v.x[] = gradient_distacew (s[-1],s[],Tsat00,1,rightd)/Delta;// (s[] - s[-1])/Delta;
                        }else{

                        }
                    }
                }else if(is_solid_l && (!is_solid_r)){// solid not_s not_s
                      // if(fabs(topo_mask[])<3 && cc>=0.5){
                        if(((fabs(topo_mask[])<3) || (topo_mask_s[]<=0 && level==level_interface)) && cc>=0.5){
                             if (cr >= cmin) {
                                if (cl >= cmin){
                                            //3.1
                                            v.x[] = (s[1]-s[])/Delta;
                                }else{
                                    double leftd_f=1.0;
                                    if(phase_flag3==1){
                                            leftd_f  = modphase1.x[];
                                    }else{
                                            leftd_f  = modphase0.x[];
                                    }
                                    double leftd_s=1.0;
                                    leftd_s = modphase_s_0.x[];
                                    if(leftd_f<leftd_s){
                                        if(leftd_f<1e-2){
                                            leftd_f=1e-2;
                                        }
                                        //3.2
                                        v.x[] = gradient_distacew (Tsat00,s[],s[1],leftd_f,1)/Delta;//(s[1] - s[])/Delta;
                                    }else{
                                        //3.3
                                        v.x[] = (s[1]-s[])/Delta;
                                    }
                                }
                        }else if (cl >= cmin){
                                 double rightd_f=1.0;
                                    if(phase_flag3==1){
                                            rightd_f  = modphase1.x[1];
                                    }else{
                                            rightd_f  = modphase0.x[1];
                                    }
                                 if(rightd_f < 1.0){
                                        if(rightd_f<1e-2){
                                            rightd_f = 1e-2;
                                        }
                                        //3.4
                                        v.x[] = (Tsat00-s[])/(rightd_f*Delta); 
                                 }else{
                                        //3.5
                                        v.x[] = (s[1]-s[])/Delta;
                                 }
                        }else{ 
                                  double leftd_f=1.0;
                                    if(phase_flag3==1){
                                            leftd_f  = modphase1.x[];
                                    }else{
                                            leftd_f  = modphase0.x[];
                                    }
                                    double leftd_s=1.0;
                                    leftd_s = modphase_s_0.x[];
                                  double rightd_f=1.0;
                                    if(phase_flag3==1){
                                            rightd_f  = modphase1.x[1];
                                    }else{
                                            rightd_f  = modphase0.x[1];
                                    }
                                    double rightd_s = 1.0;
                                    // rightd_s = modphase_s_0.x[1];
                                bool left_flag,right_flag;
                                left_flag = leftd_f < leftd_s;
                                right_flag=rightd_f < rightd_s;
                                if(leftd_f<1e-2){
                                    leftd_f = 1e-2;
                                }
                                if(rightd_f<1e-2){
                                    rightd_f = 1e-2;
                                }
                                if(left_flag && right_flag){
                                       //3.6
                                        v.x[] = gradient_distacew (Tsat00,s[],Tsat00,leftd_f,rightd_f)/Delta;//(s[1] - s[])/Delta;
                                }else if(left_flag && (!right_flag)){
                                       //3.7
                                        v.x[] = gradient_distacew (Tsat00,s[],s[1],leftd_f,1)/Delta;//(s[1] - s[])/Delta;
                                }else if((!left_flag) && right_flag){
                                       //3.8
                                       v.x[] = (Tsat00-s[])/(rightd_f*Delta);
                                }else{
                                       //3.9
                                       v.x[] = (s[1]-s[])/Delta;
                                }
                        }
                      }   
                }else if((!is_solid_l) && (is_solid_r)){ //not_s not_s solid
                      // if(fabs(topo_mask[])<3 && cc>=0.5){
                        if(((fabs(topo_mask[])<3) || (topo_mask_s[]<=0 && level==level_interface)) && cc>=0.5){
                             if (cr >= cmin) {
                                if (cl >= cmin){
                                    //3.10
                                     v.x[] = (s[]-s[-1])/Delta;  
                                }else{
                                     double leftd_f=1.0;
                                    if(phase_flag3==1){
                                            leftd_f  = modphase1.x[];
                                    }else{
                                            leftd_f  = modphase0.x[];
                                    }
                                    if(leftd_f < 1.0){
                                            if(leftd_f<1e-2){
                                                leftd_f = 1e-2;
                                            }
                                            //3.13
                                            v.x[] = (s[]-Tsat00)/(leftd_f*Delta); 
                                    }else{
                                            //3.14
                                            v.x[] = (s[]-s[-1])/Delta;
                                    }
                                     
                                }
                        }else if (cl >= cmin){
                                 double rightd_f=1.0;
                                    if(phase_flag3==1){
                                            rightd_f  = modphase1.x[1];
                                    }else{
                                            rightd_f  = modphase0.x[1];
                                    }
                                    double rightd_s = 1.0;
                                    rightd_s = modphase_s_0.x[1];
                                    if(rightd_f < rightd_s){
                                        if(rightd_f<1e-2){
                                            rightd_f = 1e-2;
                                        }
                                        //3.11
                                         v.x[] = gradient_distacew (s[-1],s[],Tsat00,1,rightd_f)/Delta;
                                    }else{
                                        //3.12
                                         v.x[] = (s[]-s[-1])/Delta;
                                    }  
                        }else{ 
                                double leftd_f=1.0;
                                    if(phase_flag3==1){
                                            leftd_f  = modphase1.x[];
                                    }else{
                                            leftd_f  = modphase0.x[];
                                    }
                                    double leftd_s=1.0;
                                    // leftd_s = modphase_s_0.x[];
                                  double rightd_f=1.0;
                                    if(phase_flag3==1){
                                            rightd_f  = modphase1.x[1];
                                    }else{
                                            rightd_f  = modphase0.x[1];
                                    }
                                    double rightd_s = 1.0;
                                    rightd_s = modphase_s_0.x[1];
                                    bool left_flag,right_flag;
                                left_flag = leftd_f < leftd_s;
                                right_flag=rightd_f < rightd_s;
                                if(leftd_f<1e-2){
                                    leftd_f = 1e-2;
                                }
                                if(rightd_f<1e-2){
                                    rightd_f = 1e-2;
                                }
                                if(left_flag && right_flag){
                                       //3.16
                                        v.x[] = gradient_distacew (Tsat00,s[],Tsat00,leftd_f,rightd_f)/Delta;//(s[1] - s[])/Delta;
                                }else if(left_flag && (!right_flag)){
                                       //3.17
                                        v.x[] = (s[]-Tsat00)/(leftd_f*Delta);
                                }else if((!left_flag) && right_flag){
                                       //3.18
                                       v.x[] = gradient_distacew (s[-1],s[],Tsat00,1,rightd_f)/Delta;
                                }else{
                                       //3.19
                                       v.x[] = (s[]-s[-1])/Delta;
                                }
                        }
                      }   

                }else if(is_solid_l && is_solid_r){ //solid not_s solid
                      // if(fabs(topo_mask[])<3 && cc>=0.5){
                       if(((fabs(topo_mask[])<3) || (topo_mask_s[]<=0 && level==level_interface)) && cc>=0.5){
                             if (cr >= cmin) {
                                if (cl >= cmin){
                                          //3.20
                                          //v.x[] = 0.0; 
                                }else{
                                    double leftd_f=1.0;
                                    if(phase_flag3==1){
                                            leftd_f  = modphase1.x[];
                                    }else{
                                            leftd_f  = modphase0.x[];
                                    }
                                    double leftd_s=1.0;
                                    leftd_s = modphase_s_0.x[];
                                    if(leftd_f<leftd_s){
                                            if(leftd_f<1e-2){
                                                leftd_f=1e-2;
                                            }
                                            //3.21
                                            v.x[] = (s[]-Tsat00)/(leftd_f*Delta);
                                    }else{
                                            //3.22
                                            //v.x[] = 0.0;
                                    }
                                }
                        }else if (cl >= cmin){
                                  double rightd_f=1.0;
                                    if(phase_flag3==1){
                                            rightd_f  = modphase1.x[1];
                                    }else{
                                            rightd_f  = modphase0.x[1];
                                    }
                                    double rightd_s = 1.0;
                                    rightd_s = modphase_s_0.x[1];
                                if(rightd_f<rightd_s){
                                        if(rightd_f<1e-2){
                                             //3.23
                                             v.x[] = (Tsat00 - s[])/(rightd_f*Delta);
                                        }
                                }else{
                                         //3.24
                                         // v.x[] = 0.0;
                                }
                        }else{ 
                                 double leftd_f=1.0;
                                    if(phase_flag3==1){
                                            leftd_f  = modphase1.x[];
                                    }else{
                                            leftd_f  = modphase0.x[];
                                    }
                                    double leftd_s=1.0;
                                    leftd_s = modphase_s_0.x[];
                                     double rightd_f=1.0;
                                    if(phase_flag3==1){
                                            rightd_f  = modphase1.x[1];
                                    }else{
                                            rightd_f  = modphase0.x[1];
                                    }
                                    double rightd_s = 1.0;
                                    rightd_s = modphase_s_0.x[1];
                                    bool left_flag,right_flag;
                                left_flag = leftd_f < leftd_s;
                                right_flag=rightd_f < rightd_s;
                                if(leftd_f<1e-2){
                                    leftd_f = 1e-2;
                                }
                                if(rightd_f<1e-2){
                                    rightd_f = 1e-2;
                                }
                                 if(left_flag && right_flag){
                                       //3.25
                                        v.x[] = gradient_distacew (Tsat00,s[],Tsat00,leftd_f,rightd_f)/Delta;//(s[1] - s[])/Delta;
                                }else if(left_flag && (!right_flag)){
                                       //3.26
                                        v.x[] = (s[]-Tsat00)/(leftd_f*Delta);
                                }else if((!left_flag) && right_flag){
                                       //3.27
                                       v.x[] = (Tsat00-s[])/(rightd_f*Delta);
                                }else{
                                       //3.28
                                       // v.x[] =0;
                                }
                        }
                      }   

                }
            }
        }
    }
  }
}

// static double vof_concentration_gradient_f2 (Point point, scalar c, scalar t,bool inverse)
// {
//   //static const double cmin = 0.5;
//   static const double cmin = 0.5;
//   double cl = c[-1], cc = c[], cr = c[1];
//   if (inverse)
//     cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
//   if (cc >= cmin) {
//     if(t.gradient != zero){
//           // printf("t.gradient !=zero\n");
//           if (cr >= cmin) {
//             if (cl >= cmin)
//         return t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
//             else
//         return (t[1]/cr - t[]/cc)/Delta;
//           }
//           else if (cl >= cmin)
//             return (t[]/cc - t[-1]/cl)/Delta;
//     }else{
//         // printf("t.gradient ==zero\n");
//     }
//   }
//   return 0.;
// }

void direction_ff(vector T_g, scalar f, scalar f2, scalar topo_mask, scalar topo_mask_s, vector direction,bool flag4){
    foreach(){
        foreach_dimension(){
            direction.x[] = 0.0;
        }
    }

    foreach(){
        if((f[]>0.0 && f[]<1.0) || ((topo_mask_s[]<=0 && topo_mask_s[]>=-2) &&
         ((topo_mask[]>=0 && (!flag4)) || (topo_mask[]<=0 && (flag4)) ) )){
            coord nn;
            nn.x = T_g.x[];
            nn.y = T_g.y[];
            double value;
            value = sqrt(sq(nn.x)+sq(nn.y));
            bool flag=false;
            if(fabs(value)>1e-20){ //check normal is not 0
                normalize(&nn);
                direction.x[] = nn.x;
                direction.y[] = nn.y;
                flag=true;
            }
            if(!flag){
                Point me = point;
                bool flag2=false;
                double distance1=HUGE;
                coord direction_temp;
                coord temp1;
                temp1.x = x;
                temp1.y = y;
                foreach_neighbor(2){ //this is not good based on distance; direction should based on temperature and distance
                    if( (!(me.i == point.i && me.j == point.j)) && (!is_boundary(cell)) && (((!flag4) && (f[]>=0.5)) || ((flag4) && (f[]<0.5))) && (f2[]<0.5)){
                        coord temp2;
                        temp2.x = x;
                        temp2.y = y;
                        double value2 = sqrt(sq(temp1.x-temp2.x)+sq(temp1.y-temp2.y));
                        if(value2<distance1){
                            coord temp3;
                            temp3.x = T_g.x[];
                            temp3.y = T_g.y[];
                            double value3 =  sqrt(sq(temp3.x)+sq(temp3.y));
                            if(value3>1e-20){ //check normal is not 0
                                  normalize(&temp3);
                                  distance1 = value2;
                                  direction_temp.x = temp3.x;
                                  direction_temp.y = temp3.y;
                                  flag2=true;
                            }
                        }
                    }
                }
                if(flag2){
                    direction.x[] = direction_temp.x;
                    direction.y[] = direction_temp.y;
                }
            }
        }
    }
}

void direction_s(vector T_g, scalar f, vector direction,bool flag4){
    foreach(){
        foreach_dimension(){
            direction.x[] = 0.0;
        }
    }

    foreach(){
        if(f[]>0.0 && f[]<1.0){
            coord nn;
            nn.x = T_g.x[];
            nn.y = T_g.y[];
            double value;
            value = sqrt(sq(nn.x)+sq(nn.y));
            bool flag=false;
            if(fabs(value)>1e-20){ //check normal is not 0
                normalize(&nn);
                direction.x[] = nn.x;
                direction.y[] = nn.y;
                flag=true;
            }
            if(!flag){
                Point me = point;
                bool flag2=false;
                double distance1=HUGE;
                coord direction_temp;
                coord temp1;
                temp1.x = x;
                temp1.y = y;
                foreach_neighbor(2){
                    if((!(me.i==point.i && me.j==point.j)) && (!is_boundary(cell)) && (((!flag4) && (f[]>=0.5)) || ((flag4) && (f[]<0.5)))){
                        coord temp2;
                        temp2.x = x;
                        temp2.y = y;
                        double value2 = sqrt(sq(temp1.x-temp2.x)+sq(temp1.y-temp2.y));
                        if(value2<distance1){
                            coord temp3;
                            temp3.x = T_g.x[];
                            temp3.y = T_g.y[];
                            double value3 =  sqrt(sq(temp3.x)+sq(temp3.y));
                            if(value3>1e-20){ //check normal is not 0
                                  normalize(&temp3);
                                  distance1 = value2;
                                  direction_temp.x = temp3.x;
                                  direction_temp.y = temp3.y;
                                  flag2=true;
                            }
                        }
                    }
                }
                if(flag2){
                    direction.x[] = direction_temp.x;
                    direction.y[] = direction_temp.y;
                }
            }
        }
    }
}
extern scalar topo_mask,topo_mask_s;
extern int level_interface;
extern vector Tlff_g,Tgff_g,T_solid_g;
extern scalar ff,css_test3,Tl,Tg,Ts; 
extern scalar css_test,css_test2;
extern vector direction_Tfg, direction_Tfl, direction_Ts;
extern scalar aiml,aimg, T_modl,T_modg;
extern scalar flux_l,flux_g;
extern scalar areasg,areasl,arealg;


extern face vector fss_test3, fss_test3_n;
extern scalar css_test3_n;

// foreach_dimension()
// static inline double* gradient_conjugates_x (Point point, scalar s, scalar f, bool inverse_flag,
// 			   coord n, coord p, int * number)
foreach_dimension()
bool gradient_conjugates_x (Point point, scalar s, scalar f, bool inverse_flag,
			   coord  n, coord p, 
         double * T1, double * T2, double *a, double* b, int* number, int* sign_flag)
{
  // coord n,p;
  // foreach_dimension(){
  //   n.x=n_temp->x;
  //   p.x=p_temp->x;
  // }
  
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fss_test3.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      *sign_flag = sign(n.x);
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fss_test3.x[i + (i < 0),j] && fss_test3.y[i,j] && fss_test3.y[i,j+1] &&
	  css_test3[i,j-1] && css_test3[i,j] && css_test3[i,j+1]){
	        v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
      }
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fss_test3.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fss_test3.y[i,j,k+m] || !fss_test3.y[i,j+1,k+m] ||
	    !fss_test3.z[i,j+m,k] || !fss_test3.z[i,j+m,k+1] ||
	    !css_test3[i,j+m,k-1] || !css_test3[i,j+m,k] || !css_test3[i,j+m,k+1])
	  defined = false;
      if (defined){
            // bi-quadratic interpolation
            v[l] =
            quadratic (z,
                    quadratic (y1,
                        (s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
                    quadratic (y1,
                        (s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
                    quadratic (y1,
                        (s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
      }
#endif // dimension == 3
      else{
	        break;
      }
    }
  if (v[0] == nodata) {
    //   if(((!(inverse_flag)) && f[]>=0.5) || (inverse_flag && f[]<0.5)){
         if(1==1){
          /**
          This is a degenerate case, we use the boundary value and the
          cell-center value to define the gradient. */
        
          // d[0] = max(1e-1, fabs(p.x/n.x));
        //   d[0] = max(1e-1, fabs(p.x/n.x));
        //   d[0] = max(1e-1,sqrt(sq(p.x)+sq(p.y)));
        // d[0] = fabs(p.x*n.x + p.y*n.y);
        d[0] = max(1e-1,fabs(p.x*n.x + p.y*n.y));
          // *coef = - 1./(d[0]*Delta);
          // return bc/(d[0]*Delta);
          //  double *data_s = (double*) malloc(2 * sizeof(double));
          *T1=s[];//data_s[0]=s[];
          *a=fabs(d[0]);//data_s[1]=d[0];
          *number=1;
          return true;
          //  return data_s;
      }else{
          return false;
      }
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
//   *coef = 0.;
  if (v[1] != nodata){ // third-order gradient
        // if(!Rcc_flag){
         if(((!Rcc_flag_l) && (!inverse_flag)) || (inverse_flag)){
            // return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
            // double *data_s = (double*) malloc(4 * sizeof(double));
            //  *T1=v[0];//data_s[0]=v[0];
            //  *T2=v[1];//data_s[1]=v[1];
            //  *a=d[0];//data_s[2]=d[0];
            //  *b=d[1];//data_s[3]=d[1];
            *T1=v[1];//data_s[0]=v[0];
            *T2=v[0];//data_s[1]=v[1];
            *b=fabs(d[0]);//data_s[3]=d[1];
            *a=fabs(d[1]) - fabs(d[0]);//data_s[2]=d[0];
            *number=2;
        }else{ //currently, even I got 2 points, I still use 2nd order, not 3rd order using 2 points
            *T1=v[0];//data_s[0]=v[0];
            *a=fabs(d[0]);//d[0];//data_s[1]=d[0];
            *number=1;
        }
     return true;
    //  return data_s;
  }
//   return (bc - v[0])/(d[0]*Delta); // second-order gradient
  //  double *data_s = (double*) malloc(2 * sizeof(double));
     *T1=v[0];//data_s[0]=v[0];
     *a=fabs(d[0]);//d[0];//data_s[1]=d[0];
     *number=1;
     return true;
    //  return data_s;
  
}



// foreach_dimension()
// static inline double* gradient_conjugatef_x (Point point, scalar s, scalar f, bool inverse_flag,
// 			   coord n, coord p, int * number)
foreach_dimension()
bool gradient_conjugatef_x (Point point, scalar s, scalar f, bool inverse_flag,
			   coord n, coord p, double * T3, double * T4, double *c, double* e, 
         int* number, int* sign_flag)
{
  bool consider_dirichlet_interface=true;
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fss_test3_n.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      *sign_flag = sign(n.x);
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
         if(!consider_dirichlet_interface){
            bool flag_f1 = (((!inverse_flag) && f[i,j-1]>=0.5) || (inverse_flag && f[i,j-1]<0.5)  );
            bool flag_f2 = (((!inverse_flag) && f[i,j]>=0.5) || (inverse_flag && f[i,j]<0.5)  );
            bool flag_f3 = (((!inverse_flag) && f[i,j+1]>=0.5) || (inverse_flag && f[i,j+1]<0.5)  );
            if (fss_test3_n.x[i + (i < 0),j] && fss_test3_n.y[i,j] && fss_test3_n.y[i,j+1] &&
          css_test3_n[i,j-1] && css_test3_n[i,j] && css_test3_n[i,j+1] && (flag_f1) && (flag_f2) && (flag_f3)){
                v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
            }else{
                break;
            }
          }else{
            bool flag_f1 = (css_test3_n[i,j-1] && (((!inverse_flag) && f[i,j-1]>0.0) || (inverse_flag && f[i,j-1]<1))  );
            bool flag_f2 = (css_test3_n[i,j] && (((!inverse_flag) && f[i,j]>0) || (inverse_flag && f[i,j]<1))  );
            bool flag_f3 = (css_test3_n[i,j+1] && (((!inverse_flag) && f[i,j+1]>0) || (inverse_flag && f[i,j+1]<1))  );
            double value1;
            if(flag_f1){
                  if((!inverse_flag && f[i,j-1]<0.5) && (!inverse_flag && f[i,j-1]>0.5) ){
                      value1 = Tsat00;
                  }else{
                      value1 = s[i,j-1];
                  }
            }
            double value2;
            if(flag_f2){
                  if((!inverse_flag && f[i,j]<0.5) && (!inverse_flag && f[i,j]>0.5) ){
                      value2 = Tsat00;
                  }else{
                      value2 = s[i,j];
                  }
            }
            double value3;
            if(flag_f3){
                  if((!inverse_flag && f[i,j+1]<0.5) && (!inverse_flag && f[i,j+1]>0.5) ){
                      value3 = Tsat00;
                  }else{
                      value3 = s[i,j+1];
                  }
            }
            if (fss_test3_n.x[i + (i < 0),j] && fss_test3_n.y[i,j] && fss_test3_n.y[i,j+1] &&
          (flag_f1) && (flag_f2) && (flag_f3)){
                v[l] = quadratic (y1, (value1), (value2), (value3));
             }else{
                  break;
              }
          }
#else // dimension == 3  // note: dimension==3 has not beed modified
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fss_test3_n.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fss_test3_n.y[i,j,k+m] || !fss_test3_n.y[i,j+1,k+m] ||
	    !fss_test3_n.z[i,j+m,k] || !fss_test3_n.z[i,j+m,k+1] ||
	    !css_test3_n[i,j+m,k-1] || !css_test3_n[i,j+m,k] || !css_test3_n[i,j+m,k+1])
	  defined = false;
      if (defined){
            // bi-quadratic interpolation
            v[l] =
            quadratic (z,
                    quadratic (y1,
                        (s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
                    quadratic (y1,
                        (s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
                    quadratic (y1,
                        (s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
      }else{
	        break;
      }
#endif // dimension == 3
      
    }
  if (v[0] == nodata) {
     if(css_test3_n[]>0.5){
        if(((!(inverse_flag)) && f[]>=0.5) || (inverse_flag && f[]<0.5)){
            /**
            This is a degenerate case, we use the boundary value and the
            cell-center value to define the gradient. */
          
            // d[0] = max(1e-1, fabs(p.x/n.x));
            // d[0] = max(1e-1,sqrt(sq(p.x)+sq(p.y)));
            // d[0] = fabs(p.x*n.x + p.y*n.y);
            d[0] = max(1e-2,fabs(p.x*n.x + p.y*n.y));

            // *coef = - 1./(d[0]*Delta);
            // return bc/(d[0]*Delta);
            //  double *data_f = (double*) malloc(2 * sizeof(double));
            *T3=s[];//data_f[0]=s[];
            *c=fabs(d[0]);//d[0];//data_f[1]=d[0];
            *number=1;
            // printf("nu=11111111111\n");


            // double d1=0.0;
            // if(f[]<1.0 && f[]>0.0){
            //     coord n_f = mycs (point, f);
            //     double alpha = plane_alpha (f[], n_f);
            //     coord p_f;
            //     double area = plane_area_center (n_f, alpha, &p_f);
            //     d1 = fabs((p_f.x-p.x)*n.x + (p_f.y-p.y)*n.y);//sqrt(sq(p_f.x-p.x)+sq(p_f.y-p.y));
            // }
            // if(d1>d[0]+1e-2){
            //     *T4=Tsat00;//data_f[1]=v[1];
            //     d[1] = d1;// - d[0];
            //     *e=fabs(d[1]) - fabs(d[0]);//data_f[3]=d[1];
            //     *number=2;
            // }
            return true;
            //  return data_f;
        }else {
            if(1==0){
                    double d0=0.0;
                    if(f[]<1.0 && f[]>0.0){
                        coord n_f = mycs (point, f);
                        double alpha = plane_alpha (f[], n_f);
                        coord p_f;
                        double area = plane_area_center (n_f, alpha, &p_f);
                        d0 = fabs((p_f.x-p.x)*n.x + (p_f.y-p.y)*n.y);//sqrt(sq(p_f.x-p.x)+sq(p_f.y-p.y));
                    }
                    if(d0>1e-2){
                        *T3=Tsat00;//data_f[0]=s[];
                        *c=d0;//d[0];//data_f[1]=d[0];
                        *number=1;
                    }else{
                        *T3=Tsat00;//data_f[0]=s[];
                        *c=1e-2;//d[0];//data_f[1]=d[0];
                        *number=1;
                    }
            }else{
                 d[0] = max(1e-2,fabs(p.x*n.x + p.y*n.y));
                *T3=T[];//data_f[0]=s[];
                *c=fabs(d[0]);//d[0];//data_f[1]=d[0];
                *number=1;
            }
         return true;
            //  return false;
        }
     }else{
         return false;
     }

     //just in order to get an value, find the nearest value;
     bool flag_find=false;
    // // //  if(!(f[]>0.0 && f[]<1.0)){  // method for non-triple point, with css_test3_n<0.5 >0.5: biggest conlinearity
    // // //         Point me = point;
    // // //         double distance1=HUGE;
    // // //         coord temp1;
    // // //         temp1.x=x;
    // // //         temp1.y=y;
    // // //         coord p_real;
    // // //         p_real.x=p.x*Delta+x;
    // // //         p_real.y=p.y*Delta+y;
    // // //         double colinearity=0.0;
    // // //         foreach_neighbor(1){
    // // //             if(!(me.i==point.i && me.j==point.j)){
    // // //                   if(css_test3_n[]>0.5 && (((!(inverse_flag)) && f[]>=0.5) || (inverse_flag && f[]<0.5))){
    // // //                       coord temp2;
    // // //                       temp2.x=x;
    // // //                       temp2.y=y;
    // // //                       double distance_temp=sqrt(sq(temp2.x-p_real.x) + sq(temp2.y-p_real.y));
    // // //                       coord temp2_dir;
    // // //                       temp2_dir.x = temp2.x/distance_temp;
    // // //                       temp2_dir.y = temp2.y/distance_temp;
    // // //                       double colinearity_temp = fabs(temp2_dir.x*n.x+temp2_dir.y*n.y);
    // // //                       if(colinearity_temp>colinearity){
    // // //                           colinearity=colinearity_temp;
    // // //                           distance1 = distance_temp;
    // // //                           *T3=s[];
    // // //                           *c=distance1;
    // // //                           *number=1;
    // // //                           flag_find=true;
    // // //                       }
    // // //                   }
    // // //             }
    // // //         }
    // // //   }else{
    // // //       //for triple point ???????????? : nearest point
    // // //       // // *T3=Tsat00;
    // // //       // // *c=distance1; //from center of solid surface to center of fluid surface
    // // //       // // *number=1;
    // // //       // // flag_find=true;

    // // //        Point me = point;
    // // //             bool flag2=false;
    // // //             double distance1=HUGE;
    // // //             coord direction_temp;
    // // //             coord temp1;
    // // //             temp1.x = x;
    // // //             temp1.y = y;
    // // //             coord p_real;
    // // //             p_real.x=p.x*Delta+x;
    // // //             p_real.y=p.y*Delta+y;
    // // //             foreach_neighbor(1){
    // // //                 if((!(me.i==point.i && me.j==point.j)) && (!is_boundary(cell)) &&  
    // // //                 (((!(inverse_flag)) && f[]>=0.5) || (inverse_flag && f[]<0.5))){
    // // //                     coord temp2;
    // // //                     temp2.x = x;
    // // //                     temp2.y = y;
    // // //                     double value2 = sqrt(sq(p_real.x-temp2.x)+sq(p_real.y-temp2.y));
    // // //                     if(value2<distance1){
    // // //                         flag2 = true;
    // // //                         distance1 = value2;
    // // //                         *c = value2;
    // // //                         *T3 = s[];
    // // //                         *number=1;
    // // //                     }
    // // //                 }
    // // //             }
    // // //             if(flag2){
    // // //                 flag_find=true;
    // // //             }
    // // //   }
      if(flag_find){
          return true;
      }else{
        return false;
      }
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
//   *coef = 0.;
 if (v[1] != nodata){ // third-order gradient
        // if(!Rcc_flag){
        if(((!Rcc_flag_l) && (!inverse_flag)) || (inverse_flag)){
            // return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
            // double *data_f = (double*) malloc(4 * sizeof(double));
            //  *T3=v[0];//data_f[0]=v[0];
            //  *T4=v[1];//data_f[1]=v[1];
            //  *c=d[0];//data_f[2]=d[0];
            //  *e=d[1];//data_f[3]=d[1];
            *T3=v[0];//data_f[0]=v[0];
            *T4=v[1];//data_f[1]=v[1];
            *c=fabs(d[0]);//d[0];//data_f[2]=d[0];
            *e=fabs(d[1]) - fabs(d[0]);//data_f[3]=d[1];
            *number=2;
        }else{ //currently, even I got 2 points, I still use 2nd order, not 3rd order using 2 points
            *T3=v[0];//data_f[0]=v[0];
            *c=fabs(d[0]);//;d[0];//data_f[1]=d[0];
            *number=1;
        }
     return true;
    //  return data_f;
  }
//   return (bc - v[0])/(d[0]*Delta); // second-order gradient
  //  double *data_f = (double*) malloc(2 * sizeof(double));
     *T3=v[0];//data_f[0]=v[0];
     *c=fabs(d[0]);//;d[0];//data_f[1]=d[0];
     *number=1;
    //   printf("nu=2222222\n");
     return true;
    //  return data_f;
  
}


static inline double five_points(double T1,double T2,double T3,double T4
        ,double a,double b,double c,double d,double kl, double kr, double* flux_temp){
        double aT1,aT2,aTf,aT3,aT4;
        //aTf*Tf = aT1*T1 + aT2*T2 + aT3*T3 + aT4*T4;
        double Tf; 
        aTf = (a+2*b)/(a+b)/b*kl + (d+2*c)/(d+c)/c*kr;
        aT1 = kl*(-(a+2*b)/(a+b)/a + 1.0/a);
        aT2 = kl*((a+2*b)/(a+b)*(1.0/a+1.0/b)-1.0/a);
        aT3 = kr*((d+2*c)/(d+c)*(1.0/d+1.0/c)-1.0/d);
        aT4 = kr*(-(d+2*c)/(d+c)/d + 1.0/d);
        if(fabs(aTf)>1e-12){
             Tf = (aT1*T1 + aT2*T2 + aT3*T3 + aT4*T4)/aTf;
             *flux_temp = kl*((T2-T1)/a + ((Tf-T2)/b-(T2-T1)/a)/(a+b)*(a+2*b));
        }else{
             Tf = -1;
        }
        return Tf;
}

static inline double one_two_points(double T1,double T3,double T4
        ,double a,double c,double d,double kl, double kr, double* flux_temp){
        double aT1,aTf,aT3,aT4;
        //aTf*Tf = aT1*T1 + aT3*T3 + aT4*T4;
        double Tf; 
        aTf = 1.0/a*kl + (d+2*c)/(d+c)/c*kr;
        aT1 = kl*(1.0/a);
        aT3 = kr*((d+2*c)/(d+c)*(1.0/d+1.0/c)-1.0/d);
        aT4 = kr*(-(d+2*c)/(d+c)/d + 1.0/d);
        if(fabs(aTf)>1e-12){
             Tf = (aT1*T1 + aT3*T3 + aT4*T4)/aTf;
             *flux_temp = kl*((Tf-T1)/a);
        }else{
             Tf = -1;
        }
        return Tf;
}

static inline double two_one_points(double T1,double T2,double T3
        ,double a,double b,double c,double kl, double kr, double* flux_temp){
        double aT1,aT2,aTf,aT3;
        //aTf*Tf = aT1*T1 + aT2*T2 + aT3*T3 + aT4*T4;
        double Tf; 
        aTf = (a+2*b)/(a+b)/b*kl + 1.0/c*kr;
        aT1 = kl*(-(a+2*b)/(a+b)/a + 1.0/a);
        aT2 = kl*((a+2*b)/(a+b)*(1.0/a+1.0/b)-1.0/a);
        aT3 = kr*(1.0/c);
        if(fabs(aTf)>1e-12){
             Tf = (aT1*T1 + aT2*T2 + aT3*T3)/aTf;
             *flux_temp = kr*((T3-Tf)/c);
        }else{
             Tf = -1;
        }
        return Tf;
}

static inline double one_one_points (double T1,double T3,
        double a,double c,double kl, double kr, double* flux_temp){
        double aT1,aTf,aT3;
        //aTf*Tf = aT1*T1 + aT2*T2 + aT3*T3 + aT4*T4;
        double Tf; 
        aTf = 1.0/a*kl + 1.0/c*kr;
        aT1 = kl*1.0/a;
        aT3 = kr*1.0/c;
        if(fabs(aTf)>1e-12){
             Tf = (aT1*T1 + aT3*T3)/aTf;
             *flux_temp = kl*((Tf-T1)/a);
        }else{
             Tf = -1;
        }
        return Tf;
}
void one_one_points_Rcc (double T1,double T3,
        double a,double c,double kl, double kr, double* solid_temp){
        // double* solid_temp = (double*)malloc(2 * sizeof(double));
        double Rs = a/kl, Rl = c/kr;
        double Rcc_v = Rcc/delta_min;
        double aTf = (Rs+Rl)*Rcc_v+Rcc_v*Rcc_v;
        // printf("a=%g c=%g \n",a,c);
        // printf("Rs=%g Rl=%g  Rcc_v=%g \n",Rs,Rl,Rcc_v);
        if(fabs(aTf)>1e-12){
            solid_temp[0] = ((Rs+Rcc_v)*Rcc_v*T3+Rl*Rcc_v*T1)/aTf; //T2_l
            solid_temp[1] = ((Rl+Rcc_v)*Rcc_v*T1+Rs*Rcc_v*T3)/aTf; //T2_s
            solid_temp[2] = (solid_temp[0] - solid_temp[1])/Rcc_v;
        }else{
            solid_temp[0] = -1;
            solid_temp[1] = -1;
        }
        
        // return solid_temp;
        // double aT1,aTf,aT3;
        // //aTf*Tf = aT1*T1 + aT2*T2 + aT3*T3 + aT4*T4;
        // double Tf; 
        // aTf = 1.0/a*kl + 1.0/c*kr;
        // aT1 = kl*1.0/a;
        // aT3 = kr*1.0/c;
        // if(fabs(aTf)>1e-12){
        //      Tf = (aT1*T1 + aT3*T3)/aTf;
        //      *flux_temp = kl*((Tf-T1)/a);
        // }else{
        //      Tf = -1;
        // }
        // return Tf;

}

double equal_of_flux(double* data_s, int number_s, double* data_f, int number_f, int * mod_temp, bool inverse_flag, double* flux_temp){
   double result=-1;
   double k1,k2;
  //  double Tks,Tkl,Tkg;
   if(!inverse_flag){
      k1 = Tks;
      k2 = Tkl; 
   }else{
      k1 = Tks;
      k2 = Tkg;
   }
   if((number_s==2) && (number_f==2)){
        double T1=data_s[0],T2=data_s[1],a=data_s[2];
        double b=data_s[3];
        double T3=data_f[0],T4=data_f[1],c=data_f[2];
        double d=data_f[3];
        *mod_temp=1;
        result=five_points(T1,T2,T3,T4,a,b,c,d,k1,k2,flux_temp);
   }else if((number_s==1) && (number_f==2)){
        double T1=data_s[0],a=data_s[1];
        double T3=data_f[0],T4=data_f[1],c=data_f[2];
        double d=data_f[3];
        *mod_temp=2;
        result=one_two_points(T1,T3,T4,a,c,d,k1,k2,flux_temp);
   }else if((number_s==2) && (number_f==1)){
        double T1=data_s[0],T2=data_s[1],a=data_s[2];
        double b=data_s[3];
        double T3=data_f[0], c=data_f[1];
        *mod_temp=3;
        result=two_one_points(T1,T2,T3,a,b,c,k1,k2,flux_temp);
   }else if((number_s==1) && (number_f==1)){
        double T1=data_s[0],a=data_s[1];
        double T3=data_f[0],c=data_f[1];
        *mod_temp=4;
        result=one_one_points(T1,T3,a,c,k1,k2,flux_temp);
   }else{
    //   result = 0.0;
    //   printf("error: number_s=%d, number_f=%d\n",number_s,number_f);
   }
      return result;
}

void equal_of_flux_Rcc(double* data_s, int number_s, double* data_f, int number_f, int * mod_temp, 
                      bool inverse_flag, double* result){
//    double result=-1;
    // double* result= NULL;;
   double k1,k2;
  //  double Tks,Tkl,Tkg;
   if(!inverse_flag){
      k1 = Tks;
      k2 = Tkl; 
   }else{
      k1 = Tks;
      k2 = Tkg;
   }
   
        if((number_s==1) && (number_f==1)){
                // printf("number_s=%d, number_f=%d\n",number_s,number_f);
                double T1=data_s[0],a=data_s[1];
                double T3=data_f[0],c=data_f[1];
                *mod_temp=4;
                one_one_points_Rcc(T1,T3,a,c,k1,k2,result);
        }else{
            //   result = 0.0;
            //   printf("error: number_s=%d, number_f=%d\n",number_s,number_f);
        }

    //   return result;
}

bool max_colinearity_f (Point point, scalar s, scalar f, bool inverse_flag,
			   coord n, coord p, double * T3, double * T4, double *c, double* e, 
         int* number)
{
    foreach_dimension()
    n.x = - n.x;
     bool flag_find=false;
     //if(!(ff[]>0.0 && ff[]<1.0)){  // method for non-triple point, with css_test3_n<0.5 >0.5: biggest conlinearity
            Point me = point;
            double distance1=HUGE;
            coord temp1;
            temp1.x=x;
            temp1.y=y;
            coord p_real;
            p_real.x=p.x*Delta+x;
            p_real.y=p.y*Delta+y;
            double colinearity=0.0;
            foreach_neighbor(1){
                if(!(me.i==point.i && me.j==point.j)){
                      if(css_test3_n[]>0.5 && (((!(inverse_flag)) && ff[]>=0.5) || (inverse_flag && ff[]<0.5))){
                          coord temp2;
                          temp2.x=x;
                          temp2.y=y;
                          double distance_temp=sqrt(sq(temp2.x-p_real.x) + sq(temp2.y-p_real.y));
                          coord temp2_dir;
                          temp2_dir.x = temp2.x/distance_temp;
                          temp2_dir.y = temp2.y/distance_temp;
                          double colinearity_temp = fabs(temp2_dir.x*n.x+temp2_dir.y*n.y);
                          if(colinearity_temp>colinearity){
                              colinearity=colinearity_temp;
                              distance1 = distance_temp;
                              *T3=s[];
                              *c=distance1;
                              *number=1;
                              flag_find=true;
                          }
                      }
                }
            }
   //   }
      if(flag_find){
          return true;
      }else{
        return false;
      }
}

bool min_distance_f (Point point, scalar s, scalar f, bool inverse_flag,
			   coord n, coord p, double * T3, double * T4, double *c, double* e, 
         int* number)
{
      foreach_dimension()
        n.x = - n.x;
      bool flag_find=false;
     
           Point me = point;
                bool flag2=false;
                double distance1=HUGE;
                coord direction_temp;
                coord temp1;
                temp1.x = x;
                temp1.y = y;
                coord p_real;
                p_real.x=p.x*Delta+x;
                p_real.y=p.y*Delta+y;
                foreach_neighbor(1){
                    if((!(me.i==point.i && me.j==point.j)) && (!is_boundary(cell)) &&  
                    (((!(inverse_flag)) && ff[]>=0.5) || (inverse_flag && ff[]<0.5)) && 
                        css_test3_n[]>0.5){
                        coord temp2;
                        temp2.x = x;
                        temp2.y = y;
                        double value2 = sqrt(sq(p_real.x-temp2.x)+sq(p_real.y-temp2.y));
                        if(value2<distance1){
                            flag2 = true;
                            distance1 = value2;
                            *c = value2;
                            *T3 = s[];
                            *number=1;
                        }
                    }
                }
                if(flag2){
                    flag_find=true;
                }
      
      if(flag_find){
          return true;
      }else{
        return false;
      }
}

bool max_colinearity_s (Point point, scalar s, scalar f, bool inverse_flag,
			   coord n, coord p, double * T1, double * T2, double *a, double* b, 
         int* number)
{
    foreach_dimension()
    n.x = - n.x;
     bool flag_find=false;
     //if(!(ff[]>0.0 && ff[]<1.0)){  // method for non-triple point, with css_test3_n<0.5 >0.5: biggest conlinearity
            Point me = point;
            double distance1=HUGE;
            coord temp1;
            temp1.x=x;
            temp1.y=y;
            coord p_real;
            p_real.x=p.x*Delta+x;
            p_real.y=p.y*Delta+y;
            double colinearity=0.0;
            foreach_neighbor(1){
                if(!(me.i==point.i && me.j==point.j)){
                      if((((!(inverse_flag)) && f[]>=0.5) || (inverse_flag && f[]<0.5))){
                          coord temp2;
                          temp2.x=x;
                          temp2.y=y;
                          double distance_temp=sqrt(sq(temp2.x-p_real.x) + sq(temp2.y-p_real.y));
                          coord temp2_dir;
                          temp2_dir.x = temp2.x/distance_temp;
                          temp2_dir.y = temp2.y/distance_temp;
                          double colinearity_temp = fabs(temp2_dir.x*n.x+temp2_dir.y*n.y);
                          if(colinearity_temp>colinearity){
                              colinearity=colinearity_temp;
                              distance1 = distance_temp;
                              *T1=s[];
                              *a=distance1;
                              *number=1;
                              flag_find=true;
                          }
                      }
                }
            }
   //   }
      if(flag_find){
          return true;
      }else{
        return false;
      }
}

double gradient_conjugate (Point point, scalar s, scalar ss,coord * nn1_temp, coord * nn2_temp, coord * p_temp,
             int* mod_temp, bool is_liquid, double* flux_temp)
{

// ///////////////////////////////////////////data_s

bool inverse_flag;
inverse_flag=false;
coord nn1,nn2,p;
foreach_dimension(){
  nn1.x=nn1_temp->x;
  nn2.x=nn2_temp->x;
  p.x=p_temp->x;
}
double result=0;
// double data_s[4];
double T1, T2, a, b;
// double data_f[4];
double T3,T4,c,d;
int sign_flag_s;
int sign_flag_f;
double data_s[4];
double data_f[4];
//near triple point
int number_s =0;
int number_f=0;
bool near_triple_flag=false;
// foreach_neighbor(2){
//     if(css_test3[]>0 && css_test3[]<1.0 && ff[]>0.0 && ff[]<1.0){
//         near_triple_flag = true;
//     }
// }
if(!near_triple_flag){
    if(is_liquid){
    inverse_flag=false;
    }else{
    inverse_flag=true;
    }
    #if dimension == 2
    bool flag1=false;
    if(1==0){
        foreach_dimension()
            if (fabs(nn1.x) >= fabs(nn1.y) && (!flag1)){
            flag1 = gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p, 
                &T1,&T2,&a,&b, &number_s,&sign_flag_s);
            // flag1 = true;
            }
    }else{
        if(fabs(nn1.x) >= fabs(nn1.y)){
            flag1 = gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p, 
                &T1,&T2,&a,&b, &number_s,&sign_flag_s);
            //  number_s_f.x[] = 1;
        }else{
            flag1 = gradient_conjugates_y (point, s, css_test3, inverse_flag, nn1, p, 
                &T1,&T2,&a,&b, &number_s,&sign_flag_s);
            // number_s_f.x[] = 2;
        }
    }

    #else // dimension == 3
    if (fabs(nn1.x) >= fabs(nn1.y)) {
        if (fabs(nn1.x) >= fabs(nn1.z))
        // *data_s = gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p,  &number_s);
        gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p, 
            &T1,&T2,&a,&b, &number_s,&sign_flag_s);
    }
    else if (fabs(nn1.y) >= fabs(nn1.z)){
        // *data_s = gradient_conjugates_y (point, s, css_test3, inverse_flag, nn1, p,  &number_s);
        gradient_conjugates_y (point, s, css_test3, inverse_flag, nn1, p, 
            &T1,&T2,&a,&b, &number_s,&sign_flag_s);
    }else{
        // *data_s =  gradient_conjugates_z (point, s,  css_test3, inverse_flag, nn1,  p, &number_s);
        gradient_conjugates_z (point, s, css_test3, inverse_flag, nn1, p, 
            &T1,&T2,&a,&b, &number_s,&sign_flag_s);
    }
    #endif // dimension == 3
    // if(number_s==0){
    //    max_colinearity_s (point, s, css_test3, inverse_flag, nn1, p, 
    //       &T1,&T2,&a,&b, &number_s);
    // }
    // //    return data_s;

    // ///////////////////////////////////////////data_f
    if(is_liquid){
    inverse_flag=false;
    }else{
    inverse_flag=true;
    }
    // number_f =0;
    #if dimension == 2
    bool flag2=false;
    foreach_dimension()
        if (fabs(nn2.x) >= fabs(nn2.y) && (!flag2)){
        // *data_f;// = gradient_conjugatef_x (point, ss,  ff, inverse_flag, nn2, p, &number_f);
        flag2 = gradient_conjugatef_x (point, ss, ff, inverse_flag, nn2, p, 
            &T3,&T4,&c,&d, &number_f,&sign_flag_f);
        // flag2 = true;
        }
    #else // dimension == 3
    if (fabs(nn2.x) >= fabs(nn2.y)) {
        if (fabs(nn2.x) >= fabs(nn2.z))
        // *data_f;// = gradient_conjugatef_x (point, ss, ff, inverse_flag, nn2, p,  &number_f);
        gradient_conjugatef_x (point, ss, ff, inverse_flag, nn2, p, 
            &T3,&T4,&c,&d, &number_f,&sign_flag_f);

    }reach_neighbo
    else if (fabs(nn2.y) >= fabs(nn2.z)){
        // *data_f;// = gradient_conjugatef_y (point, ss, ff, inverse_flag, nn2, p,  &number_f);
        gradient_conjugatef_y (point, ss, ff, inverse_flag, nn2, p, 
            &T3,&T4,&c,&d, &number_f,&sign_flag_f);
    }else{
        // *data_f;// =  gradient_conjugatef_z (point, ss,  ff, inverse_flag, nn2,  p, &number_f);
        gradient_conjugatef_z (point, ss, ff, inverse_flag, nn2, p, 
            &T3,&T4,&c,&d, &number_f,&sign_flag_f);
    }
    #endif // dimension == 3
}else{ //if near triple point
    if(css_test3[]>=0.5){
        number_s = 1;
        T1 = s[];
        // data_s[1] = max(1e-1,sqrt(sq(p.x)+sq(p.y)));
        a = max(1e-1,fabs(p.x*nn1.x + p.y*nn1.y));
        int ii,jj;
        if(is_liquid){
            ii = round(merge_to_me_l_position.x[]);
            jj = round(merge_to_me_l_position.y[]);
        }else{
            ii = round(merge_to_me_g_position.x[]);
            jj = round(merge_to_me_g_position.y[]);
        }
        number_f = 1;
        T3 = ss[ii,jj];
        c = max(1e-1,fabs((ii-p.x)*nn2.x + (jj-p.y)*nn2.y));
    }else{ //css_test3[] can not == 0.5, so in initial tuning 1e-20, if ==0.5
        number_s = 1;
        int ii = round(merge_to_me_s_position.x[]);
        int jj = round(merge_to_me_s_position.y[]);
        T1 = s[ii,jj];
        a = max(1e-1,fabs((ii-p.x)*nn1.x + (jj-p.y)*nn1.y));
        number_f = 1;
        T3 = ss[];
        // data_f[1] = max(1e-1,sqrt(sq(p.x)+sq(p.y)));
        c = max(1e-1,fabs(p.x*nn2.x + p.y*nn2.y));
    }
}
  // if(number_f==0){
  //   if(!(ff[]>0.0 && ff[]<1.0)){
  //     max_colinearity_f (point, ss, ff, inverse_flag, nn2, p, 
  //       &T3,&T4,&c,&d, &number_f);
  //   }else{
  //      min_distance_f (point, ss, ff, inverse_flag, nn2, p, 
  //        &T3,&T4,&c,&d, &number_f);
  //   }
  // }


   if(number_s>=1 && number_f>=1){
        if(number_s==1){
              data_s[0]=T1;
              data_s[1]=a;
            //   printf("number_s=1 a=%g T1=%g ",a,T1);
        }else if(number_s==2){
              data_s[0]=T1;
              data_s[1]=T2;
              data_s[2]=a;
              data_s[3]=b;
            //   printf("number_s=2 a=%g b=%g T1=%g T2=%g ",a,b,T1,T2);
        }
        if(number_f==1){
              data_f[0]=T3;
              data_f[1]=c;
            //   printf("number_f=1 c=%g T3=%g ",c,T3);
        }else if(number_f==2){
              data_f[0]=T3;
              data_f[1]=T4;
              data_f[2]=c;
              data_f[3]=d;
            //   printf("number_f=2 c=%g d=%g T3=%g T4=%g ",c,d,T3,T4);
        }
        // result = sign_flag_s*equal_of_flux(data_s, number_s, data_f, number_f, mod_temp, inverse_flag, flux_temp);
        result = equal_of_flux(data_s, number_s, data_f, number_f, mod_temp, inverse_flag, flux_temp);
        // printf("Tf=%g \n",result);
   }else{
        // printf("no equal_of_flux, mod_temp=%d Tsurface_get_from average\n",*mod_temp);
       
   }
   
// if(is_liquid){
//      number_s_f.x[] = number_s;
//      number_s_f.y[] = number_f;
// }
  //  free(data_f);
  //  free(data_s);
   return result;

}

void gradient_conjugate_Rcc (Point point, scalar s, scalar ss,coord * nn1_temp, coord * nn2_temp, coord * p_temp,
             int* mod_temp, bool is_liquid, double* result)
{

// ///////////////////////////////////////////data_s

bool inverse_flag;
inverse_flag=false;
coord nn1,nn2,p;
foreach_dimension(){
  nn1.x=nn1_temp->x;
  nn2.x=nn2_temp->x;
  p.x=p_temp->x;
}
// double result=0;
// double data_s[4];
double T1, T2, a, b;
// double data_f[4];
double T3,T4,c,d;
int sign_flag_s;
int sign_flag_f;
double data_s[4];
double data_f[4];
//near triple point
int number_s =0;
int number_f=0;
bool near_triple_flag=false;
// foreach_neighbor(2){
//     if(css_test3[]>0 && css_test3[]<1.0 && ff[]>0.0 && ff[]<1.0){
//         near_triple_flag = true;
//     }
// }
if(!near_triple_flag){
    if(is_liquid){
    inverse_flag=false;
    }else{
    inverse_flag=true;
    }
    #if dimension == 2
    bool flag1=false;
    foreach_dimension()
        if (fabs(nn1.x) >= fabs(nn1.y) && (!flag1)){
        flag1 = gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p, 
            &T1,&T2,&a,&b, &number_s,&sign_flag_s);
        // flag1 = true;
        }
    #else // dimension == 3
    if (fabs(nn1.x) >= fabs(nn1.y)) {
        if (fabs(nn1.x) >= fabs(nn1.z))
        // *data_s = gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p,  &number_s);
        gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p, 
            &T1,&T2,&a,&b, &number_s,&sign_flag_s);
    }
    else if (fabs(nn1.y) >= fabs(nn1.z)){
        // *data_s = gradient_conjugates_y (point, s, css_test3, inverse_flag, nn1, p,  &number_s);
        gradient_conjugates_y (point, s, css_test3, inverse_flag, nn1, p, 
            &T1,&T2,&a,&b, &number_s,&sign_flag_s);
    }else{
        // *data_s =  gradient_conjugates_z (point, s,  css_test3, inverse_flag, nn1,  p, &number_s);
        gradient_conjugates_z (point, s, css_test3, inverse_flag, nn1, p, 
            &T1,&T2,&a,&b, &number_s,&sign_flag_s);
    }
    #endif // dimension == 3
    // if(number_s==0){
    //    max_colinearity_s (point, s, css_test3, inverse_flag, nn1, p, 
    //       &T1,&T2,&a,&b, &number_s);
    // }
    // //    return data_s;

    // ///////////////////////////////////////////data_f
    if(is_liquid){
    inverse_flag=false;
    }else{
    inverse_flag=true;
    }
    // number_f =0;
    #if dimension == 2
    bool flag2=false;
    foreach_dimension()
        if (fabs(nn2.x) >= fabs(nn2.y) && (!flag2)){
        // *data_f;// = gradient_conjugatef_x (point, ss,  ff, inverse_flag, nn2, p, &number_f);
        flag2 = gradient_conjugatef_x (point, ss, ff, inverse_flag, nn2, p, 
            &T3,&T4,&c,&d, &number_f,&sign_flag_f);
        // flag2 = true;
        }
    #else // dimension == 3
    if (fabs(nn2.x) >= fabs(nn2.y)) {
        if (fabs(nn2.x) >= fabs(nn2.z))
        // *data_f;// = gradient_conjugatef_x (point, ss, ff, inverse_flag, nn2, p,  &number_f);
        gradient_conjugatef_x (point, ss, ff, inverse_flag, nn2, p, 
            &T3,&T4,&c,&d, &number_f,&sign_flag_f);

    }reach_neighbo
    else if (fabs(nn2.y) >= fabs(nn2.z)){
        // *data_f;// = gradient_conjugatef_y (point, ss, ff, inverse_flag, nn2, p,  &number_f);
        gradient_conjugatef_y (point, ss, ff, inverse_flag, nn2, p, 
            &T3,&T4,&c,&d, &number_f,&sign_flag_f);
    }else{
        // *data_f;// =  gradient_conjugatef_z (point, ss,  ff, inverse_flag, nn2,  p, &number_f);
        gradient_conjugatef_z (point, ss, ff, inverse_flag, nn2, p, 
            &T3,&T4,&c,&d, &number_f,&sign_flag_f);
    }
    #endif // dimension == 3
}else{ //if near triple point
    if(css_test3[]>=0.5){
        number_s = 1;
        T1 = s[];
        // data_s[1] = max(1e-1,sqrt(sq(p.x)+sq(p.y)));
        a = max(1e-1,fabs(p.x*nn1.x + p.y*nn1.y));
        int ii,jj;
        if(is_liquid){
            ii = round(merge_to_me_l_position.x[]);
            jj = round(merge_to_me_l_position.y[]);
        }else{
            ii = round(merge_to_me_g_position.x[]);
            jj = round(merge_to_me_g_position.y[]);
        }
        number_f = 1;
        T3 = ss[ii,jj];
        c = max(1e-1,fabs((ii-p.x)*nn2.x + (jj-p.y)*nn2.y));
    }else{ //css_test3[] can not == 0.5, so in initial tuning 1e-20, if ==0.5
        number_s = 1;
        int ii = round(merge_to_me_s_position.x[]);
        int jj = round(merge_to_me_s_position.y[]);
        T1 = s[ii,jj];
        a = max(1e-1,fabs((ii-p.x)*nn1.x + (jj-p.y)*nn1.y));
        number_f = 1;
        T3 = ss[];
        // data_f[1] = max(1e-1,sqrt(sq(p.x)+sq(p.y)));
        c = max(1e-1,fabs(p.x*nn2.x + p.y*nn2.y));
    }
}
  // if(number_f==0){
  //   if(!(ff[]>0.0 && ff[]<1.0)){
  //     max_colinearity_f (point, ss, ff, inverse_flag, nn2, p, 
  //       &T3,&T4,&c,&d, &number_f);
  //   }else{
  //      min_distance_f (point, ss, ff, inverse_flag, nn2, p, 
  //        &T3,&T4,&c,&d, &number_f);
  //   }
  // }

//   double* result = NULL;;

   if(number_s==1 && number_f==1){ //only 1-1 point allowed for Rcc
        if(number_s==1){
              data_s[0]=T1;
              data_s[1]=a;
            //   printf("number_s=1 a=%g T1=%g ",a,T1);
        }
        if(number_f==1){
              data_f[0]=T3;
              data_f[1]=c;
            //   printf("number_f=1 c=%g T3=%g ",c,T3);
        }
        equal_of_flux_Rcc(data_s, number_s, data_f, number_f, mod_temp, inverse_flag, result);
   }else{
        // printf("no equal_of_flux, mod_temp=%d Tsurface_get_from average\n",*mod_temp);
   }
   

  //  free(data_f);
  //  free(data_s);
//    return result;

}


void gradient_conjugate_for_solid (Point point, scalar s, scalar ss, coord * nn1_temp, coord * p_temp, 
                     int* mod_temp, bool is_liquid, double* flux_temp)
{
// ///////////////////////////////////////////data_s
int number_s =0;
bool inverse_flag;
inverse_flag=false;
coord nn1,p;
foreach_dimension(){
  nn1.x=nn1_temp->x;
//   nn2.x=nn2_temp->x;
  p.x=p_temp->x;
}
double result=0;
// double data_s[4];
double T1, T2, a, b;
int sign_flag_s;
#if dimension == 2
 bool flag1=false;
 if(1==0){
  foreach_dimension()
    if (fabs(nn1.x) >= fabs(nn1.y) && (!flag1)){
      flag1 = gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p, 
        &T1,&T2,&a,&b, &number_s,&sign_flag_s);
      // flag1 = true;
    }
 }else{
    if (fabs(nn1.x) >= fabs(nn1.y)){
      flag1 = gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p, 
        &T1,&T2,&a,&b, &number_s,&sign_flag_s);
      // flag1 = true;
    //    number_s_f.x[] = 1;
    }else{
        flag1 = gradient_conjugates_y (point, s, css_test3, inverse_flag, nn1, p, 
        &T1,&T2,&a,&b, &number_s,&sign_flag_s);
        // number_s_f.x[] = 2;
    }
 }
#else // dimension == 3
  if (fabs(nn1.x) >= fabs(nn1.y)) {
    if (fabs(nn1.x) >= fabs(nn1.z))
      // *data_s = gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p,  &number_s);
      gradient_conjugates_x (point, s, css_test3, inverse_flag, nn1, p, 
        &T1,&T2,&a,&b, &number_s,&sign_flag_s);
  }
  else if (fabs(nn1.y) >= fabs(nn1.z)){
    // *data_s = gradient_conjugates_y (point, s, css_test3, inverse_flag, nn1, p,  &number_s);
    gradient_conjugates_y (point, s, css_test3, inverse_flag, nn1, p, 
        &T1,&T2,&a,&b, &number_s,&sign_flag_s);
  }else{
    // *data_s =  gradient_conjugates_z (point, s,  css_test3, inverse_flag, nn1,  p, &number_s);
    gradient_conjugates_z (point, s, css_test3, inverse_flag, nn1, p, 
          &T1,&T2,&a,&b, &number_s,&sign_flag_s);
  }
#endif // dimension == 3

if(is_liquid){
  inverse_flag=false;
}else{
  inverse_flag=true;
}
// double data_s[4];
   if(number_s>=1){
        if(number_s==1){
            //   data_s[0]=T1;
            //   data_s[1]=a;
            //   *flux_temp = sign_flag_s*((ss[] - T1)/a);
            *flux_temp = Tks*((ss[] - T1)/a);
            *mod_temp = 6;
        }else if(number_s==2){
            //   data_s[0]=T1;
            //   data_s[1]=T2;
            //   data_s[2]=a;
            //   data_s[3]=b;
            //    *flux_temp = sign_flag_s*((T2-T1)/a + ((ss[]-T2)/b-(T2-T1)/a)/(a+b)*(a+2*b));
            *flux_temp = Tks*((T2-T1)/a + ((ss[]-T2)/b-(T2-T1)/a)/(a+b)*(a+2*b));
            *mod_temp = 6;
        }else{
            printf("no gradient_conjugate_for_solid 2\n");
        }
   }else{
        printf("no gradient_conjugate_for_solid 2\n");
   }
   

//   //  free(data_f);
//   //  free(data_s);
//    return result;

}


























// extern face vector fss_test3;
// double embed_flux_conjugate (scalar aim, scalar T_mod, scalar s, scalar ss, face vector mu)
//s is solid, ss is liquid, sss is gas;
void embed_flux_conjugate (scalar s, scalar ss, scalar sss, scalar aiml, scalar aimg, scalar T_modl, scalar T_modg, scalar flux_l, scalar flux_g)
{
 
foreach(){
            /**
         If the cell does not contain a fragment of embedded boundary, the
        flux is zero. */
        T_modl[] = 0;
        T_modg[] = 0;
        aiml[]=0;
        aimg[]=0;
        flux_l[] = 0.0;
        flux_g[] =0.0;
      if (!(css_test3[] >= 1. || css_test3[] <= 0.)){

                 
              /**
               We compute the normal, area and barycenter of the fragment of embedded
              boundary contained within the cell. */

              coord n = facet_normal (point, css_test3, fss_test3);
              coord p;
              double alpha = plane_alpha (css_test3[], n);
              double area = plane_area_center (n, alpha, &p);
              if (metric_embed_factor)  
                  area *= metric_embed_factor (point, p);
              

              /**
               If the boundary condition is Dirichlet, we need to compute the
              normal gradient. */

              // double coef = 0.;

              
              //   if (dirichlet) {
                  normalize(&n);

                //   number_s_f.x[] = n.x;
                //  number_s_f.y[] = n.y;

                  coord nn1;
                  if(!structure_normal_flag){
                    nn1.x = direction_Ts.x[];
                    nn1.y = direction_Ts.y[];
                  }else{
                    nn1.x = n.x;
                    nn1.y = n.y;
                  }
                  double value_temp1= sqrt(sq(nn1.x)+sq(nn1.y));
                  //nn1 is the direction of temperature higher, now transfer it to direction of the phase. only in solidface cell
                 

                  coord nn2;
                 if(!structure_normal_flag){
                    nn2.x = direction_Tfl.x[];
                    nn2.y = direction_Tfl.y[];
                 }else{
                    nn2.x = -n.x;
                    nn2.y = -n.y;
                 }
                  double value_temp2= sqrt(sq(nn2.x)+sq(nn2.y));
                 

                  coord nn3;
                if(!structure_normal_flag){
                    nn3.x = direction_Tfg.x[];
                    nn3.y = direction_Tfg.y[];
                 }else{
                    nn3.x = -n.x;
                    nn3.y = -n.y;
                 }
                  double value_temp3= sqrt(sq(nn3.x)+sq(nn3.y));
                  

                  // int interface_num=1; //interface_num==1, interface between, solid and liquid
                  int mod_templ=-1;
                  T_modl[] = mod_templ;
                  int mod_tempg=-1;
                  T_modg[] = mod_tempg;
                //n is the normal out of solid
                //nn1 changes it's direction to the side same with out of solid
                //nn2 is the normal direction out of fluid ?????????????????????
                if((value_temp1>1e-20) && (value_temp2>1e-20) && ff[]>0.0){
                        // number_s_f.x[] = n.x;
                        //  number_s_f.y[] = n.y;
                         double same_direction = n.x*nn1.x + n.y*nn1.y; 
                        if(same_direction<0.0){
                            nn1.x = -nn1.x;
                            nn1.y = -nn1.y;
                        }
                         same_direction = (-n.x)*nn2.x + (-n.y)*nn2.y; // -n.x -n.y is the direction out fluid
                        if(same_direction<0.0){
                            nn2.x = -nn2.x;
                            nn2.y = -nn2.y;
                        }
                        // double grad = gradient_conjugate (point,s,ss,nn1,nn2,p,&mod_temp);
                        bool is_liquid=true;
                        coord pp;
                        foreach_dimension(){
                          pp.x=p.x;
                        }
                        // if(!Rcc_flag){
                        if(!Rcc_flag_l){
                                double gradl;
                                double flux1=0.0;
                                gradl=gradient_conjugate (point, s, ss, &nn1, &nn2, &pp, &mod_templ,is_liquid,&flux1);
                                //       
                                if(mod_templ!=-1 && fabs(gradl+1)>1e-20){                          //Ts,Tfl
                                    aiml[] = gradl   ;
                                    T_modl[] = mod_templ;
                                    flux_l[] = flux1;
                                }
                        }else{
                                double gradl[3];;
                                // double flux1=0.0;
                                gradl[2] = 0;
                 gradient_conjugate_Rcc (point, s, ss, &nn1, &nn2, &pp, &mod_templ,is_liquid,gradl);
                                //  
                                if(gradl!=NULL){     
                                    // if(mod_templ!=-1 && fabs(gradl+1)>1e-20){  
                                    if(mod_templ!=-1 && (fabs(gradl[0]+1)>1e-20)){                          //Ts,Tfl
                                        aiml[] = gradl[0] ;
                                        aiml_s[] = gradl[1] ; 
                                        T_modl[] = mod_templ;
                                        flux_l[] = gradl[2] ;
                                    }
                                    // free(gradl);
                                }
                        }
                        // flux_l[] = flux1;
                }

                if((value_temp1>1e-20) && (value_temp3>1e-20) && ff[]<1.0){
                        double same_direction = n.x*nn1.x + n.y*nn1.y;
                        if(same_direction<0.0){
                            nn1.x = -nn1.x;
                            nn1.y = -nn1.y;
                        }
                        same_direction = (-n.x)*nn3.x + (-n.y)*nn3.y; // -n.x -n.y is the direction out fluid
                        if(same_direction<0.0){
                            nn3.x = -nn3.x;
                            nn3.y = -nn3.y;
                        }
                        // double grad = gradient_conjugate (point,s,ss,nn1,nn2,p,&mod_temp);
                        bool is_liquid=false;
                        coord pp;
                        foreach_dimension(){
                          pp.x=p.x;
                        }
                        // if(!Rcc_flag){
                        if(!Rcc_flag_g){
                                double gradg;
                                double flux1=0.0;
                                gradg=gradient_conjugate (point, s, sss, &nn1, &nn3, &pp, &mod_tempg, is_liquid, &flux1);
                                //                                 //Ts,Tfg
                                if(mod_tempg!=-1 && fabs(gradg+1)>1e-20){ 
                                    aimg[] = gradg   ;
                                    T_modg[] = mod_tempg;
                                    flux_g[] = flux1;
                                }
                        }else{
                                double gradg[3];
                                // double flux1=0.0;
                                gradg[2]=0.0;
                                gradient_conjugate_Rcc (point, s, sss, &nn1, &nn3, &pp, &mod_tempg, is_liquid,gradg);
                                //                                 //Ts,Tfg
                                // if(mod_tempg!=-1 && fabs(gradg+1)>1e-20){ 
                                if(gradg!=NULL){          
                                        if(mod_tempg!=-1 && fabs(gradg[0]+1)>1e-20){ 
                                            aimg[] = gradg[0]   ;
                                            aimg_s[] = gradg[1] ;
                                            T_modg[] = mod_tempg;
                                            flux_g[] = gradg[2] ;
                                        }
                                        // free(gradg);
                                }          
                         }
                        // flux_g[] = flux1;
                }
                                                  //Ts,                    , Tf
              //   }

              
              // // We retrieve the (average) value of $\mu$ without the metric. */
              
              // double mua = 0., fa = 0.;
              // foreach_dimension() {
              //     mua += mu.x[] + mu.x[1];
              //     fa  += fm.x[] + fm.x[1];
              // }
              // return - mua/(fa + SEPS)*grad*area/Delta;
      }
  
    }
/////////////////////////////////////////
//get from known solid surface temperature
////////////////////////////////////////////
    foreach(){
       if ((!(css_test3[] >= 1. || css_test3[] <= 0.)) && (intersect_true[]==0)){
              coord n = facet_normal (point, css_test3, fss_test3);
              coord p;
              double alpha = plane_alpha (css_test3[], n);
              double area = plane_area_center (n, alpha, &p);
              if (metric_embed_factor)  
                  area *= metric_embed_factor (point, p);
               double distance=HUGE;
               if((ff[]<1.0 && ff[]>0.0)){
                     coord nl = mycs (point, ff);
                      coord pl;
                      double alphal = plane_alpha (ff[], nl);
                      double area = plane_area_center (nl, alphal, &pl);
                      if (metric_embed_factor)  
                            area *= metric_embed_factor (point, pl);
                      distance =  sqrt(sq(p.x-pl.x)+sq(p.y-pl.y));
                }
                 Point me = point;
                 double xx = p.x*Delta+x;
                 double yy = p.y*Delta+y;
          //!first for T_modl==-1
          if(fabs(T_modl[]+1)<1e-12 && (ff[]>0.0)){
                double total=0.0;
                double weight=0.0;
                // //  if((ff[]<1.0 && ff[]>0.0)){
                // //     weight += 1.0/distance;
                // //     total += weight*Tsat00;
                // //     // printf("weight_self=%g Tsat00=%g\n",weight,Tsat00);
                // //  }
               // Point me = point;
                // foreach_neighbor(1){
                foreach_neighbor(2){
                    if(!(me.i==point.i && me.j==point.j)){
                        // if((ff[]<=1.0 && ff[]>=0.5) && T_modl[]>0.5){
                        if(T_modl[]>0.5){
                            double weight_local;
                            weight_local = 1.0/(sqrt(sq(x-xx) + sq(y - yy))/Delta);
                            weight += weight_local;
                            total += weight_local*aiml[];
                            // printf("weight_local=%g aiml=%g\n",weight_local,aiml[]);
                        // }else if((ff[]<0.5 && ff[]>0.0) && (fabs(T_modl[]+1)<1e-12)){
                        // }else if((intersect_true[]==1) && (fabs(T_modl[]+1)<1e-12)){
                     }else if((intersect_true[]==1)){
                            double weight_local;
                            weight_local = 1.0/(sqrt(sq(x-xx) + sq(y - yy))/Delta);
                            weight += weight_local;
                            total += weight_local*Tsat00;
                            // printf("weight_local=%g Tsat00=%g\n",weight_local,Tsat00);
                        }
                    }
                }
                if(weight>1e-12){
                    aiml[] = total/weight;
                    // flux_l[] = ;
                    T_modl[]=6;
                }
          }

          //then for T_modg==-1
          if(fabs(T_modg[]+1)<1e-12 && (ff[]<1.0)){
                double total=0.0;
                double weight=0.0;
                // //  if((ff[]<1.0 && ff[]>0.0)){
                // //     weight = 1.0/distance;
                // //     total = weight*Tsat00;
                // //  }
                //Point me = point;
                // foreach_neighbor(1){
                foreach_neighbor(2){
                    if(!(me.i==point.i && me.j==point.j)){
                        // if(((ff[]<=0.5 && ff[]>=0.0)) && T_modg[]>0.5){
                        if(T_modg[]>0.5){
                            double weight_local;
                            // weight_local = sqrt(sq(point.i - p.x) + sq(point.j - p.y));
                            weight_local = 1.0/(sqrt(sq(x-xx) + sq(y - yy))/Delta);
                            weight += weight_local;
                            total += weight_local*aimg[];
                        // }else if((ff[]>0.5 && ff[]<1.0) && (fabs(T_modg[]+1)<1e-12)){
                        // }else if((intersect_true[]==1) && (fabs(T_modg[]+1)<1e-12)){
                         }else if((intersect_true[]==1)){
                            double weight_local;
                            weight_local = 1.0/(sqrt(sq(x-xx) + sq(y - yy))/Delta);
                            weight += 1.0/weight_local;
                            total += weight_local*Tsat00;
                        }
                    }
                }
                 if(weight>1e-12){
                    aimg[] = total/weight;
                    // flux_g[] = ;
                    T_modg[]=6;
                }
          }//if
          
      }//if
    }//foreach


     foreach(){
       if ((!(css_test3[] >= 1. || css_test3[] <= 0.)) && (intersect_true[]==1)){
              coord n = facet_normal (point, css_test3, fss_test3);
              coord p;
              double alpha = plane_alpha (css_test3[], n);
              double area = plane_area_center (n, alpha, &p);
              if (metric_embed_factor)  
                  area *= metric_embed_factor (point, p);
               double distance=HUGE;
               if((ff[]<1.0 && ff[]>0.0)){
                     coord nl = mycs (point, ff);
                      coord pl;
                      double alphal = plane_alpha (ff[], nl);
                      double area = plane_area_center (nl, alphal, &pl);
                      if (metric_embed_factor)  
                            area *= metric_embed_factor (point, pl);
                      distance =  sqrt(sq(p.x-pl.x)+sq(p.y-pl.y));
                }
                 Point me = point;
                 double xx = p.x*Delta+x;
                 double yy = p.y*Delta+y;
          //!first for T_modl==-1
          if(fabs(T_modl[]+1)<1e-12 && (ff[]>0.0)){
                double total=0.0;
                double weight=0.0;
                 if((ff[]<1.0 && ff[]>0.0)){
                    weight += 1.0/distance;
                    total += weight*Tsat00;
                    // printf("weight_self=%g Tsat00=%g\n",weight,Tsat00);
                 }
               // Point me = point;
                // foreach_neighbor(1){
                foreach_neighbor(2){
                    if(!(me.i==point.i && me.j==point.j)){
                        // if((ff[]<=1.0 && ff[]>=0.5) && T_modl[]>0.5){
                        if(T_modl[]>0.5){
                            double weight_local;
                            weight_local = 1.0/(sqrt(sq(x-xx) + sq(y - yy))/Delta);
                            weight += weight_local;
                            total += weight_local*aiml[];
                            // printf("weight_local=%g aiml=%g\n",weight_local,aiml[]);
                        }
                        // else if((ff[]<0.5 && ff[]>0.0) && (fabs(T_modl[]+1)<1e-12)){
                        //     double weight_local;
                        //     weight_local = 1.0/(sqrt(sq(x-xx) + sq(y - yy))/Delta);
                        //     weight += weight_local;
                        //     total += weight_local*Tsat00;
                        //     // printf("weight_local=%g Tsat00=%g\n",weight_local,Tsat00);
                        // }
                    }
                }
                if(weight>1e-12){
                    aiml[] = total/weight;
                    // flux_l[] = ;
                    T_modl[]=6;
                }
          }

          //then for T_modg==-1
          if(fabs(T_modg[]+1)<1e-12 && (ff[]<1.0)){
                double total=0.0;
                double weight=0.0;
                 if((ff[]<1.0 && ff[]>0.0)){
                    weight = 1.0/distance;
                    total = weight*Tsat00;
                 }
                //Point me = point;
                // foreach_neighbor(1){
                foreach_neighbor(2){
                    if(!(me.i==point.i && me.j==point.j)){
                        // if(((ff[]<=0.5 && ff[]>=0.0)) && T_modg[]>0.5){
                        if(T_modg[]>0.5){
                            double weight_local;
                            // weight_local = sqrt(sq(point.i - p.x) + sq(point.j - p.y));
                            weight_local = 1.0/(sqrt(sq(x-xx) + sq(y - yy))/Delta);
                            weight += weight_local;
                            total += weight_local*aimg[];
                        }
                        // else if((ff[]>0.5 && ff[]<1.0) && (fabs(T_modg[]+1)<1e-12)){
                        //     double weight_local;
                        //     weight_local = 1.0/(sqrt(sq(x-xx) + sq(y - yy))/Delta);
                        //     weight += 1.0/weight_local;
                        //     total += weight_local*Tsat00;
                        // }
                    }
                }
                 if(weight>1e-12){
                    aimg[] = total/weight;
                    // flux_g[] = ;
                    T_modg[]=6;
                }
          }//if
          
      }//if
    }//foreach
////////////////////////////////////////
// get flux for -1
//////////////////////////////////////
if(1==0){
foreach(){
     
      if (!(css_test3[] >= 1. || css_test3[] <= 0.)){
              coord n = facet_normal (point, css_test3, fss_test3);
              coord p;
              double alpha = plane_alpha (css_test3[], n);
              double area = plane_area_center (n, alpha, &p);
              if (metric_embed_factor)  
                  area *= metric_embed_factor (point, p);

                  normalize(&n);

                  coord nn1;
                  if(!structure_normal_flag){
                    nn1.x = direction_Ts.x[];
                    nn1.y = direction_Ts.y[];
                  }else{
                    nn1.x = n.x;
                    nn1.y = n.y;
                  }
                  double value_temp1= sqrt(sq(nn1.x)+sq(nn1.y));
                  //nn1 is the direction of temperature higher, now transfer it to direction of the phase. only in solidface cell
                 

                  coord nn2;
                  if(!structure_normal_flag){
                    nn2.x = direction_Tfl.x[];
                    nn2.y = direction_Tfl.y[];
                 }else{
                    nn2.x = -n.x;
                    nn2.y = -n.y;
                 }
                  double value_temp2= sqrt(sq(nn2.x)+sq(nn2.y));
                 

                  coord nn3;
                  if(!structure_normal_flag){
                    nn3.x = direction_Tfg.x[];
                    nn3.y = direction_Tfg.y[];
                 }else{
                    nn3.x = -n.x;
                    nn3.y = -n.y;
                 }
                  double value_temp3= sqrt(sq(nn3.x)+sq(nn3.y));
                  
                //n is the normal out of solid
                //nn1 changes it's direction to the side same with out of solid
                //nn2 is the normal direction out of fluid ?????????????????????
                if((fabs(T_modl[]+1)<1e-12 && fabs(aiml[])>1e-12) && ((value_temp1>1e-20) && (value_temp2>1e-20))){
                         double same_direction = n.x*nn1.x + n.y*nn1.y; 
                        if(same_direction<0.0){
                            nn1.x = -nn1.x;
                            nn1.y = -nn1.y;
                        }
                         same_direction = (-n.x)*nn2.x + (-n.y)*nn2.y; // -n.x -n.y is the direction out fluid
                        if(same_direction<0.0){
                            nn2.x = -nn2.x;
                            nn2.y = -nn2.y;
                        }
                        // double grad = gradient_conjugate (point,s,ss,nn1,nn2,p,&mod_temp);
                        bool is_liquid=true;
                        double flux1=0.0;
                        coord pp;
                        foreach_dimension(){
                          pp.x=p.x;
                        }
                        int mod_temp=-1;
                        // gradl=gradient_conjugate (point, s, ss, &nn1, &nn2, &pp, &mod_templ,is_liquid,&flux1);
                        gradient_conjugate_for_solid (point, s, aiml, &nn1, &pp, &mod_temp,is_liquid, &flux1);
                        if(mod_temp==6){
                            flux_l[] = flux1;//Tks*flux1;
                            T_modl[] = mod_temp;
                        }
                }

                if((fabs(T_modg[]+1)<1e-12 && fabs(aimg[])>1e-12) &&  ((value_temp1>1e-20) && (value_temp3>1e-20))){
                        double same_direction = n.x*nn1.x + n.y*nn1.y;
                        if(same_direction<0.0){
                            nn1.x = -nn1.x;
                            nn1.y = -nn1.y;
                        }
                        same_direction = (-n.x)*nn3.x + (-n.y)*nn3.y; // -n.x -n.y is the direction out fluid
                        if(same_direction<0.0){
                            nn3.x = -nn3.x;
                            nn3.y = -nn3.y;
                        }
                        // double grad = gradient_conjugate (point,s,ss,nn1,nn2,p,&mod_temp);
                        bool is_liquid=false;
                        double flux1=0.0;
                        coord pp;
                        foreach_dimension(){
                          pp.x=p.x;
                        }
                         int mod_temp=-1;
                        // gradg=gradient_conjugate (point, s, sss, &nn1, &nn3, &pp, &mod_tempg, is_liquid, &flux1);
                        gradient_conjugate_for_solid (point, s, aimg, &nn1, &pp, &mod_temp,is_liquid, &flux1);
                         if(mod_temp==6){
                            flux_g[] = flux1;//Tks*flux1;
                            T_modg[] = 6;
                         }
                }
      }
  
    }
}

}//function


// void get_solidsurfaceT(scalar topo_mask, scalar topo_mask_s,int level_interface){
 void  get_solidsurfaceT(){
    //1.seach interface cell
    get_topomask_soilid(topo_mask_s,css_test3,level_interface);
    //2.gradient of the field in region topomask_s=0,+-1.+-2
    // vector Tlff_l[],Tgff_l[],T_solid_g[];
    bool flag;
    flag=false; //fluid
    gradient_ff({Tl},{Tlff_g},ff,css_test3,topo_mask,topo_mask_s,flag);
    flag=true; //gas
    gradient_ff({Tg},{Tgff_g},ff,css_test3,topo_mask,topo_mask_s,flag);
    flag=false; //solid 
    gradient_s({Ts},{T_solid_g},css_test3,topo_mask_s,flag);

        bool flag4;
        flag4=true;
        direction_ff(Tgff_g,ff,css_test3,topo_mask,topo_mask_s,direction_Tfg,flag4);
        flag4=false;
        direction_ff(Tlff_g,ff,css_test3,topo_mask,topo_mask_s,direction_Tfl,flag4);
        flag4=false;
        direction_s(T_solid_g,css_test3,direction_Ts,flag4);

     embed_flux_conjugate(Ts, Tl, Tg, aiml, aimg, T_modl,T_modg, flux_l, flux_g);
}


////////////////////////////////////////////////////////////
//// merge_part
///////////////////////////////////////////////////////////

// extern scalar merge_to_me_s_c,merge_to_me_s_energy;
// extern vector merge_to_me_s_position;

// extern scalar merge_to_me_l_c,merge_to_me_l_energy;
// extern vector merge_to_me_l_position;

// extern scalar merge_to_me_g_c,merge_to_me_g_energy;
// extern vector merge_to_me_g_position;
//solid surface
//1.each surface<0.5 get merged to fraction>1

#define IS_CONDITION_MET(ii, jj) ((ii==0 && jj==-1) || (ii==1 && jj==0) || (ii==0 && jj==1) || (ii==-1 && jj==0))

bool is_f_small_part(double f_value) {
    if ( (f_value > 0.0 && f_value < 0.5)) {
        return true;
    }
    return false;
}

bool is_f_big_part(double f_value) {
    if ((f_value >= 0.5 && f_value<=1.0)) {
        return true;
    }
    return false;
}
// input inverse and ff5 get the merged cell
void merge_get_c_position(bool inverse, scalar ff5, scalar merge_c,vector merge_position){
        foreach(){
            merge_c[] = 0;
            merge_position.x[]=0;
            merge_position.y[]=0;
            // if(level==level_interface && is_f_small_part(ff5[],inverse) &&  (!is_boundary(cell))){
                if(level==level_interface && (is_f_small_part(ff5[]))){
                coord mixed;
                mixed.x=x,mixed.y=y;//,mixed.z=z;
                // coord nn=interface_normal7(point,ff5,hh2);;
              coord nn =  mycs (point,ff5); //2023050603
              coord pp;
              double alpha2 = plane_alpha (ff5[], nn);
              double area = plane_area_center (nn, alpha2, &pp);
              if (metric_embed_factor)  
                  area *= metric_embed_factor (point, pp);
              pp.x = pp.x*Delta + mixed.x;
              pp.y = pp.y*Delta + mixed.y;
              normalize(&nn);

                int pppx,pppy;
                double colinearity=0.0;
                foreach_neighbor(1){
                 if(level==level_interface && (is_f_big_part(ff5[]))){   
                    coord pure;
                    pure.x = x, pure.y = y;//, pure.z = z;
                    int ii=round((pure.x-mixed.x)/Delta);
                    int jj=round((pure.y-mixed.y)/Delta);
                    bool flag2;
                    flag2 = IS_CONDITION_MET(ii,jj);
                        if(flag2){      
                            coord puretomixed;
                            puretomixed.x = pure.x - pp.x;
                            puretomixed.y = pure.y - pp.y;
                            double leng1 = sqrt(sq(puretomixed.x)+sq(puretomixed.y));//+sq(puretomixed.z));
                            double colinearity_temp;
                            colinearity_temp = fabs(nn.x*puretomixed.x/leng1+nn.y*puretomixed.y/leng1); // before 20230614
                            if(colinearity<colinearity_temp){
                                colinearity = colinearity_temp;
                                pppx = ii;
                                pppy = jj;
                            }
                                        
                        }
                    }
                } //foreach_neighbor
                if(colinearity>0){
                    merge_position.x[] = pppx;
                    merge_position.y[] = pppy;
                    // merge_c[] = inverse ? (1.0 - ff[]) : ff[];
                    merge_c[] = ff5[];
                }
                
            }//topo_mask
    }

} 


//assgn temperature
// outside in .c: Ts,Tl,Tg

/// not residual, not relax; but after final relax, adding this new solver
#include "./poissonsolver.h"
// #include "./poissonsolver-face-slg.h"
void solver_new(scalar poisson_s,double percent_s,double percent_l,double percent_g){

      //get merge_c and merge_energy for Ts,Tl,Tg
      bool flag=false;
      merge_get_c_position(flag, css_test3,merge_to_me_s_c,merge_to_me_s_position);
      merge_get_c_position(flag, css_test,merge_to_me_l_c,merge_to_me_l_position);
    //   merge_get_c_position(true, css_test2,merge_to_me_g_c,merge_to_me_g_position);
      merge_get_c_position(flag, css_test2,merge_to_me_g_c,merge_to_me_g_position);
      //get temperature for cells with level==level_interface && 0<f[]<0.5 = temperature of its merge aim 

         //check get_solidsurfaceT

//    foreach(){
//     //   Tl[] = 0.0; //Tsat00;
//     //   Tg[] = 0.0; //Tsat00;
//     //   Ts[] = 0.0; //Tsat00;
//       Tl[] = Tsat00;
//       Tg[] = Tsat00;
//       Ts[] = Tsat00;
//       if(css_test3[]>=0.5){
//             Ts[] = T[];
//       }else{
//           if(css_test[]>=0.5){
//               Tl[] = T[];
//           }
//           if(css_test2[]>=0.5){
//               Tg[] = T[]; 
//           }
//       }
//    }

foreach(){
   if(fabs(merge_to_me_s_c[])>1e-7){
       int ii,jj;
       ii = merge_to_me_s_position.x[];
       jj = merge_to_me_s_position.y[];
       Ts[] = Ts[ii,jj];
   }
}

foreach(){
   if(fabs(merge_to_me_l_c[])>1e-7){
       int ii,jj;
       ii = merge_to_me_l_position.x[];
       jj = merge_to_me_l_position.y[];
       Tl[] = Tl[ii,jj];
   }
}

foreach(){
   if(fabs(merge_to_me_g_c[])>1e-7){
       int ii,jj;
       ii = merge_to_me_g_position.x[];
       jj = merge_to_me_g_position.y[];
       Tg[] = Tg[ii,jj];
   }
}

// assgn restriction function
Ts.ff6 = css_test3;
Ts.restriction = restriction_zero;
Tl.ff6 = css_test;
Tl.restriction = restriction_zero;
Tg.ff6 = css_test2;
Tg.restriction = restriction_zero;

Ts.refine = Ts.prolongation = refine_embed_linear_css_test3;
Tl.refine = Tl.prolongation = refine_embed_linear_css_test;
Tg.refine = Tg.prolongation = refine_embed_linear_css_test2;


restriction({Ts,Tl,Tg});


poisson_solver(poisson_s,percent_s,percent_l,percent_g);

}



//////////////////////////////////////////////
////////////////////////////////////////////
////////////get-css-fss-area
//////////////////////////////////////////
/////////////////////////////////////////////
//areasl,areasg,arealg
//css_test,css_test2,css_test3
//fss_test,fss_test2,fss_test3
#include "./line9-4-4-for-basilisk-css-test.h"

 extern   scalar temp27;
  extern    scalar temp28;
  extern    scalar temp211;
   extern   scalar temp212;

   extern   scalar temp25;
   extern   scalar temp26;
    extern  scalar temp29;
    extern  scalar temp210;
void get_css_fss_areaslg_triple_point(){
    // scalar temp27[];
    // scalar temp28[];
    // scalar temp211[];
    // scalar temp212[];

    // scalar temp25[];
    // scalar temp26[];
    // scalar temp29[];
    // scalar temp210[];
    foreach(){
        areasl[]=0.0;
        areasg[]=0.0;
        arealg[]=0.0;

        scalar temp27[]=0.0;
        scalar temp28[]=0.0;
        scalar temp211[]=0.0;
        scalar temp212[]=0.0;

        scalar temp25[]=0.0;
        scalar temp26[]=0.0;
        scalar temp29[]=0.0;
        scalar temp210[]=0.0;
    }
    //initial 
    foreach(){
        //use ff not css_test and nor css_test2
        if(css_test3[]<1.0){
            if(ff[]>0.0 && ff[]<1.0){
                coord nlg = mycs (point, ff);//interface_normal(point, ff);//interface_normal7(point,ff,hhh); //
                double alphalg = plane_alpha (ff[], nlg);
                coord plg;
                double arealg_temp = plane_area_center (nlg, alphalg, &plg);
                if (metric_embed_factor)  
                    arealg_temp *= metric_embed_factor (point, plg);
                arealg[] = arealg_temp;
            }
        }
        if(css_test3[]>0.0 && css_test3[]<1.0){
            coord nsf = mycs (point, css_test3);//interface_normal(point, ff);//interface_normal7(point,ff,hhh); //
            double alphasf = plane_alpha (css_test3[], nsf);
            coord psf;
            double areasf_temp = plane_area_center (nsf, alphasf, &psf);
            if (metric_embed_factor)  
                areasf_temp *= metric_embed_factor (point, psf);
                
            //liquid existing
            if(ff[]>0.0){
                areasl[] = areasf_temp;
            }
            if(ff[]<1.0){
                    areasg[] = areasf_temp;
            }
        }
    }

    foreach(){
        if(css_test3[]>0.0 && css_test3[]<1.0 && ff[]>0.0 && ff[]<1.0){
            coord nsf = mycs (point, css_test3);//interface_normal(point, ff);//interface_normal7(point,ff,hhh); //
            double alphasf = plane_alpha (css_test3[], nsf);
            coord psf;
            double areasf_temp = plane_area_center (nsf, alphasf, &psf);
            if (metric_embed_factor)  
                areasf_temp *= metric_embed_factor (point, psf);

            coord nlg = mycs (point, ff);//interface_normal(point, ff);//interface_normal7(point,ff,hhh); //
            double alphalg = plane_alpha (ff[], nlg);
            coord plg;
            double arealg_temp = plane_area_center (nlg, alphalg, &plg);
            if (metric_embed_factor)  
                arealg_temp *= metric_embed_factor (point, plg);
            //in basilisk
            //solid: nsf.x*x + nsf.y*y = alphasf
            //fg surface: nlg.x*x + nlg.y*y = alphalg
            //origin in the center of the square

            //but in my code line9
            //in line-9-4-for-basilisk: n.x*x + n.y*y + alpha = 0
            //so there is some difference between alphas;
            //origin in the left bottom corner of the square

            //coordinate transfer
            //x1 = x + 0.5, y1 = y+0.5
            //x=x1-0.5; y=y1-0.5
            //nlg.x*x + nlg.y*y = alphalg --> nlg.x*(x1-0.5) + nlg.y*(y1-0.5) = alphalg
            // -->nlg.x*x1 + nlg.y*y1-0.5*(nlg.x+nlg.y)-alphalg = 0
            double a1,b1,c1,a2,b2,c2;
            a1 = nsf.x,b1=nsf.y,c1=-0.5*(nsf.x+nsf.y)-alphasf;
            a2 = nlg.x,b2=nlg.y,c2=-0.5*(nlg.x+nlg.y)-alphalg; 
            //a1*x+b1*y+c1=0,a2*x+b2*y+c2=0;
            // printf("a1=%g,b1=%g,c1=%g,a2=%g,b2=%g,c2=%g\n",a1,b1,c1,a2,b2,c2);
            double data2[13];
            // cut_line_test_in_basilisk(a1,b1,c1,a2,b2,c2,data2);
            cut_line_test_in_basilisk((coord){a1,b1},c1,(coord){a2,b2},c2,data2);
            css_test2[] = data2[0]; //=css_test2_value;
            css_test[]  = data2[1]; //=css_test_value;
            arealg[]    = data2[2]; //=arealg_value;
            areasl[] = data2[3]; //=areasl_value;
            areasg[] = data2[4]; //=areasg_value;
            // fss_test.x[]  = data2[5]; //=fss_test_small_x;
            // fss_test.y[]  = data2[6]; //=fss_test_small_y;
            // fss_test.x[1] = data2[7]; //=fss_test_big_x;
            // fss_test.y[1] = data2[8]; //=fss_test_big_y;
            // fss_test2.x[] = data2[9]; //=fss_test2_small_x;
            // fss_test2.y[] = data2[10]; //=fss_test2_small_y;
            // fss_test2.x[1] = data2[11]; //=fss_test2_big_x;
            // fss_test2.y[1] =data2[12]; //=fss_test2_big_y;
            // Point me = point;
            // coord me_coord;
            // me_coord.x = x;
            // me_coord.y = y;
            // double data27=data2[7];
            // double data28=data2[8];
            // double data211=data2[11];
            // double data212=data2[12];
            temp27[]=data2[7];
            temp28[]=data2[8];
            temp211[]=data2[11];
            temp212[]=data2[12];

            temp25[]=data2[5];
            temp26[]=data2[6];
            temp29[]=data2[9];
            temp210[]=data2[10];
            // foreach_neighbor(1){
            //     int delta_x = round((x - me_coord.x)/Delta);
            //     int delta_y = round((y - me_coord.y)/Delta);
            //     if(delta_x==1 && delta_y==0 && (!is_boundary(cell))){
            //         fss_test.x[] = 1;//data27; //data2[7];
            //         fss_test2.x[] = 1;//data211; //data2[11];
            //     }
            //     if(delta_x==0 && delta_y==1 && (!is_boundary(cell)) ){
            //         fss_test.y[] = 1;//data28; //data2[8];
            //         fss_test2.y[] = 1;//data212; //data2[12];
            //     }
            // }
            if(fabs(arealg[])>1e-30){
                if (metric_embed_factor){  
                     arealg[] *= y;//metric_embed_factor (point, plg);
                }
            }
            if(fabs(areasg[])>1e-30){
                if (metric_embed_factor){  
                     areasg[] *= y;//metric_embed_factor (point, plg);
                }
            }
             if(fabs(areasl[])>1e-30){
                if (metric_embed_factor){  
                     areasl[] *= y;//metric_embed_factor (point, plg);
                }
            }



        }
    }

    foreach_face(x){
        if(level==level_interface){
            if(fabs(temp27[-1,0])>1e-20){
                fss_test.x[] = temp27[-1,0];
            }
            // if(fabs(temp28[0,-1])>1e-20){
            //     fss_test.y[] = temp28[0,-1];
            // }
            if(fabs(temp211[-1,0])>1e-20){
                fss_test2.x[] = temp211[-1,0];
            }
            // if(fabs(temp212[0,-1])>1e-20){
            //     fss_test2.y[] = temp212[0,-1];
            // }

            if(fabs(temp25[])>1e-20){
                fss_test.x[] = temp25[];
            }
            // if(fabs(temp26[])>1e-20){
            //     fss_test.y[] = temp26[];
            // }
            if(fabs(temp29[])>1e-20){
                fss_test2.x[] = temp29[];
            }
            // if(fabs(temp210[])>1e-20){
            //     fss_test2.y[] = temp210[];
            // }
        }
    }

    foreach_face(y){
        if(level==level_interface){
            // if(fabs(temp27[-1,0])>1e-20){
            //     fss_test.x[] = temp27[-1,0];
            // }
            if(fabs(temp28[0,-1])>1e-20){
                fss_test.y[] = temp28[0,-1];
            }
            // if(fabs(temp211[-1,0])>1e-20){
            //     fss_test2.x[] = temp211[-1,0];
            // }
            if(fabs(temp212[0,-1])>1e-20){
                fss_test2.y[] = temp212[0,-1];
            }

            // if(fabs(temp25[])>1e-20){
            //     fss_test.x[] = temp25[];
            // }
            if(fabs(temp26[])>1e-20){
                fss_test.y[] = temp26[];
            }
            // if(fabs(temp29[])>1e-20){
            //     fss_test2.x[] = temp29[];
            // }
            if(fabs(temp210[])>1e-20){
                fss_test2.y[] = temp210[];
            }
        }
    }
}