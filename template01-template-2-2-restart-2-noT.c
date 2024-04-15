//#include "grid/quadtree.h"
//#include "embed.h"
#include "./embed-update-2.h"
//#include "axi.h"
//#include "embed-update.h"
#include "./axi-official-update.h"
//#include "embed-update.h"
//#include "embed.h"
//#include "utils.h"
//#include "run.h"
//#include "navier-stokes/centered.h"
#include "./axi-centered.h"
//#include "two-phase.h"
#include "axi-two-phase.h"
#include "./embed-solid.h"
#include "./linear2-tree-1-2.h"
//#include "navier-stokes/conserving.h"
//#include "tension.h"
//#include "reduced.h"
//#include "henry.h"
// #include "./axi-diffusion1-2.h"
#include "./axi-diffusion1.h"
#include "./diffusion5.h"
#include "./my-diffusion.h"  //for project_source
#include "./my-tension.h"
#include "./reduced.h"
#include "./poisson-ps-usf2-2.h"     //new
#include "./poisson-ps-usf2g.h"     //new
// #include "./my-contact-embed.h"
#include "./my-contact-embed-triple.h"
#include "./axi-css_test.h"
#include "./axi-css_test2.h"
#include "./axi-css_test3.h"
#include "./getsolid.h"

// #include "tag.h"
// #include "./my-tag.h"

#ifndef M_PI
  #define M_PI 3.141592653589793238
#endif

//#include "run.h"
// #include "view.h"

// #define BASE_TIME 0.0892672
// #define BASE_TIME 0
// #define BASE_TIME 0.0943463
// #define BASE_TIME 0.105438
#define BASE_TIME 0.0158541

clock_t start_time,end_time;
double during_time;
double originx=0.0,originy=0.0;
double L0_pysical;
int case_number = 3; //casenumber==2 : tanguy
int sink_number=0;
// int maxl= 11, minl = 4;
int maxl= MESHN , minl = 4;
bool restartsymbol = false;//true;

int globali,outstep=1;
double out_interval=5e-7;//0.01;//1e-6;//0.01;//1e-6;//2e-8;//1e-6;//1e-3;//0.1;
// double out_interval2=5e-6;//2e-8;//5e-6;
double out_interval2=5e-7;//5e-6;//0.01;//1e-7;;//0.01;//1e-5;//0.01;//1e-6;//1e-7;//5e-6;//1e-7;//2e-8;//5e-6;
double tend = 0.00022 + BASE_TIME;
// double tend = 0.15 + BASE_TIME;
// double tend = 0.0006 + 0.0892672;;//0.0892684;//1.5e-5 + 0.0892672;
//;//0.0006 + 0.0892672;//0.090481;//0.1;//0.102172; //0.5;//0.062169;//0.5;//1.0;
//double Tkg,Tkl,Tcpg,Tcpl,hfg,Trhol,Trhog;
double Tsat0,Tsat00,T_inf,T_sup;
double bubble_radius;
double centerx,centery;
double Tsub,Twall_init;
double Tmax=0.0;
bool tune_cs_flag2=true;
double tune2_value=0.0;
double Length_heat=0.004;
double Rwater;// = 461.52277955158064;
double omega;// = 0.0345l
double Rcc;
bool Rcc_flag_l = false; //default
bool Rcc_flag_g = false; //default
scalar aiml_s[],aimg_s[];

scalar T[];
scalar Tlff[],Tgff[],T_solid[];
face vector modphase0[],modphase1[];
face vector modphase_s_1[],modphase_s_0[];
vector hhh[];

scalar corner_ff[];

vector direction_Tfg[],direction_Tfl[],direction_Ts[];
//vector hhh_ff_oppo[];
#ifndef rho_func
# define rho_func(ff) (clamp(ff,0.,1.)*(Trhol - Trhog) + Trhog)
#endif

scalar topo_mask[];
scalar topo_mask_g[];

scalar topo_mask_s[];

scalar masstr[], source_pc[], vtr[], source_pc2[];

//scalar css_test3_n[];
scalar poisson_source2[];
//when add navier stokes
scalar ps[],psg[];
face vector usf[],usfg[],ulf[],ugf[];
//gradient of Temperature
 scalar phase0Tgrad[],phase1Tgrad[];
//  scalar tempgg[],tempgl[];
 scalar flux_lg_l[],flux_lg_g[];

//  vector number_s_f[];

double source_total,total_area;
// double dt;
int phase_flag3;
int level_interface;
face vector fss_test[];
scalar css_test[];
face vector fss_test2[];
scalar css_test2[];
//for solid
scalar fs_solid[];
scalar css_test3[];
face vector fss_test3[];
//for not solid
scalar css_test3_n[];
face vector fss_test3_n[];
scalar deltac[];
double CFL_evap = 0.1;
bool dump_each_event=false;//true;
int dump_each_event_interval=1;
bool point_trace=false;
double tracex,tracey;
double delta_min;
bool out_flag;
bool energy_advecting_flag;

bool flag_average_source=false;//false;//true;//false;
bool flag_cant_smaller_than_half=false;//false;//true;//false;//true;
bool flag_get_u_ghost=false;//true;//false;

scalar T_oold[];
char Tini_file[80];
 scalar smf[],ssmf[],sssmf[],ssssmf[];
//  scalar css_test3_n[];
 double thickness;
 double theta0=  CONTACTANGLE ;

 bool restart_Tsat=false;
 int restart_Tsat_i=0;

// scalar temp_vof[];
 vector smallmodl[],bigmodl[];
 vector smallmodg[],bigmodg[];
 scalar resl[],resg[];
 scalar ress[];

 vector flux_show[];
//  scalar flux_x[],flux_y[];
scalar intersect_true[];

 bool poisson_check=true;//false;
 bool surface_heat=false;//true;
 bool surface_heat_restart=false;//true;
//T[left] = dirichlet(T_inf);
//T[left] = dirichlet(Tsat00);
//T[bottom] = T_inf; //neumann(0.);
bool use_Tslg=true;


scalar res_ps[];


//temp scalar
vector Tlff_g[],Tgff_g[],T_solid_g[];
scalar Tl[],Tg[],Ts[];
scalar aiml[], aimg[],T_modl[], T_modg[];
scalar flux_l[],flux_g[]; //flux +: out solid;
scalar areasl[],areasg[],arealg[];

scalar merge_to_me_s_c[];
vector merge_to_me_s_position[];
scalar merge_to_me_s_energy[];

scalar merge_to_me_l_c[];
vector merge_to_me_l_position[];
scalar merge_to_me_l_energy[];

scalar merge_to_me_g_c[];
vector merge_to_me_g_position[];
scalar merge_to_me_g_energy[];

vector ulf_v[];
scalar modify_near_region[];
    scalar temp27[];
    scalar temp28[];
    scalar temp211[];
    scalar temp212[];

    scalar temp25[];
    scalar temp26[];
    scalar temp29[];
    scalar temp210[];

     scalar ff_old[],ff_old2[];

     scalar div_numerical[];
     scalar div_numerical2[];

     scalar poisson_source2_f[]; 
/////////case_number == 1 // scriven leon test
    //  u.n[right] =neumann(0.);
    // p[right]   = dirichlet(0.);
    // pf[right]  = dirichlet(0.);

    //       u.n[top] = neumann(0.);
    // p[top]   = dirichlet(0.);
    // pf[top]  = dirichlet(0.);

//////////////////case_number == 2
          //      u.n[embed] = dirichlet (0.);
          //   u.t[embed] = dirichlet (0.);
          // ulf.n[bottom] = 0.;
          // ulf.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which


          //   //u.r[embed] = dirichlet (0.);

          //   T[left] = dirichlet(T_inf);
          //   T[right] = dirichlet(Tsub);

          //      //be careful, in tanguy's test, this should be wall boundary condition.
          //     //  u.n[top] =neumann(0.);
          //     // p[top]   = dirichlet(0.);
          //     // pf[top]  = dirichlet(0.);

          //      u.n[right] =neumann(0.);
          //     p[right]   = dirichlet(0.);
          //     pf[right]  = dirichlet(0.);
//////////////////case_number == 3 -- boundary condition
   //left
        //heat flux 425000 = -k_s \times \frac{\partial{T}}{\partial{n}}
          double heat_flux;
          double range_heat_flux;

bool flag_topos_advect_uf=false;

  double thickbottom;
   double groove_h, groove_w_bottom,groove_w_top;

  //  int mesu_case_number=178;
  // topo_mask_s==0 using u advect not ulf+mov; vof.h, mov_interface, poisson_us_2.h

  
        //step heat flux'  -- tempalte01-template2.c
          // T[left] = neumann(y<range_heat_flux?(heat_flux/Tks):0.0); //step heat flux, heat_flux = 425000
        //linear heat flux -- tempalte01-template2-2.c
         // T[left] = neumann(max(heat_flux*(1.0-y/0.004)/Tks,0.0)); //linear heat flux, heat_flux = 481000
        //default u and p boundary condition
   //right
        //zero heat flux = defalut symetric
        // outflow boundary for u and p
          u.n[right] =neumann(0.);
          p[right]   = dirichlet(0.);
          pf[right]  = dirichlet(0.);
   //bottom
        //bottom axi-symmetric axis = default
   //top
        //default symmetric boundary condition for T
        //outflow for u and p
        u.n[top] =neumann(0.);
        p[top]   = dirichlet(0.);
        pf[top]  = dirichlet(0.);

        u.n[embed] = dirichlet (0.);
        u.t[embed] = dirichlet (0.);
        // u.t[embed] = neumann(0.);
        // // p[embed] = neumann (0.);

        ps[embed] = neumann (0.);
        p[embed] = neumann (0.);



int main (int argc, char **argv)
{
   start_time = clock();
end_time = start_time;

  if(case_number==1){
    init_grid(1<<minl);
      DT = 0.005;
      double L0_pysical = 0.5;
      L0 = L0_pysical;
      Tkg =  0.007;  //thermal conductivity
      Tkl = 0.07;
      Trhog = 0.25; //density
      Trhol = 2.5;
      Tcpg = 1.0;  //specific heat capicity
      Tcpl = 2.5;
      hfg = 100.0;
      Tsat0 = 1.0;//3.0; //3.0;
      Tsat00 = Tsat0;
      T_sup = 2.0;
      T_inf =  Tsat0 + T_sup;
      ff.sigma = 0.001;

      bubble_radius = 0.12;
       sprintf(Tini_file,"leontest-3-5000-0.dat");
      centerx=0.0;
      centery=0.0;
  }else if(case_number==2){ //tangug copy from case_number == 8
          L0_pysical = 0.0045;//1.0; //1.0; //8.0; // 8.0;
          L0 = L0_pysical;
          printf("before init_grid,pid=%d\n",pid());

          originx = -0.0015-L0_pysical/((1<<maxl)*2.0);
          originy=0.0;
          origin ( originx, originy); //the domain is now bigger //origin (-0.00225, -0.00225, -0.0015);
          init_grid(1<<minl);
          printf("after init,pid=%d\n",pid());
          DT = 1e-7;//2e-7;//5e-7;//5e-8;//1e-7;//5e-7 ;//1e-7;//5e-7;//5e-8;// 5e-8;// 1e-7;//5e-8;//1e-7; //1e-7; //0.0005 ; //0.0005 for rhol/rhog=100;//0.005;

          ///////table 12.4 of Leon's thesis
          Tkg =  0.024; //0.025; //0.007;  //thermal conductivity
          Tkl = 0.677; //0.6;  //0.07;
          Trhog = 0.5974;//0.6; //0.0025; //density
          Trhol = 958.0; //958.0; //2.5;
          Tcpg = 2034.0;//2080.0; //1.0;  //specific heat capicity
          Tcpl = 4216.0;//4216.0; //2.5;
          Tsat0 = 373.12;//373.15; //1.0;//3.0; //3.0;
          Tsat00 = Tsat0;
          T_sup = 2.0;//10.0; //2.0;
          Tsub= Tsat00 - 10.0;
          //T_inf = ;//1.0; //1.0; //= Tsat0 + T_sup;
          Twall_init =  Tsat0 + T_sup;

          hfg = 2.256e6;//2.256e6; //100.0;
          Trhos= 1.4;//8910.0; //10.0;// 10.0;
          Tcps=1500.0;//4940.0; //10.0;//10.0;
          Tks=1.0; //90.9; //10.0; //10.0;//10.0;
          //Tkl*(Twall_init-Tsub)/thickness_of_liquid = Tks*(T_inf-Twall_init)/thickness_of_solid
          T_inf = Tkl/Tks*(Twall_init-Tsub)*(0.5/1.0) + Twall_init;
          ff.sigma = 0.058; //0.059;
          mu1 = 2.82*0.0001;//2.82*0.0001;
          mu2 = 1.228*0.00001;//1.23*0.00001;
          G.x = -9.8;

          k1 = Tkl/(Trhol*Tcpl);//10.0; //k1=kl/rhol_cpl
          k0 = Tkg/(Trhog*Tcpg);//10.0; //k0=kg/rhoc_cpg

           bubble_radius = 2e-4; //6e-5; //0.24;
          //position of bubble
          centerx=bubble_radius*cos(50.0*M_PI/180.0);;
          centery=0.0;
          //centerz=bubble_radius*cos(50.0*M_PI/180.0);//L0/5.0; //0.2 ;//L0/2.0+0.07;


          thickness = 0.001;//centerz;//0.5;//0.257;
          theta0=50.0;

          delta_min = L0_pysical/((1<<maxl)*1.0);

}else if(case_number==3){ // case_number==3 lubomir bucci's experiment
          // L0_pysical = 0.0014; //1.4mm
          // L0_pysical = 0.0029;//0.003;//0.0029;
          // L0 = L0_pysical;
           sink_number = 2; //not 2 but represent multiple

          tracex = 2.14648e-5;
          tracey = 2.03711e-5;


          // L0_pysical = 0.00598;//0.003;//0.0029;
          // thickness = 0.001;
          L0_pysical = 0.00098;//0.003;//0.0029;
          thickness = 0.0001;
          L0 = L0_pysical;

          Length_heat=0.004;
          bool tune_cs_flag=false;
          if(tune_cs_flag){
              // double p=0.001/L0_pysical*(1<<maxl);
              // double pp = p - round(p);
              // if(pp>=0.0){
              //     if(pp<0.2){
              //        pp=0.2;
              //     }
              // }else{
              //     if(fabs(pp)<0.2){
              //         pp=-0.2;
              //     }
              // }
              // originx = -(round(p)+pp)/(1<<maxl)*L0_pysical;
          }else{
            // originx = -0.001;
            originx = -thickness;
          }

          energy_advecting_flag = false;

//method2:tune_cs_flag2 enables the volume fraction between (cs_min_fluid<cs<1-cs_min_fluid)
          // if(tune_cs_flag2){
          //     double p=0.001/L0_pysical*(1<<maxl);
          //     double pp = p - round(p);
          //     if(pp>=0.0){
          //         if(pp<cs_min_fliud){
          //            tune2_value = (cs_min_fliud-pp)/(1<<maxl)*L0_pysical;
          //         }
          //     }else{
          //         if(fabs(pp)<cs_min_fliud){
          //             tune2_value = -(cs_min_fliud-fabs(pp))/(1<<maxl)*L0_pysical;
          //         }
          //     }

          // }
///method3:tune_cs_flag2 enables the volume fraction between (cs>cs_min_fluid)
          double cs_min_fliud = 0.9;
          if(tune_cs_flag2){
              // double p=0.001/L0_pysical*(1<<maxl);
              double p=thickness/L0_pysical*(1<<maxl);
              double pp = p - round(p);
              if(pp>0.0){
                  if((1.0-pp)<cs_min_fliud){
                     tune2_value = ((1.0-pp)-cs_min_fliud)/(1<<maxl)*L0_pysical;
                  }
              }else{
                  if(-pp<cs_min_fliud){
                      tune2_value = (-pp-cs_min_fliud)/(1<<maxl)*L0_pysical;
                  }
              }

          }

          delta_min = L0_pysical/((1<<maxl)*1.0);

          // originx = -0.0001;//-0.0005;//-0.0001; //thickness of solid is 0.1mm
          originy = 0.0;
          origin ( originx, originy); //the domain is now bigger //origin (-0.00225, -0.00225, -0.0015);
      if(surface_heat){
          init_grid(1<<minl);
      }
     if(!surface_heat){
          if(maxl == 10){
              DT = 1e-7; //need to be check
              //0.00298--3e-7
          }else if(maxl == 11){
              DT = 1e-8;//1e-9;//5e-9; //need to be check
          }else if(maxl == 12){
              DT = 2e-9;
          }else if(maxl == 8){
              DT = 1e-6;
              //0.00298-5e-6
          }else if(maxl == 9){
              DT = 1e-6;
          }else if(maxl == 13){
              DT = 1e-8;
          }else if(maxl ==14){
              DT = 2e-9;
          }   
     }else{
          if(maxl == 10){
              DT = 1e-5; //need to be check
          }else if(maxl == 11){
               DT = 1e-6; //need to be check
          }else if(maxl == 12){
              // DT = 2e-8;
          }else if(maxl == 8){
              DT = 1e-5;
          }else if(maxl == 9){
              DT = 1e-6;
          }else if(maxl == 13){
              DT = 1e-8;
          }else if(maxl ==14){
              DT = 2e-9;
          }    
     }
          thickbottom = 5e-5; //distance from heating surface to solid bottom
          // groove_h=1e-4, groove_w_bottom=1e-4, groove_w_top=1e-4;
          // groove_h=4e-5, groove_w_bottom=4e-5, groove_w_top=4e-5; //40
          // groove_h=5e-4, groove_w_bottom=5e-4, groove_w_top=5e-4; //500
          groove_h=2e-5, groove_w_bottom=2e-5, groove_w_top=2e-5; //20

          ///////table 12.4 of Leon's thesis
          Tkg =  0.0246; //thermal conductivity
          Tkl = 0.677;
          Trhog = 0.598; //density
          Trhol = 958.0;
          Tcpg = 2080.0; //specific heat capicity
          Tcpl = 4220.0;
          Tsat0 = 373.15; //need to be checked;
          Tsat00 = Tsat0;
          T_sup = 12.55;

          hfg = 2.26e6; //latent heat

          // Trhos = 4510.0; //titanium
          // Tcps = 544.0;
          // Tks = 17.0;

          Trhos = 3980.0; //sapphire
          Tcps = 929.0;
          Tks = 25.1;
          //Tkl*(Twall_init-Tsub)/thickness_of_liquid = Tks*(T_inf-Twall_init)/thickness_of_solid
          T_inf = Tsat00 + T_sup;//
         ff.sigma = 0.0589; //0.059;
          // ff.sigma = 0.0; //mesu7
          mu1 = 2.82*0.0001;//2.82*0.0001;
          mu2 = 1.22*0.00001;//1.23*0.00001;
           G.x = -9.8;

          k1 = Tkl/(Trhol*Tcpl);//10.0; //k1=kl/rhol_cpl
          k0 = Tkg/(Trhog*Tcpg);//10.0; //k0=kg/rhoc_cpg

          Rcc_flag_l = true;
          Rcc_flag_g = false;
          Rwater = 461.52277955158064;
          omega = 0.0345;
          Rcc = sqrt(Tsat0*Tsat0*Tsat0*2.0*pi*Rwater)/omega/hfg/hfg/Trhog; //eq 6.8 lubo thesis
          printf("value of pi:%g\n",pi);
          printf("value of Rcc:%g\n",Rcc);

          //  bubble_radius = 1e-4;//1e-5;//5e-4;//1e-5; // R0=0.00001
           bubble_radius = 5e-6;//1e-5;//1e-5;//5e-4;//1e-5; // R0=0.00001
          //position of bubble

          centery=0.0;
          //centerz=bubble_radius*cos(50.0*M_PI/180.0);//L0/5.0; //0.2 ;//L0/2.0+0.07;
          //thickness = 0.001;//centerz;//0.5;//0.257;
          theta0= CONTACTANGLE ;//60;//90.0; //65.65;
          centerx=bubble_radius*cos(theta0*M_PI/180.0);;
        //step heat flux'  -- tempalte01-template2.c
        // heat_flux = 0; //mesu3-6
              // surface_heat=false;//false;
              if(surface_heat){
                heat_flux = 481000.0;//425000.0;
              }else{
                heat_flux = 425000.0;//0.0;//425000.0;
              }
        //linear heat flux'  -- tempalte01-template2-2.c
        //  heat_flux = 481000.0;//425000.0;
          range_heat_flux = (1+sqrt(2.0))/2*0.0015;

}

  run();
}

double center1 (double s0, double s1, double s2) {
  return (s2-s0)/2.0;
}

event defaults (i=0) {
  ff.tracers = list_copy ({Tlff,Tgff});
  for (scalar s in ff.tracers){
    //  s.gradient = superbee;//generic_limiter;//superbee;//minmod2;
    //  s.gradient = sweby;//generic_limiter;//superbee;//minmod2;
     s.gradient =  center1;
  }
  for (scalar s in {Tgff}){
    s.inverse = true;
    s.khaki = true;
  }
  for (scalar s in {Tlff}){
    s.inverse = false;
    s.khaki = false;
 }

 for (scalar s in {Tg}){
    s.inverse = true;
    s.khaki = true;
  }
  for (scalar s in {Ts,Tl}){
    s.inverse = false;
    s.khaki = false;
 }

    #if TREE
     ulf.x.refine = refine_face;
     usf.x.refine = refine_face;
     foreach_dimension(){
        ulf.x.prolongation = refine_embed_face_x;
        usf.x.prolongation = refine_embed_face_x;

        // ugf.x.prolongation = refine_embed_face_x;
        // usfg.x.prolongation = refine_embed_face_x;
     }
    #endif


}



// #define EPS32 0.0
// struct Topo_m2 {
//     scalar topo_mask;
// };
// void get_topomask(struct Topo_m2 q) {
//      scalar topo_mask = q.topo_mask;

//      foreach(){
//         int phase_sign = (ff[]>=0.5-EPS32) ? 1 : -1;
// 	topo_mask[] = 3*phase_sign;
//         if(ff[]<1.0-EPS32 && ff[]>EPS32){
//         // if(f[]<1.0-EPS32 && f[]>EPS32 && (css_test3_n[]>0.0)){
//  		        topo_mask[] = 0;
//         }
//      }
//      boundary({topo_mask});
//      foreach(){
//         if(level==level_interface){
//        //  if(level==level_interface && (css_test3_n[]>0.0)){
//             bool is1= false;
//             foreach_dimension(){
//                int temp1=topo_mask[1];
//                int temp2=topo_mask[-1];
//             if(temp1==0 || temp2==0){
//                is1 = true;
//                }
//             }
//             int temp3=topo_mask[];
//             if (is1 && temp3!=0){
//               // topo_mask[] = (ff[]>=0.5-EPS32) ? 1 : -1 ;
// 	      if(temp3==3){
// 	          topo_mask[] = 1;
// 	      }else if(topo_mask[]==-3){
// 	          topo_mask[] = -1;
// 	      }
//             }
//         }

//      }
//   for(int phase=0;phase<=1;phase++){
//      foreach(){
//          if(level==level_interface){
//         //  if(level==level_interface && (css_test3_n[]>0.0)){
//                int phase_sign = ff[]>=0.5-EPS32 ? 1 : -1;
//             bool is1= false;
//                foreach_dimension(){
//                   int temp1=topo_mask[1];
//                   int temp2=topo_mask[-1];
//                if( temp1==(2*phase-1) || temp2==(2*phase-1)){
//                   is1 = true;
//                   }
//                }
//                int temp3 = topo_mask[];
//                if (is1 && temp3==3*phase_sign){
//                   topo_mask[] = 2*(2*phase-1) ;
//                }
//          }
//      }
//   }



// }



#include "distance.h" //for input_xy
// // Function to calculate signed distance from point (x4, y4) to a line defined by ax + by = c
// double signedDistanceToLine(double a, double b, double c, double x4, double y4) {
//     return (a*x4 + b*y4 - c) / sqrt(a*a + b*b);
// }
// // Adjusted function to use signed distance for precise position determination
// double positionRelativeToPolyline(double x4, double y4, double delta, double base, double groove_h, double groove_w) {
//     // Define original lines based on the provided equations and parameters
//     double a1 = 0.0, b1 = 1.0, c1 = base + groove_h; // Top horizontal line
//     double a2 = 1.0, b2 = 0.0, c2 = groove_w / 2.0;  // Vertical line
//     double a3 = 0.0, b3 = 1.0, c3 = base;            // Bottom horizontal line

//     // Calculate signed distances to the three lines
//     double distToLine1 = signedDistanceToLine(a1, b1, c1, x4, y4);
//     double distToLine2 = signedDistanceToLine(a2, b2, c2, x4, y4);
//     double distToLine3 = signedDistanceToLine(a3, b3, c3, x4, y4);

//     double mindist = HUGE;
//     if(fabs(distToLine1)<fabs(distToLine2)){
//         if(fabs(distToLine1)<fabs(distToLine3)){
//             mindist = distToLine1;
//         }else{
//             mindist = distToLine3;
//         }
//     }else{
//         if(fabs(distToLine2)<fabs(distToLine3)){
//             mindist = distToLine2;
//         }else{
//             mindist = distToLine3;
//         }
//     }
// }
// double solid_phi(double base, vertex scalar phi){
//     // line1 a1x+b1y=c1, line2 a2x+b2y=c2, line3 a3x+b3y=c3. normal toward normal of soid, 
//     // mod of normal = 1
//     // double base=0.0;
   
//     // double a1 = 0.0, b1 = 1.0, c1 = base + groove_h;
//     // double a2 = 1.0, b2 = 0.0, c2 = groove_w/2.0;
//     // double a3 = 0.0, b3 = 1.0, c3 = base;
//     double min_delta = L0/(1<<maxl);

//     double below_base = base-5*min_delta;
//     double above_groove = groove_h + base + 5*min_delta;
//     double scale = 100.0;
//     foreach_vertex(){
//         // if(x<below_base){
//         //       phi[] = -1.0/scale;
//         // }else if(x>above_groove){
//         //       phi[] = 1.0/scale;
//         // }else{
//         //       if((y>groove_w+5*min_delta) && (x<groove_h + base-5*min_delta)){
//         //             phi[] = -1.0/scale;
//         //       }else if((y<groove_w-5*min_delta) && (x>base+5*min_delta)){
//         //             phi[] = 1.0/scale;
//         //       }else{
//                     phi[] = positionRelativeToPolyline(x, y, min_delta, base, groove_h, groove_w)/scale;
//         //       }
//         // }
//     }

// }

// Define a new struct named PointSolid
typedef struct {
    double x, y;
} PointSolid;

double crossProduct(double x1, double y1, double x2, double y2) {
    return x1 * y2 - y1 * x2;
}

void noless_than_csfluid(double* coordxy,double originxy,double x1,double y1,double x2,double y2,bool flag_modifyx){ //css_test3
    double deltax = x2-x1, deltay = y2-y1;
        double A,B;
        // if(fabs(deltax)<1e-20){
        if(flag_modifyx){
              A = 1;
              B = 0;
        }else{
              A = 0;
              B = 1;
        }
        double cp = -crossProduct(deltax, deltay, A, B);

    double cs_min_fliud = 0.9;
    if(cp<0){
        cs_min_fliud = 1.0 -cs_min_fliud;
    }
    double delta_min = L0_pysical/((1<<maxl)*1.0);
    double tunevalue=0.0;
    double p=(*coordxy-originxy)/delta_min;
    double pp = p - round(p);
    if(pp>0.0){
        if(cp > 0){
          // if((1.0-pp)<cs_min_fliud){
              tunevalue = ((1.0-pp)-cs_min_fliud)*delta_min;
          // }
        }else{
          // if((1.0-pp)>cs_min_fliud){
              tunevalue = ((1.0-pp)-cs_min_fliud)*delta_min;
          // }
        }
    }else{
        if(cp > 0){
            // if(-pp<cs_min_fliud){
              tunevalue = (-pp-cs_min_fliud)*delta_min;
            // }
        }else{
            // if(-pp<cs_min_fliud){
              tunevalue = (-pp-cs_min_fliud)*delta_min;
            // }
        }
        
    }
    *coordxy = *coordxy + tunevalue;

}

double signedDistance(double x1, double y1, double x2, double y2, double x0, double y0) {
    double A = x0 - x1;
    double B = y0 - y1;
    double C = x2 - x1;
    double D = y2 - y1;

    double dot = A * C + B * D;
    double len_sq = C * C + D * D;
    double param = -1;
    if (len_sq != 0) {
        param = dot / len_sq;
    }

    double xx, yy;

    if (param < 0) {
        xx = x1;
        yy = y1;
    } else if (param > 1) {
        xx = x2;
        yy = y2;
    } else {
        xx = x1 + param * C;
        yy = y1 + param * D;
    }

    double dx = x0 - xx;
    double dy = y0 - yy;
    double dist = sqrt(dx * dx + dy * dy);

    double cp = crossProduct(C, D, A, B);
    if (cp < 0) {
        dist = -dist;
    }

    return dist;
}

double signedDistanceToSegmentedLine(double x0, double y0, double points[][2], int numPoints) {
    double min_d = signedDistance(points[0][0], points[0][1], points[1][0], points[1][1], x0, y0);
    double d;

    for (int i = 1; i < numPoints - 1; i++) {
        d = signedDistance(points[i][0], points[i][1], points[i + 1][0], points[i + 1][1], x0, y0);
        if (fabs(d) < fabs(min_d)) {
            min_d = d;
        }
    }

    return min_d;
}


PointSolid* addPoint(PointSolid* points, int* size, double x, double y) {
    PointSolid* temp = realloc(points, (*size + 1) * sizeof(PointSolid));
    if (!temp) {
        perror("Failed to allocate memory");
        free(points);
        exit(EXIT_FAILURE);
    }
    temp[*size].x = x;
    temp[*size].y = y;
    (*size)++;
    return temp;
}

void modifyPointAtIndex(PointSolid* points, int size, int index, double newX, double newY) {
    if (index >= 0 && index < size) {
        points[index].x = newX;
        points[index].y = newY;
    } else {
        printf("Invalid index. Index must be between 0 and %d.\n", size - 1);
    }
}

void modifylasttwopoints (PointSolid* points, int size) {
    bool flag_modifyx;
    double x1,y1,x2,y2;
    x1 = points[size-2].x;
    y1 = points[size-2].y;
    x2 = points[size-1].x;
    y2 = points[size-1].y;
    double origin_value;
    double deltax = x2-x1, deltay = y2-y1;
        double A,B;
        if(fabs(deltax)<1e-20){
              flag_modifyx = true;
        }else{
              flag_modifyx = false;
        }
    if(flag_modifyx){
        origin_value = originx;
        noless_than_csfluid(&x1,origin_value,x1,y1,x2,y2,flag_modifyx);
        x2 = x1;
    }else{
        origin_value = originy;
        noless_than_csfluid(&y1,origin_value,x1,y1,x2,y2,flag_modifyx); 
        y2 = y1;
    }
    int indexToModify = size - 2;
    modifyPointAtIndex(points, size, indexToModify, x1, y1);
    indexToModify = size - 1;
    modifyPointAtIndex(points, size, indexToModify, x2, y2);
    
}

PointSolid* generatePoints(double base, double groove_h, double groove_w_bottom, double groove_w_top, double Length, int* size) {
    *size = 0;
    PointSolid* points = NULL;
    double x = base, y = 0;

    points = addPoint(points, size, x, y);
    points = addPoint(points, size, x, y += groove_w_bottom / 2.0);
    modifylasttwopoints(points,*size); //modify x of last two points

    points = addPoint(points, size, x += groove_h, y);
    modifylasttwopoints(points,*size); //modify y of last two points

    points = addPoint(points, size, x, y += groove_w_top);
    modifylasttwopoints(points,*size); //modify x of last two points

    points = addPoint(points, size, x -= groove_h, y);
    modifylasttwopoints(points,*size); //modify y of last two points

    while (y + groove_w_bottom + groove_w_top < Length) {
        points = addPoint(points, size, x, y += groove_w_bottom);
        modifylasttwopoints(points,*size); //modify x of last two points

        points = addPoint(points, size, x += groove_h, y);
        modifylasttwopoints(points,*size); //modify y of last two points

        points = addPoint(points, size, x, y += groove_w_top);
        modifylasttwopoints(points,*size); //modify x of last two points


        points = addPoint(points, size, x -= groove_h, y);
        modifylasttwopoints(points,*size); //modify y of last two points
    }

    points = addPoint(points, size, x, y += L0);
    return points;
}

// double signedDistanceToSegmentedLine(double x0, double y0) {
// void solid_phi(double base, vertex scalar phi){
//     double x1 = base + groove_h, y1 = 2*L0;
//     double x2 = base + groove_h, y2 = groove_w/2.0;
//     double x3 = base, y3 = groove_w/2.0;
//     double x4 = base, y4 = 0;

//     foreach_vertex(){
//         double x0 = x, y0 = y;

//         double d1 = signedDistance(x1, y1, x2, y2, x0, y0);
//         double d2 = signedDistance(x2, y2, x3, y3, x0, y0);
//         double d3 = signedDistance(x3, y3, x4, y4, x0, y0);

//         double min_d = d1;
//         if (fabs(d2) < fabs(min_d)) min_d = d2;
//         if (fabs(d3) < fabs(min_d)) min_d = d3;

//         phi[] = min_d;
//     }
// }
void solid_phi(double base, vertex scalar phi){
    // double points[][2] = {
    //     {base + groove_h, 2 * L0},
    //     {base + groove_h, groove_w / 2.0},
    //     {base, groove_w / 2.0},
    //     {base, 0}
    // };
    // int numPoints = sizeof(points) / sizeof(points[0]);
    int size;
    double length_limit;
    length_limit = min(Length_heat,L0_pysical-2*groove_w_bottom);
    PointSolid* points = generatePoints(base, groove_h, groove_w_bottom, groove_w_top, length_limit, &size);
    foreach_vertex(){
        double x0 = x, y0 = y;

        // double min_d = signedDistanceToSegmentedLine(x0, y0, points, numPoints);
        double min_d = -signedDistanceToSegmentedLine(x0, y0, points, size);

        phi[] = min_d;
    }

    free(points);
}



// void smooth_for_arm(struct smooth1 q);
event init (t = 0)
{
  // CFL = 0.5;//0.01; // mesu3-14
  CFL = 0.4;//0.2;//0.02; //mesu3-15
  if(poisson_check){
      if(pid()==0){
           char name93[80];
            sprintf(name93,"poisson_check.dat");
            FILE * fp93 = fopen(name93,"w");
            fclose(fp93);
      }
  }
  if(point_trace){
      if(pid()==0){
           char name93[80];
            sprintf(name93,"point_trace.dat");
            FILE * fp93 = fopen(name93,"w");
            fclose(fp93);
      }
  }
  T_oold.nodump = true;
  printf("hello init\n");
  // if((!restartsymbol) && (!restore ("restart"))){
  if((!restartsymbol)){
     printf("hello init2\n");
      // if(restore ("filerestart/level10-dump-final8-i35763-t0.102166",list = all)){
      // if(restore ("transfer8-2/level10-dump-final8-i35763-t0.102166",list = all)){
       if(!surface_heat){
          // if(restore ("transfer8-4/level10-dump-final8-i35763-t0.102166",list = all)){
          // if(restore ("../pre-lub-temperature/dump-final-i52349-t0.0896353-maxl8-box3",list = all)){
          if(sink_number==0){
              // if(restore ("transfer8-4/level10-dump-final8-i35763-t0.102166",list = all)){
              // if(restore ("../pre-lub-temperature/dump-final-i52349-t0.0896353-maxl8-box3",list = all)){
              if(restore ("../pre-lub-temperature/level9-dump-final-i51777-t0.0892672",list = all)){
                printf("hello init 3\n");
                foreach_face(){
                  uf.x[]=0.0;
                  ulf.x[]=0.0;
                  ugf.x[]=0.0;
                }
                restart_Tsat = true;
                  // int maxl_temp=0;
                  // foreach(reduction(max:maxl_temp)){
                  //     if(level>maxl_temp){
                  //         maxl_temp = level;
                  //     }
                  // }
                  // if(maxl_temp<maxl){
                  //     refine (level<maxl);
                  //     char dumpname[80];
                  //     sprintf(dumpname,"outfacets/transfer-from%d-to%d",maxl_temp,maxl);
                  //     dump(dumpname);
                  // }

              }
            }else if(sink_number==1){
                 if(restore ("../pre-one-sink/level10-dump-final-i252970-t0.109052",list = all)){
                  printf("hello init 3\n");
                  foreach_face(){
                    uf.x[]=0.0;
                    ulf.x[]=0.0;
                    ugf.x[]=0.0;
                  }
                  restart_Tsat = true;
                }
            }else if(sink_number==2){ //dump-final-i244586-t0.105438
                // if(restore ("../pre-multi-sinks/level10-dump-final-i218856-t0.0943463",list = all)){
                 if(restore ("../pre-multi-sinks-181-smallregion1mm/level9-dump-final-i83212-t0.0158541",list = all)){
                  printf("hello init 3\n");
                  foreach_face(){
                    uf.x[]=0.0;
                    ulf.x[]=0.0;
                    ugf.x[]=0.0;
                  }
                  restart_Tsat = true;
                }
                // char name[80];
                // if(mesu_case_number==178){
                //     sprintf(name,"../pre-multi-sinks-groove40/level10-dump-final-");
                // }else if(mesu_case_number==179){
                //     sprintf(name,"../pre-multi-sinks-groove500/level10-dump-final-");
                // }
                // if(restore (name,list = all)){
                //       printf("hello init 5\n");
                //       foreach_face(){
                //         uf.x[]=0.0;
                //         ulf.x[]=0.0;
                //         ugf.x[]=0.0;
                //       }
                //       restart_Tsat = true;
                // }
            }
       }else{
          if(sink_number==1 && surface_heat_restart){
             if(restore ("../pre-one-sink/level10-dump-i231971-t0.1",list = all)){
                 printf("hello init 4\n");
                  foreach_face(){
                    uf.x[]=0.0;
                    ulf.x[]=0.0;
                    ugf.x[]=0.0;
                  }
                  // surface_heat_restart = true;
             }
          }else if(sink_number==2 && surface_heat_restart){
            //  if(restore ("../pre-one-sink/level10-dump-i231971-t0.1",list = all)){
            //      printf("hello init 4\n");
            //       foreach_face(){
            //         uf.x[]=0.0;
            //         ulf.x[]=0.0;
            //         ugf.x[]=0.0;
            //       }
            //       // surface_heat_restart = true;
            //  }
          }
       }

      if(case_number==1){
          #if TREE
          refine (sq(2.*bubble_radius) - sq(x) - sq(y) > 0 &&
            level < maxl);
            level_interface = maxl;
          #endif

    #if EMBED
          foreach(){
            cs[] =1.0;
          }
          foreach_face(){
            fs.x[] =1.0;
          }
          boundary({cs,fs.x,fs.y});
          cm_update (cm, cs, fs);
          fm_update (fm, cs, fs);
          restriction ({cm, fm, cs, fs});

            // scalar cmv = cm;
            // foreach()
            //   cmv[] = y;
            // cm[top] = dirichlet(y);
            // cm[bottom] = dirichlet(y);

            //   face vector fmv = fm;
            // foreach_face()
            //   fmv.x[] = max(y, 1e-20);
            // fm.t[top] = dirichlet(y);
            // fm.t[bottom] = dirichlet(y);

              for (scalar s in {u}){
                    s.third = true;
              }
              #if EMBED
                // foreach_dimension(){
                //   //u.x.gradient = embed_face_gradient_x ;
                //   u.x.gradient = minmod;// gradients ;
                // }
              #endif
              // u.n[embed] = dirichlet(0);
              // u.t[embed] = dirichlet(0);
    #endif


          foreach(){
            css_test3_n[] = 1.0;
          }


          fraction (ff, - (sq(bubble_radius) - sq(x-tune2_value) - sq(y)));
        // fraction (ff, x - 0.5);
          foreach(){
              css_test3_n[] = 1.0;
              ff[] = clamp (ff[], 0., 1.);
              masstr[] = 0.0;
              source_pc[] = 0.0;
              vtr[] = 0.0;
          }

              // foreach(){
              //   T[] = T_inf; //Tsat00;
              // }
          ////////////////////////////////////////////////////////
          coord * pp = input_xy (fopen (Tini_file, "r"));
            coord * ii = pp;
            scalar T_input_flag[];
            foreach(){
              T_input_flag[]=0;
            }
            int num1=0;
            double aa,bb;
            double aa_old,bb_old;
            while (ii->x != nodata) {
            // i->x += amp*noise(), i->y += amp*noise();
              printf ("%g %g\n", ii->x, ii->y);
              aa=ii->x, bb=ii->y;
            // bool check_whole = false;
                foreach(){
                    double rr = sqrt(pow(x-centerx,2.0)+pow(y-centery,2.0));//+pow(z-centerz,2.0));
                  if(num1==0){
                        if(rr<aa && T_input_flag[]==0){
                            T_input_flag[] = 1;
                            T[] = bb;
                        }
              //           check_whole = true;
                  }else{
                  if(rr<aa && T_input_flag[]==0){
                      T_input_flag[]=1;
                      T[] = (rr-aa_old)/(aa-aa_old)*bb + (aa-rr)/(aa-aa_old)*bb_old;
              //                   check_whole = true;
                        }
                  }
                }
            // if(!check_whole){
            //    break;
            // }
              aa_old=aa,bb_old=bb;
              ii++;
              num1++;
            }
            int check_flag=0;//false;
            foreach(){
              double rr = sqrt(pow(x-centerx,2.0)+pow(y-centery,2.0));//+pow(z-centerz,2.0));
              if(T_input_flag[]==0 && rr>=aa){
                T[] = bb;
                T_input_flag[] = 1;
              }
            }

            foreach(reduction(max:check_flag)){
                if(T_input_flag[]==0 && (check_flag==0)){
                      check_flag=1;
                }
            }
            if(check_flag==1){
                // fprintf(stderr,"input error\n");
            }

            foreach(){
              if(ff[]<0.5){
                    T[] = Tsat00;
              }
            }

            char name93[80];
            sprintf(name93,"check-initialTlff-pid%d.dat",pid());
            FILE * fp93 = fopen(name93,"w");
            int num2=0;
            foreach(){
              //fprintf (fp93,"%g %g %g %g %g %g\n",  x, y, z, T[], Tgff[],Tlff[]);
              fprintf (fp93,"%g %g %g\n",  x, y, T[]);
            }
            fclose(fp93);
            MPI_Barrier (MPI_COMM_WORLD);
            if(pid()==0){
                      char command1[150];
                      sprintf(command1, "LC_ALL=C cat check-initialTlff-pid*.dat > outfacets/check-initialTlff-%g",t);
                      system(command1);

                      char command7[150];
                      sprintf(command7, "LC_ALL=C rm -rf check-initialTlff-pid*.dat");
                      system(command7);
                  }
        // }else if(case_number==2 || (case_number==3) && !restore (file = "restart_Tsat0")){// tanguy, same with original case_number == 8
        }else if((case_number==2 || (case_number==3)) && (!restart_Tsat)){// tanguy, same with original case_number == 8

            if(case_number==2){
                    refine (sq(3.*bubble_radius) - sq(x-centerx) - sq(y-centery) > 0 && level < maxl);
                    level_interface = maxl;
                    fraction (ff, (- (sq(bubble_radius) - sq(x - tune2_value - centerx) - sq(y - centery))  ));
                    foreach(){
                            ff[] = clamp (ff[], 0., 1.);
                    }
                    boundary({ff});
                    const scalar c[] = theta0*pi/180.;
                    contact_angle = c;
                    foreach(){
                          ff_oppo[] = 1.0 - ff[];
                    }
            // }else if(case_number==3 && (!restore (file = "restart_Tsat0"))){
               }else if(case_number==3 && (!restart_Tsat)){
                  //  refine (sq(3.*bubble_radius) - sq(x-centerx) - sq(y-centery) > 0 && level < maxl);
                  //   level_interface = maxl;
                  //   fraction (ff, (- (sq(bubble_radius) - sq(x - centerx) - sq(y - centery))  ));
                  //   foreach(){
                  //           ff[] = clamp (ff[], 0., 1.);
                  //   }
                  //   boundary({ff});
                  //   const scalar c[] = theta0*pi/180.;
                  //   contact_angle = c;
                  //   foreach(){
                  //         ff_oppo[] = 1.0 - ff[];
                  //   }
                  if(1==0){
                     refine (x < 0.0002 && level < maxl);
                  }else{
                       for(int ii=1;ii<=10;ii++){
                          vertex scalar phi_temp3[];
                          foreach_vertex(){
                              phi_temp3[] =  -(x-(-thickbottom+tune2_value));
                          }
                          fractions (phi_temp3, poisson_source2_f);


                            vertex scalar phi_temp[];
                            solid_phi(tune2_value, phi_temp);
                            foreach_vertex(){
                              phi_temp[] = -phi_temp[];
                            }

                            fractions (phi_temp, css_test3, fss_test3);

                          scalar df2[],df3[];
                          smooth_for_arm(df2,css_test3);
                          smooth_for_arm(df3,poisson_source2_f);
                          adapt_wavelet ({df2,df3},(double[]){0.001,0.001}, maxlevel = maxl, minlevel = minl);
                        }
                  }





                    level_interface = maxl;
                    foreach(){
                            ff[] = 1.0;
                    }
                    boundary({ff});
                    const scalar c[] = theta0*pi/180.;
                    contact_angle = c;
                    foreach(){
                          ff_oppo[] = 1.0 - ff[];
                    }
            }
              ff_oppo.refine=ff_oppo.prolongation = fraction_refine;
              #if TREE
                  for (scalar c in interfaces) {
                      c.refine = c.prolongation = fraction_refine;
                      c.dirty = true;
                  }
              #endif
               vertex scalar phi_temp3[];
               foreach_vertex(){
                  phi_temp3[] =  -(x-(-thickbottom+tune2_value));
               }
              fractions (phi_temp3, poisson_source2_f);




               vertex scalar phi_temp[];
              //  foreach_vertex(){
              //     phi_temp[] =  -(x-tune2_value);//+Delta;
              //     // phi_temp[] =  -x + L0_pysical/(1<<maxl)/2.0;
              //  }
              solid_phi(tune2_value, phi_temp);
              foreach_vertex(){
                phi_temp[] = -phi_temp[];
              }

               fractions (phi_temp, css_test3, fss_test3);
              foreach(){
                    fs_solid[] = css_test3[];
                    css_test3_n[] = 1.0 - css_test3[];
               }
               foreach_face(){
                      fss_test3_n.x[] = 1.0 - fss_test3.x[];
              }
              boundary({css_test3,fss_test3});
              boundary({fs_solid});
              boundary({css_test3_n,fss_test3_n});
              foreach(){
                    cs[] = css_test3_n[];
              }
              foreach_face(){
                    fs.x[] = fss_test3_n.x[];
              }
              boundary({cs,fs});
               #if EMBED
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});

              #endif
               MPI_Barrier(MPI_COMM_WORLD);
              if(case_number==2){
                  foreach(){
                        if(x<0.0){
                            T[] = T_inf - (T_inf - Twall_init)*(x-originx)/fabs(originx);
                        }else{
                            T[] = Twall_init - (Twall_init-Tsub)*(x-0.0)/fabs(L0_pysical-fabs(originx));
                        }
                    }
              }else if(case_number==3){
                if(!surface_heat_restart){
                  foreach(){
                        T[] = Tsat0;
                    }
                }
              }
                boundary({T});
                fs_solid.prolongation = fs_solid.refine = fraction_refine;
                css_test3.refine = embed_fraction_refine_s;
                css_test3.prolongation = fraction_refine; //embed_fraction_refine_s; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3.x.prolongation = embed_face_fraction_refine_s_x;
                css_test3_n.refine = embed_fraction_refine_s_n;
                css_test3_n.prolongation = fraction_refine; //embed_fraction_refine_s_n; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3_n.x.prolongation = embed_face_fraction_refine_s_n_x;

  //new fss_test fss_test2 css_test css_test2 -- 20230821

               foreach(){
                  if(css_test3_n[]<=0){
                      css_test[]=0;
                      css_test2[]=0;
                  }else if(css_test3_n[]>=1){
                      css_test[]=ff[];
                      css_test2[]=1.0-ff[];
                  }else{
                      //?????????????????????
                      css_test[]=ff[]*css_test3_n[];
                      css_test2[]=(1.0-ff[])*css_test3_n[];
                  }
               }
      face vector fs_temp[];
       if(1==1){
               vertex scalar phi21[];
               vof2dist(css_test,phi21);
               scalar ff_temp[];
               fractions(phi21, ff_temp, fss_test);
               vof2dist(css_test2,phi21);
               fractions(phi21, ff_temp, fss_test2);
               foreach_face(){
                  fs_temp.x[] = fss_test.x[];
               }


       }else{
                vertex scalar phi21[];
                scalar ff_temp[];
                vof2dist(ff, phi21);
                fractions (phi21, ff_temp, fs_temp);

                foreach(){
                  css_test[] = 0.0;
                  css_test2[] = 0.0;
                }
                foreach_face(){
                  fss_test.x[] = 0.0;
                  fss_test2.x[]=0.0;
                }
                foreach(){
                  if(css_test3_n[]>0.0){
                    css_test[] = ff_temp[];
                    css_test2[] = 1.0 - ff_temp[];
                  }
                }
                foreach_face(){
                  if(fss_test3_n.x[]>0.0){
                      fss_test.x[] = fs_temp.x[];
                      fss_test2.x[] = 1.0 - fs_temp.x[];
                  }
                }
       }
                 css_test.refine = embed_fraction_refine_css_test;
                css_test.prolongation = embed_fraction_refine_css_test;//fraction_refine;
                foreach_dimension(){
                          fss_test.x.refine = fss_test.x.prolongation = embed_face_fraction_fss_test_refine_x;
                          //fss_test.x.dirty = false;
                }
                css_test2.refine = embed_fraction_refine_css_test2;
                css_test2.prolongation = embed_fraction_refine_css_test2;//fraction_refine;
                foreach_dimension()
                          fss_test2.x.prolongation = embed_face_fraction_fss_test2_refine_x;

                  //add 20230915
              //  get_css_fss_areaslg_triple_point();

                restriction({css_test,fss_test});
                restriction({css_test2,fss_test2});
                restriction({css_test3,fss_test3});
                restriction({css_test3_n,fss_test3_n});

                heights(ff,hhh);
              //  heights(ff_oppo,hhh_ff_oppo);
                 printf("event init 482,pid=%d\n",pid());
                      foreach(){
                        source_pc[] = 0.0;
                        masstr[]=0.0;
                        vtr[] = 0.0;
                }

        //  }else if(restore (file = "restart_Tsat0") && case_number==3){
          }else if(restart_Tsat && case_number==3){
              int maxl_temp=0;
              foreach(reduction(max:maxl_temp)){
                  if(level>maxl_temp){
                      maxl_temp = level;
                  }
              }
              if(maxl_temp<maxl){

                    scalar ff_temp2[];
                    double max_yy=0.0;
                    double temp2=L0/(1<<maxl);

                    max_yy = 5.0*bubble_radius;
                    srand(time(0));
                    foreach(){
                      
                      if(css_test3_n[]>0.0){
                        double temp = sqrt(sq(max_yy) - (sq(x) + sq(y)));
                        if(temp>0.0){
                            
                            // Generate a random floating-point number between 0 and 1
                            double randomValue = ((double)rand() / RAND_MAX-0.5);
                          ff_temp2[] = randomValue*0.001;//temp+randomValue*temp2*8 ;
                        }else{
                          ff_temp2[] = 0.0;
                        }
                      }else{
                        ff_temp2[] = 0.0;
                      }

                    }


                  for(int ii=1;ii<=10;ii++){
                    vertex scalar phi_temp3[];
                    foreach_vertex(){
                        phi_temp3[] =  -(x-(-thickbottom+tune2_value));
                    }
                    fractions (phi_temp3, poisson_source2_f);

                    scalar df1[],df2[],df3[];
                    smooth_for_arm(df1,ff);
                    smooth_for_arm(df2,css_test3_n);
                    smooth_for_arm(df3,poisson_source2_f);
                    adapt_wavelet ({df1,df2,T,ff_temp2,df3},(double[]){0.001,0.001,0.01,0}, maxlevel = maxl, minlevel = minl);
                  }
              }

              // restart_Tsat = true;
              fprintf(stderr,"init restart_Tsat && case_number==3\n");
              foreach(){
                T_oold[] = T[];
              }
                level_interface = maxl;
                    fraction (ff, (- (sq(bubble_radius) - sq(x - tune2_value - centerx) - sq(y - centery))  ));
                    foreach(){
                            ff[] = clamp (ff[], 0., 1.);
                    }
                    boundary({ff});
                    const scalar c[] = theta0*pi/180.;
                    contact_angle = c;
                    foreach(){
                          ff_oppo[] = 1.0 - ff[];
                    }
              ff_oppo.refine=ff_oppo.prolongation = fraction_refine;
              foreach(){
                  if(cs[]>=1){
                      if(ff[]<=0){
                          T[] = Tsat00;
                      }
                  }
              }
              #if TREE
                  for (scalar c in interfaces) {
                      c.refine = c.prolongation = fraction_refine;
                      c.dirty = true;
                  }
              #endif
               vertex scalar phi_temp[];
              //  foreach_vertex(){
              //   // phi_temp[] =  -x;
              //   phi_temp[] =  -(x-tune2_value);
              //     // phi_temp[] =  -x + L0_pysical/(1<<maxl)/2.0; //;//+Delta;
              //  }
              solid_phi(tune2_value, phi_temp);
              foreach_vertex(){
                phi_temp[] = -phi_temp[];
              }

              fractions (phi_temp, css_test3, fss_test3);
              foreach(){
                    fs_solid[] = css_test3[];
                    css_test3_n[] = 1.0 - css_test3[];
               }
               foreach_face(){
                      fss_test3_n.x[] = 1.0 - fss_test3.x[];
              }
              boundary({css_test3,fss_test3});
              boundary({fs_solid});
              boundary({css_test3_n,fss_test3_n});
              foreach(){
                    cs[] = css_test3_n[];
              }
              foreach_face(){
                    fs.x[] = fss_test3_n.x[];
              }
              boundary({cs,fs});
               #if EMBED
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});

              #endif
               MPI_Barrier(MPI_COMM_WORLD);
                boundary({T});
                fs_solid.prolongation = fs_solid.refine = fraction_refine;
                css_test3.refine = embed_fraction_refine_s;
                css_test3.prolongation = fraction_refine; //embed_fraction_refine_s; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3.x.prolongation = embed_face_fraction_refine_s_x;
                css_test3_n.refine = embed_fraction_refine_s_n;
                css_test3_n.prolongation = fraction_refine; //embed_fraction_refine_s_n; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3_n.x.prolongation = embed_face_fraction_refine_s_n_x;

      foreach(){
                  if(css_test3_n[]<=0){
                      css_test[]=0;
                      css_test2[]=0;
                  }else if(css_test3_n[]>=1){
                      css_test[]=ff[];
                      css_test2[]=1.0-ff[];
                  }else{
                      //?????????????????????
                      css_test[]=ff[]*css_test3_n[];
                      css_test2[]=(1.0-ff[])*css_test3_n[];
                  }
               }
      face vector fs_temp[];
       if(1==1){
               vertex scalar phi21[];
               vof2dist(css_test,phi21);
               scalar ff_temp[];
               fractions(phi21, ff_temp, fss_test);
               vof2dist(css_test2,phi21);
               fractions(phi21, ff_temp, fss_test2);
               foreach_face(){
                  fs_temp.x[] = fss_test.x[];
               }


       }else{

                vertex scalar phi21[];
                scalar ff_temp[];
                vof2dist(ff, phi21);
                fractions (phi21, ff_temp, fs_temp);

                foreach(){
                  css_test[] = 0.0;
                  css_test2[] = 0.0;
                }
                foreach_face(){
                  fss_test.x[] = 0.0;
                  fss_test2.x[]=0.0;
                }
                foreach(){
                  if(css_test3_n[]>0.0){
                    css_test[] = ff_temp[];
                    css_test2[] = 1.0 - ff_temp[];
                  }
                }
                foreach_face(){
                  if(fss_test3_n.x[]>0.0){
                      fss_test.x[] = fs_temp.x[];
                      fss_test2.x[] = 1.0 - fs_temp.x[];
                  }
                }
          }
                 css_test.refine = embed_fraction_refine_css_test;
                css_test.prolongation = embed_fraction_refine_css_test;//fraction_refine;
                foreach_dimension(){
                          fss_test.x.refine = fss_test.x.prolongation = embed_face_fraction_fss_test_refine_x;
                          //fss_test.x.dirty = false;
                }
                css_test2.refine = embed_fraction_refine_css_test2;
                css_test2.prolongation = embed_fraction_refine_css_test2;//fraction_refine;
                foreach_dimension()
                          fss_test2.x.prolongation = embed_face_fraction_fss_test2_refine_x;

                 //add 20230915
              //  get_css_fss_areaslg_triple_point();

                restriction({css_test,fss_test});
                restriction({css_test2,fss_test2});
                restriction({css_test3,fss_test3});
                restriction({css_test3_n,fss_test3_n});

                heights(ff,hhh);
              //  heights(ff_oppo,hhh_ff_oppo);
                 printf("event init 482,pid=%d\n",pid());
                foreach(){
                        source_pc[] = 0.0;
                        masstr[]=0.0;
                        vtr[] = 0.0;
                }
         }
          ///////////////////////////////////////////////////////
        //  T.refine  = refine_bilinear;
        //  T.restriction = restriction_volume_average;//restriction_average;//restriction_volume_average;
        //  T.dirty = true;
          globali = 0;
            foreach_face(){
              ulf.x[] =0.0;
              usf.x[] =0.0;
              uf.x[] = 0.0;

              ugf.x[] =0.0;
              usfg.x[] =0.0;
            }

           boundary({uf,ulf,usf});
           foreach(){
              p[] = 0.;
           }
           foreach(){
              foreach_dimension(){
                 g.x[] = 0.0;
              }
           }
             for (scalar s in {u}){
                    s.third = true;
              }

          aiml.refine = refine_aim;
          aiml.restriction = restriction_aim;
          aimg.refine = refine_aim;
          aimg.restriction = restriction_aim;
        //  int maxl_temp=0;
        //       foreach(reduction(max:maxl_temp)){
        //           if(level>maxl_temp){
        //               maxl_temp = level;
        //           }
        //       }
              // if(maxl_temp<maxl){
              //     // refine (level<maxl);
              //     refine (level<maxl && (x>-0.0002 && x<0.0002));
              //     char dumpname[80];
              //     if(maxl<13){
              //         // sprintf(dumpname,"outfacets/transfer-from%d-to%d",maxl_temp,maxl);
              //     }
              //     // dump(dumpname);
              // }

    }else if(restartsymbol && case_number==3){
        if(restore("filerestart/restart-n",list=all)){
                    level_interface = maxl;
                    const scalar c[] = theta0*pi/180.;
                    contact_angle = c;
                    foreach(){
                          ff_oppo[] = 1.0 - ff[];
                    }
              ff_oppo.refine=ff_oppo.prolongation = fraction_refine;
              #if TREE
                  for (scalar c in interfaces) {
                      c.refine = c.prolongation = fraction_refine;
                      c.dirty = true;
                  }
              #endif
               vertex scalar phi_temp[];
              //  foreach_vertex(){
              //     // phi_temp[] =  -x;//+Delta;
              //     phi_temp[] =  -(x-tune2_value);
              //     // phi_temp[] =  -x + L0_pysical/(1<<maxl)/2.0;
              //  }
               solid_phi(tune2_value, phi_temp);
              foreach_vertex(){
                phi_temp[] = -phi_temp[];
              }


               fractions (phi_temp, css_test3, fss_test3);
              foreach(){
                    fs_solid[] = css_test3[];
                    css_test3_n[] = 1.0 - css_test3[];
               }
               foreach_face(){
                      fss_test3_n.x[] = 1.0 - fss_test3.x[];
              }
              boundary({css_test3,fss_test3});
              boundary({fs_solid});
              boundary({css_test3_n,fss_test3_n});
              foreach(){
                    cs[] = css_test3_n[];
              }
              foreach_face(){
                    fs.x[] = fss_test3_n.x[];
              }
              boundary({cs,fs});
               #if EMBED
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});

              #endif
               MPI_Barrier(MPI_COMM_WORLD);
                fs_solid.prolongation = fs_solid.refine = fraction_refine;
                css_test3.refine = embed_fraction_refine_s;
                css_test3.prolongation = fraction_refine; //embed_fraction_refine_s; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3.x.prolongation = embed_face_fraction_refine_s_x;
                css_test3_n.refine = embed_fraction_refine_s_n;
                css_test3_n.prolongation = fraction_refine; //embed_fraction_refine_s_n; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3_n.x.prolongation = embed_face_fraction_refine_s_n_x;

       foreach(){
                  if(css_test3_n[]<=0){
                      css_test[]=0;
                      css_test2[]=0;
                  }else if(css_test3_n[]>=1){
                      css_test[]=ff[];
                      css_test2[]=1.0-ff[];
                  }else{
                      //?????????????????????
                      css_test[]=ff[]*css_test3_n[];
                      css_test2[]=(1.0-ff[])*css_test3_n[];
                  }
               }
      face vector fs_temp[];
       if(1==1){
               vertex scalar phi21[];
               vof2dist(css_test,phi21);
               scalar ff_temp[];
               fractions(phi21, ff_temp, fss_test);
               vof2dist(css_test2,phi21);
               fractions(phi21, ff_temp, fss_test2);
               foreach_face(){
                  fs_temp.x[] = fss_test.x[];
               }


       }else{
                vertex scalar phi21[];
                scalar ff_temp[];
                vof2dist(ff, phi21);
                fractions (phi21, ff_temp, fs_temp);

                foreach(){
                  css_test[] = 0.0;
                  css_test2[] = 0.0;
                }
                foreach_face(){
                  fss_test.x[] = 0.0;
                  fss_test2.x[]=0.0;
                }
                foreach(){
                  if(css_test3_n[]>0.0){
                    css_test[] = ff_temp[];
                    css_test2[] = 1.0 - ff_temp[];
                  }
                }
                foreach_face(){
                  if(fss_test3_n.x[]>0.0){
                      fss_test.x[] = fs_temp.x[];
                      fss_test2.x[] = 1.0 - fs_temp.x[];
                  }
                }
        }
                 css_test.refine = embed_fraction_refine_css_test;
                css_test.prolongation = embed_fraction_refine_css_test;//fraction_refine;
                foreach_dimension(){
                          fss_test.x.refine = fss_test.x.prolongation = embed_face_fraction_fss_test_refine_x;
                          //fss_test.x.dirty = false;
                }
                css_test2.refine = embed_fraction_refine_css_test2;
                css_test2.prolongation = embed_fraction_refine_css_test2;//fraction_refine;
                foreach_dimension()
                          fss_test2.x.prolongation = embed_face_fraction_fss_test2_refine_x;

                //add 20230915
              //  get_css_fss_areaslg_triple_point();

                restriction({css_test,fss_test});
                restriction({css_test2,fss_test2});
                restriction({css_test3,fss_test3});
                restriction({css_test3_n,fss_test3_n});

                heights(ff,hhh);
        }else{
             printf("restartsymbol error\n");
        }
    }else{
       printf("restart - istep=%d:pid=%d\n",i,pid());
    }


    ///initial cm_css_test
    metric_css_test();
    metric_css_test2();
    metric_css_test3();
  }

// event diffusion_T (i++){
//   face vector D[];
//   scalar rhocp[];
//   dt = dtnext(DT);
//   fprintf(stderr,"dt=%g\n",dt);
//   foreach(){
//     rhocp[] = Trhol*Tcpl*cm[];
//   }
//   foreach_face(){
//     D.x[] = Tkl*fm.x[];
//   }
//   diffusion1(T, dt, D, theta=rhocp);

// }

// event stability(i++,last){
//  // if(!surface_heat){
//       double vtr_max=0.0;
//       double CFL_evap = 0.1;
//        double Delta_min = HUGE;
//       foreach(reduction(max:vtr_max) reduction(min:Delta_min)){
//            vtr_max = max(vtr_max,fabs(vtr[]));
//            if(Delta_min>Delta){
//                Delta_min = Delta;
//            }
//       }
//       //double Delta_min = Delta;
//       if(fabs(vtr_max)>EPS){
//          double dtmax_evap = CFL_evap*Delta_min/vtr_max;
//         if(dtmax>dtmax_evap){
//              dtmax = dtmax_evap;
//            //  fprintf(stderr,"dtmax_evap>dtmax dtmax_evap=%g\n",dtmax_evap);
//         }
//       }
//   //}
// }

event before_vof (i++,last){

   get_topomask(topo_mask);// add 20230103
        // get aiml,aimg. for Tgrad_leon
            // foreach(){
            //     Tl[] = Tsat00;
            //     Tg[] = Tsat00;
            //     Ts[] = Tsat00;
            //     if(css_test3[]>=0.5){
            //           Ts[] = T[];
            //     }else{
            //         if(css_test[]>=0.5){
            //             Tl[] = T[];
            //         }
            //         if(css_test2[]>=0.5){
            //             Tg[] = T[]; 
            //         }
            //     }
            // }
            // boundary({Ts,Tl,Tg});
            // embed_flux_conjugate(Ts, Tl, Tg, aiml, aimg, T_modl,T_modg, flux_l, flux_g);
   mass_transfer_rate();
  //  printf("after  mass_transfer \n");
 }

event vof (i++,last){

    foreach(){
    ff_old[] = ff[];
    ff_old2[] = ff[];
  }
  boundary({ff_old});




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////get css_test
 fprintf(stderr,"sweep2 vof 1378\n");
                vertex scalar phi[];
                // fraction(css_test3_n,x-tune2_value);
                // foreach(){
                //   css_test3_n[] =  clamp (css_test3_n[], 0., 1.);
                // }
                // vof2dist(css_test3_n,phi);
                 solid_phi(tune2_value, phi);
                boundary ({phi});
                fractions (phi, css_test3_n, fss_test3_n);
                              foreach(){
                                  css_test3[] = 1.0 - css_test3_n[];
                                  fs_solid[] = css_test3[];
                                }
                                foreach_face(){
                                    fss_test3.x[] = 1.0 - fss_test3_n.x[]; 
                                }

                
                  foreach(){
                    cs[] = css_test3_n[];
                  }
                  foreach_face(){
                    fs.x[] = fss_test3_n.x[]; 
                  }
                  boundary({cs,fs});

   fprintf(stderr,"sweep2 vof 1404\n");
               #if EMBED     
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});
              #endif
    fprintf(stderr,"sweep2 vof 1411\n");
                
                fs_solid.prolongation = fs_solid.refine = fraction_refine;
                css_test3.refine = embed_fraction_refine_s;
                css_test3.prolongation = fraction_refine; //embed_fraction_refine_s; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3.x.prolongation = embed_face_fraction_refine_s_x;
                css_test3_n.refine = embed_fraction_refine_s_n;
                css_test3_n.prolongation = fraction_refine; //embed_fraction_refine_s_n; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3_n.x.prolongation = embed_face_fraction_refine_s_n_x;
       foreach(){
                  if(css_test3_n[]<=0){
                      css_test[]=0;
                      css_test2[]=0;
                  }else if(css_test3_n[]>=1){
                      css_test[]=ff[];
                      css_test2[]=1.0-ff[];
                  }else{
                      //?????????????????????
                      css_test[]=ff[]*css_test3_n[];
                      css_test2[]=(1.0-ff[])*css_test3_n[];
                  }
               } 
      
      face vector fs_temp[];
       if(1==1){        
               vertex scalar phi21[];
               vof2dist(css_test,phi21);
               scalar ff_temp[];
               fractions(phi21, ff_temp, fss_test);
               vof2dist(css_test2,phi21);
               fractions(phi21, ff_temp, fss_test2);
               foreach_face(){
                  fs_temp.x[] = fss_test.x[];
               }

              
       }
                 css_test.refine = embed_fraction_refine_css_test;
                css_test.prolongation = embed_fraction_refine_css_test;//fraction_refine;
                foreach_dimension(){
                          fss_test.x.refine = fss_test.x.prolongation = embed_face_fraction_fss_test_refine_x;
                          //fss_test.x.dirty = false;
                }
                css_test2.refine = embed_fraction_refine_css_test2;
                css_test2.prolongation = embed_fraction_refine_css_test2;//fraction_refine;
                foreach_dimension()
                          fss_test2.x.prolongation = embed_face_fraction_fss_test2_refine_x;

                 //add 20230915
   fprintf(stderr,"sweep2 vof 1462\n");
               get_css_fss_areaslg_triple_point(); 
    fprintf(stderr,"sweep2 vof 1464\n");
////////////////////////////////////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////


if(dump_each_event && (i%dump_each_event_interval==0)){
        char dumpname[80];
    // topo_mask_s.nodump = true;
      foreach_dimension(){
        hhh.x.nodump = true;
        //  hhh_ff_oppo.x.nodump = true;
      }
      sprintf(dumpname,"outfacets/dumpbeforevof-i%d-t%g",i,t);
      dump(dumpname);


}
  if(point_trace){
      double dd=-1;
      foreach(reduction(max:dd)){
        if(fabs(x-tracex)<1e-6 && fabs(y-tracey)<1e-6){
           dd = T[];
        }
      }
      if(pid()==0){
           printf("you yong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
           char name93[80];
            sprintf(name93,"point_trace.dat");
            FILE * fp93 = fopen(name93,"a");
            fprintf(fp93,"%.8f %d %g\n",t,1,dd); //1 represent before vof
            fclose(fp93);
      }
  }

    foreach_face(){
        ugf.x[] = 0.0;
    }

  //  foreach(){
  //     psg[]=-1;
  //  }
  //  boundary(flux_out);
///////////// test
    // foreach_face(){
    //     uf.x[] = 0.01*fm.x[];
    //     ulf.x[] = 0.01*fm.x[];
    //     ugf.x[] = 0.01*fm.x[];
    // }

//////////////
// boundary({T,Tlff,Tgff});

     get_topomask(topo_mask);
     event ("properties");
  //    foreach(){
  //    // Tlff[] = Tsat00*(ff[])*Trhol*Tcpl;
  //    // Tgff[] = Tsat00*(1.0-ff[])*Trhog*Tcpg;
  //    // T_solid[] = Tsat00*Trhos*Tcps;
  //     if(css_test3_n[]>0.0){

  //       // Tlff[] = T[]*ff[]*Trhol*Tcpl;
  //       // Tgff[] = T[]*(1.0-ff[])*Trhog*Tcpg;

  //       if(ff[]>=0.5){
  //           Tlff[] = T[]*ff[]*Trhol*Tcpl;
  //           Tgff[] = Tsat00*(1.0-ff[])*Trhog*Tcpg;
  //       }else{
  //           Tlff[] = Tsat00*ff[]*Trhol*Tcpl;
  //           Tgff[] = T[]*(1.0-ff[])*Trhog*Tcpg;
  //       }
  //     }else{ //inside solid, temperature???(extrapolate?)
  //       Tlff[] = T[]*ff[]*Trhol*Tcpl;
  //       Tgff[] = T[]*(1.0-ff[])*Trhog*Tcpg;
  //     }

  //     if(css_test3_n[]<1.0){ //css_test3[]>0.0
  //         T_solid[]=T[]*Trhos*Tcps;
  //     }
  //  }
   if(energy_advecting_flag){
     foreach(){
      Tlff[] = Tsat00*(ff[])*Trhol*Tcpl;
      Tgff[] = Tsat00*(1.0-ff[])*Trhog*Tcpg;
      T_solid[] = Tsat00*Trhos*Tcps;
      if(css_test3_n[]>0.0){
        Tlff[] = T[]*ff[]*Trhol*Tcpl;
        Tgff[] = T[]*(1.0-ff[])*Trhog*Tcpg;
      }else{ //inside solid, temperature???(extrapolate?)
        Tlff[] = T[]*ff[]*Trhol*Tcpl;
        Tgff[] = T[]*(1.0-ff[])*Trhog*Tcpg;
      }

      if(css_test3_n[]<1.0){ //css_test3[]>0.0
          T_solid[]=T[]*Trhos*Tcps;
      }
   }
}else{ //ff*T advecting
    if((!use_Tslg) || (i<=83212)){
        foreach(){
          Tlff[] = Tsat00*(ff[]);
          Tgff[] = Tsat00*(1.0-ff[]);
          T_solid[] = Tsat00;
          if(css_test3_n[]>0.5){
            if(ff[]>=0.5){
                Tlff[] = T[]*ff[];
            }
            if(ff[]<0.5){
                Tgff[] = T[]*(1.0-ff[]);
            }
          }else{ //inside solid, temperature???(extrapolate?)
            Tlff[] = T[]*ff[];
            Tgff[] = T[]*(1.0-ff[]);
          }

          if(css_test3_n[]<=0.5){ //css_test3[]>0.0
              T_solid[]=T[];
          }else{
              T_solid[]=T[];
          }
      }
    }else{
         foreach(){
            Tlff[] = Tl[]*(ff[]);
            Tgff[] = Tg[]*(1.0-ff[]);
            T_solid[] = Ts[];
        }
    }
 }

   Tlff.restriction = restriction_Tlff;
   Tgff.restriction = restriction_Tgff;
  Tlff.refine = Tlff.prolongation = conservative_refine;
  Tgff.refine = Tgff.prolongation = conservative_refine;
  //  Tlff.refine = bilinear;
  //  Tgff.refine = bilinear;
   boundary({Tlff,Tgff});



  //  get_u_ghost(u);

}




event after_vof (i++){




////////////////////////////////////////
    // foreach_face(){
    //     uf.x[] = 0.01*0;
    //     ulf.x[] = 0.01*0;
    //     ugf.x[] = 0.01*0;
    // }
///////////////////////////////////////
if(energy_advecting_flag){
    foreach(){
      double val_tot=0.0;
      double wei=0.0;
      double val1=0.0;
      if(css_test3_n[]>0.0){
         if(ff[]>0.0){
            val_tot += Tlff[];
            wei += ff[]*Trhol*Tcpl;
         }
         if(1.0-ff[]>0.0){
            val_tot += Tgff[];
            wei += (1.0-ff[])*Trhog*Tcpg;
         }
          val1 = val_tot/wei;
      }
      //  if(topo_mask[]==0){//|| (val1<Tsat00 && (topo_mask[]==-1 || topo_mask[]==-2))){
      //      val1 = Tsat00;
      //  }
      if(css_test3_n[]<1.0){ //css_test3[]>0.0
          T[] = val1*css_test3_n[] + T_solid[]/(Trhos*Tcps)*(1.0-css_test3_n[]);
      }else{
          T[] = val1;
      }
   }
   }else{
      foreach(){
          double val_tot=0.0;
          double wei=0.0;
          double val1=0.0;
          if(css_test3_n[]<=0.5){
              T[] = T_solid[];
          }else{
              if(ff[]>=0.5){
                  T[] = Tlff[]/ff[];
              }else{
                  T[] = Tgff[]/(1.0-ff[]);
              }
          }
      }
      if(use_Tslg){
          foreach(){
            if(css_test3_n[]>0.0){
                if(ff[]>0.0){
                    Tl[] = Tlff[]/ff[];
                    if(Tl[]<Tsat00){
                        Tl[]=Tsat00;
                    }
                }else{
                  Tl[]=Tsat00;
                }
                if(ff[]<1.0){
                    Tg[] = Tgff[]/(1.0-ff[]);
                    if(Tg[]<Tsat00){
                        Tg[]=Tsat00;
                    }
                }else{
                  Tg[]=Tsat00;
                }
            }else{
                Tl[]=Tsat00;
                Tg[]=Tsat00;
            }
          }
      }
}
 across_interface(ff_old,ff,T);
if(dump_each_event && (i%dump_each_event_interval==0)){
            char dumpname[80];
      // topo_mask_s.nodump = true;
        foreach_dimension(){
          hhh.x.nodump = true;
          //  hhh_ff_oppo.x.nodump = true;
        }
        // foreach(){
        //   flux_x[] = flux_show.x[];
        //   flux_y[] = flux_show.y[];
        // }
        sprintf(dumpname,"outfacets/dumpaftervof-i%d-t%g",i,t);

        dump(dumpname);
}

  if(point_trace){
     double dd=-1;
      foreach(reduction(max:dd)){
        if(fabs(x-tracex)<1e-6 && fabs(y-tracey)<1e-6){
           dd = T[];
        }
      }
      if(pid()==0){
           char name93[80];
            sprintf(name93,"point_trace.dat");
            FILE * fp93 = fopen(name93,"a");
            fprintf(fp93,"%.8f %d %g\n",t,2,dd); //2 represent after vof
            fclose(fp93);
      }
  }

  foreach(){
    ff_old[] = ff[];
  }
  boundary({ff_old});

      boundary({ff});
   // boundary({T});
      event ("properties");
      get_topomask(topo_mask);  //add this 26-3
      LevelSetShift2VOFChange(dt);
      // flag_topos_advect_uf=true;

      double over_half=1e-2;
      double base=0.6;//0.5;
      mov_interface_dc(flag_topos_advect_uf,flag_cant_smaller_than_half,base,over_half);

//////////////////////////////////////////////////////////
///////////////////  topo_mask_s[]==0 ff[]>0.5
       foreach(){
          ff[] = clamp (ff[], 0., 1.);
      }
      boundary({ff});
      get_topomask_soilid(topo_mask_s,css_test3,level_interface);
      if(flag_cant_smaller_than_half){ //too much cell <0.5 , pressure oscilation
          foreach(){
              bool flag=false;
              // foreach_neighbor(){
              foreach_neighbor(2){
                  if(intersect_true[]==1){
                      flag=true;
                  }
              }
              if(topo_mask_s[]==0 && (!flag)){
                  if(ff[]<base){
                      if(ff_old2[]>=base+over_half){
                        ff[] = ff_old2[];
                      }else{
                        if(ff_old2[]<=1e-6){
                            ff[] = 0.0;
                        }else{ //ff_old2 between EPS and 0.5
                            ff[] = base+over_half;
                        }
                      }
                  }
              }
          }
      }


      foreach(){
          ff[] = clamp (ff[], 0., 1.);
      }
      boundary({ff});
       across_interface(ff_old,ff,T);
    event ("properties");
    get_topomask(topo_mask);

if(dump_each_event && (i%dump_each_event_interval==0)){
          char dumpname2[80];
      // topo_mask_s.nodump = true;
        foreach_dimension(){
          hhh.x.nodump = true;
          //  hhh_ff_oppo.x.nodump = true;
        }
        sprintf(dumpname2,"outfacets/dumpafterlevel-i%d-t%g",i,t);
        dump(dumpname2);
}
  if(point_trace){
      double dd=-1;
      foreach(reduction(max:dd)){
        if(fabs(x-tracex)<1e-6 && fabs(y-tracey)<1e-6){
           dd = T[];
        }
      }
      if(pid()==0){
           char name93[80];
            sprintf(name93,"point_trace.dat");
            FILE * fp93 = fopen(name93,"a");
            fprintf(fp93,"%.8f %d %g\n",t,3,dd); //3 represent after levelset
            fclose(fp93);
      }
  }
  // if(i%3==0){
    scalar ff_remove[];
    foreach(){
      ff_remove[]=0.0;
      if(cs[]>0.0){
        ff_remove[] = ff[];
      }
    }
    remove_droplets (ff_remove,3,1e-10, false);
    foreach(){
       if(cs[]>0.0){
          ff[] = ff_remove[];
       }
    }
  // }

    // remove_droplets (ff,3,1e-10, false);
}

event diffusionT_one (i++){
  ////////////////////////////////////////copy from event init
              // vertex scalar phi_temp[];
              //  foreach_vertex(){
              //     phi_temp[] = -x;
              //  }
              //  fractions (phi_temp, css_test3, fss_test3);
              // foreach(){
              //       fs_solid[] = css_test3[];
              //       css_test3_n[] = 1.0 - css_test3[];
              //  }
              //  foreach_face(){
              //         fss_test3_n.x[] = 1.0 - fss_test3.x[];
              // }
              // boundary({css_test3,fss_test3});
              // boundary({fs_solid});
              // boundary({css_test3_n,fss_test3_n});
              // foreach(){
              //       cs[] = css_test3_n[];
              // }
              // foreach_face(){
              //       fs.x[] = fss_test3_n.x[];
              // }
              // boundary({cs,fs});
              //  #if EMBED
              //     restriction ({cs, fs});
              //     cm_update (cm, cs, fs);
              //     fm_update (fm, cs, fs);
              //     restriction ({cm, fm});
              // #endif

   foreach(){
      cs[] = 1.0; //css_test3_n[];
    }
    foreach_face(){
      fs.x[] = 1.0;//fss_test3_n.x[];
    }
    boundary({cs,fs});
               #if EMBED
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});
              #endif
  fs_solid.prolongation = fs_solid.refine = fraction_refine;
 // fs_solid.dirty = true; // boundary conditions need to be updated
   vertex scalar phi30[];
   boundary({fs_solid});
   vof2dist(fs_solid, phi30);
   fractions (phi30, css_test3, fss_test3);
   boundary({css_test3,fss_test3});
   scalar tempfs[];
   foreach(){
       tempfs[] = fs_solid[];
       fs_solid[] = 1.0 - tempfs[];
   }
   vertex scalar phi31[];
   boundary({fs_solid});
   vof2dist(fs_solid, phi31);
   fractions (phi31, css_test3_n, fss_test3_n);
   boundary({css_test3_n,fss_test3_n});
   foreach(){
       fs_solid[] = tempfs[];
   }
             //   fs_solid.prolongation = fs_solid.refine = fraction_refine;
                css_test3.refine = embed_fraction_refine_s;
                css_test3.prolongation = fraction_refine; //embed_fraction_refine_s; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3.x.prolongation = embed_face_fraction_refine_s_x;
                css_test3_n.refine = embed_fraction_refine_s_n;
                css_test3_n.prolongation = fraction_refine; //embed_fraction_refine_s_n; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3_n.x.prolongation = embed_face_fraction_refine_s_n_x;

      foreach(){
                  if(css_test3_n[]<=0){
                      css_test[]=0;
                      css_test2[]=0;
                  }else if(css_test3_n[]>=1){
                      css_test[]=ff[];
                      css_test2[]=1.0-ff[];
                  }else{
                      //?????????????????????
                      css_test[]=ff[]*css_test3_n[];
                      css_test2[]=(1.0-ff[])*css_test3_n[];
                  }
               }
      face vector fs_temp[];
       if(1==1){
               vertex scalar phi21[];
               vof2dist(css_test,phi21);
               scalar ff_temp[];
               fractions(phi21, ff_temp, fss_test);
               vof2dist(css_test2,phi21);
               fractions(phi21, ff_temp, fss_test2);
               foreach_face(){
                  fs_temp.x[] = fss_test.x[];
               }


       }else{
                vertex scalar phi21[];
                scalar ff_temp[];
                vof2dist(ff, phi21);
                fractions (phi21, ff_temp, fs_temp);

                foreach(){
                  css_test[] = 0.0;
                  css_test2[] = 0.0;
                }
                foreach_face(){
                  fss_test.x[] = 0.0;
                  fss_test2.x[]=0.0;
                }
                foreach(){
                  if(css_test3_n[]>0.0){
                    css_test[] = ff_temp[];
                    css_test2[] = 1.0 - ff_temp[];
                  }
                }
                foreach_face(){
                  if(fss_test3_n.x[]>0.0){
                      fss_test.x[] = fs_temp.x[];
                      fss_test2.x[] = 1.0 - fs_temp.x[];
                  }
                }

          }
                 css_test.refine = embed_fraction_refine_css_test;
                css_test.prolongation = embed_fraction_refine_css_test;//fraction_refine;
                foreach_dimension(){
                          fss_test.x.refine = fss_test.x.prolongation = embed_face_fraction_fss_test_refine_x;
                          //fss_test.x.dirty = false;
                }
                css_test2.refine = embed_fraction_refine_css_test2;
                css_test2.prolongation = embed_fraction_refine_css_test2;//fraction_refine;
                foreach_dimension()
                          fss_test2.x.prolongation = embed_face_fraction_fss_test2_refine_x;

                 //add 20230915
              //  get_css_fss_areaslg_triple_point();

                restriction({css_test,fss_test});
                restriction({css_test2,fss_test2});
                restriction({css_test3,fss_test3});
                restriction({css_test3_n,fss_test3_n});

                heights(ff,hhh);
              //  heights(ff_oppo,hhh_ff_oppo);


  ////////////////////////////////////////

  // face vector D[];
  // scalar rhocp[];
  //  dt = dtnext(DT);
  // printf("dt=%g\n",dt);
  // fprintf(stderr,"dt=%g\n",dt);
  // fprintf(stderr,"T_inf=%g\n",T_inf);
  // foreach(){
  //   rhocp[] = Trhol*Tcpl * cm[];
  // }
  // foreach_face(){
  //   D.x[] = Tkl * fm.x[];
  // }
  // diffusion1(T, dt, D, theta=rhocp);
  // boundary({T});
  // heights(ff,hhh);
  get_modphase01_2(ff,modphase1,modphase0);
  get_css_fss_areaslg_triple_point();

  #if EMBED
    restriction ({cm_css_test, fss_test});
    cm_css_test_update (cm_css_test, css_test, fss_test);
    fm_fss_test_update (fm_fss_test, css_test, fss_test);
    restriction ({cm_css_test, fss_test});
  #endif

  #if EMBED
    restriction ({cm_css_test2, fss_test2});
    cm_css_test2_update (cm_css_test2, css_test2, fss_test2);
    fm_fss_test2_update (fm_fss_test2, css_test2, fss_test2);
    restriction ({cm_css_test2, fss_test2});
  #endif

  #if EMBED
    restriction ({cm_css_test3, fss_test3});
    cm_css_test3_update (cm_css_test3, css_test3, fss_test3);
    fm_fss_test3_update (fm_fss_test3, css_test3, fss_test3);
    restriction ({cm_css_test3, fss_test3});
  #endif

  foreach(){
    foreach_dimension(){
      smallmodl.x[]=-10.0;
      bigmodl.x[]=-10.0;
      smallmodg.x[]=-10.0;
      bigmodg.x[]=-10.0;

    }
}

//check modphase

  face vector D[];
  //scalar rhocp[];
  scalar tag_phase[];
  scalar T_old[];

  // foreach(){
  //   T_old[] = T[];
  //   ff_old[] = ff[];
  //   if(ff[]>=0.5){
  //       tag_phase[] = 1;
  //   }else{
  //       tag_phase[] = 0;
  //   }
  // }
  get_topomask(topo_mask);
  restriction({ff});

 event ("properties") ;
#if !EMBED
  // for(phase_flag3=0;phase_flag3<=1;phase_flag3++){
  //     foreach(){
  //             if(tag_phase[] == phase_flag3){
  //                 T[] = T_old[];
  //             }else{
  //                         //set boundary for , topo_mask can be used
  //                 T[] = Tsat0;
  //             }
  //     }
  //     if(phase_flag3==0){

  //             foreach_face(){
  //             // #if EMBED
  //             //      double temp1,temp2;
  //             //      temp1 = fs.x[];
  //             //       if(fs.x[]<EPS_FS){
  //             //          temp1 = EPS_FS;
  //             //       }
  //             //       temp2 = max(y,1e-20);
  //             //      D.x[] = Tkg*temp1*temp2;//Tkg*fss_test.x[];
  //             // #else
  //                   D.x[] = Tkg*fm.x[];//Tkg*fss_test.x[];
  //             // #endif
  //             }
  //             foreach(){
  //                  if(css_test3_n[]>=1){ // add this on 20230614
  //                        rhocp[] = Trhog*Tcpg*cm[];
  //                  }
  //             }

  //     }else{
  //             foreach_face(){
  //               // #if EMBED
  //               //    double temp1,temp2;
  //               //    temp1 = fs.x[];
  //               //     if(fs.x[]<EPS_FS){
  //               //        temp1 = EPS_FS;
  //               //     }
  //               //     temp2 = max(y,1e-20);
  //               //    D.x[] = Tkl*temp1*temp2;//Tkg*fss_test.x[];
  //               //  #else
  //                   D.x[] = Tkl*fm.x[];//Tkg*fss_test.x[];
  //               //  #endif
  //             }
  //             foreach(){
  //                 if(css_test3_n[]>=1){ // add this on 20230614
  //                        rhocp[] = Trhol*Tcpl*cm[];
  //                 }

  //              }
  //     }
  //     event ("properties");
  //     // diffusion1(T, dt, D, theta=rhocp); //20230918
  //    // event ("properties");
  //     foreach(){
  //         if(tag_phase[] == phase_flag3){
  //             T_old[] = T[];
  //         }
  //         boundary({T});
  //     }
  // }
  // foreach(){
  //     T[] = T_old[];
  // }
  // foreach(){
  //         //if(css_test3_n[]>=1.0){
  //            if(css_test3_n[]>=1.0){
  //                 if(fabs(ff[]-0.5)<1e-6){
  //                     T[] = Tsat00;
  //                  }
  //         }
  // }
#else
  //  //  face vector D[];
  //         foreach_face(){
  //            if(fss_test3_n.x[]>=1){
  //                  // D.x[] = fss_test.x[]*Tkl + Tkg*(1. - fss_test.x[]);
  //                   D.x[] = (fs_temp.x[]*Tkl + Tkg*(1. - fs_temp.x[]))*fm.x[]; //*max(y,1e-20);//*fm.x[]; //max(y,1e-20);
  //            }else if(fss_test3_n.x[]<=0){
  //                   D.x[] = Tks*fm.x[]; //*max(y,1e-20);//*fm.x[];;
  //                   //printf("D.x = %g\n",D.x[]);
  //            }else{
  //                   double val_fluid,val_solid;
  //                  // val_fluid = fss_test.x[]*Tkl + Tkg*(1. - fss_test.x[]);
  //                   val_fluid = fs_temp.x[]*Tkl + Tkg*(1. - fs_temp.x[]);
  //                   val_solid = Tks;
  //                   D.x[] = (fss_test3_n.x[]*val_fluid + (1.0-fss_test3_n.x[])*val_solid)*fm.x[]; //*max(y,1e-20);//*fm.x[];//max(y,1e-20);
  //            }
  //         }
    for(phase_flag3=0;phase_flag3<=1;phase_flag3++){ //for fliud and gas
          //     foreach(){
          //       //if(css_test3_n[]>0.0){
          //       if(css_test3_n[]>=1.0){
          //             if(tag_phase[] == phase_flag3){
          //                 T[] = T_old[];
          //             }else{
          //                 //set boundary for , topo_mask can be used
          //                 T[] = Tsat00;
          //             }
          //       }else{  //T inside solid
          //             T[] = T_old[];
          //       }
          //     }
          //     boundary({T});
          // //    T.flag_Tsat0_0 = false;
          // //    T.phase_flag33 = phase_flag3;
          // //    T.ff2 = ff;
          //     //T.restriction = restriction_cprho_average23_like2;//for liquid ff, for gas 1-ff
          //    //event ("properties");
          // event ("properties");
          //         if(poisson_check){
          //         if(pid()==0){
          //             char name93[80];
          //               sprintf(name93,"poisson_check.dat");
          //               FILE * fp93 = fopen(name93,"a");
          //               fprintf(fp93,"\n\n");
          //               fprintf(fp93,"i=%d maxres(phase_flag3=%d):",i,phase_flag3);
          //               fclose(fp93);
          //         }
          //     }
          //     foreach(){
          //       if(css_test3_n[]>=1){
          //           if(phase_flag3==0){
          //             rhocp[] = Trhog*Tcpg*cm[];
          //           }else{
          //             rhocp[] = Trhol*Tcpl*cm[];
          //           }
          //       }
          //     }
          // foreach_face(){
          //    if(fss_test3_n.x[]>=1){
          //          // D.x[] = fss_test.x[]*Tkl + Tkg*(1. - fss_test.x[]);
          //           // D.x[] = (fs_temp.x[]*Tkl + Tkg*(1. - fs_temp.x[]))*fm.x[]; //*max(y,1e-20);//*fm.x[]; //max(y,1e-20);
          //           if(phase_flag3==1){
          //               D.x[] = Tkl*fm.x[];
          //           }else{
          //               D.x[] = Tkg*fm.x[];
          //           }
          //    }else if(fss_test3_n.x[]<=0){
          //           D.x[] = Tks*fm.x[]; //*max(y,1e-20);//*fm.x[];;
          //           //printf("D.x = %g\n",D.x[]);
          //    }else{
          //           double val_fluid,val_solid;
          //          // val_fluid = fss_test.x[]*Tkl + Tkg*(1. - fss_test.x[]);
          //           val_fluid = fs_temp.x[]*Tkl + Tkg*(1. - fs_temp.x[]);
          //           val_solid = Tks;
          //           D.x[] = (fss_test3_n.x[]*val_fluid + (1.0-fss_test3_n.x[])*val_solid)*fm.x[]; //*max(y,1e-20);//*fm.x[];//max(y,1e-20);
          //    }
          // }
                // diffusion1(T,dt,D,theta=rhocp); //20230918
          // event ("properties");
          //     foreach(){
          //         if(css_test3_n[]>=1.0){
          //             if(tag_phase[] == phase_flag3){
          //                 T_old[] = T[];
          //             }
          //         }else{ // css_test3_n<1.0
          //               T_old[] = T[];  //for solid res is divided into three part, one in gas, one in liquid and the other in solid
          //         }
          //     }
          //     boundary({T});
          }
          // foreach(){
          //     T[] = T_old[];
          // }
          // foreach_face(){
          //    if(fss_test3_n.x[]>=1){
          //          // D.x[] = fss_test.x[]*Tkl + Tkg*(1. - fss_test.x[]);
          //           D.x[] = (fs_temp.x[]*Tkl + Tkg*(1. - fs_temp.x[]))*fm.x[]; //*max(y,1e-20);//*fm.x[]; //max(y,1e-20);
          //           // if(phase_flag3==1){
          //           //     D.x[] = Tkl*fm.x[];
          //           // }else{
          //           //     D.x[] = Tkg*fm.x[];
          //           // }
          //    }else if(fss_test3_n.x[]<=0){
          //           D.x[] = Tks*fm.x[]; //*max(y,1e-20);//*fm.x[];;
          //           //printf("D.x = %g\n",D.x[]);
          //    }else{
          //           double val_fluid,val_solid;
          //          // val_fluid = fss_test.x[]*Tkl + Tkg*(1. - fss_test.x[]);
          //           val_fluid = fs_temp.x[]*Tkl + Tkg*(1. - fs_temp.x[]);
          //           val_solid = Tks;
          //           D.x[] = (fss_test3_n.x[]*val_fluid + (1.0-fss_test3_n.x[])*val_solid)*fm.x[]; //*max(y,1e-20);//*fm.x[];//max(y,1e-20);
          //    }
          // }


//         // if(i%dump_each_event_interval0==0){
// if(dump_each_event){
//           char dumpname2[80];
//         // topo_mask_s.nodump = true;
//           foreach_dimension(){
//             hhh.x.nodump = true;
//             //  hhh_ff_oppo.x.nodump = true;
//           }
//           sprintf(dumpname2,"outfacets/dumpaterdiff1-i%d-t%g",i,t);
//           dump(dumpname2);
//         // }
// }

  // if(point_trace){
  //     double dd=-1;
  //     foreach(reduction(max:dd)){
  //       if(fabs(x-tracex)<1e-6 && fabs(y-tracey)<1e-6){
  //          dd = T[];
  //       }
  //     }
  //     if(pid()==0){
  //          char name93[80];
  //           sprintf(name93,"point_trace.dat");
  //           FILE * fp93 = fopen(name93,"a");
  //           fprintf(fp93,"%.8f %d %g\n",t,4,dd); //1 represent after diff1
  //           fclose(fp93);
  //     }
  // }

      //  printf("finish Tlff Tgff\n");
       ////////recalculate - rhocp
       //event ("properties");
       /////D.x is unchanged
       ////////////////// here, solve the diffusion equation for T_solid as one-fluid method
       // source_term for triple point
      //    scalar tempggg[],tempgll[];//,poisson_source2[];
      // //  // Tgrad_leon(ff,phase0Tgrad,phase1Tgrad);
      //  Tgrad_leon(ff,tempggg,tempgll);

      // // // // // // // // // // //  Tgrad_leon(ff,phase0Tgrad,phase1Tgrad); //comment 20231013
      // if(i>=1){

      // Tgrad_leon(ff,phase0Tgrad,phase1Tgrad);
      // }
      //  Tgrad_leon(ff,phase0Tgrad,phase1Tgrad);

      //poisson_source2 is the heater soulid face

        vertex scalar phi_temp2[];
      
               foreach_vertex(){
                  // phi_temp[] =  -x;//+Delta;
                  phi_temp2[] =  -(x-(-thickbottom+tune2_value));
                  // phi_temp[] =  -x + L0_pysical/(1<<maxl)/2.0;
               }
        face vector poisson_source2_f_fs[];
               fractions (phi_temp2, poisson_source2_f, poisson_source2_f_fs);
              foreach(){
                    poisson_source2_f[] = clamp (poisson_source2_f[], 0., 1.);
               }

       foreach(){
           poisson_source2[] = 0.0;
       }
       double percent_s=1.0,percent_l=0.0,percent_g=0.0;
       if(surface_heat){
          foreach(){
            if(poisson_source2_f[]>0 && poisson_source2_f[]<1.0){
                  coord n = interface_normal(point, poisson_source2_f);//interface_normal7(point,ff,hhh); //
                  double alphaa = plane_alpha (poisson_source2_f[], n);
                  coord pp;
                  double area1 = plane_area_center (n, alphaa, &pp);
                  if (metric_embed_factor){
                      area1 *= metric_embed_factor (point, pp);
                      // printf("metric_embed_factor is true!!!!!!!!!!!!!!!!!!!!!!!!\n");
                  }

                    poisson_source2[] = max(heat_flux*(1.0-y/Length_heat),0.0)*area1/Delta;
            }
          }
       }else{
          foreach(){
            if(poisson_source2_f[]>0 && poisson_source2_f[]<1.0){
                  coord n = interface_normal(point, poisson_source2_f);//interface_normal7(point,ff,hhh); //
                  double alphaa = plane_alpha (poisson_source2_f[], n);
                  coord pp;
                  double area1 = plane_area_center (n, alphaa, &pp);
                  if (metric_embed_factor){
                      area1 *= metric_embed_factor (point, pp);
                      // printf("metric_embed_factor is true!!!!!!!!!!!!!!!!!!!!!!!!\n");
                  }

                    // poisson_source2[] = max(heat_flux*(1.0-y/0.004),0.0)*area1/Delta;
                    poisson_source2[] = ((y<range_heat_flux?(heat_flux):0.0)*area1/Delta);
            }
          }
       }
// if(!restartsymbol){
//     if(case_number==3 && (!restart_Tsat) && surface_heat){
//       foreach(){
//            if(css_test3[]>1e-6 && css_test3[]<1.0-1e-6){
//                 // poisson_source2[] += ((y<range_heat_flux?(heat_flux):0.0)*area1/Delta);
//                 coord n = interface_normal(point, css_test3);//interface_normal7(point,ff,hhh); //
//                 double alphaa = plane_alpha (css_test3[], n);
//                 coord pp;
//                 double area1 = plane_area_center (n, alphaa, &pp);

//                 // poisson_source2[] += (max(heat_flux*(1.0-y/0.004)/Tks,0.0)*area1/Delta)*cm[];
//               if(surface_heat){
//                 poisson_source2[] += (max(heat_flux*(1.0-y/0.004),0.0)*area1/Delta)*cm[];
//               }else{
//                 poisson_source2[] += (heat_flux*area1/Delta)*cm[];
//               }

//                 // fprintf(stderr,"area1=%g heat_flux=%g range=%g ps2[]=%g\n",area1,heat_flux,range_heat_flux,poisson_source2[]);
//            }
//       }
//     }
// }
//        foreach(){
//         if((css_test3_n[]<1.0-1e-6 && css_test3_n[]>1e-6) && (ff[]>1e-6 && ff[]<1.0-1e-6)){
//           coord n = interface_normal(point, ff);//interface_normal7(point,ff,hhh); //
//           double alphaa = plane_alpha (ff[], n);
//           coord pp;
//           double area1 = plane_area_center (n, alphaa, &pp);
//           poisson_source2[] += (Tkl*tempgl[]-Tkg*tempgg[])*area1/Delta*cm[]; //20230603 -- 20230626
//           //poisson_source2[] = 0.0;
//         }
//        }

        // T.restriction = restriction_three_phase1;
        // T.refine = refine_three_phase1;
        // T.prolongation = refine_three_phase1; // this should be consistent with function used in posson3.h;
        // T.dirty = true;
      // #if EMBED
      //     int flag4=0;
      //     foreach(reduction(+:flag4)){
      //         if(css_test3_n[]<1.0){
      //             flag4+=1;
      //         }
      //     }
      //     if(flag4>=1){
      //         printf("flag4=%d\n",flag4);
      //         event ("properties");
      //         scalar poisson_source3[];
      //         foreach(){
      //           poisson_source3[] = poisson_source2[];
      //         }
      //         //diffusion5(T, dt, D, theta=rhocp);
      //         // diffusion5(T, dt, D, r=poisson_source3, theta=rhocp); //20230603 //20200918
      //         event ("properties");
      //     }

      // #endif

       foreach(){
          //if(css_test3_n[]>=1.0){
             if(css_test3_n[]>=1.0){
                  if(fabs(ff[]-0.5)<1e-6){
                      T[] = Tsat00;
                   }
          }
      }


// if(dump_each_event){
//             // // if(i%dump_each_event_interval0==0){
//                       char dumpname3[80];
//               // topo_mask_s.nodump = true;
//                 foreach_dimension(){
//                   hhh.x.nodump = true;
//                   //  hhh_ff_oppo.x.nodump = true;
//                 }
//                 sprintf(dumpname3,"outfacets/dumpaterdiff5-i%d-t%g",i,t);
//                 dump(dumpname3);
//             // // }
// }

//   if(point_trace){
//       double dd=-1;
//       foreach(reduction(max:dd)){
//         if(fabs(x-tracex)<1e-6 && fabs(y-tracey)<1e-6){
//            dd = T[];
//         }
//       }
//       if(pid()==0){
//            char name93[80];
//             sprintf(name93,"point_trace.dat");
//             FILE * fp93 = fopen(name93,"a");
//             fprintf(fp93,"%.8f %d %g\n",t,5,dd); //5 represent after diff5
//             fclose(fp93);
//       }
//   }
#endif



    foreach(){
      cs[] = css_test3_n[];
    }
    foreach_face(){
      fs.x[] = fss_test3_n.x[];
    }
    boundary({cs,fs});
               #if EMBED
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});

              #endif
     event ("properties") ;
                // restriction({css_test,fss_test});
                // restriction({css_test2,fss_test2});
                // restriction({css_test3,fss_test3});
                // restriction({css_test3_n,fss_test3_n});

   get_modphase01_2(css_test3,modphase_s_1,modphase_s_0);
  //  get_solidsurfaceT();
   get_css_fss_areaslg_triple_point(); //
  #if EMBED
    restriction ({cm_css_test, fss_test});
    cm_css_test_update (cm_css_test, css_test, fss_test);
    fm_fss_test_update (fm_fss_test, css_test, fss_test);
    restriction ({cm_css_test, fss_test});
  #endif
  #if EMBED
    restriction ({cm_css_test2, fss_test2});
    cm_css_test2_update (cm_css_test2, css_test2, fss_test2);
    fm_fss_test2_update (fm_fss_test2, css_test2, fss_test2);
    restriction ({cm_css_test2, fss_test2});
  #endif
  #if EMBED
    restriction ({cm_css_test3, fss_test3});
    cm_css_test3_update (cm_css_test3, css_test3, fss_test3);
    fm_fss_test3_update (fm_fss_test3, css_test3, fss_test3);
    restriction ({cm_css_test3, fss_test3});
  #endif

  if(!use_Tslg){
      foreach(){
          //   Tl[] = 0.0; //Tsat00;
          //   Tg[] = 0.0; //Tsat00;
          //   Ts[] = 0.0; //Tsat00;
            Tl[] = Tsat00;
            Tg[] = Tsat00;
            Ts[] = Tsat00;
            if(css_test3[]>=0.5){
                  Ts[] = T[];
            }else{
                if(css_test[]>=0.5){
                    Tl[] = T[];
                }
                if(css_test2[]>=0.5){
                    Tg[] = T[];
                }
            }
        }
        boundary({Ts,Tl,Tg});
  }else{

  }
   solver_new(poisson_source2,percent_s,percent_l,percent_g);
   foreach(){
      if(css_test3[]>=0.5){
          T[] = Ts[];
      }else{
          if(css_test[]>=css_test2[]){
              T[] = Tl[];
          }else{
              T[] = Tg[];
          }
      }
   }
   boundary({T});

}




event poisson_ps (i++,last){ // laplacian ps = m_dot


    //source term is non-zero only in cells with maximum depth() near interface
    //so here we only need to solve the poisson equation in that region
 get_topomask(topo_mask);
 event ("properties");

 //poisson_ps_usf1(alpha_rho,ps,usf,topo_mask,level_interface,dt);


 foreach_face(){
    usf.x[] = 0.0;
    usfg.x[] = 0.0;
 }

if(i>=0){
   if(!restartsymbol){
        // if(!(case_number==3 && (!restore (file = "restart_Tsat0")))){
        if(!(case_number==3 && (!restart_Tsat))){
             if(1==1){
                    // scalar div_numerical[];
                    foreach() {
                      div_numerical[] = 0.;
                      if(cs[]>0.0 && topo_mask[]==0){
                          foreach_dimension() 
                            div_numerical[] += uf.x[1] - uf.x[];
                          // div[] /= dt*Delta;
                          div_numerical[] /= Delta;
                      }
                    }
                      // add phase change source term
                      //div[] = div[] - sv[]*cm[]/dt; //\nabola uf - source_pc
                      // div[] = div[] - sv[]/dt;
                      //div[] = div[] - sv[]/dt; //\nabola uf - source_pc
                      poisson_ps_usf1(alpha,ps,usf,topo_mask,div_numerical,level_interface,dt);       /////////////////////checked for 3D
                    foreach() {
                      div_numerical2[] = 0.;
                      if(cs[]>0.0 && topo_mask[]==0){
                          foreach_dimension() 
                            div_numerical2[] += usf.x[1] - usf.x[];
                          // div[] /= dt*Delta;
                          div_numerical2[] /= Delta;
                      }
                    }
                    
            }else{
              poisson_ps_usf1(alpha,ps,usf,topo_mask,source_pc2,level_interface,dt);       /////////////////////checked for 3D

            }// poisson_ps_usf1(alpha,ps,usf,topo_mask,source_pc,level_interface,dt);
        }else{
          foreach(){
            ps[] = 0.0;
            source_pc[] = 0.0;
            source_pc2[] = 0.0;
          }
          foreach_face(){
            usf.x[] = 0.0;
            usfg.x[] = 0.0;
          }
        }
   }else{
        poisson_ps_usf1(alpha,ps,usf,topo_mask,source_pc2,level_interface,dt);       /////////////////////checked for 3D
   }
}

 foreach(){
    topo_mask_g[] = -topo_mask[];
 }
//  poisson_ps_usf1g (alpha,psg,usfg,topo_mask_g,level_interface,dt);


 //poisson_ps_usf1();

 /*
    //output ps and usf topo_mask and level to check
    char name115[80];
    sprintf(name115,"poisson-check%g.dat", t);
    FILE * fp115 = fopen(name115,"w");
    foreach(){
        if(fabs(source_pc[])>1e-15){
           fprintf(fp115,"%g %g %g %g %g %g %d\n",x,y,source_pc[],masstr[],ps[],topo_mask[],level);
        }
    }
    fclose(fp115);
 */


}




event get_ulf(i++,last){
  //  ulf_function(); // checked for 3D
  if(i==51840){
      dump("i51840");
  }
  if(i>=1){
    ulf_ugf_function(topo_mask,topo_mask_g);
  }
}






#if TREE
event update_Tl_Tg(i++,last) {
  // css_test3.restriction = restriction_average;
   // css_test3_n.restriction = restriction_average;

}




double smff(double cc[3][3][3]){
    double value;
	    // sf[] = (8.*ff[] +
		  //   4.*(ff[-1] + ff[1] + ff[0,1] + ff[0,-1] + ff[0,0,1] + ff[0,0,-1]) +
		  //   2.*(ff[-1,1] + ff[-1,0,1] + ff[-1,0,-1] + ff[-1,-1] +
			// ff[0,1,1] + ff[0,1,-1] + ff[0,-1,1] + ff[0,-1,-1] +
			// ff[1,1] + ff[1,0,1] + ff[1,-1] + ff[1,0,-1]) +
		  //   ff[1,-1,1] + ff[-1,1,1] + ff[-1,1,-1] + ff[1,1,1] +
		  //   ff[1,1,-1] + ff[-1,-1,-1] + ff[1,-1,-1] + ff[-1,-1,1])/64.;
	    value = (8.*cc[1][1][1] +
		    4.*(cc[0][1][1] + cc[2][1][1] + cc[1][2][1] + cc[1][0][1] + cc[1][1][2] + cc[1][1][0]) +
		    2.*(cc[0][2][1] + cc[0][1][2] + cc[0][1][0] + cc[0][0][1] +
			cc[1][2][2] + cc[1][2][0] + cc[1][0][2] + cc[1][0][0] +
			cc[2][2][1] + cc[2][1][2] + cc[2][0][1] + cc[2][1][0]) +
		    cc[2][0][2] + cc[0][2][2] + cc[0][2][0] + cc[2][2][2] +
		    cc[2][2][0] + cc[0][0][0] + cc[2][0][0] + cc[0][0][2])/64.;

        return value;
}


// event extract1 (t+=out_interval)
event extract1 (t=BASE_TIME; t+=out_interval)
// event extract1 (i++)
// event extract1 (i+=10)
{
  char name33[80];
  sprintf(name33,"mass_record2-pid%d.dat",pid());
  FILE * fp33 = fopen(name33,"w");
  foreach(){
    //fprintf(fp33,"%g %g %g %g\n",x,y,T[],cs[]);
    fprintf(fp33,"%g %g %g %g %g %g %g %g %g %g\n",x,y,ff[],T[],source_pc[],topo_mask[],phase0Tgrad[],phase1Tgrad[],aiml[],aimg[]);
  }

  fclose(fp33);
   MPI_Barrier(MPI_COMM_WORLD);
   if(pid()==0){
                  char command1[150];
                  sprintf(command1, "LC_ALL=C cat mass_record2-pid*.dat > outfacets/mass_record2-%g",t);
                  system(command1);

                  char command7[150];
                  sprintf(command7, "LC_ALL=C rm -rf mass_record2-pid*.dat");
                  system(command7);
    }


        char names[80];
	sprintf(names, "interface%d", pid());
	FILE * fp = fopen (names, "w");
	output_facets (ff,fp);
	fclose(fp);

  MPI_Barrier(MPI_COMM_WORLD);
  if(pid()==0){
    char command[80];
    sprintf(command, "LC_ALL=C  cat interface* > outfacets/interface_%d_%g.dat",i,t);
    system(command);
  }


  during_time = (clock() - end_time)/(double)CLOCKS_PER_SEC;
end_time = clock();
      if(pid()==0){
              char name117[80];
              sprintf(name117,"time-performance.dat");
              FILE * fp117 = fopen(name117,"a");
              fprintf(fp117,"%d %g %g %g\n",i,t,during_time,(end_time-start_time)/(double)CLOCKS_PER_SEC);
              fclose(fp117);
      }

}

// event pictures (i++)
// event pictures (t+=out_interval2)
event pictures (t=BASE_TIME; t+=out_interval2)
// event pictures (i+=10)
// event pictures (i+=1)
{
  foreach(){
      ulf_v.x[] = (ulf.x[]/y+ulf.x[1,0]/y)/2.0;
      ulf_v.y[] = (ulf.y[]/max(1e-20,y-Delta/2.0)+ulf.y[0,1]/max(1e-20,y+Delta/2.0))/2.0;
  }

  printf("t=%g\n",t);
    char dumpname[80];
 // topo_mask_s.nodump = true;
  foreach_dimension(){
     hhh.x.nodump = true;
    //  hhh_ff_oppo.x.nodump = true;
  }

  // get_u_ghost(u);

  sprintf(dumpname,"dump-i%d-t%g",i,t);
  dump(dumpname);
//   char name[80];
// #if 1
//   foreach_dimension(){
//      hhh.x.nodump = true;
//     // hhh_ff_oppo.x.nodump = true;
//   }
//   sprintf (name, "outfacets/dump-%g", t);
//   dump (name);
// #endif
//   clear();
//  // double ty[] = { - 0.6, - 0.5, - 0.5, - 0.5};
//   view (fov = 19, quat = {0.707,0.707,0,0},
// 	tx = 0.0, ty = -L0/2.0,
// 	width = 800, height = 400);
//   squares ("T", spread = -1, linear = true, map = cool_warm);
//  // draw_vof ("f");
//   mirror ({0,1}) {
//  //   squares ("u.x", spread = -1, linear = true, map = cool_warm);
//       draw_vof ("ff");
//       squares ("T", spread = -1, linear = true, map = cool_warm);
//       //isoline("T",val=2.5);
//       //isoline("T",val=2);
//   }
//   //sprintf (name, "final-%d.png", dcase);
//   //isoline("T",val=2.5);
//   //isoline("T",val=2);
//   cells();
//   sprintf (name, "outfacets/image-%05d-%g.png", i, t);
//   save (name);
}


event adapt (i++) {

if(dump_each_event && (i%dump_each_event_interval==0)){
        char dumpname[80];
    // topo_mask_s.nodump = true;
      foreach_dimension(){
        hhh.x.nodump = true;
        //  hhh_ff_oppo.x.nodump = true;
      }
      sprintf(dumpname,"outfacets/dumpbeforeadapt-i%d-t%g",i,t);
      dump(dumpname);


}

  // T.refine = refine_bilinear;
  //adapt_wavelet ({T}, (double[]){0.01}, maxlevel = maxl, minlevel= minl);
  scalar df1[],df2[],df3[];
 vertex scalar phi_temp3[];
        foreach_vertex(){
            phi_temp3[] =  -(x-(-thickbottom+tune2_value));
        }
        fractions (phi_temp3, poisson_source2_f);

  // foreach(){
  //  df[] = max(ff[],css_test3[]);
  // }
// printf("10-1 adapt_wavelet pid%d\n",pid());
boundary({ff});
  smooth_for_arm(df1,ff);
  smooth_for_arm(df2,css_test3_n);
  smooth_for_arm(df3,poisson_source2_f);
 // adapt_wavelet ({T,ff}, (double[]){0.01,0.001}, maxlevel = maxl, minlevel= minl);
  ff.refine = fraction_refine;
if(!restartsymbol){
      if(case_number == 1){
        adapt_wavelet ({df1},(double[]){0.001}, maxlevel = maxl, minlevel = minl);
      }else if(case_number==2){
        adapt_wavelet ({df1,df2},(double[]){0.001,0.001}, maxlevel = maxl, minlevel = minl);
      // }else if(case_number == 3 && !restore (file = "restart_Tsat0")){
        }else if(case_number == 3 && (!restart_Tsat)){ //surface_heat
      // adapt_wavelet ({df2,T},(double[]){0.001,0.01}, maxlevel = maxl, minlevel = maxl);
      // }else if(case_number == 3 && restore (file = "restart_Tsat0")){
                 scalar T_adapt[];
                foreach(){
                  T_adapt[] = T[];
                }
                T_adapt.refine  = bilinear_no_cs2;
                T_adapt.restriction = restriction_average;
        adapt_wavelet ({df1,df2,T_adapt,df3},(double[]){0.001,0.001,0.1,0.001}, maxlevel = maxl, minlevel = minl);


      }else if((case_number == 3 && restart_Tsat)){
           scalar T_adapt[];
        foreach(){
          T_adapt[] = T[];
        }
        T_adapt.refine  = bilinear_no_cs2;
        T_adapt.restriction = restriction_average;
        // adapt_wavelet ({df1,df2,T_adapt},(double[]){0.001,0.001,0.01}, maxlevel = maxl, minlevel = minl);

         scalar ff_temp2[];
        double max_yy=0.0;
        double temp2=L0/(1<<maxl);
        foreach(reduction(max:max_yy)){
            if(css_test3_n[]>0.0){
                if(ff[]<1.0 && ff[]>0.0 && (y>max_yy)){
                    max_yy = y;
                }
            }
        }
        max_yy = 1.1*max_yy;
        // srand(time(0));
        // foreach(){
          
        //   if(css_test3_n[]>0.0){
        //     double temp = sqrt(sq(max_yy) - (sq(x) + sq(y)));
        //     if(temp>0.0){
                
        //         // Generate a random floating-point number between 0 and 1
        //         double randomValue = ((double)rand() / RAND_MAX-0.5);
        //       ff_temp2[] = randomValue*0.001;//temp+randomValue*temp2*8 ;
        //     }else{
        //       ff_temp2[] = 0.0;
        //     }
        //   }else{
        //     ff_temp2[] = 0.0;
        //   }

        // }
        // adapt_wavelet ({df1,df2,T_adapt,ff_temp2},(double[]){0.001,0.001,0.001,0}, maxlevel = maxl, minlevel = minl);
        // adapt_wavelet ({df1,df2,T_adapt},(double[]){0.001,0.001,0.1}, maxlevel = maxl, minlevel = minl);

       


         adapt_wavelet ({df1,df2,T_adapt,df3},(double[]){0.001,0.001,0.1,0.001}, maxlevel = maxl, minlevel = minl);






        // adapt_wavelet ({df1,df2,T_adapt},(double[]){0.001,0.001,0.001}, maxlevel = maxl, minlevel = minl);
          // adapt_wavelet ({df1,df2},(double[]){0.001,0.001}, maxlevel = maxl, minlevel = minl);
      }else{
        printf("case error\n");
        exit(1);
      }
}else{
  // adapt_wavelet ({df1,df2},(double[]){0.001,0.001}, maxlevel = maxl, minlevel = minl);
   scalar T_adapt[];
        foreach(){
          T_adapt[] = T[];
        }
        T_adapt.refine  = bilinear_no_cs2;
        T_adapt.restriction = restriction_average;
        // adapt_wavelet ({df1,df2,T_adapt},(double[]){0.001,0.001,0.01}, maxlevel = maxl, minlevel = minl);



        scalar ff_temp2[];
        double max_yy=0.0;
        double temp2=L0/(1<<maxl);
        foreach(reduction(max:max_yy)){
            if(css_test3_n[]>0.0){
                if(ff[]<1.0 && ff[]>0.0 && (y>max_yy)){
                    max_yy = y;
                }
            }
        }
        max_yy = 1.1*max_yy;
        // srand(time(0));
        // foreach(){
        //   if(css_test3_n[]>0.0){
        //     double temp = sqrt(sq(max_yy) - (sq(x) + sq(y)));
        //     if(temp>0.0){
        //         // Generate a random floating-point number between 0 and 1
        //         double randomValue = ((double)rand() / RAND_MAX-0.5);
        //       ff_temp2[] = randomValue*0.001;//temp+randomValue*temp2*8 ;
        //     }else{
        //       ff_temp2[] = 0.0;
        //     }
        //   }else{
        //     ff_temp2[] = 0.0;
        //   }

        // }
        // adapt_wavelet ({df1,df2,T_adapt,ff_temp2},(double[]){0.001,0.001,0.001,0}, maxlevel = maxl, minlevel = minl);
        // adapt_wavelet ({df1,df2,T_adapt},(double[]){0.001,0.001,0.001}, maxlevel = maxl, minlevel = minl);
        // adapt_wavelet ({df1,df2,T_adapt},(double[]){0.001,0.001,0.1}, maxlevel = maxl, minlevel = minl);
         adapt_wavelet ({df1,df2,T_adapt,df3},(double[]){0.001,0.001,0.1,0.001}, maxlevel = maxl, minlevel = minl);
}
  if(!restartsymbol){
    if(case_number==2 && i==0){
      for(int ii=1;ii<=5;ii++){
          fraction (ff, (- (sq(bubble_radius) - sq(x - tune2_value- centerx) - sq(y - centery))));
              foreach(){
                  ff[] = clamp (ff[], 0., 1.);
              }
              boundary({ff});
          smooth_for_arm(df1,ff);
                vertex scalar phi_temp[];
                foreach_vertex(){
                      // phi_temp[] = 0.0; //- z + thickness;
                      // phi_temp[] = 0.25 - z;// z;
                      // phi_temp[] = -x;// z;
                      phi_temp[] =  -(x-tune2_value);
                  }
                  fractions (phi_temp, css_test3, fss_test3);
                  foreach(){
                    fs_solid[] = css_test3[];
                    css_test3_n[] = 1.0 - css_test3[];
                  }
                  foreach_face(){
                      fss_test3_n.x[] = 1.0 - fss_test3.x[];
                  }
          smooth_for_arm(df2,css_test3_n);
         // adapt_wavelet ({df1,df2,T},(double[]){0.001,0.001,0.001}, maxlevel = maxl, minlevel = minl);
          adapt_wavelet ({df1,df2},(double[]){0.001,0.001}, maxlevel = maxl, minlevel = minl);

   }
    // tanguy, same with original case_number == 8
               level_interface = maxl;
               fraction (ff, (- (sq(bubble_radius) - sq(x - tune2_value - centerx) - sq(y - centery))  ));
               foreach(){
                      ff[] = clamp (ff[], 0., 1.);
               }
               boundary({ff});
               const scalar c[] = theta0*pi/180.;
               contact_angle = c;
               foreach(){
                    ff_oppo[] = 1.0 - ff[];
               }

              ff_oppo.refine=ff_oppo.prolongation = fraction_refine;
              #if TREE
                  for (scalar c in interfaces) {
                      c.refine = c.prolongation = fraction_refine;
                      c.dirty = true;
                  }
              #endif
              //  vertex scalar phi_temp[];

              //  foreach_vertex(){
              //     phi_temp[] = 1;//-x;
              //  }
              //  fractions (phi_temp, css_test3, fss_test3);
              // foreach(){
              //       fs_solid[] = css_test3[];
              //       css_test3_n[] = 1.0 - css_test3[];
              //  }
              //  foreach_face(){
              //         fss_test3_n.x[] = 1.0 - fss_test3.x[];
              // }
              // boundary({css_test3,fss_test3});
              // boundary({fs_solid});
              // boundary({css_test3_n,fss_test3_n});
              // foreach(){
              //       cs[] = css_test3_n[];
              // }
              // foreach_face(){
              //       fs.x[] = fss_test3_n.x[];
              // }
              // boundary({cs,fs});
                vertex scalar phi[];
                fraction(css_test3_n,x-tune2_value);
                foreach(){
                  css_test3_n[] =  clamp (css_test3_n[], 0., 1.);
                }
                vof2dist(css_test3_n,phi);
                boundary ({phi});
                fractions (phi, css_test3_n, fss_test3_n);
                              foreach(){
                                  css_test3[] = 1.0 - css_test3_n[];
                                  fs_solid[] = css_test3[];
                                }
                                foreach_face(){
                                    fss_test3.x[] = 1.0 - fss_test3_n.x[];
                                }


                  foreach(){
                    cs[] = css_test3_n[];
                  }
                  foreach_face(){
                    fs.x[] = fss_test3_n.x[];
                  }
                  boundary({cs,fs});
               #if EMBED
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});
              #endif
              //  foreach(){
              //       if(x<0.0){
              //           T[] = T_inf - (T_inf - Twall_init)*(x-(-0.0015))/0.0015;
              //       }else{
              //           T[] = Twall_init - (Twall_init-Tsub)*(x-0.0)/0.003;
              //       }
              //   }
               if(case_number==2){
                  foreach(){
                        if(x<0.0){
                            T[] = T_inf - (T_inf - Twall_init)*(x-originx)/fabs(originx);
                        }else{
                            T[] = Twall_init - (Twall_init-Tsub)*(x-0.0)/fabs(L0_pysical-fabs(originx));
                        }
                    }
              }else if(case_number==3){
                  foreach(){
                        T[] = Tsat00;
                    }
              }
                boundary({T});
                fs_solid.prolongation = fs_solid.refine = fraction_refine;
                css_test3.refine = embed_fraction_refine_s;
                css_test3.prolongation = fraction_refine; //embed_fraction_refine_s; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3.x.prolongation = embed_face_fraction_refine_s_x;
                css_test3_n.refine = embed_fraction_refine_s_n;
                css_test3_n.prolongation = fraction_refine; //embed_fraction_refine_s_n; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3_n.x.prolongation = embed_face_fraction_refine_s_n_x;

        foreach(){
                  if(css_test3_n[]<=0){
                      css_test[]=0;
                      css_test2[]=0;
                  }else if(css_test3_n[]>=1){
                      css_test[]=ff[];
                      css_test2[]=1.0-ff[];
                  }else{
                      //?????????????????????
                      css_test[]=ff[]*css_test3_n[];
                      css_test2[]=(1.0-ff[])*css_test3_n[];
                  }
               }

        face vector fs_temp[];
       if(1==1){
               vertex scalar phi21[];
               vof2dist(css_test,phi21);
               scalar ff_temp[];
               fractions(phi21, ff_temp, fss_test);
               vof2dist(css_test2,phi21);
               fractions(phi21, ff_temp, fss_test2);
               foreach_face(){
                  fs_temp.x[] = fss_test.x[];
               }


       }else{

                vertex scalar phi21[];
                scalar ff_temp[];
                vof2dist(ff, phi21);
                fractions (phi21, ff_temp, fs_temp);

                foreach(){
                  css_test[] = 0.0;
                  css_test2[] = 0.0;
                }
                foreach_face(){
                  fss_test.x[] = 0.0;
                  fss_test2.x[]=0.0;
                }
                foreach(){
                  if(css_test3_n[]>0.0){
                    css_test[] = ff_temp[];
                    css_test2[] = 1.0 - ff_temp[];
                  }
                }
                foreach_face(){
                  if(fss_test3_n.x[]>0.0){
                      fss_test.x[] = fs_temp.x[];
                      fss_test2.x[] = 1.0 - fs_temp.x[];
                  }
                }
    }
                 css_test.refine = embed_fraction_refine_css_test;
                css_test.prolongation = embed_fraction_refine_css_test;//fraction_refine;
                foreach_dimension(){
                          fss_test.x.refine = fss_test.x.prolongation = embed_face_fraction_fss_test_refine_x;
                          //fss_test.x.dirty = false;
                }
                css_test2.refine = embed_fraction_refine_css_test2;
                css_test2.prolongation = embed_fraction_refine_css_test2;//fraction_refine;
                foreach_dimension()
                          fss_test2.x.prolongation = embed_face_fraction_fss_test2_refine_x;



                restriction({css_test,fss_test});
                restriction({css_test2,fss_test2});
                restriction({css_test3,fss_test3});
                restriction({css_test3_n,fss_test3_n});

              //add 20230915
               get_css_fss_areaslg_triple_point();

  #if EMBED
    restriction ({cm_css_test, fss_test});
    cm_css_test_update (cm_css_test, css_test, fss_test);
    fm_fss_test_update (fm_fss_test, css_test, fss_test);
    restriction ({cm_css_test, fss_test});
  #endif
  #if EMBED
    restriction ({cm_css_test2, fss_test2});
    cm_css_test2_update (cm_css_test2, css_test2, fss_test2);
    fm_fss_test2_update (fm_fss_test2, css_test2, fss_test2);
    restriction ({cm_css_test2, fss_test2});
  #endif
  #if EMBED
    restriction ({cm_css_test3, fss_test3});
    cm_css_test3_update (cm_css_test3, css_test3, fss_test3);
    fm_fss_test3_update (fm_fss_test3, css_test3, fss_test3);
    restriction ({cm_css_test3, fss_test3});
  #endif



                heights(ff,hhh);
              //  heights(ff_oppo,hhh_ff_oppo);
                 printf("event init 482,pid=%d\n",pid());
                      foreach(){
                        source_pc[] = 0.0;
                        masstr[]=0.0;
                        vtr[] = 0.0;
                      }

          ///////////////////////////////////////////////////////
        //  T.refine  = refine_bilinear;
        //  T.restriction = restriction_volume_average;//restriction_average;//restriction_volume_average;
        //  T.dirty = true;
          globali = 0;
            foreach_face(){
              ulf.x[] =0.0;
              usf.x[] =0.0;
              uf.x[] = 0.0;

              ugf.x[] =0.0;
              usfg.x[] =0.0;

            }
            foreach(){
              p[]=0.;
              pf[]=0.0;
              foreach_dimension(){
                u.x[]=0.0;
                g.x[]=0.0;
              }
            }
            boundary({u.x,u.y});
 // }else if(case_number==3 && i==0 && restore (file = "restart_Tsat0")){
  }else if(case_number==3 && restart_Tsat_i==0 && restart_Tsat){
    restart_Tsat_i=1;
    fprintf(stderr,"event adapt=%d t=%g restart_Tsat_i=%d\n",i,t,restart_Tsat_i);
    for(int ii=1;ii<=0;ii++){
          // fraction (ff, (- (sq(bubble_radius) - sq(x - centerx) - sq(y - centery))));
          fraction (ff, (- (sq(bubble_radius) - sq(x - tune2_value - centerx) - sq(y - centery))  ));
              foreach(){
                  ff[] = clamp (ff[], 0., 1.);
              }
              boundary({ff});
          smooth_for_arm(df1,ff);
                vertex scalar phi_temp[];
                // foreach_vertex(){
                //       // phi_temp[] = 0.0; //- z + thickness;
                //       // phi_temp[] = 0.25 - z;// z;
                //       //  phi_temp[] = -x ;
                //       phi_temp[] =  -(x-tune2_value);
                //       // phi_temp[] = -x + L0_pysical/(1<<maxl)/2.0;// z;
                //   }
                 solid_phi(tune2_value, phi_temp);
              foreach_vertex(){
                phi_temp[] = -phi_temp[];
              }
                  fractions (phi_temp, css_test3, fss_test3);
                  foreach(){
                    fs_solid[] = css_test3[];
                    css_test3_n[] = 1.0 - css_test3[];
                  }
                  foreach_face(){
                      fss_test3_n.x[] = 1.0 - fss_test3.x[];
                  }
          smooth_for_arm(df2,css_test3_n);

           vertex scalar phi_temp3[];
               foreach_vertex(){
                  phi_temp3[] =  -(x-(-thickbottom+tune2_value));
               }
              fractions (phi_temp3, poisson_source2_f);
               smooth_for_arm(df3,poisson_source2_f);



         // adapt_wavelet ({df1,df2,T},(double[]){0.001,0.001,0.001}, maxlevel = maxl, minlevel = minl);
          adapt_wavelet ({df1,df2,df3},(double[]){0.001,0.001,0.001}, maxlevel = maxl, minlevel = minl);

   }
    // tanguy, same with original case_number == 8
               level_interface = maxl;
              //  fraction (ff, (- (sq(bubble_radius) - sq(x - centerx) - sq(y - centery))  ));
              fraction (ff, (- (sq(bubble_radius) - sq(x - tune2_value - centerx) - sq(y - centery))  ));
               foreach(){
                      ff[] = clamp (ff[], 0., 1.);
               }
               boundary({ff});
               const scalar c[] = theta0*pi/180.;
               contact_angle = c;
               foreach(){
                    ff_oppo[] = 1.0 - ff[];
               }

              ff_oppo.refine=ff_oppo.prolongation = fraction_refine;
              #if TREE
                  for (scalar c in interfaces) {
                      c.refine = c.prolongation = fraction_refine;
                      c.dirty = true;
                  }
              #endif
                // // vertex scalar phi[];
                // // fraction(css_test3_n,x);
                // // foreach(){
                // //   css_test3_n[] =  clamp (css_test3_n[], 0., 1.);
                // // }
                // // vof2dist(css_test3_n,phi);
                // // boundary ({phi});
                // // fractions (phi, css_test3_n, fss_test3_n);
                // //               foreach(){
                // //                   css_test3[] = 1.0 - css_test3_n[];
                // //                   fs_solid[] = css_test3[];
                // //                 }
                // //                 foreach_face(){
                // //                     fss_test3.x[] = 1.0 - fss_test3_n.x[];
                // //                 }

                //  vertex scalar phi_temp[];
                // foreach_vertex(){
                //       // phi_temp[] = 0.0; //- z + thickness;
                //       // phi_temp[] = 0.25 - z;// z;
                //       //  phi_temp[] = -x ;
                //       phi_temp[] =  -(x-tune2_value);
                //       // phi_temp[] = -x + L0_pysical/(1<<maxl)/2.0;// z;
                //   }
                //   fractions (phi_temp, css_test3, fss_test3);
                //   foreach(){
                //     fs_solid[] = css_test3[];
                //     css_test3_n[] = 1.0 - css_test3[];
                //   }
                //   foreach_face(){
                //       fss_test3_n.x[] = 1.0 - fss_test3.x[];
                //   }
vertex scalar phi[];
                // fraction(css_test3_n,x - tune2_value);
                // foreach(){
                //   css_test3_n[] =  clamp (css_test3_n[], 0., 1.);
                // }
                // vof2dist(css_test3_n,phi);
                solid_phi(tune2_value, phi);
                boundary ({phi});
                fractions (phi, css_test3_n, fss_test3_n);
                              foreach(){
                                  css_test3[] = 1.0 - css_test3_n[];
                                  fs_solid[] = css_test3[];
                                }
                                foreach_face(){
                                    fss_test3.x[] = 1.0 - fss_test3_n.x[]; 
                                }
                  //new

                  foreach(){
                    cs[] = css_test3_n[];
                  }
                  foreach_face(){
                    fs.x[] = fss_test3_n.x[];
                  }
                  boundary({cs,fs});
               #if EMBED
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});
              #endif
                  foreach(){
                        T[]=T_oold[];
                    }

                boundary({T});
                fs_solid.prolongation = fs_solid.refine = fraction_refine;
                css_test3.refine = embed_fraction_refine_s;
                css_test3.prolongation = fraction_refine; //embed_fraction_refine_s; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3.x.prolongation = embed_face_fraction_refine_s_x;
                css_test3_n.refine = embed_fraction_refine_s_n;
                css_test3_n.prolongation = fraction_refine; //embed_fraction_refine_s_n; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3_n.x.prolongation = embed_face_fraction_refine_s_n_x;

       foreach(){
                  if(css_test3_n[]<=0){
                      css_test[]=0;
                      css_test2[]=0;
                  }else if(css_test3_n[]>=1){
                      css_test[]=ff[];
                      css_test2[]=1.0-ff[];
                  }else{
                      //?????????????????????
                      css_test[]=ff[]*css_test3_n[];
                      css_test2[]=(1.0-ff[])*css_test3_n[];
                  }
               }
      face vector fs_temp[];
       if(1==1){
               vertex scalar phi21[];
               vof2dist(css_test,phi21);
               scalar ff_temp[];
               fractions(phi21, ff_temp, fss_test);
               vof2dist(css_test2,phi21);
               fractions(phi21, ff_temp, fss_test2);
               foreach_face(){
                  fs_temp.x[] = fss_test.x[];
               }


       }else{
                vertex scalar phi21[];
                scalar ff_temp[];
                vof2dist(ff, phi21);
                fractions (phi21, ff_temp, fs_temp);

                foreach(){
                  css_test[] = 0.0;
                  css_test2[] = 0.0;
                }
                foreach_face(){
                  fss_test.x[] = 0.0;
                  fss_test2.x[]=0.0;
                }
                foreach(){
                  if(css_test3_n[]>0.0){
                    css_test[] = ff_temp[];
                    css_test2[] = 1.0 - ff_temp[];
                  }
                }
                foreach_face(){
                  if(fss_test3_n.x[]>0.0){
                      fss_test.x[] = fs_temp.x[];
                      fss_test2.x[] = 1.0 - fs_temp.x[];
                  }
                }
            }
                 css_test.refine = embed_fraction_refine_css_test;
                css_test.prolongation = embed_fraction_refine_css_test;//fraction_refine;
                foreach_dimension(){
                          fss_test.x.refine = fss_test.x.prolongation = embed_face_fraction_fss_test_refine_x;
                          //fss_test.x.dirty = false;
                }
                css_test2.refine = embed_fraction_refine_css_test2;
                css_test2.prolongation = embed_fraction_refine_css_test2;//fraction_refine;
                foreach_dimension()
                          fss_test2.x.prolongation = embed_face_fraction_fss_test2_refine_x;

                 //add 20230915
               get_css_fss_areaslg_triple_point();
                restriction({css_test,fss_test});
                restriction({css_test2,fss_test2});
                restriction({css_test3,fss_test3});
                restriction({css_test3_n,fss_test3_n});


  #if EMBED
    restriction ({cm_css_test, fss_test});
    cm_css_test_update (cm_css_test, css_test, fss_test);
    fm_fss_test_update (fm_fss_test, css_test, fss_test);
    restriction ({cm_css_test, fss_test});
  #endif
  #if EMBED
    restriction ({cm_css_test2, fss_test2});
    cm_css_test2_update (cm_css_test2, css_test2, fss_test2);
    fm_fss_test2_update (fm_fss_test2, css_test2, fss_test2);
    restriction ({cm_css_test2, fss_test2});
  #endif
  #if EMBED
    restriction ({cm_css_test3, fss_test3});
    cm_css_test3_update (cm_css_test3, css_test3, fss_test3);
    fm_fss_test3_update (fm_fss_test3, css_test3, fss_test3);
    restriction ({cm_css_test3, fss_test3});
  #endif



                heights(ff,hhh);
              //  heights(ff_oppo,hhh_ff_oppo);
                 printf("event init 482,pid=%d\n",pid());
                      foreach(){
                        source_pc[] = 0.0;
                        masstr[]=0.0;
                        vtr[] = 0.0;
                      }

          ///////////////////////////////////////////////////////
        //  T.refine  = refine_bilinear;
        //  T.restriction = restriction_volume_average;//restriction_average;//restriction_volume_average;
        //  T.dirty = true;
          globali = 0;
            // foreach_face(){
            //   ulf.x[] =0.0;
            //   usf.x[] =0.0;
            //   uf.x[] = 0.0;

            //   ugf.x[] =0.0;
            //   usfg.x[] =0.0;

            // }
            // foreach(){
            //   p[]=0.;
            //   pf[]=0.0;
            //   foreach_dimension(){
            //     u.x[]=0.0;
            //     g.x[]=0.0;
            //   }
            // }
            boundary({u.x,u.y});
  // }else if(case_number==3 && i==0 && (!restore (file = "restart_Tsat0"))){
    }else if(case_number==3 && i==0 && (!restart_Tsat)){
         for(int ii=1;ii<=5;ii++){
              foreach(){
                  ff[] = 1.0;
              }
              boundary({ff});
          smooth_for_arm(df1,ff);
                vertex scalar phi_temp[];
                // foreach_vertex(){
                //       // phi_temp[] = 0.0; //- z + thickness;
                //       // phi_temp[] = 0.25 - z;// z;
                //       // phi_temp[] = -x;// z;
                //       phi_temp[] =  -(x-tune2_value);
                //       // phi_temp[] =  -x + L0_pysical/(1<<maxl)/2.0;
                //   }
                solid_phi(tune2_value, phi_temp);
              foreach_vertex(){
                phi_temp[] = -phi_temp[];
              }
                  fractions (phi_temp, css_test3, fss_test3);
                  foreach(){
                    fs_solid[] = css_test3[];
                    css_test3_n[] = 1.0 - css_test3[];
                  }
                  foreach_face(){
                      fss_test3_n.x[] = 1.0 - fss_test3.x[];
                  }
          smooth_for_arm(df2,css_test3_n);
         // adapt_wavelet ({df1,df2,T},(double[]){0.001,0.001,0.001}, maxlevel = maxl, minlevel = minl);
        //  adapt_wavelet ({df2},(double[]){0.001}, maxlevel = maxl, minlevel = minl);

   }
    // tanguy, same with original case_number == 8
               level_interface = maxl;
               foreach(){
                      ff[] = 1.0;
               }
               boundary({ff});
               const scalar c[] = theta0*pi/180.;
               contact_angle = c;
               foreach(){
                    ff_oppo[] = 1.0 - ff[];
               }

              ff_oppo.refine=ff_oppo.prolongation = fraction_refine;
              #if TREE
                  for (scalar c in interfaces) {
                      c.refine = c.prolongation = fraction_refine;
                      c.dirty = true;
                  }
              #endif
                vertex scalar phi[];
                fraction(css_test3_n,(x-tune2_value));
                foreach(){
                  css_test3_n[] =  clamp (css_test3_n[], 0., 1.);
                }
                vof2dist(css_test3_n,phi);
                boundary ({phi});
                fractions (phi, css_test3_n, fss_test3_n);
                              foreach(){
                                  css_test3[] = 1.0 - css_test3_n[];
                                  fs_solid[] = css_test3[];
                                }
                                foreach_face(){
                                    fss_test3.x[] = 1.0 - fss_test3_n.x[];
                                }


                  foreach(){
                    cs[] = css_test3_n[];
                  }
                  foreach_face(){
                    fs.x[] = fss_test3_n.x[];
                  }
                  boundary({cs,fs});
               #if EMBED
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});
              #endif
              if(!(surface_heat && surface_heat_restart)){
                  foreach(){
                        T[]=Tsat00;
                    }
              }

                boundary({T});
                fs_solid.prolongation = fs_solid.refine = fraction_refine;
                css_test3.refine = embed_fraction_refine_s;
                css_test3.prolongation = fraction_refine; //embed_fraction_refine_s; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3.x.prolongation = embed_face_fraction_refine_s_x;
                css_test3_n.refine = embed_fraction_refine_s_n;
                css_test3_n.prolongation = fraction_refine; //embed_fraction_refine_s_n; //fraction_refine; // this is in fraction
                foreach_dimension()
                  fss_test3_n.x.prolongation = embed_face_fraction_refine_s_n_x;
       foreach(){
                  if(css_test3_n[]<=0){
                      css_test[]=0;
                      css_test2[]=0;
                  }else if(css_test3_n[]>=1){
                      css_test[]=ff[];
                      css_test2[]=1.0-ff[];
                  }else{
                      //?????????????????????
                      css_test[]=ff[]*css_test3_n[];
                      css_test2[]=(1.0-ff[])*css_test3_n[];
                  }
               }

      face vector fs_temp[];
       if(1==1){
               vertex scalar phi21[];
               vof2dist(css_test,phi21);
               scalar ff_temp[];
               fractions(phi21, ff_temp, fss_test);
               vof2dist(css_test2,phi21);
               fractions(phi21, ff_temp, fss_test2);
               foreach_face(){
                  fs_temp.x[] = fss_test.x[];
               }


       }else{
                vertex scalar phi21[];
                scalar ff_temp[];
                vof2dist(ff, phi21);
                fractions (phi21, ff_temp, fs_temp);

                foreach(){
                  css_test[] = 0.0;
                  css_test2[] = 0.0;
                }
                foreach_face(){
                  fss_test.x[] = 0.0;
                  fss_test2.x[]=0.0;
                }
                foreach(){
                  if(css_test3_n[]>0.0){
                    css_test[] = ff_temp[];
                    css_test2[] = 1.0 - ff_temp[];
                  }
                }
                foreach_face(){
                  if(fss_test3_n.x[]>0.0){
                      fss_test.x[] = fs_temp.x[];
                      fss_test2.x[] = 1.0 - fs_temp.x[];
                  }
                }
          }
                 css_test.refine = embed_fraction_refine_css_test;
                css_test.prolongation = embed_fraction_refine_css_test;//fraction_refine;
                foreach_dimension(){
                          fss_test.x.refine = fss_test.x.prolongation = embed_face_fraction_fss_test_refine_x;
                          //fss_test.x.dirty = false;
                }
                css_test2.refine = embed_fraction_refine_css_test2;
                css_test2.prolongation = embed_fraction_refine_css_test2;//fraction_refine;
                foreach_dimension()
                          fss_test2.x.prolongation = embed_face_fraction_fss_test2_refine_x;

                 //add 20230915
               get_css_fss_areaslg_triple_point();

                restriction({css_test,fss_test});
                restriction({css_test2,fss_test2});
                restriction({css_test3,fss_test3});
                restriction({css_test3_n,fss_test3_n});
  #if EMBED
    restriction ({cm_css_test, fss_test});
    cm_css_test_update (cm_css_test, css_test, fss_test);
    fm_fss_test_update (fm_fss_test, css_test, fss_test);
    restriction ({cm_css_test, fss_test});
  #endif
  #if EMBED
    restriction ({cm_css_test2, fss_test2});
    cm_css_test2_update (cm_css_test2, css_test2, fss_test2);
    fm_fss_test2_update (fm_fss_test2, css_test2, fss_test2);
    restriction ({cm_css_test2, fss_test2});
  #endif
  #if EMBED
    restriction ({cm_css_test3, fss_test3});
    cm_css_test3_update (cm_css_test3, css_test3, fss_test3);
    fm_fss_test3_update (fm_fss_test3, css_test3, fss_test3);
    restriction ({cm_css_test3, fss_test3});
  #endif


                heights(ff,hhh);
              //  heights(ff_oppo,hhh_ff_oppo);
                 printf("event init 482,pid=%d\n",pid());
                      foreach(){
                        source_pc[] = 0.0;
                        masstr[]=0.0;
                        vtr[] = 0.0;
                      }

          ///////////////////////////////////////////////////////
        //  T.refine  = refine_bilinear;
        //  T.restriction = restriction_volume_average;//restriction_average;//restriction_volume_average;
        //  T.dirty = true;
          globali = 0;
            foreach_face(){
              ulf.x[] =0.0;
              usf.x[] =0.0;
              uf.x[] = 0.0;

              ugf.x[] =0.0;
              usfg.x[] =0.0;

            }
            foreach(){
              p[]=0.;
              pf[]=0.0;
              foreach_dimension(){
                u.x[]=0.0;
                g.x[]=0.0;
              }
            }
            boundary({u.x,u.y});
  }
}


     level_interface = maxl;



     ///reset cs after adapt //20230614
             vertex scalar phi_temp[];
              //  foreach_vertex(){
              //     // phi_temp[] =  -x;//+Delta;
              //     phi_temp[] =  -(x-tune2_value);
              //     // phi_temp[] =  -x + L0_pysical/(1<<maxl)/2.0;
              //  }
              solid_phi(tune2_value, phi_temp);
              foreach_vertex(){
                phi_temp[] = -phi_temp[];
              }
               fractions (phi_temp, css_test3, fss_test3);
              foreach(){
                    fs_solid[] = css_test3[];
                    css_test3_n[] = 1.0 - css_test3[];
               }
               foreach_face(){
                      fss_test3_n.x[] = 1.0 - fss_test3.x[];
              }
              boundary({css_test3,fss_test3});
              boundary({fs_solid});
              boundary({css_test3_n,fss_test3_n});
              foreach(){
                    cs[] = css_test3_n[];
              }
              foreach_face(){
                    fs.x[] = fss_test3_n.x[];
              }
              boundary({cs,fs});
               #if EMBED
                  restriction ({cs, fs});
                  cm_update (cm, cs, fs);
                  fm_update (fm, cs, fs);
                  restriction ({cm, fm});
              #endif


     if(case_number==3 && (restart_Tsat || restartsymbol)){
           double Delta_min = L0_pysical/(1<<maxl);
           foreach(){
                if(cs[]<=0 && fabs(ff[])>1e-6 && x<-3*Delta_min) {
                      // ff[] = 0.0;
                      ff[] = 1.0;
                }
           }

           foreach(){
              ff[] = clamp (ff[], 0., 1.);
              if(fabs(ff[])<1e-12){
                    ff[]=0.0;
              }else if(fabs(1-ff[])<1e-12){
                    ff[] = 1.0;
              }
          }
          foreach(){
                if(fabs(cs[])>=0){
                     if(ff[]<1e-10){
                        ff[]=0;
                     }else if(ff[]>1-1e-10){
                        ff[]=1.0;
                     }
                }
           }
           foreach(){
                if((level != level_interface) && (ff[]<1.0 && ff[]>0.0)){
                        ff[] = (ff[]>0.5?1.0:0.0);
                }
           }
          //  foreach(){
          //       if(ff[]>0 && ff[]<1){
          //             int num4=0;
          //             double num2=0.0;
          //            // double xx=x, yy =y;
          //            int iii=point.i;
          //            int jjj=point.j;
          //             foreach_neighbor(1){

          //                 if((ff[]>0.0) && (ff[]<1.0) && !(point.i==iii && point.j==jjj)){
          //                     num4=num4+1;
          //                     num2=num2+ff[];
          //                 }
          //             }
          //             if(num4==0){
          //                 ff[]=((num2/(9.0)<0.5)?0.0:1.0);
          //             }
          //       }
          //  }


           //mesu3-15 only valid for cs adapt: eliminate cs inside bubble
            //mesu3-15 only valid for cs adapt: eliminate cs inside bubble
          //  foreach(){
          //       if(cs[]>0 && cs[]<1 && ff[]>0 && ff[]<1){
          //             int ii = point.i;
          //             int jj = point.j;
          //             int num=0;
          //             int num0=0;
          //             int num1=0;
          //             foreach_neighbor(2){
          //                 bool judge = (ff[]>0 && ff[]<1);
          //                 if(level==level_interface){
          //                     if(cs[]>=1 && judge){
          //                         num++;
          //                     }
          //                     if(cs[]>=1 && ff[]>=1 && (fabs(point.i-ii)<=1) && (fabs(point.j-jj)<=1)){
          //                         num1++;
          //                     }
          //                     if(cs[]>=1 && ff[]<=0 && (fabs(point.i-ii)<=1) && (fabs(point.j-jj)<=1)){
          //                         num0++;
          //                     }
          //                 }
          //             }
          //             if(num==0 && (intersect_true[]==1)){
          //                 ff[]=((num1>num0)?1:0);
          //             }
          //       }
          //  }

          // foreach(){
          //   if(cs[]>0.0 && cs[]<1.0){
          //       if(css_test[]<1e-3){
          //           css_test[]=0.0;
          //           css_test2[]=cs[]-css_test[];
          //           ff[]=0.0;
          //       }
          //   }
          // }

          foreach(){
            if(ff[]<1.0 && ff[]>0.0 && cs[]>0){
                if(T[]<Tsat00){
                  T[]=Tsat00;
                }
            }
          }
          boundary({T});

           //-new
          //  foreach(){
          //       if(cs[]>0 && cs[]<1 && ff[]>0 && ff[]<1){
          //           double value = css_test[] + css_test2[];
          //           bool flag = false;
          //            if(value>0){
          //               flag = (css_test[]/value)<0.1;
          //            }
          //           if(flag){
          //             int ii = point.i;
          //             int jj = point.j;
          //             int num=0;
          //             int num0=0;
          //             int num1=0;
          //             foreach_neighbor(2){
          //                 if(!(point.i==ii && point.j==jj)){
          //                     if(cs[]>0 && cs[]<1 && ff[]<1.0 && ff[]>0.0){
          //                         num++;
          //                         double value2 = css_test[] + css_test2[];
          //                         bool flag2 = false;
          //                         if(value2>0){
          //                             flag2 = (css_test[]/value)<0.1;
          //                         }
          //                         if(flag2){
          //                             num0++;
          //                         }

          //                     }
          //                 }
          //             }
          //             if(num>0 && (num==num0)){
          //                 ff[] = 0;
          //                 css_test[] = 0;
          //                 css_test2[] = 1.0;
          //             }
                      
          //         }//if(css_test[]/(css_test[]+css_test2[])<0.1)
          //       }
          //  }
     }

}
#endif




event logfile (i++,last){
    double sb = 0.;
    double sb2 = 0.0;
    double area=0.;
  int number2=0;
  //double radius2D = bubble_radius;
  double radius3D ;//= bubble_radius;
  double dd=-1;
  foreach(reduction(+:sb) reduction(+:number2) reduction(+:area) reduction(max:dd)) {
    double dv = (1. - ff[])*dv();
    sb += dv; //dv = rdrdz.   dv= rdrdz*2pi
    number2++;

    if (ff[] > 1e-6 && ff[] < 1. - 1e-6) {
      coord n = interface_normal (point, ff), p;
      double alpha = plane_alpha (ff[], n);
      // area of the bubble interface
      area += y*pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    }

    if(fabs(x - tracex)<1e-6 && fabs(y-tracey)<1e-6){
        dd=T[];
    }
  }
  //double vv=0.0;
  sb2 = 2.0*M_PI*sb*2.0;
  radius3D = pow(3.0/(4.0*M_PI)*sb2,1.0/3.0);
  //radius3D = pow(3.0/(4.0*M_PI)*sb,1.0/3.0);
  // double dd=0;
  if(pid()==0){
       char name82[80];
       sprintf(name82,"volume.dat");
       FILE * fp82 = fopen(name82,"a");
      //radius2D = pow(sb/M_PI,0.5);
      fprintf (fp82,
        "%.8f %.8f %.8f %.8f %d %g %g\n",
        t, sb, sb2, radius3D, number2, area, dd);

        fclose(fp82);
  }


  // if(!restore (file = "restart_Tsat0") && case_number==3){
 if(!restartsymbol){
      if((!restart_Tsat) && case_number==3){
          int flag_stop = 0;
          double Delta_min = L0_pysical/(1<<maxl);
          Tmax = 0.0;
          foreach(reduction(max:flag_stop) reduction(max:Tmax)){
              // if(y<Delta_min && x<=Delta_min && x>0){
              if(topo_mask_s[]>=0){
                // if(T[]>=T_inf){
                if(aiml[]>=T_inf){
                    flag_stop = 1;
                    Tmax = T[];
                    printf("flag_stop=1 t=end=%g origin_T=%g\n",t,Tmax);
                }
              }
          }
          if(flag_stop==1){
                    char dumpname[80];
                    foreach_dimension(){
                      hhh.x.nodump = true;
                      //  hhh_ff_oppo.x.nodump = true;
                    }
                    sprintf(dumpname,"dump-final-i%d-t%g",i,t);
                    dump(dumpname);


                // return 1;
                exit(1);
          }
      }
 }

if(dump_each_event && (i%dump_each_event_interval==0)){
        char dumpname[80];
    // topo_mask_s.nodump = true;
      foreach_dimension(){
        hhh.x.nodump = true;
        //  hhh_ff_oppo.x.nodump = true;
      }
      sprintf(dumpname,"outfacets/dumpafteradapt-i%d-t%g",i,t);
      dump(dumpname);


}
}

// event logfile2(i++)
// {
//     char names[80];
// 	sprintf(names, "interface%d", pid());
// 	FILE * fp = fopen (names, "w");
// 	output_facets (ff,fp);
// 	fclose(fp);

//   MPI_Barrier(MPI_COMM_WORLD);
//   if(pid()==0){
//     char command[80];
//     sprintf(command, "LC_ALL=C  cat interface* > outfacets/interface_%d_%g.dat",i,t);
//     system(command);
//   }
// }
//event pictures (t = end)
//event pictures (t += out_interval)
event end (t = end ) {
  // if(case_number==3 && !restore (file = "restart_Tsat0")){
  if(case_number==3 && (!restart_Tsat)){
   printf("t=end=%g\n",t);
  }
}

event stop (t = tend) {
  return 1;
}

