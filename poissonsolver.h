extern scalar flux_lg_l,flux_lg_g;
extern scalar phase0Tgrad,phase1Tgrad;
extern scalar resl,resg,ress;
// extern scalar cm_css_test,cm_css_test2,cm_css_test3;
// extern face vector fm_fss_test,fm_fss_test2,fm_fss_test3;
extern scalar ff,ff_oppo;

extern scalar T_modl,T_modg,aiml,aimg;
extern int maxl;

void get_flux_lg(scalar flux_lg_l1, scalar flux_lg_g1, scalar phase0Tgrad, scalar phase1Tgrad){
     foreach(){
        flux_lg_l1[]=0.0;
        flux_lg_g1[]=0.0;
        if(level==level_interface){
            if(css_test3[]<1.0){
                    if(css_test[]>0 && css_test[]<1.0){
                            // embed dirichlet flux
                            double c1, e=0.0;
                            double c1_2,e_2=0.0;
                            double arealg=0.0;
                        if(1==1){
                            double mua = 0., fa = 0.;
                            foreach_dimension() {
                            // mua += alpha.x[] + alpha.x[1];
                            fa  += fss_test.x[] + fss_test.x[1];
                            }
                            // coord nl = mycs (point, css_test);//interface_normal(point, ff);//interface_normal7(point,ff,hhh); //
                            // double alphal = plane_alpha (css_test[], nl);
                            coord nl = mycs (point, ff);//interface_normal(point, ff);//interface_normal7(point,ff,hhh); //
                            double alphal = plane_alpha (ff[], nl);
                            coord ppp;
                            arealg = plane_area_center (nl, alphal, &ppp);
                            if (metric_embed_factor)  
                                    arealg *= metric_embed_factor (point, ppp);
                            double gradl2 = phase1Tgrad[];// for fluid this is value for from water to out direction
                            mua = Tkl; //Tkl*max(y,1e-20);
                            fa =1.0;
                            // c1 =  - mua/(fa + SEPS)*gradl2*arealg/Delta;
                            c1 =  mua/(fa + SEPS)*gradl2*arealg/Delta;
                            e = 0.0;
                  //////////////////////////set flux limit for flux_lg_l            
                            double limit;
                             coord coord_pa;
                             coord_pa.x=x,coord_pa.y=y;//,mixed.z=z;
                             double merge_plus_self_volume=css_test[];
                            if(css_test[]<0.5){
                                int ii = merge_to_me_l_position.x[];
                                int jj = merge_to_me_l_position.y[];
                                if(!(ii==0 && jj==0)){
                                    merge_plus_self_volume += css_test[ii,jj];
                                }
                            }else{
                                
                                foreach_neighbor(1){
                                    coord coord_child;
                                    coord_child.x = x, coord_child.y = y;//, pure.z = z;
                                    int ii_real=round((coord_pa.x-coord_child.x)/Delta);
                                    int jj_real=round((coord_pa.y-coord_child.y)/Delta);
                                    int ii_aim=merge_to_me_l_position.x[];
                                    int jj_aim=merge_to_me_l_position.y[];
                                    if((!(ii_aim==0 && jj_aim==0)) && (ii_real==ii_aim && jj_real==jj_aim)){
                                        merge_plus_self_volume+=css_test[];
                                    }
                                }
                            }
                            limit = Trhol*Tcpl*(Tsat00-Tl[])*merge_plus_self_volume*max(1e-20,y)/dt;
                            // limit = mua/(fa + SEPS)*(Tsat0-Tl[])/Delta/Delta*merge_plus_self_volume;
                            if(sign(c1)==sign(limit) && true){
                                flux_lg_l1[] = fabs(c1)<fabs(limit)?c1:limit;
                            }else{
                                flux_lg_l1[] = c1;
                            }
                        }

                        double arealg2=0.0;
                        if(1==1){
                            double mua = 0., fa = 0.;
                            foreach_dimension() {
                                //mua += alpha.x[] + alpha.x[1];
                                fa  += fss_test2.x[] + fss_test2.x[1];
                            }
                            //mua= k0;//Tkg;
                            // coord ng = interface_normal(point, css_test2);//ff_oppo);// interface_normal7(point,ff_oppo,hhh_ff_oppo); // interface_normal7(point,css_test,hhh); //interface_normal7
                            coord ng = mycs(point, ff_oppo);//ff_oppo);// interface_normal7(point,ff_oppo,hhh_ff_oppo); // 
                            coord ppp2;
                            // double alphag = plane_alpha (css_test2[], ng);
                            double alphag = plane_alpha (ff_oppo[], ng);
                            arealg2 = plane_area_center (ng, alphag, &ppp2);
                            if (metric_embed_factor)  
                                    arealg2 *= metric_embed_factor (point, ppp2);
    //direction of phase1Tg and phase0Tg is from fluid to out; for liquid it is out, for gas it is in; so add one -, finally it is +
                            double gradg2 = -phase0Tgrad[]; 
                            mua = Tkg; //Tkg*max(y,1e-20);
                            fa= 1.0;
                            // c1_2 =  - mua/(fa + SEPS)*gradg2*arealg2/Delta;
                            c1_2 =  mua/(fa + SEPS)*gradg2*arealg2/Delta;
                        e_2 = 0.0;
        //////////////////////////set flux limit for flux_lg_g                
                        double limit;
                             coord coord_pa;
                             coord_pa.x=x,coord_pa.y=y;//,mixed.z=z;
                              double merge_plus_self_volume=css_test2[];
                            if(css_test2[]<0.5){
                                int ii = merge_to_me_g_position.x[];
                                int jj = merge_to_me_g_position.y[];
                                if(!(ii==0 && jj==0)){
                                    merge_plus_self_volume += css_test2[ii,jj];
                                }
                            }else{
                                foreach_neighbor(1){
                                    coord coord_child;
                                    coord_child.x = x, coord_child.y = y;//, pure.z = z;
                                    int ii_real=round((coord_pa.x-coord_child.x)/Delta);
                                    int jj_real=round((coord_pa.y-coord_child.y)/Delta);
                                    int ii_aim=merge_to_me_g_position.x[];
                                    int jj_aim=merge_to_me_g_position.y[];
                                    if((!(ii_aim==0 && jj_aim==0)) && (ii_real==ii_aim && jj_real==jj_aim)){
                                        merge_plus_self_volume+=css_test2[];
                                    }
                                }
                            }
                            // limit = mua/(fa + SEPS)*(Tsat0-Tg[])/Delta/Delta*merge_plus_self_volume;
                             limit = Trhog*Tcpg*(Tsat00-Tg[])*merge_plus_self_volume*max(1e-20,y)/dt;
                             if(sign(c1_2)==sign(limit) && true){
                                flux_lg_g1[] = fabs(c1_2)<fabs(limit)?c1_2:limit;
                            }else{
                                flux_lg_g1[] = c1_2;
                            }
                        }

                }
        }
        }
    }
}

// Function to return the maximum value among three integers
double max_function(double x, double y, double z) {
    if (x > y) {
        if (x > z) {
            return x;
        } else {
            return z;
        }
    } else {
        if (y > z) {
            return y;
        } else {
            return z;
        }
    }
}

void poisson_solver(scalar poisson_s, double percent_s,double percent_l,double percent_g){

int maxitt=2000,beta1=1.2;
int itt;
double T_tolerance = 1e-6;
// if(maxl==14){
//     T_tolerance = 1e-5;
// }
if(maxl==11){
    T_tolerance = 1e-8;
}
foreach(){
    ff_oppo[] = 1.0 - ff[];
}
get_solidsurfaceT(); //get aiml aimg T_modl T_modg flux_l flux_g
Tgrad_leon(ff,phase0Tgrad,phase1Tgrad,aiml,aimg);
get_flux_lg(flux_lg_l, flux_lg_g, phase0Tgrad, phase1Tgrad);


scalar flux_s_6[],flux_l_6[],flux_g_6[];
foreach(){
    flux_s_6[] = 0.0;
    flux_l_6[] = 0.0;
    flux_g_6[] = 0.0;
}

foreach(){
    if(T_modl[]==6){
        // double limit_s = Trhos*Tcps*(aiml[]-Ts[])*css_test3[]*max(1e-20,y)/dt;
        double limit_s = Trhos*Tcps*(aiml[]-Ts[])*1.0*max(1e-20,y)/dt;
        //if css_test3<0.5,Ts is actrually is the temperature of its merged cell's temperature
        double value_s=0.0;
        value_s = (aiml[] - Ts[])/Delta*areasl[]/Delta ;
        flux_s_6[] += (fabs(value_s)<fabs(limit_s)?value_s:limit_s) ;//distance use Delta for stability

        //  double limit_l = Trhol*Tcpl*(aiml[]-Tl[])*css_test[]*max(1e-20,y)/dt;
         double limit_l = Trhol*Tcpl*(aiml[]-Tl[])*1.0*max(1e-20,y)/dt;
        //if css_test<0.5,Ts is actrually is the temperature of its merged cell's temperature
        double value_l=0.0;
        value_l = (aiml[] - Tl[])/Delta*areasl[]/Delta;
        flux_l_6[] +=  (fabs(value_l)<fabs(limit_l)?value_l:limit_l);//distance use Delta for stability
    }
    if(T_modg[]==6){
        // double limit_s = Trhos*Tcps*(aimg[]-Ts[])*css_test3[]*max(1e-20,y)/dt;
        double limit_s = Trhos*Tcps*(aimg[]-Ts[])*1.0*max(1e-20,y)/dt;
        //if css_test3<0.5,Ts is actrually is the temperature of its merged cell's temperature
        double value_s=0.0;
        value_s = (aimg[] - Ts[])/Delta*areasg[]/Delta ;
        flux_s_6[] +=  (fabs(value_s)<fabs(limit_s)?value_s:limit_s) ;//distance use Delta for stability

        //  double limit_g = Trhog*Tcpg*(aimg[]-Tg[])*css_test2[]*max(1e-20,y)/dt;
         double limit_g = Trhog*Tcpg*(aimg[]-Tg[])*1.0*max(1e-20,y)/dt;
        //if css_test2<0.5,Ts is actrually is the temperature of its merged cell's temperature
        double value_g=0.0;
        value_g =(aimg[] - Tg[])/Delta*areasg[]/Delta;
        flux_g_6[] +=  (fabs(value_g)<fabs(limit_g)?value_g:limit_g);//distance use Delta for stability
    }

}
///////////////////
//let the units of flux_l, flux_g, flux_lg_l, flux_lg_g, flux_s_6 be the same
foreach(){
    if(level==level_interface){ //flux accross solid surface
        if(css_test3[]>0.0 && css_test3[]<1.0){
                        double flux_i=0.0;
                        // flux_i = flux_l[]/Delta*areasl[]/Delta + flux_g[]/Delta*areasg[]/Delta + percent_s*poisson_s[] + flux_s_6[];
                        flux_l[] = flux_l[]/Delta*areasl[]/Delta;
                        flux_g[] = flux_g[]/Delta*areasg[]/Delta;
                        flux_s_6[] = 0; //T_modl or T_modg  !=6
                        flux_l_6[] = 0;
                        flux_g_6[] = 0;
        }else{
            flux_l[] = 0;
            flux_g[] = 0;
            flux_s_6[] = 0;
            flux_l_6[] = 0;
            flux_g_6[] = 0;
        }

        if(cs[]>0.0){
            if(ff[]<1.0 && ff[]>0.0){
                flux_lg_l[] = flux_lg_l[];
                flux_lg_g[] = flux_lg_g[];
            }
        }else{
            flux_lg_l[] = 0;
            flux_lg_g[] = 0;
        }
        poisson_s[] = poisson_s[];
    }else{
        flux_l[] = 0;
        flux_g[] = 0;
        flux_lg_l[] = 0;
        flux_lg_g[] = 0;
        poisson_s[] = 0;
        flux_s_6[] = 0;
        flux_l_6[] = 0;
        flux_g_6[] = 0;
    }
}
flux_l.restriction = restriction_flux_sum;
flux_g.restriction = restriction_flux_sum;
flux_s_6.restriction = restriction_flux_sum;
flux_l_6.restriction = restriction_flux_sum;
flux_g_6.restriction = restriction_flux_sum;
flux_lg_l.restriction = restriction_flux_sum;
flux_lg_g.restriction = restriction_flux_sum;
poisson_s.restriction = restriction_flux_sum;

restriction({flux_l,flux_g,flux_lg_l,flux_lg_g,poisson_s});
restriction({flux_s_6,flux_l_6,flux_g_6});

scalar Ts_old[],Tl_old[],Tg_old[];
///flux_l is the flux across solid and fluid surface
///flux_g is the flux across solid and gas surface
//flux_lg_l is the flux across liquid and gas for liquid 
//flux_lg_g is the flux across liquid and gas for gas
foreach(){
    Ts_old[] = Ts[];
    Tl_old[] = Tl[];
    Tg_old[] = Tg[];
}
Ts_old.ff6 = css_test3;
Ts_old.restriction = restriction_zero;
Tl_old.ff6 = css_test;
Tl_old.restriction = restriction_zero;
Tg_old.ff6 = css_test2;
Tg_old.restriction = restriction_zero;

Ts_old.refine = Ts_old.prolongation = refine_embed_linear_css_test3;
Tl_old.refine = Tl_old.prolongation = refine_embed_linear_css_test;
Tg_old.refine = Tg_old.prolongation = refine_embed_linear_css_test2;
boundary({Ts_old,Tl_old,Tg_old});
restriction({Ts_old,Tl_old,Tg_old});
double above_value_merge=above_value;
scalar caculatecell_s[],caculatecell_l[],caculatecell_g[];
foreach(){
    caculatecell_s[]=0;
    caculatecell_l[]=0;
    caculatecell_g[]=0;
    if(css_test3[]>0){
        if(level==level_interface){
            if(css_test3[]>=above_value_merge){
                caculatecell_s[]=1;
            }
        }else{
            if(css_test3[]>0){
                caculatecell_s[]=1;
            }
        }
    }
    if(css_test3[]<1.0){
        if(level==level_interface){
            if(css_test[]>=above_value_merge){
                caculatecell_l[]=1;
            }
        }else{
            if(css_test3[]<1 && css_test[]>0){
                caculatecell_l[]=1;
            }
        }
    }
    if(css_test3[]<1.0){
        if(level==level_interface){
            if(css_test2[]>=above_value_merge){
                caculatecell_g[]=1;
            }
        }else{
            if(css_test3[]<1.0 && css_test2[]>0){
                caculatecell_g[]=1;
            }
        }
    }
}
for(itt=1;itt<=maxitt;itt++){
        
    foreach(){
        // double cm_s_value = css_test3[]*max(y,1e-20)
        // double fs_s_l_value = fss_test3[-1,0]*max(y,1e-20);
        // double fs_s_r_value = fss_test3[1,0]*max(y,1e-20);
        // double fs_s_b_value = fss_test3[0,-1]*max(y-Delta/2.0,1e-20);
        // double fs_s_t_value = fss_test3[0,1]*max(y+Delta/2.0,1e-20);

        if(css_test3[]>0){
            ///////////////////////////////solid
            bool flag=false;
            double total= percent_s*poisson_s[];
            double A3=Trhos*Tcps/dt*cm_css_test3[];
            double A4=Trhos*Tcps/dt*cm_css_test3[];
            if(level==level_interface){
                if(css_test3[]>=above_value_merge){
                    flag = true;
                    if(css_test3[]>0.0 && css_test3[]<1.0){ //flux accross solid surface
                        double flux_i=0.0;
                        // flux_i = flux_l[]/Delta*areasl[]/Delta + flux_g[]/Delta*areasg[]/Delta + percent_s*poisson_s[] + flux_s_6[];
                        flux_i = flux_l[] + flux_g[]  + flux_s_6[];
                        total += flux_i;
                    }
                        //left
                        // if(css_test3[-1,0]>=0.5){
                        if(css_test3[-1,0]>=above_value_merge){
                            double Al;
                            Al = Tks*fm_fss_test3.x[]*1.0/(Delta*Delta);
                            total +=   Al*Ts[-1,0];
                            A3 +=    Al; 
                        }else if(css_test3[-1,0]>0.0){
                            //merge judge
                            if(fabs(merge_to_me_s_c[-1,0])>1e-12){
                                int ii = merge_to_me_s_position.x[-1,0];
                                int jj = merge_to_me_s_position.y[-1,0];
                                if((ii==1) && (jj==0)){ //now merge time
                                    //flux across surface
                                    double flux_i=0.0;
                                    flux_i = flux_l[-1,0] + flux_g[-1,0] + percent_s*poisson_s[-1,0] + flux_s_6[-1,0];
                                    total += flux_i;
                                    A3 += Trhos*Tcps/dt*cm_css_test3[-1,0];
                                    A4 += Trhos*Tcps/dt*cm_css_test3[-1,0];
                                }
                            }
                        } // css_test3[-1,0]

                        //right
                       // if(css_test3[1,0]>=0.5){
                         if(css_test3[1,0]>=above_value_merge){
                            double Ar;
                            Ar = Tks*fm_fss_test3.x[1,0]*1.0/(Delta*Delta);
                            total +=   Ar*Ts[1,0];
                            A3 +=    Ar; 
                        }else if(css_test3[1,0]>0.0){
                            //merge judge
                            if(fabs(merge_to_me_s_c[1,0])>1e-12){
                                int ii = merge_to_me_s_position.x[1,0];
                                int jj = merge_to_me_s_position.y[1,0];
                                if((ii==-1) && (jj==0)){
                                    //flux across surface
                                    double flux_i=0.0;
                                    flux_i = flux_l[1,0] + flux_g[1,0] + percent_s*poisson_s[1,0]+flux_s_6[1,0];
                                    total += flux_i;
                                    A3 += Trhos*Tcps/dt*cm_css_test3[1,0];
                                    A4 += Trhos*Tcps/dt*cm_css_test3[1,0];
                                }
                            }
                        }// css_test3[1,0]

                         //bottom
                        // if(css_test3[0,-1]>=0.5){
                        if(css_test3[0,-1]>=above_value_merge){
                            double Al;
                            Al = Tks*fm_fss_test3.y[]*1.0/(Delta*Delta);
                            total +=   Al*Ts[0,-1];
                            A3 +=    Al; 
                        }else if(css_test3[0,-1]>0.0){
                            //merge judge
                            if(fabs(merge_to_me_s_c[0,-1])>1e-12){
                                int ii = merge_to_me_s_position.x[0,-1];
                                int jj = merge_to_me_s_position.y[0,-1];
                                if((ii==0) && (jj==1)){
                                    //flux across surface
                                    double flux_i=0.0;
                                    flux_i = flux_l[0,-1] + flux_g[0,-1] + percent_s*poisson_s[0,-1] + flux_s_6[0,-1];
                                    total += flux_i;
                                    A3 += Trhos*Tcps/dt*cm_css_test3[0,-1];
                                    A4 += Trhos*Tcps/dt*cm_css_test3[0,-1];
                                }
                            }
                        } // css_test3[0,-1]

                         //top
                        // if(css_test3[0,1]>=0.5){
                        if(css_test3[0,1]>=above_value_merge){
                            double Ar;
                            Ar = Tks*fm_fss_test3.y[0,1]*1.0/(Delta*Delta);
                            total +=   Ar*Ts[0,1];
                            A3 +=    Ar; 
                        }else if(css_test3[0,1]>0.0){
                            //merge judge
                            if(fabs(merge_to_me_s_c[0,1])>1e-12){
                                int ii = merge_to_me_s_position.x[0,1];
                                int jj = merge_to_me_s_position.y[0,1];
                                if((ii==0) && (jj==-1)){
                                    //flux across surface
                                    double flux_i=0.0;
                                    flux_i = flux_l[0,1] + flux_g[0,1] + percent_s*poisson_s[0,1]+ flux_s_6[0,1];
                                    total += flux_i;
                                    A3 += Trhos*Tcps/dt*cm_css_test3[0,1];
                                    A4 += Trhos*Tcps/dt*cm_css_test3[0,1];
                                }
                            }
                        } // css_test3[-1,0]
                }
            }else{ //level<level_interface
                // if(css_test3[]>=0.5){
                if(css_test3[]>0){
                    flag=true;
                    if(css_test3[]>0.0 && css_test3[]<1.0){ //flux accross solid surface
                        double flux_i=0.0;
                        // flux_i = flux_l[]/Delta*areasl[]/Delta + flux_g[]/Delta*areasg[]/Delta + percent_s*poisson_s[] + flux_s_6[];
                        flux_i = flux_l[] + flux_g[]  + flux_s_6[];
                        total += flux_i;
                    }
                    //4faces
                    foreach_dimension(){
                        if(css_test3[-1]>0.0){
                           double Al;
                            Al = Tks*fm_fss_test3.x[]*1.0/(Delta*Delta);
                            total +=   Al*Ts[-1];
                            A3 +=    Al; 
                        }

                        if(css_test3[1]>0.0){
                            double Ar;
                            Ar = Tks*fm_fss_test3.x[1]*1.0/(Delta*Delta);
                            total +=   Ar*Ts[1];
                            A3 +=    Ar; 
                        }
                    }
                }//if
            }
            if(flag){
                Ts[] = (total+A4*Ts_old[])*beta1/A3 + (1.0-beta1)*Ts[];
            }
        } 



        //////////////
        ///////////////////////////////fluid
        //consider flux_l and flux across solid
      //    if(css_test3[]<0.5 && css_test[]>=0.5){
        // if(css_test[]>=above_value_merge){
        if(css_test3[]<1.0){
            bool flag=false;
            double total_l= percent_l*poisson_s[];
            double A3_l=Trhol*Tcpl/dt*cm_css_test[];
            double A4_l=Trhol*Tcpl/dt*cm_css_test[];       
             if(level==level_interface){
                if(css_test[]>=above_value_merge){
                    flag=true;
                //if there is solid interface in the cell
                if(css_test3[]>0.0 && css_test3[]<1.0){
                    double flux_i=0.0;
                    // flux_i = (-flux_l[])/Delta*areasl[]/Delta + percent_l*poisson_s[] + flux_l_6[];
                    flux_i = (-flux_l[])  + flux_l_6[];
                    total_l += flux_i;
                }
                //if there is fluid interface in the cell
                if(css_test3[]<1.0){
                    if(css_test[]>0 && css_test[]<1.0){
                        // flux across the interface between gas and liquid
                        // - mua/(fa + SEPS)*gradl2*arealg/Delta;
                        total_l += flux_lg_l[];
                    }
                }
                    //left
                    if(css_test3[-1,0]<1.0){
                        // if(css_test3[-1,0]<0.5 && css_test[-1,0]>=0.5){
                        if(css_test[-1,0]>=above_value_merge){
                                    double Al;
                                    Al = Tkl*fm_fss_test.x[]*1.0/(Delta*Delta);
                                    total_l +=   Al*Tl[-1,0];
                                    A3_l +=    Al; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_l_c[-1,0])>1e-12){
                                 int ii = merge_to_me_l_position.x[-1,0];
                                 int jj = merge_to_me_l_position.y[-1,0];
                                 if((ii==1) && (jj==0)){
                                        // if(css_test3[-1,0]<1.0 && css_test3[-1,0]>0.5){ //0.5< <1.0
                                        if(css_test3[-1,0]<1.0 && css_test3[-1,0]>0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_l[-1,0])  + percent_l*poisson_s[-1,0] + flux_l_6[-1,0];
                                            total_l += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test[-1,0]>0.0 && css_test[-1,0]<0.5){
                                        if(css_test[-1,0]>0.0 && css_test[-1,0]<above_value_merge){
                                                //lg interface
                                                total_l += flux_lg_l[-1,0];
                                        }
                                        A3_l += Trhol*Tcpl/dt*cm_css_test[-1,0];
                                        A4_l += Trhol*Tcpl/dt*cm_css_test[-1,0];
                                 }
                            }
                        }//merge
                    }

                   //right
                   if(css_test3[1,0]<1.0){
                        // if(css_test3[1,0]<0.5 && css_test[1,0]>=0.5){
                        if(css_test[1,0]>=above_value_merge){
                                    double Ar;
                                    Ar = Tkl*fm_fss_test.x[1,0]*1.0/(Delta*Delta);
                                    total_l +=   Ar*Tl[1,0];
                                    A3_l +=    Ar; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_l_c[1,0])>1e-12){
                                 int ii = merge_to_me_l_position.x[1,0];
                                 int jj = merge_to_me_l_position.y[1,0];
                                 if((ii==-1) && (jj==0)){
                                        // if(css_test3[1,0]<1.0 && css_test3[1,0]>0.5){ //0.5< <1.0
                                        if(css_test3[1,0]<1.0 && css_test3[1,0]>0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_l[1,0])  + percent_l*poisson_s[1,0]+ flux_l_6[1,0];
                                            total_l += flux_i;
                                        }
                                        //merge_energy_from fluid
                                       // if(css_test[1,0]>0.0 && css_test[1,0]<0.5){
                                        if(css_test[1,0]>0.0 && css_test[1,0]<above_value_merge){
                                                //merge judge
                                                total_l += flux_lg_l[1,0];
                                        }
                                        A3_l += Trhol*Tcpl/dt*cm_css_test[1,0];
                                        A4_l += Trhol*Tcpl/dt*cm_css_test[1,0];
                                 }
                            }
                        }//merge
                    }

                     //bottom
                   if(css_test3[0,-1]<1.0){
                        // if(css_test3[0,-1]<0.5 && css_test[0,-1]>=0.5){
                            if(css_test[0,-1]>=above_value_merge){
                                    double Al;
                                    Al = Tkl*fm_fss_test.y[]*1.0/(Delta*Delta);
                                    total_l +=   Al*Tl[0,-1];
                                    A3_l +=    Al; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_l_c[0,-1])>1e-12){
                                 int ii = merge_to_me_l_position.x[0,-1];
                                 int jj = merge_to_me_l_position.y[0,-1];
                                 if((ii==0) && (jj==1)){
                                        // if(css_test3[0,-1]<1.0 && css_test3[0,-1]>0.5){ //0.5< <1.0
                                        if(css_test3[0,-1]<1.0 && css_test3[0,-1]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_l[0,-1]) + percent_l*poisson_s[0,-1]+ flux_l_6[0,-1];
                                            total_l += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test[0,-1]>0.0 && css_test[0,-1]<0.5){
                                        if(css_test[0,-1]>0.0 && css_test[0,-1]<above_value_merge){
                                                //merge judge
                                                total_l += flux_lg_l[0,-1];
                                        }
                                        A3_l += Trhol*Tcpl/dt*cm_css_test[0,-1];
                                        A4_l += Trhol*Tcpl/dt*cm_css_test[0,-1];
                                 }
                            }
                        }//merge
                    }

                   //top
                   if(css_test3[0,1]<1.0){
                        // if(css_test3[0,1]<0.5 && css_test[0,1]>=0.5){
                       if(css_test[0,1]>=above_value_merge){ 
                                    double Ar;
                                    Ar = Tkl*fm_fss_test.y[0,1]*1.0/(Delta*Delta);
                                    total_l +=   Ar*Tl[0,1];
                                    A3_l +=    Ar; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_l_c[0,1])>1e-12){
                                 int ii = merge_to_me_l_position.x[0,1];
                                 int jj = merge_to_me_l_position.y[0,1];
                                 if((ii==0) && (jj==-1)){
                                        // if(css_test3[0,1]<1.0 && css_test3[0,1]>0.5){ //0.5< <1.0
                                        if(css_test3[0,1]<1.0 && css_test3[0,1]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_l[0,1])  + percent_l*poisson_s[0,1]+ flux_l_6[0,1];
                                            total_l += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test[0,1]>0.0 && css_test[0,1]<0.5){
                                        if(css_test[0,1]>0.0 && css_test[0,1]<above_value_merge){
                                                //merge judge
                                                total_l += flux_lg_l[0,1];
                                        }
                                        A3_l += Trhol*Tcpl/dt*cm_css_test[0,1];
                                        A4_l += Trhol*Tcpl/dt*cm_css_test[0,1];
                                 }
                            }
                        }//merge
                    }
                }
            }else{ 
                //level<level_interface
               // if(css_test3[]<0.5 && css_test[]>=0.5){
                if(css_test3[]<1 && css_test[]>0){
                    flag = true;
                     if(css_test3[]>0.0 && css_test3[]<1.0){
                        double flux_i=0.0;
                        // flux_i = (-flux_l[])/Delta*areasl[]/Delta + percent_l*poisson_s[] + flux_l_6[];
                        flux_i = (-flux_l[])  + flux_l_6[];
                        total_l += flux_i;
                    }
                    //if there is fluid interface in the cell
                    if(css_test3[]<1.0){
                        if(css_test[]>0 && css_test[]<1.0){
                            // flux across the interface between gas and liquid
                            // - mua/(fa + SEPS)*gradl2*arealg/Delta;
                            total_l += flux_lg_l[];
                        }
                    }
                    foreach_dimension(){
                        // if(css_test3[-1]<0.5 && css_test[-1]>0.0){
                        if(css_test[-1]>0.0){
                           double Al;
                            Al = Tkl*fm_fss_test.x[]*1.0/(Delta*Delta);
                            total_l +=   Al*Tl[-1];
                            A3_l +=    Al; 
                        }

                        // if(css_test3[1]<0.5 && css_test[1]>0.0){
                        if(css_test[1]>0.0){
                            double Ar;
                            Ar = Tkl*fm_fss_test.x[1]*1.0/(Delta*Delta);
                            total_l +=   Ar*Tl[1];
                            A3_l +=    Ar; 
                        }
                    }
                }//if
            }
           if(flag){
                Tl[] = (total_l+A4_l*Tl_old[])*beta1/A3_l + (1.0-beta1)*Tl[]; 
            }
    }


     //////////////
        ///////////////////////////////gas
        //consider flux_l and flux across solid
       //    if(css_test3[]<0.5 && css_test2[]>=0.5){
         if(css_test3[]<1.0){
            bool flag=false;
            double total_g=percent_g*poisson_s[];
            double A3_g=Trhog*Tcpg/dt*cm_css_test2[];
            double A4_g=Trhog*Tcpg/dt*cm_css_test2[];       
            if(level==level_interface){
                if(css_test2[]>=above_value_merge){
                    //if there is solid interface in the cell
                    flag=true;
                //if there is solid interface in the cell
                    if(css_test3[]>0.0 && css_test3[]<1.0){
                        double flux_i=0.0;
                        // flux_i = (-flux_g[])/Delta*areasg[]/Delta  + percent_g*poisson_s[] + flux_g_6[];
                        flux_i = (-flux_g[])   + flux_g_6[];
                        total_g += flux_i;
                    }
                //if there is fluid interface in the cell
                if(css_test3[]<1.0){
                    if(css_test2[]>0 && css_test2[]<1.0){
                        // flux across the interface between gas and liquid
                        total_g += flux_lg_g[];
                    }
                }
                    //left
                    if(css_test3[-1,0]<1.0){
                        // if(css_test3[-1,0]<0.5 && css_test2[-1,0]>=0.5){
                        if(css_test2[-1,0]>=above_value_merge){
                                    double Al;
                                    Al = Tkg*fm_fss_test2.x[]*1.0/(Delta*Delta);
                                    total_g +=   Al*Tg[-1,0];
                                    A3_g +=    Al; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_g_c[-1,0])>1e-12){
                                 int ii = merge_to_me_g_position.x[-1,0];
                                 int jj = merge_to_me_g_position.y[-1,0];
                                 if((ii==1) && (jj==0)){
                                         // if(css_test3[-1,0]<1.0 && css_test3[-1,0]>0.5){ //0.5< <1.0
                                        if(css_test3[-1,0]<1.0 && css_test3[-1,0]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_g[-1,0]) + percent_g*poisson_s[-1,0]+ flux_g_6[-1,0];
                                            total_g += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test2[-1,0]>0.0 && css_test2[-1,0]<0.5){
                                        if(css_test2[-1,0]>0.0 && css_test2[-1,0]<above_value_merge){
                                                //lg interface
                                                total_g += flux_lg_g[-1,0];
                                        }
                                        A3_g += Trhog*Tcpg/dt*cm_css_test2[-1,0];
                                        A4_g += Trhog*Tcpg/dt*cm_css_test2[-1,0];
                                 }
                            }
                        }//merge
                    }

                   //right
                   if(css_test3[1,0]<1.0){
                        // if(css_test3[1,0]<0.5 && css_test2[1,0]>=0.5){
                        if(css_test2[1,0]>=above_value_merge){
                                    double Ar;
                                    Ar = Tkg*fm_fss_test2.x[1,0]*1.0/(Delta*Delta);
                                    total_g +=   Ar*Tg[1,0];
                                    A3_g +=    Ar; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_g_c[1,0])>1e-12){
                                 int ii = merge_to_me_g_position.x[1,0];
                                 int jj = merge_to_me_g_position.y[1,0];
                                 if((ii==-1) && (jj==0)){
                                        // if(css_test3[1,0]<1.0 && css_test3[1,0]>0.5){ //0.5< <1.0
                                        if(css_test3[1,0]<1.0 && css_test3[1,0]>0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_g[1,0])+ percent_g*poisson_s[1,0]+ flux_g_6[1,0];
                                            total_g += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test2[1,0]>0.0 && css_test2[1,0]<0.5){
                                        if(css_test2[1,0]>0.0 && css_test2[1,0]<above_value_merge){
                                                //merge judge
                                                total_g += flux_lg_g[1,0];
                                        }
                                         A3_g += Trhog*Tcpg/dt*cm_css_test2[1,0];
                                         A4_g += Trhog*Tcpg/dt*cm_css_test2[1,0];
                                 }
                            }
                        }//merge
                    }

                     //bottom
                   if(css_test3[0,-1]<1.0){
                        // if(css_test3[0,-1]<0.5 && css_test2[0,-1]>=0.5){
                        if(css_test2[0,-1]>=above_value_merge){
                                    double Al;
                                    Al = Tkg*fm_fss_test2.y[]*1.0/(Delta*Delta);
                                    total_g +=   Al*Tg[0,-1];
                                    A3_g +=    Al; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_g_c[0,-1])>1e-12){
                                 int ii = merge_to_me_g_position.x[0,-1];
                                 int jj = merge_to_me_g_position.y[0,-1];
                                 if((ii==0) && (jj==1)){
                                        // if(css_test3[0,-1]<1.0 && css_test3[0,-1]>0.5){ //0.5< <1.0
                                        if(css_test3[0,-1]<1.0 && css_test3[0,-1]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_g[0,-1])  + percent_g*poisson_s[0,-1]+ flux_g_6[0,-1];
                                            total_g += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test2[0,-1]>0.0 && css_test2[0,-1]<0.5){
                                        if(css_test2[0,-1]>0.0 && css_test2[0,-1]<above_value_merge){
                                                //merge judge
                                                total_g += flux_lg_g[0,-1];
                                        }
                                         A3_g += Trhog*Tcpg/dt*cm_css_test2[0,-1];
                                         A4_g += Trhog*Tcpg/dt*cm_css_test2[0,-1];
                                 }
                            }
                        }//merge
                    }

                   //top
                   if(css_test3[0,1]<1.0){
                        // if(css_test3[0,1]<0.5 && css_test2[0,1]>=0.5){
                        if(css_test2[0,1]>=above_value_merge){
                                    double Ar;
                                    Ar = Tkg*fm_fss_test2.y[0,1]*1.0/(Delta*Delta);
                                    total_g +=   Ar*Tg[0,1];
                                    A3_g +=    Ar; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_g_c[0,1])>1e-12){
                                 int ii = merge_to_me_g_position.x[0,1];
                                 int jj = merge_to_me_g_position.y[0,1];
                                 if((ii==0) && (jj==-1)){
                                       // if(css_test3[0,1]<1.0 && css_test3[0,1]>0.5){ //0.5< <1.0
                                        if(css_test3[0,1]<1.0 && css_test3[0,1]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_g[0,1]) + percent_g*poisson_s[0,1]+ flux_g_6[0,1];
                                            total_g += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test2[0,1]>0.0 && css_test2[0,1]<0.5){
                                        if(css_test2[0,1]>0.0 && css_test2[0,1]<above_value_merge){
                                                //merge judge
                                                total_g += flux_lg_g[0,1];
                                        }
                                        A3_g += Trhog*Tcpg/dt*cm_css_test2[0,1];
                                        A4_g += Trhog*Tcpg/dt*cm_css_test2[0,1];
                                 }
                            }
                        }//merge
                    }
                }
            }else{ 
                //level<level_interface
                if(css_test3[]<1.0 && css_test2[]>0){
                    flag=true;
                    //4faces
                    if(css_test3[]>0.0 && css_test3[]<1.0){
                            double flux_i=0.0;
                            // flux_i = (-flux_g[])/Delta*areasg[]/Delta  + percent_g*poisson_s[] + flux_g_6[];
                            flux_i = (-flux_g[])   + flux_g_6[];
                            total_g += flux_i;
                        }
                    //if there is fluid interface in the cell
                    if(css_test3[]<1.0){
                        if(css_test2[]>0 && css_test2[]<1.0){
                            // flux across the interface between gas and liquid
                            total_g += flux_lg_g[];
                        }
                    }
                    foreach_dimension(){
                        // if(css_test3[-1]<0.5 && css_test2[-1]>0.0){
                        if(css_test2[-1]>0.0){
                           double Al;
                            Al = Tkg*fm_fss_test2.x[]*1.0/(Delta*Delta);
                            total_g +=   Al*Tg[-1];
                            A3_g +=    Al; 
                        }

                        // if(css_test3[1]<0.5 && css_test2[1]>0.0){
                        if(css_test2[1]>0.0){
                            double Ar;
                            Ar = Tkg*fm_fss_test2.x[1]*1.0/(Delta*Delta);
                            total_g +=   Ar*Tg[1];
                            A3_g +=    Ar; 
                        }
                    }
                }//if
            }
           if(flag){
                Tg[] = (total_g+A4_g*Tg_old[])*beta1/A3_g + (1.0-beta1)*Tg[]; 
            }
      }
      
} //foreach

boundary({Tg,Tl,Ts});
restriction({Tg,Tl,Ts});
//////////////////////////////////////////////////
////////////residual
/////////////////////////////////////////////////
if(itt%1==0){
    double maxres_total=0.0;
    foreach(reduction(max:maxres_total)){
        double maxres_s=0.0;
        double maxres_l=0.0;
        double maxres_g=0.0;
       // if(css_test3[]>=0.5){
        bool flag_s=false;
        bool flag_l=false;
        bool flag_g=false;
        if(css_test3[]>0){
            ///////////////////////////////solid
            double total= percent_s*poisson_s[];
            double A3=Trhos*Tcps/dt*cm_css_test3[];
            double A4=Trhos*Tcps/dt*cm_css_test3[];
            if(level==level_interface){
                flag_s=true;
                if(css_test3[]>=above_value_merge){
                    if(css_test3[]>0.0 && css_test3[]<1.0){ //flux accross solid surface
                        double flux_i=0.0;
                        // flux_i = flux_l[]/Delta*areasl[]/Delta + flux_g[]/Delta*areasg[]/Delta  + percent_s*poisson_s[]+ flux_s_6[];
                        flux_i = flux_l[] + flux_g[]  + flux_s_6[];
                        total += flux_i;
                    }
                        //left
                       // if(css_test3[-1,0]>=0.5){
                        if(css_test3[-1,0]>=above_value_merge){
                            double Al;
                            Al = Tks*fm_fss_test3.x[]*1.0/(Delta*Delta);
                            total +=   Al*Ts[-1,0];
                            A3 +=    Al; 
                        }else if(css_test3[-1,0]>0.0){
                            //merge judge
                            if(fabs(merge_to_me_s_c[-1,0])>1e-12){
                                int ii = merge_to_me_s_position.x[-1,0];
                                int jj = merge_to_me_s_position.y[-1,0];
                                if((ii==1) && (jj==0)){
                                    //flux across surface
                                    double flux_i=0.0;
                                    flux_i = flux_l[-1,0] + flux_g[-1,0] + percent_s*poisson_s[-1,0]+ flux_s_6[-1,0];
                                    total += flux_i;
                                    A3 += Trhos*Tcps/dt*cm_css_test3[-1,0];
                                    A4 += Trhos*Tcps/dt*cm_css_test3[-1,0];
                                }
                                
                            }
                        } // css_test3[-1,0]

                        //right
                        // if(css_test3[1,0]>=0.5){
                        if(css_test3[1,0]>=above_value_merge){
                            double Ar;
                            Ar = Tks*fm_fss_test3.x[1,0]*1.0/(Delta*Delta);
                            total +=   Ar*Ts[1,0];
                            A3 +=    Ar; 
                        }else if(css_test3[1,0]>0.0){
                            //merge judge
                            if(fabs(merge_to_me_s_c[1,0])>1e-12){
                                int ii = merge_to_me_s_position.x[1,0];
                                int jj = merge_to_me_s_position.y[1,0];
                                if((ii==-1) && (jj==0)){
                                    //flux across surface
                                    double flux_i=0.0;
                                    flux_i = flux_l[1,0] + flux_g[1,0] + percent_s*poisson_s[1,0]+ flux_s_6[1,0];
                                    total += flux_i;
                                    A3 += Trhos*Tcps/dt*cm_css_test3[1,0];
                                    A4 += Trhos*Tcps/dt*cm_css_test3[1,0];
                                }
                            }
                        }// css_test3[1,0]

                         //bottom
                        // if(css_test3[0,-1]>=0.5){
                        if(css_test3[0,-1]>=above_value_merge){
                            double Al;
                            Al = Tks*fm_fss_test3.y[]*1.0/(Delta*Delta);
                            total +=   Al*Ts[0,-1];
                            A3 +=    Al; 
                        }else if(css_test3[0,-1]>0.0){
                            //merge judge
                            if(fabs(merge_to_me_s_c[0,-1])>1e-12){
                                int ii = merge_to_me_s_position.x[0,-1];
                                int jj = merge_to_me_s_position.y[0,-1];
                                if((ii==0) && (jj==1)){
                                    //flux across surface
                                    double flux_i=0.0;
                                    flux_i = flux_l[0,-1] + flux_g[0,-1] + percent_s*poisson_s[0,-1] + flux_s_6[0,-1];
                                    total += flux_i;
                                    A3 += Trhos*Tcps/dt*cm_css_test3[0,-1];
                                    A4 += Trhos*Tcps/dt*cm_css_test3[0,-1];
                                }
                            }
                        } // css_test3[0,-1]

                         //top
                        // if(css_test3[0,1]>=0.5){
                        if(css_test3[0,1]>=above_value_merge){
                            double Ar;
                            Ar = Tks*fm_fss_test3.y[0,1]*1.0/(Delta*Delta);
                            total +=   Ar*Ts[0,1];
                            A3 +=    Ar; 
                        }else if(css_test3[0,1]>0.0){
                            //merge judge
                            if(fabs(merge_to_me_s_c[0,1])>1e-12){
                                int ii = merge_to_me_s_position.x[0,1];
                                int jj = merge_to_me_s_position.y[0,1];
                                if((ii==0) && (jj==-1)){
                                    //flux across surface
                                    double flux_i=0.0;
                                    flux_i = flux_l[0,1] + flux_g[0,1]  + percent_s*poisson_s[0,1]+ flux_s_6[0,1];
                                    total += flux_i;
                                     A3 += Trhos*Tcps/dt*cm_css_test3[0,1];
                                     A4 += Trhos*Tcps/dt*cm_css_test3[0,1];
                                }
                            }
                        } // css_test3[-1,0]
                }

            }else{ //level<level_interface
                 // if(css_test3[]>=0.5){
                if(css_test3[]>0){
                    flag_s = true;
                     if(css_test3[]>0.0 && css_test3[]<1.0){ //flux accross solid surface
                        double flux_i=0.0;
                        // flux_i = flux_l[]/Delta*areasl[]/Delta + flux_g[]/Delta*areasg[]/Delta  + percent_s*poisson_s[]+ flux_s_6[];
                        flux_i = flux_l[] + flux_g[]  + flux_s_6[];
                        total += flux_i;
                    }
                    //4faces
                    foreach_dimension(){
                        if(css_test3[-1]>0.0){
                           double Al;
                            Al = Tks*fm_fss_test3.x[]*1.0/(Delta*Delta);
                            total +=   Al*Ts[-1];
                            A3 +=    Al; 
                        }

                        if(css_test3[1]>0.0){
                            double Ar;
                            Ar = Tks*fm_fss_test3.x[1]*1.0/(Delta*Delta);
                            total +=   Ar*Ts[1];
                            A3 +=    Ar; 
                        }
                    }
                }//if
            }
             if(flag_s){
                maxres_s = fabs(total+A4*Ts_old[]-A3*Ts[]);
            }
        } 



        //////////////
        ///////////////////////////////fluid
        //consider flux_l and flux across solid
      //    if(css_test3[]<0.5 && css_test[]>=0.5){
        if(css_test3[]<1.0){
            double total_l= percent_l*poisson_s[];
            double A3_l=Trhol*Tcpl/dt*cm_css_test[];
            double A4_l=Trhol*Tcpl/dt*cm_css_test[];       
            if(level==level_interface){
                if(css_test[]>=above_value_merge){
                    flag_l = true;
                //if there is solid interface in the cell
                    if(css_test3[]>0.0 && css_test3[]<1.0){
                        double flux_i=0.0;
                        // flux_i = (-flux_l[])/Delta*areasl[]/Delta + percent_l*poisson_s[]+ flux_l_6[];
                        flux_i = (-flux_l[]) + flux_l_6[];
                        total_l += flux_i;
                    }
                //if there is fluid interface in the cell
                if(css_test3[]<1.0){
                    if(css_test[]>0 && css_test[]<1.0){
                        // flux across the interface between gas and liquid
                        total_l += flux_lg_l[];
                    }
                }
                    //left
                    if(css_test3[-1,0]<1.0){
                        // if(css_test3[-1,0]<0.5 && css_test[-1,0]>=0.5){
                        if(css_test[-1,0]>=above_value_merge){
                                    double Al;
                                    Al = Tkl*fm_fss_test.x[]*1.0/(Delta*Delta);
                                    total_l +=   Al*Tl[-1,0];
                                    A3_l +=    Al; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_l_c[-1,0])>1e-12){
                                 int ii = merge_to_me_l_position.x[-1,0];
                                 int jj = merge_to_me_l_position.y[-1,0];
                                 if((ii==1) && (jj==0)){
                                        // if(css_test3[-1,0]<1.0 && css_test3[-1,0]>0.5){ //0.5< <1.0
                                        if(css_test3[-1,0]<1.0 && css_test3[-1,0]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_l[-1,0])  + percent_l*poisson_s[-1,0]+ flux_l_6[-1,0];
                                            total_l += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test[-1,0]>0.0 && css_test[-1,0]<0.5){
                                        if(css_test[-1,0]>0.0 && css_test[-1,0]<above_value_merge){
                                                //lg interface
                                                total_l += flux_lg_l[-1,0];
                                        }
                                        A3_l += Trhol*Tcpl/dt*cm_css_test[-1,0];
                                        A4_l += Trhol*Tcpl/dt*cm_css_test[-1,0];
                                 }
                            }
                        }//merge
                    }

                   //right
                   if(css_test3[1,0]<1.0){
                        // if(css_test3[1,0]<0.5 && css_test[1,0]>=0.5){
                        if(css_test[1,0]>=above_value_merge){
                                    double Ar;
                                    Ar = Tkl*fm_fss_test.x[1,0]*1.0/(Delta*Delta);
                                    total_l +=   Ar*Tl[1,0];
                                    A3_l +=    Ar; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_l_c[1,0])>1e-12){
                                 int ii = merge_to_me_l_position.x[1,0];
                                 int jj = merge_to_me_l_position.y[1,0];
                                 if((ii==-1) && (jj==0)){
                                        // if(css_test3[1,0]<1.0 && css_test3[1,0]>0.5){ //0.5< <1.0
                                        if(css_test3[1,0]<1.0 && css_test3[1,0]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_l[1,0]) + percent_l*poisson_s[1,0]+ flux_l_6[1,0];
                                            total_l += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test[1,0]>0.0 && css_test[1,0]<0.5){
                                        if(css_test[1,0]>0.0 && css_test[1,0]<above_value_merge){
                                                //merge judge
                                                total_l += flux_lg_l[1,0];
                                        }
                                        A3_l += Trhol*Tcpl/dt*cm_css_test[1,0];
                                        A4_l += Trhol*Tcpl/dt*cm_css_test[1,0];

                                 }
                            }
                        }//merge
                    }

                     //bottom
                   if(css_test3[0,-1]<1.0){
                        // if(css_test3[0,-1]<0.5 && css_test[0,-1]>=0.5){
                        if(css_test[0,-1]>=above_value_merge){
                                    double Al;
                                    Al = Tkl*fm_fss_test.y[]*1.0/(Delta*Delta);
                                    total_l +=   Al*Tl[0,-1];
                                    A3_l +=    Al; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_l_c[0,-1])>1e-12){
                                 int ii = merge_to_me_l_position.x[0,-1];
                                 int jj = merge_to_me_l_position.y[0,-1];
                                 if((ii==0) && (jj==1)){
                                        // if(css_test3[0,-1]<1.0 && css_test3[0,-1]>0.5){ //0.5< <1.0
                                        if(css_test3[0,-1]<1.0 && css_test3[0,-1]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_l[0,-1]) + percent_l*poisson_s[0,-1]+ flux_l_6[0,-1];
                                            total_l += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test[0,-1]>0.0 && css_test[0,-1]<0.5){
                                        if(css_test[0,-1]>0.0 && css_test[0,-1]<above_value_merge){
                                                //merge judge
                                                total_l += flux_lg_l[0,-1];
                                        }
                                        A3_l += Trhol*Tcpl/dt*cm_css_test[0,-1];
                                        A4_l += Trhol*Tcpl/dt*cm_css_test[0,-1];
                                 }
                            }
                        }//merge
                    }

                   //top
                   if(css_test3[0,1]<1.0){
                        // if(css_test3[0,1]<0.5 && css_test[0,1]>=0.5){
                        if(css_test[0,1]>=above_value_merge){
                                    double Ar;
                                    Ar = Tkl*fm_fss_test.y[0,1]*1.0/(Delta*Delta);
                                    total_l +=   Ar*Tl[0,1];
                                    A3_l +=    Ar; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_l_c[0,1])>1e-12){
                                 int ii = merge_to_me_l_position.x[0,1];
                                 int jj = merge_to_me_l_position.y[0,1];
                                 if((ii==0) && (jj==-1)){
                                         // if(css_test3[0,1]<1.0 && css_test3[0,1]>0.5){ //0.5< <1.0
                                        if(css_test3[0,1]<1.0 && css_test3[0,1]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_l[0,1])  + percent_l*poisson_s[0,1]+ flux_l_6[0,1];
                                            total_l += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test[0,1]>0.0 && css_test[0,1]<0.5){
                                        if(css_test[0,1]>0.0 && css_test[0,1]<above_value_merge){
                                                //merge judge
                                                total_l += flux_lg_l[0,1];
                                        }
                                         A3_l += Trhol*Tcpl/dt*cm_css_test[0,1];
                                         A4_l += Trhol*Tcpl/dt*cm_css_test[0,1];
                                 }
                            }
                        }//merge
                    }
                }
                
            }else{ 
                //level<level_interface
               // if(css_test3[]<0.5 && css_test[]>=0.5){
                 if(css_test3[]<1 && css_test[]>0){
                    flag_l=true;
                    //if there is solid interface in the cell
                    if(css_test3[]>0.0 && css_test3[]<1.0){
                            double flux_i=0.0;
                            // flux_i = (-flux_l[])/Delta*areasl[]/Delta + percent_l*poisson_s[]+ flux_l_6[];
                            flux_i = (-flux_l[]) + flux_l_6[];
                            total_l += flux_i;
                        }
                    //if there is fluid interface in the cell
                    if(css_test3[]<1.0){
                        if(css_test[]>0 && css_test[]<1.0){
                            // flux across the interface between gas and liquid
                            total_l += flux_lg_l[];
                        }
                    }
                    //4faces
                    foreach_dimension(){
                        // if(css_test3[-1]<0.5 && css_test[-1]>0.0){
                        if(css_test[-1]>0.0){
                           double Al;
                            Al = Tkl*fm_fss_test.x[]*1.0/(Delta*Delta);
                            total_l +=   Al*Tl[-1];
                            A3_l +=    Al; 
                        }

                        // if(css_test3[1]<0.5 && css_test[1]>0.0){
                        if(css_test[1]>0.0){
                            double Ar;
                            Ar = Tkl*fm_fss_test.x[1]*1.0/(Delta*Delta);
                            total_l +=   Ar*Tl[1];
                            A3_l +=    Ar; 
                        }
                    }
                }//if
            }
            // Tl[] = (total_l+A4_l*Tl[])*beta1/A3_l + (1.0-beta1)*Tl[]; 
            if(flag_l){
                maxres_l = fabs(total_l+A4_l*Tl_old[]-A3_l*Tl[]);
            }
    }


     //////////////
        ///////////////////////////////gas
        //consider flux_l and flux across solid
      //    if(css_test3[]<0.5 && css_test2[]>=0.5){
         if(css_test3[]<1.0){
            double total_g= percent_g*poisson_s[];
            double A3_g=Trhog*Tcpg/dt*cm_css_test2[];
            double A4_g=Trhog*Tcpg/dt*cm_css_test2[];       
            if(level==level_interface){
                if(css_test2[]>=above_value_merge){
                    flag_g = true;
                //if there is solid interface in the cell
                    if(css_test3[]>0.0 && css_test3[]<1.0){
                        double flux_i=0.0;
                        // flux_i = (-flux_g[])/Delta*areasg[]/Delta  + percent_g*poisson_s[]+ flux_g_6[];
                        flux_i = (-flux_g[]) + flux_g_6[];
                        total_g += flux_i;
                    }
                //if there is fluid interface in the cell
                if(css_test3[]<1.0){
                    if(css_test2[]>0 && css_test2[]<1.0){
                        // flux across the interface between gas and liquid
                        total_g += flux_lg_g[];
                    }
                }
                    //left
                    if(css_test3[-1,0]<1.0){
                        // if(css_test3[-1,0]<0.5 && css_test2[-1,0]>=0.5){
                        if(css_test2[-1,0]>=above_value_merge){
                                    double Al;
                                    Al = Tkg*fm_fss_test2.x[]*1.0/(Delta*Delta);
                                    total_g +=   Al*Tg[-1,0];
                                    A3_g +=    Al; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_g_c[-1,0])>1e-12){
                                 int ii = merge_to_me_g_position.x[-1,0];
                                 int jj = merge_to_me_g_position.y[-1,0];
                                 if((ii==1) && (jj==0)){
                                       // if(css_test3[-1,0]<1.0 && css_test3[-1,0]>0.5){ //0.5< <1.0
                                        if(css_test3[-1,0]<1.0 && css_test3[-1,0]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_g[-1,0]) + percent_g*poisson_s[-1,0]+ flux_g_6[-1,0];
                                            total_g += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test2[-1,0]>0.0 && css_test2[-1,0]<0.5){
                                        if(css_test2[-1,0]>0.0 && css_test2[-1,0]<above_value_merge){
                                                //lg interface
                                                total_g += flux_lg_g[-1,0];
                                        }
                                        A3_g += Trhog*Tcpg/dt*cm_css_test2[-1,0];
                                        A4_g += Trhog*Tcpg/dt*cm_css_test2[-1,0];
                                 }
                            }
                        }//merge
                    }

                   //right
                   if(css_test3[1,0]<1.0){
                        // if(css_test3[1,0]<0.5 && css_test2[1,0]>=0.5){
                        if(css_test2[1,0]>=above_value_merge){
                                    double Ar;
                                    Ar = Tkg*fm_fss_test2.x[1,0]*1.0/(Delta*Delta);
                                    total_g +=   Ar*Tg[1,0];
                                    A3_g +=    Ar; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_g_c[1,0])>1e-12){
                                 int ii = merge_to_me_g_position.x[1,0];
                                 int jj = merge_to_me_g_position.y[1,0];
                                 if((ii==-1) && (jj==0)){
                                        // if(css_test3[1,0]<1.0 && css_test3[1,0]>0.5){ //0.5< <1.0
                                        if(css_test3[1,0]<1.0 && css_test3[1,0]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_g[1,0]) + percent_g*poisson_s[1,0]+ flux_g_6[1,0];
                                            total_g += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test2[1,0]>0.0 && css_test2[1,0]<0.5){
                                        if(css_test2[1,0]>0.0 && css_test2[1,0]<above_value_merge){
                                                //merge judge
                                                total_g += flux_lg_g[1,0];
                                        }
                                         A3_g += Trhog*Tcpg/dt*cm_css_test2[1,0];
                                         A4_g += Trhog*Tcpg/dt*cm_css_test2[1,0];
                                 }
                            }
                        }//merge
                    }

                     //bottom
                   if(css_test3[0,-1]<1.0){
                        // if(css_test3[0,-1]<0.5 && css_test2[0,-1]>=0.5){
                        if(css_test2[0,-1]>=above_value_merge){
                                    double Al;
                                    Al = Tkg*fm_fss_test2.y[]*1.0/(Delta*Delta);
                                    total_g +=   Al*Tg[0,-1];
                                    A3_g +=    Al; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_g_c[0,-1])>1e-12){
                                 int ii = merge_to_me_g_position.x[0,-1];
                                 int jj = merge_to_me_g_position.y[0,-1];
                                 if((ii==0) && (jj==1)){
                                       // if(css_test3[0,-1]<1.0 && css_test3[0,-1]>0.5){ //0.5< <1.0
                                        if(css_test3[0,-1]<1.0 && css_test3[0,-1]>0.0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_g[0,-1]) + percent_g*poisson_s[0,-1]+ flux_g_6[0,-1];
                                            total_g += flux_i;
                                        }
                                        //merge_energy_from fluid
                                       // if(css_test2[0,-1]>0.0 && css_test2[0,-1]<0.5){
                                        if(css_test2[0,-1]>0.0 && css_test2[0,-1]<above_value_merge){
                                                //merge judge
                                                total_g += flux_lg_g[0,-1];
                                        }
                                        A3_g += Trhog*Tcpg/dt*cm_css_test2[0,-1];
                                        A4_g += Trhog*Tcpg/dt*cm_css_test2[0,-1];
                                 }
                            }
                        }//merge
                    }

                   //top
                   if(css_test3[0,1]<1.0){
                        // if(css_test3[0,1]<0.5 && css_test2[0,1]>=0.5){
                        if(css_test2[0,1]>=above_value_merge){
                                    double Ar;
                                    Ar = Tkg*fm_fss_test2.y[0,1]*1.0/(Delta*Delta);
                                    total_g +=   Ar*Tg[0,1];
                                    A3_g +=    Ar; 
                        }else{ //merge
                            //merge_energy from solid
                            if(fabs(merge_to_me_g_c[0,1])>1e-12){
                                 int ii = merge_to_me_g_position.x[0,1];
                                 int jj = merge_to_me_g_position.y[0,1];
                                 if((ii==0) && (jj==-1)){
                                       // if(css_test3[0,1]<1.0 && css_test3[0,1]>0.5){ //0.5< <1.0
                                        if(css_test3[0,1]<1.0 && css_test3[0,1]>0){ //0.5< <1.0
                                            //solid interface
                                            double flux_i=0.0;
                                            flux_i = (-flux_g[0,1]) + percent_g*poisson_s[0,1]+ flux_g_6[0,1];
                                            total_g += flux_i;
                                        }
                                        //merge_energy_from fluid
                                        // if(css_test2[0,1]>0.0 && css_test2[0,1]<0.5){
                                        if(css_test2[0,1]>0.0 && css_test2[0,1]<above_value_merge){
                                                //merge judge
                                                total_g += flux_lg_g[0,1];
                                        }
                                        A3_g += Trhog*Tcpg/dt*cm_css_test2[0,1];
                                        A4_g += Trhog*Tcpg/dt*cm_css_test2[0,1];
                                 }
                            }
                        }//merge
                    }
                }
            }else{ 
                //level<level_interface
                // if(css_test3[]<0.5 && css_test2[]>=0.5){
                if(css_test3[]<1.0 && css_test2[]>0){
                    flag_g=true;
                     if(css_test3[]>0.0 && css_test3[]<1.0){
                            double flux_i=0.0;
                            // flux_i = (-flux_g[])/Delta*areasg[]/Delta  + percent_g*poisson_s[]+ flux_g_6[];
                            flux_i = (-flux_g[]) + flux_g_6[];
                            total_g += flux_i;
                        }
                    //if there is fluid interface in the cell
                    if(css_test3[]<1.0){
                        if(css_test2[]>0 && css_test2[]<1.0){
                            // flux across the interface between gas and liquid
                            total_g += flux_lg_g[];
                        }
                    }
                    //4faces
                    foreach_dimension(){
                        // if(css_test3[-1]<0.5 && css_test2[-1]>0.0){
                        if(css_test2[-1]>0.0){
                           double Al;
                            Al = Tkg*fm_fss_test2.x[]*1.0/(Delta*Delta);
                            total_g +=   Al*Tg[-1];
                            A3_g +=    Al; 
                        }

                       // if(css_test3[1]<0.5 && css_test2[1]>0.0){
                        if(css_test2[1]>0.0){
                            double Ar;
                            Ar = Tkg*fm_fss_test2.x[1]*1.0/(Delta*Delta);
                            total_g +=   Ar*Tg[1];
                            A3_g +=    Ar; 
                        }
                    }
                }//if
            }
            // Tg[] = (total_g+A4_g*Tg[])*beta1/A3_g + (1.0-beta1)*Tg[];
             if(flag_g){
                maxres_g = fabs(total_g+A4_g*Tg_old[]-A3_g*Tg[]);
            }
        }
        resl[] = maxres_l;
        resg[] = maxres_g;
        ress[] = maxres_s;
        maxres_total = max_function(maxres_s/(Trhos*Tcps),maxres_l/(Trhol*Tcpl),maxres_g/(Trhog*Tcpg));
      
    }//foreach residual

if(poisson_check){
      if(pid()==0){
            char name93[80]; 
            sprintf(name93,"poisson_check.dat");
            FILE * fp93 = fopen(name93,"a");
            int num2=0;
              fprintf (fp93," %g\n",  maxres_total);
            fclose(fp93);
      }
  }
    
	  if (maxres_total<T_tolerance){
            if(pid()==0){
                printf("max_total=%g ,itt=%d step\n",maxres_total,itt);
            }
	       break;
	  }

      boundary({Ts,Tl,Tg});
}//every 4 times


    }//itt
    if(pid()==0){
        printf("3T itt=%d step\n",itt);
    }

if(poisson_check){
      if(pid()==0){
            char name93[80]; 
            sprintf(name93,"poisson_check.dat");
            FILE * fp93 = fopen(name93,"a");
            int num2=0;
              fprintf (fp93," \n \n \n");
            fclose(fp93);
      }
  }

}