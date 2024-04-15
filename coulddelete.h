  scalar flux_lg_l[],flux_lg_g[];
  foreach(){
      flux_lg_l[]=0.0;
      flux_lg_g[]=0.0;
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
                        coord nl = mycs (point, css_test);//interface_normal(point, ff);//interface_normal7(point,ff,hhh); //
                        double alphal = plane_alpha (css_test[], nl);
                        coord ppp;
                        arealg = plane_area_center (nl, alphal, &ppp);
                        if (metric_embed_factor)  
                                arealg *= metric_embed_factor (point, ppp);
                        double gradl2 = phase1Tg[];// for fluid this is value for from water to out direction
                        mua = Tkl*max(y,1e-20);
                        fa =1.0;
                        c1 =  - mua/(fa + SEPS)*gradl2*arealg/Delta;
                        e = 0.0;
                        flux_lg_l[] = c1;
                    }

                    double arealg2=0.0;
                    if(1==1){
                        double mua = 0., fa = 0.;
                        foreach_dimension() {
                            //mua += alpha.x[] + alpha.x[1];
                            fa  += fss_test2.x[] + fss_test2.x[1];
                          }
                        //mua= k0;//Tkg;
                        coord ng = interface_normal(point, css_test2);//ff_oppo);// interface_normal7(point,ff_oppo,hhh_ff_oppo); // interface_normal7(point,css_test,hhh); //interface_normal7
                          coord ppp2;
                          double alphag = plane_alpha (css_test2[], ng);
                          arealg2 = plane_area_center (n, alphag, &ppp2);
                          if (metric_embed_factor)  
                                arealg2 *= metric_embed_factor (point, ppp2);
   //direction of phase1Tg and phase0Tg is from fluid to out; for liquid it is out, for gas it is in; so add one -, finally it is +
                        double gradg2 = -phase0Tg[]; 
                        mua = Tkg*max(y,1e-20);
                        fa= 1.0;
                        c1_2 =  - mua/(fa + SEPS)*gradg2*arealg2/Delta;
                       e_2 = 0.0;
                       flux_lg_g[] = c1_2;
                    }

             }
       }
    }
}