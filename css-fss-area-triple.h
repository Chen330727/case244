//areasl,areasg,arealg
//css_test,css_test2,css_test3
//fss_test,fss_test2,fss_test3
// #include "./line9-4-4-for-basilisk-css-test.h"
void get_css_fss_areaslg_triple_point(){
    foreach(){
        areasl[]=0.0;
        areasg[]=0.0;
        arealg[]=0.0;
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
        if(css_test[]>0.0 && css_test[]<1.0 && ff[]<0.0 && ff[]>0.0){
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

            //solid: nsf.x*x + nsf.y*y = alphasf
            //fg surface: nlg.x*x + nlg.y*y = alphalg

            //in line-9-4-for-basilisk: n.x*x + n.y*y + alpha = 0
            //so there is some difference between alphas;


        }
    }
}