#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 17;
    float avdat=  0.1733793E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 18]={
        -0.17240975E-03,-0.39866837E+02, 0.65869278E+02, 0.81262314E+00,
         0.24181795E+01,-0.43841198E-01, 0.72037324E-01,-0.14697678E+01,
        -0.13378297E+01, 0.33472430E-01,-0.47330443E-01,-0.14377631E-01,
        -0.28308472E-01, 0.70746593E-01,-0.69704484E-02, 0.94212554E-02,
         0.13457537E-02,
              0.      };
    int ientry=0;
    int i;
    if (ientry==0){
        ientry=1;
        for(i=0;i<m;i++){
            if(xmin[i]==xmax[i]) continue;
            scale[i]=2./(xmax[i]-xmin[i]);
        }
    }
//  normalize variables between -1 and +1
    float x1 =1.+(x[  0]-xmax[  0])*scale[  0];
    float dx1 =dx[  0]*scale[  0];
    float x2 =1.+(x[  1]-xmax[  1])*scale[  1];
    float dx2 =dx[  1]*scale[  1];
    float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
    float dx3 =dx[  2]*scale[  2];
    float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
    float dx4 =dx[  3]*scale[  3];
//  set up monomials   functions
    float x11 = x1;
    float x12 = x11*x1;
    float x21 = x2;
    float x22 = x21*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x41 = x4;
    float x42 = x41*x4;

//                 function

    float v_target_x                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]        *x31    
        +coeff[  3]*x11*x21        
        +coeff[  4]        *x31*x41
        +coeff[  5]    *x21        
        +coeff[  6]            *x41
        +coeff[  7]    *x21*x31    
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]*x11        *x41
        +coeff[  9]*x11*x22        
        +coeff[ 10]*x11*x21    *x41
        +coeff[ 11]    *x22*x31    
        +coeff[ 12]    *x21*x31*x41
        +coeff[ 13]        *x31*x42
        +coeff[ 14]*x12    *x31    
        +coeff[ 15]*x11    *x32    
        +coeff[ 16]        *x32*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  8]*x41+coeff[  9]*x22+coeff[ 10]*x21*x41+2*coeff[ 14]*x11*x31+coeff[ 15]*x32,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]+coeff[  7]*x31+2*coeff[  9]*x21*x11+coeff[ 10]*x11*x41+2*coeff[ 11]*x21*x31+coeff[ 12]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x41+coeff[  7]*x21+coeff[ 11]*x22+coeff[ 12]*x21*x41+coeff[ 13]*x42+coeff[ 14]*x12+2*coeff[ 15]*x31*x11+2*coeff[ 16]*x31*x42,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  6]+coeff[  8]*x11+coeff[ 10]*x11*x21+coeff[ 12]*x21*x31+2*coeff[ 13]*x41*x31+2*coeff[ 16]*x41*x32,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.5485286E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.33343911E-01,-0.39247921E+02, 0.64822624E+02,-0.91186875E+00,
         0.82560766E+00, 0.24816869E+00, 0.53375959E-03, 0.23526482E-02,
        -0.73899878E-02, 0.53904797E-02,-0.89232152E-03,-0.29884049E+00,
         0.90518290E+00,-0.64353520E+00,-0.17964412E-01,-0.66046221E-02,
         0.14188768E-01, 0.11676691E-01, 0.29187324E-02,-0.13541923E-01,
              0.      };
    int ientry=0;
    int i;
    if (ientry==0){
        ientry=1;
        for(i=0;i<m;i++){
            if(xmin[i]==xmax[i]) continue;
            scale[i]=2./(xmax[i]-xmin[i]);
        }
    }
//  normalize variables between -1 and +1
    float x1 =1.+(x[  0]-xmax[  0])*scale[  0];
    float dx1 =dx[  0]*scale[  0];
    float x2 =1.+(x[  1]-xmax[  1])*scale[  1];
    float dx2 =dx[  1]*scale[  1];
    float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
    float dx3 =dx[  2]*scale[  2];
    float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
    float dx4 =dx[  3]*scale[  3];
//  set up monomials   functions
    float x11 = x1;
    float x12 = x11*x1;
    float x13 = x12*x1;
    float x14 = x13*x1;
    float x21 = x2;
    float x22 = x21*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]    *x21    *x41
        +coeff[  4]            *x42
        +coeff[  5]    *x22        
        +coeff[  6]        *x34    
        +coeff[  7]*x12*x21    *x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x11*x21*x31*x41
        +coeff[  9]    *x21*x32*x41
        +coeff[ 10]*x14    *x32    
        +coeff[ 11]*x12            
        +coeff[ 12]*x11    *x31    
        +coeff[ 13]        *x32    
        +coeff[ 14]        *x31    
        +coeff[ 15]    *x21*x32    
        +coeff[ 16]*x11            
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]*x12*x21        
        +coeff[ 18]    *x22    *x41
        +coeff[ 19]*x11    *x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  7]*x11*x21*x41+coeff[  8]*x21*x31*x41+4*coeff[ 10]*x13*x32+2*coeff[ 11]*x11+coeff[ 12]*x31+coeff[ 16]+2*coeff[ 17]*x11*x21+coeff[ 19]*x31*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  5]*x21+coeff[  7]*x12*x41+coeff[  8]*x11*x31*x41+coeff[  9]*x32*x41+coeff[ 15]*x32+coeff[ 17]*x12+2*coeff[ 18]*x21*x41,2)
        +dx3*dx3*pow(0+4*coeff[  6]*x33+coeff[  8]*x11*x21*x41+2*coeff[  9]*x31*x21*x41+2*coeff[ 10]*x31*x14+coeff[ 12]*x11+2*coeff[ 13]*x31+coeff[ 14]+2*coeff[ 15]*x31*x21+coeff[ 19]*x11*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41+coeff[  7]*x12*x21+coeff[  8]*x11*x21*x31+coeff[  9]*x21*x32+coeff[ 18]*x22+coeff[ 19]*x11*x31,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1199161E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.98469784E-04,-0.50748106E-01, 0.58753278E-01,-0.50236071E-02,
         0.16777309E-02, 0.16200874E-01,-0.12953618E-01,-0.66189650E-02,
         0.40745013E-03,-0.64397824E-03, 0.62390072E-02, 0.40231968E-03,
        -0.57370140E-03,-0.24109300E-04, 0.63340100E-04,-0.15851752E-03,
        -0.20292471E-03, 0.28266691E-03, 0.95387266E-04,-0.41776993E-05,
              0.      };
    int ientry=0;
    int i;
    if (ientry==0){
        ientry=1;
        for(i=0;i<m;i++){
            if(xmin[i]==xmax[i]) continue;
            scale[i]=2./(xmax[i]-xmin[i]);
        }
    }
//  normalize variables between -1 and +1
    float x1 =1.+(x[  0]-xmax[  0])*scale[  0];
    float dx1 =dx[  0]*scale[  0];
    float x2 =1.+(x[  1]-xmax[  1])*scale[  1];
    float dx2 =dx[  1]*scale[  1];
    float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
    float dx3 =dx[  2]*scale[  2];
    float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
    float dx4 =dx[  3]*scale[  3];
//  set up monomials   functions
    float x11 = x1;
    float x12 = x11*x1;
    float x21 = x2;
    float x22 = x21*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x41 = x4;
    float x42 = x41*x4;

//                 function

    float v_target_theta                            =avdat
        +coeff[  0]                
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]*x12            
        +coeff[  4]    *x22        
        +coeff[  5]*x11    *x31    
        +coeff[  6]        *x32    
        +coeff[  7]    *x21    *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x11            
        +coeff[  9]        *x31    
        +coeff[ 10]            *x42
        +coeff[ 11]    *x21*x32*x41
        +coeff[ 12]        *x32*x42
        +coeff[ 13]*x11*x21*x31    
        +coeff[ 14]    *x21*x32    
        +coeff[ 15]*x11    *x31*x41
        +coeff[ 16]*x12*x22        
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12*x21    *x41
        +coeff[ 18]*x12        *x41
        +coeff[ 19]*x11*x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  3]*x11+coeff[  5]*x31+coeff[  8]+coeff[ 13]*x21*x31+coeff[ 15]*x31*x41+2*coeff[ 16]*x11*x22+2*coeff[ 17]*x11*x21*x41+2*coeff[ 18]*x11*x41+coeff[ 19]*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  7]*x41+coeff[ 11]*x32*x41+coeff[ 13]*x11*x31+coeff[ 14]*x32+2*coeff[ 16]*x21*x12+coeff[ 17]*x12*x41+coeff[ 19]*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  5]*x11+2*coeff[  6]*x31+coeff[  9]+2*coeff[ 11]*x31*x21*x41+2*coeff[ 12]*x31*x42+coeff[ 13]*x11*x21+2*coeff[ 14]*x31*x21+coeff[ 15]*x11*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  7]*x21+2*coeff[ 10]*x41+coeff[ 11]*x21*x32+2*coeff[ 12]*x41*x32+coeff[ 15]*x11*x31+coeff[ 17]*x12*x21+coeff[ 18]*x12+coeff[ 19]*x11*x21,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.7719393E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.38175196E-09, 0.52933129E-02, 0.58960151E-02,-0.10589760E-01,
        -0.97927395E-02, 0.17523251E-01,-0.11333450E+01,-0.50518092E-01,
         0.38588490E+01, 0.56705896E-01,-0.45437360E+01, 0.18334973E+01,
         0.22785677E-01, 0.30585181E-01, 0.18228540E+01,-0.24658298E-01,
        -0.35290096E-01,-0.63272843E+01, 0.74116483E+01,-0.29171522E+01,
              0.      };
    int ientry=0;
    int i;
    if (ientry==0){
        ientry=1;
        for(i=0;i<m;i++){
            if(xmin[i]==xmax[i]) continue;
            scale[i]=2./(xmax[i]-xmin[i]);
        }
    }
//  normalize variables between -1 and +1
    float x1 =1.+(x[  0]-xmax[  0])*scale[  0];
    float dx1 =dx[  0]*scale[  0];
    float x2 =1.+(x[  1]-xmax[  1])*scale[  1];
    float dx2 =dx[  1]*scale[  1];
    float x3 =1.+(x[  2]-xmax[  2])*scale[  2];
    float dx3 =dx[  2]*scale[  2];
    float x4 =1.+(x[  3]-xmax[  3])*scale[  3];
    float dx4 =dx[  3]*scale[  3];
//  set up monomials   functions
    float x11 = x1;
    float x12 = x11*x1;
    float x13 = x12*x1;
    float x14 = x13*x1;
    float x21 = x2;
    float x22 = x21*x2;
    float x23 = x22*x2;
    float x24 = x23*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;

//                 function

    float v_target_phi                              =avdat
        +coeff[  0]                
        +coeff[  1]        *x31    
        +coeff[  2]*x11*x21        
        +coeff[  3]    *x21*x31    
        +coeff[  4]*x11        *x41
        +coeff[  5]        *x31*x41
        +coeff[  6]*x13            
        +coeff[  7]*x11*x22        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]*x12    *x31    
        +coeff[  9]    *x22*x31    
        +coeff[ 10]*x11    *x32    
        +coeff[ 11]        *x33    
        +coeff[ 12]*x13*x22        
        +coeff[ 13]*x11*x24        
        +coeff[ 14]*x14    *x31    
        +coeff[ 15]*x12*x22*x31    
        +coeff[ 16]    *x24*x31    
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x13    *x32    
        +coeff[ 18]*x12    *x33    
        +coeff[ 19]*x11    *x34    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  2]*x21+coeff[  4]*x41+3*coeff[  6]*x12+coeff[  7]*x22+2*coeff[  8]*x11*x31+coeff[ 10]*x32+3*coeff[ 12]*x12*x22+coeff[ 13]*x24+4*coeff[ 14]*x13*x31+2*coeff[ 15]*x11*x22*x31+3*coeff[ 17]*x12*x32+2*coeff[ 18]*x11*x33+coeff[ 19]*x34,2)
        +dx2*dx2*pow(0+coeff[  2]*x11+coeff[  3]*x31+2*coeff[  7]*x21*x11+2*coeff[  9]*x21*x31+2*coeff[ 12]*x21*x13+4*coeff[ 13]*x23*x11+2*coeff[ 15]*x21*x12*x31+4*coeff[ 16]*x23*x31,2)
        +dx3*dx3*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[  8]*x12+coeff[  9]*x22+2*coeff[ 10]*x31*x11+3*coeff[ 11]*x32+coeff[ 14]*x14+coeff[ 15]*x12*x22+coeff[ 16]*x24+2*coeff[ 17]*x31*x13+3*coeff[ 18]*x32*x12+4*coeff[ 19]*x33*x11,2)
        +dx4*dx4*pow(0+coeff[  4]*x11+coeff[  5]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
