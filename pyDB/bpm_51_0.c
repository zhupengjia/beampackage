#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 11;
    float avdat=  0.3007742E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
         0.32820772E-05,-0.39357697E+02, 0.65275352E+02, 0.43696517E+00,
        -0.72512621E+00,-0.79813731E+00, 0.13240423E+01,-0.46012402E-02,
         0.76291976E-02, 0.60632052E-02,-0.25531179E-02,
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
    float x21 = x2;
    float x22 = x21*x2;
    float x31 = x3;
    float x41 = x4;
    float x42 = x41*x4;

//                 function

    float v_target_x                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]        *x31    
        +coeff[  3]*x11*x21        
        +coeff[  4]*x11        *x41
        +coeff[  5]    *x21*x31    
        +coeff[  6]        *x31*x41
        +coeff[  7]    *x21        
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]            *x41
        +coeff[  9]        *x31*x42
        +coeff[ 10]*x11*x22        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  4]*x41+coeff[ 10]*x22,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]*x31+coeff[  7]+2*coeff[ 10]*x21*x11,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  5]*x21+coeff[  6]*x41+coeff[  9]*x42,2)
        +dx4*dx4*pow(0+coeff[  4]*x11+coeff[  6]*x31+coeff[  8]+2*coeff[  9]*x41*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 13;
    float avdat=  0.3288230E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 14]={
        -0.17452467E-01,-0.39170357E+02, 0.64952637E+02, 0.14970511E+00,
        -0.54574549E+00,-0.45423475E+00, 0.49272060E+00,-0.48037444E-03,
         0.22580313E-03,-0.21162736E+00, 0.63371348E+00, 0.49435389E-02,
        -0.67654639E-02,
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
        +coeff[  3]    *x22        
        +coeff[  4]    *x21    *x41
        +coeff[  5]        *x32    
        +coeff[  6]            *x42
        +coeff[  7]*x14            
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]        *x34    
        +coeff[  9]*x12            
        +coeff[ 10]*x11    *x31    
        +coeff[ 11]*x11*x21*x31    
        +coeff[ 12]    *x21*x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  7]*x13+2*coeff[  9]*x11+coeff[ 10]*x31+coeff[ 11]*x21*x31,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  3]*x21+coeff[  4]*x41+coeff[ 11]*x11*x31+coeff[ 12]*x32,2)
        +dx3*dx3*pow(0+2*coeff[  5]*x31+4*coeff[  8]*x33+coeff[ 10]*x11+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  4]*x21+2*coeff[  6]*x41,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.6198929E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.43627970E-04,-0.55273090E-01, 0.65954469E-01, 0.80718938E-02,
        -0.64721913E-02,-0.32495707E-02, 0.80068590E-03,-0.24970022E-02,
         0.31294769E-02,-0.55314118E-04, 0.35255613E-04,-0.10983197E-03,
         0.15633830E-03, 0.23921982E-03,-0.34575051E-03,-0.33740656E-04,
         0.25909260E-03,-0.30381646E-03,-0.11955875E-03, 0.16886563E-03,
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
        +coeff[  3]*x11    *x31    
        +coeff[  4]        *x32    
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x22        
        +coeff[  7]*x12            
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]            *x42
        +coeff[  9]        *x31    
        +coeff[ 10]*x11            
        +coeff[ 11]*x11*x21*x31    
        +coeff[ 12]    *x21*x32    
        +coeff[ 13]    *x21*x32*x41
        +coeff[ 14]        *x32*x42
        +coeff[ 15]*x12        *x41
        +coeff[ 16]*x11    *x31*x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]        *x32*x41
        +coeff[ 18]*x12*x22        
        +coeff[ 19]*x12*x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  3]*x31+2*coeff[  7]*x11+coeff[ 10]+coeff[ 11]*x21*x31+2*coeff[ 15]*x11*x41+coeff[ 16]*x31*x41+2*coeff[ 18]*x11*x22+2*coeff[ 19]*x11*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  5]*x41+2*coeff[  6]*x21+coeff[ 11]*x11*x31+coeff[ 12]*x32+coeff[ 13]*x32*x41+2*coeff[ 18]*x21*x12+coeff[ 19]*x12*x41,2)
        +dx3*dx3*pow(0+coeff[  3]*x11+2*coeff[  4]*x31+coeff[  9]+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+2*coeff[ 13]*x31*x21*x41+2*coeff[ 14]*x31*x42+coeff[ 16]*x11*x41+2*coeff[ 17]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  5]*x21+2*coeff[  8]*x41+coeff[ 13]*x21*x32+2*coeff[ 14]*x41*x32+coeff[ 15]*x12+coeff[ 16]*x11*x31+coeff[ 17]*x32+coeff[ 19]*x12*x21,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4719375E-05;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.20457591E-08,-0.57034850E-01, 0.68718828E-01, 0.28439844E-02,
        -0.51943925E-02,-0.48004114E-02, 0.87057315E-02, 0.49027632E-04,
        -0.29588313E-04,-0.52681207E-04,-0.90125446E-04, 0.25963842E-04,
         0.18586109E-03,-0.33863365E-04, 0.15506583E-03,-0.43798884E-03,
         0.31277724E-03, 0.15614320E-04,-0.15054150E-03, 0.86670858E-04,
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
    float x21 = x2;
    float x22 = x21*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;

//                 function

    float v_target_phi                              =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]        *x31    
        +coeff[  3]*x11*x21        
        +coeff[  4]    *x21*x31    
        +coeff[  5]*x11        *x41
        +coeff[  6]        *x31*x41
        +coeff[  7]            *x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]    *x21        
        +coeff[  9]*x11        *x42
        +coeff[ 10]        *x33*x41
        +coeff[ 11]*x11*x22        
        +coeff[ 12]        *x31*x42
        +coeff[ 13]*x13*x21        
        +coeff[ 14]    *x22*x31*x41
        +coeff[ 15]    *x21*x31*x42
        +coeff[ 16]        *x31*x43
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]    *x22*x31    
        +coeff[ 18]    *x21*x31*x41
        +coeff[ 19]*x12*x21*x31    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[  9]*x42+coeff[ 11]*x22+3*coeff[ 13]*x12*x21+2*coeff[ 19]*x11*x21*x31,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  8]+2*coeff[ 11]*x21*x11+coeff[ 13]*x13+2*coeff[ 14]*x21*x31*x41+coeff[ 15]*x31*x42+2*coeff[ 17]*x21*x31+coeff[ 18]*x31*x41+coeff[ 19]*x12*x31,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+3*coeff[ 10]*x32*x41+coeff[ 12]*x42+coeff[ 14]*x22*x41+coeff[ 15]*x21*x42+coeff[ 16]*x43+coeff[ 17]*x22+coeff[ 18]*x21*x41+coeff[ 19]*x12*x21,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  7]+2*coeff[  9]*x41*x11+coeff[ 10]*x33+2*coeff[ 12]*x41*x31+coeff[ 14]*x22*x31+2*coeff[ 15]*x41*x21*x31+3*coeff[ 16]*x42*x31+coeff[ 18]*x21*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
