#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  8;
    float avdat=  0.3007165E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[  9]={
         0.17103780E-05,-0.38642895E+02, 0.64414642E+02, 0.40067407E+00,
         0.12146097E+01,-0.73236096E+00,-0.66461331E+00, 0.16151741E-02,
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
    float x31 = x3;
    float x41 = x4;

//                 function

    float v_target_x                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]        *x31    
        +coeff[  3]*x11*x21        
        +coeff[  4]        *x31*x41
        +coeff[  5]    *x21*x31    
        +coeff[  6]*x11        *x41
        +coeff[  7]    *x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  6]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]*x31+coeff[  7],2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x41+coeff[  5]*x21,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  6]*x11,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 15;
    float avdat=  0.2483026E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 16]={
        -0.18005319E-01,-0.38473625E+02, 0.64121407E+02,-0.50121516E+00,
         0.45140710E+00, 0.13806960E+00, 0.29339848E-03,-0.13684933E-02,
         0.32445502E-02,-0.16181588E-02,-0.42428743E-03,-0.18069975E+00,
         0.53213513E+00,-0.37243137E+00,-0.11873957E-02,
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
        +coeff[ 14]    *x21*x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  7]*x11*x21*x41+coeff[  8]*x21*x31*x41+4*coeff[ 10]*x13*x32+2*coeff[ 11]*x11+coeff[ 12]*x31,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  5]*x21+coeff[  7]*x12*x41+coeff[  8]*x11*x31*x41+coeff[  9]*x32*x41+coeff[ 14]*x32,2)
        +dx3*dx3*pow(0+4*coeff[  6]*x33+coeff[  8]*x11*x21*x41+2*coeff[  9]*x31*x21*x41+2*coeff[ 10]*x31*x14+coeff[ 12]*x11+2*coeff[ 13]*x31+2*coeff[ 14]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41+coeff[  7]*x12*x21+coeff[  8]*x11*x21*x31+coeff[  9]*x21*x32,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.6615363E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.44696168E-04,-0.55288993E-01, 0.65973639E-01, 0.80583356E-02,
        -0.64663920E-02,-0.31267535E-02, 0.74906566E-03,-0.24917489E-02,
         0.30552268E-02,-0.55265300E-04, 0.35216919E-04,-0.89418660E-04,
         0.13287037E-03, 0.22346248E-03,-0.31502920E-03,-0.67912537E-04,
         0.31542472E-03,-0.32238354E-03,-0.11235826E-03, 0.15776804E-03,
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
        +coeff[ 18]*x12*x21    *x41
        +coeff[ 19]*x12        *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  3]*x31+2*coeff[  7]*x11+coeff[ 10]+coeff[ 11]*x21*x31+2*coeff[ 15]*x11*x41+coeff[ 16]*x31*x41+2*coeff[ 18]*x11*x21*x41+2*coeff[ 19]*x11*x42,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  5]*x41+2*coeff[  6]*x21+coeff[ 11]*x11*x31+coeff[ 12]*x32+coeff[ 13]*x32*x41+coeff[ 18]*x12*x41,2)
        +dx3*dx3*pow(0+coeff[  3]*x11+2*coeff[  4]*x31+coeff[  9]+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+2*coeff[ 13]*x31*x21*x41+2*coeff[ 14]*x31*x42+coeff[ 16]*x11*x41+2*coeff[ 17]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  5]*x21+2*coeff[  8]*x41+coeff[ 13]*x21*x32+2*coeff[ 14]*x41*x32+coeff[ 15]*x12+coeff[ 16]*x11*x31+coeff[ 17]*x32+coeff[ 18]*x12*x21+2*coeff[ 19]*x41*x12,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4717777E-05;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.25069742E-08,-0.57049830E-01, 0.68737909E-01, 0.29028426E-02,
        -0.52475529E-02,-0.48692352E-02, 0.87693902E-02, 0.48988000E-04,
        -0.29565130E-04,-0.52476229E-04,-0.13583319E-03, 0.25910118E-04,
         0.95172822E-04,-0.12358285E-03,-0.46825022E-03, 0.32557701E-03,
        -0.47143309E-04, 0.16557422E-03, 0.16988348E-03, 0.54909273E-04,
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
    float x23 = x22*x2;
    float x31 = x3;
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;

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
        +coeff[ 16]    *x22*x31    
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x12*x21*x31    
        +coeff[ 18]    *x23*x31    
        +coeff[ 19]*x13        *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[  9]*x42+coeff[ 11]*x22+3*coeff[ 13]*x12*x21+2*coeff[ 17]*x11*x21*x31+3*coeff[ 19]*x12*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  8]+2*coeff[ 11]*x21*x11+coeff[ 13]*x13+2*coeff[ 14]*x21*x31*x41+coeff[ 15]*x31*x42+2*coeff[ 16]*x21*x31+coeff[ 17]*x12*x31+3*coeff[ 18]*x22*x31,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+3*coeff[ 10]*x32*x41+coeff[ 12]*x42+coeff[ 14]*x22*x41+coeff[ 15]*x21*x42+coeff[ 16]*x22+coeff[ 17]*x12*x21+coeff[ 18]*x23,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  7]+2*coeff[  9]*x41*x11+coeff[ 10]*x33+2*coeff[ 12]*x41*x31+coeff[ 14]*x22*x31+2*coeff[ 15]*x41*x21*x31+coeff[ 19]*x13,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
