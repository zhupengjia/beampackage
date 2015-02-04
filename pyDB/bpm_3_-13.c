#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 12;
    float avdat=  0.1237288E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 13]={
        -0.64973772E-03,-0.37926899E+02, 0.37612317E+01, 0.63399105E+02,
        -0.61777749E+01,-0.36704674E-01, 0.55649925E-01, 0.66078395E-01,
        -0.10121660E+00,-0.39459579E-02, 0.58683706E-02,-0.55616377E-02,
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

//                 function

    float v_target_x                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11*x21        
        +coeff[  6]*x11        *x41
        +coeff[  7]    *x21*x31    
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]        *x31*x41
        +coeff[  9]    *x21*x32*x41
        +coeff[ 10]*x12*x22        
        +coeff[ 11]*x12*x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41+2*coeff[ 10]*x11*x22+2*coeff[ 11]*x11*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  7]*x31+coeff[  9]*x32*x41+2*coeff[ 10]*x21*x12+coeff[ 11]*x12*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x21+coeff[  8]*x41+2*coeff[  9]*x31*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  8]*x31+coeff[  9]*x21*x32+coeff[ 11]*x12*x21,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 11;
    float avdat= -0.1388697E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
        -0.68234401E-02,-0.37446914E+01,-0.37926212E+02, 0.61507764E+01,
         0.63397957E+02,-0.14309901E+00, 0.12391151E+00, 0.54696472E-02,
        -0.23103161E-01, 0.41335378E-01, 0.23623554E-01,
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

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11    *x31    
        +coeff[  6]        *x32    
        +coeff[  7]    *x22        
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]    *x21    *x41
        +coeff[  9]*x12            
        +coeff[ 10]            *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x31+2*coeff[  9]*x11,2)
        +dx2*dx2*pow(0+coeff[  2]+2*coeff[  7]*x21+coeff[  8]*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x11+2*coeff[  6]*x31,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  8]*x21+2*coeff[ 10]*x41,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1789744E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.25880231E-05,-0.21781974E-01,-0.53585429E-01, 0.33663522E-01,
         0.63088670E-01,-0.73855968E-04, 0.11193584E-03,-0.17975667E-03,
        -0.14751505E-03, 0.11667734E-03, 0.10845865E-03, 0.68372763E-04,
         0.19833979E-04, 0.48111880E-04, 0.57281322E-04,-0.85413803E-05,
        -0.11287080E-03,-0.34575032E-05,-0.43495206E-05, 0.48359789E-05,
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
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;

//                 function

    float v_target_theta                            =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11*x21        
        +coeff[  6]*x11        *x41
        +coeff[  7]        *x31*x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x11    *x31    
        +coeff[  9]    *x21*x31    
        +coeff[ 10]        *x32    
        +coeff[ 11]        *x33    
        +coeff[ 12]        *x31*x42
        +coeff[ 13]*x12            
        +coeff[ 14]*x13            
        +coeff[ 15]*x11*x22        
        +coeff[ 16]*x12    *x31    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x21*x32    
        +coeff[ 18]            *x43
        +coeff[ 19]    *x22*x34    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41+coeff[  8]*x31+2*coeff[ 13]*x11+3*coeff[ 14]*x12+coeff[ 15]*x22+2*coeff[ 16]*x11*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  9]*x31+2*coeff[ 15]*x21*x11+coeff[ 17]*x32+2*coeff[ 19]*x21*x34,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x41+coeff[  8]*x11+coeff[  9]*x21+2*coeff[ 10]*x31+3*coeff[ 11]*x32+coeff[ 12]*x42+coeff[ 16]*x12+2*coeff[ 17]*x31*x21+4*coeff[ 19]*x33*x22,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  7]*x31+2*coeff[ 12]*x41*x31+3*coeff[ 18]*x42,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.9168636E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.71171430E-05,-0.53578906E-01, 0.22040226E-01, 0.63054502E-01,
        -0.34107450E-01, 0.18342525E-03,-0.15897595E-03, 0.14934510E-03,
         0.21824631E-04,-0.29373976E-04,-0.77709601E-05, 0.18234106E-05,
        -0.64094371E-06,-0.31095728E-05,-0.18853916E-04,-0.52884468E-04,
        -0.13802407E-03, 0.24421263E-03,-0.27390200E-03, 0.78496250E-05,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;

//                 function

    float v_target_phi                              =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11    *x31    
        +coeff[  6]        *x32    
        +coeff[  7]*x11        *x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]    *x21    *x41
        +coeff[  9]            *x42
        +coeff[ 10]        *x32*x41
        +coeff[ 11]*x13*x21        
        +coeff[ 12]*x11*x23        
        +coeff[ 13]*x13        *x41
        +coeff[ 14]        *x32*x43
        +coeff[ 15]*x12            
        +coeff[ 16]*x11*x21        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]    *x21*x31    
        +coeff[ 18]        *x31*x41
        +coeff[ 19]*x12*x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x31+coeff[  7]*x41+3*coeff[ 11]*x12*x21+coeff[ 12]*x23+3*coeff[ 13]*x12*x41+2*coeff[ 15]*x11+coeff[ 16]*x21+2*coeff[ 19]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  8]*x41+coeff[ 11]*x13+3*coeff[ 12]*x22*x11+coeff[ 16]*x11+coeff[ 17]*x31+coeff[ 19]*x12,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x11+2*coeff[  6]*x31+2*coeff[ 10]*x31*x41+2*coeff[ 14]*x31*x43+coeff[ 17]*x21+coeff[ 18]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+coeff[  8]*x21+2*coeff[  9]*x41+coeff[ 10]*x32+coeff[ 13]*x13+3*coeff[ 14]*x42*x32+coeff[ 18]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
