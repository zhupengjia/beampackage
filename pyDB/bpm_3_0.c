#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 13;
    float avdat=  0.1223386E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 14]={
        -0.17245908E-03,-0.38655281E+02, 0.40660124E+01, 0.64255653E+02,
        -0.66478233E+01,-0.38616095E-01, 0.57689130E-01, 0.69426775E-01,
        -0.10496221E+00,-0.53779124E-02,-0.42122817E-02, 0.35900676E-02,
         0.23259381E-02,
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
    float x41 = x4;
    float x42 = x41*x4;

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
        +coeff[  9]*x11    *x31    
        +coeff[ 10]            *x42
        +coeff[ 11]*x12            
        +coeff[ 12]    *x22        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41+coeff[  9]*x31+2*coeff[ 11]*x11,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  7]*x31+2*coeff[ 12]*x21,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x21+coeff[  8]*x41+coeff[  9]*x11,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  8]*x31+2*coeff[ 10]*x41,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 10;
    float avdat= -0.1153302E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 11]={
        -0.65843132E-02,-0.40461173E+01,-0.38654888E+02, 0.66150737E+01,
         0.64255150E+02,-0.14485626E+00, 0.12536602E+00,-0.10695244E-01,
         0.16007189E-01, 0.41816548E-01,
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
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x31+2*coeff[  9]*x11,2)
        +dx2*dx2*pow(0+coeff[  2]+2*coeff[  7]*x21+coeff[  8]*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x11+2*coeff[  6]*x31,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  8]*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1664320E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.27542110E-05,-0.22410274E-01,-0.53318199E-01, 0.34402438E-01,
         0.62674925E-01,-0.75840151E-04, 0.11382624E-03,-0.18420913E-03,
        -0.14473253E-03, 0.12075407E-03, 0.10923740E-03,-0.16198048E-03,
         0.14790794E-03,-0.10685354E-04,-0.79409974E-05, 0.46493591E-04,
         0.30255385E-04,-0.34543315E-04,-0.35140465E-05, 0.67342073E-04,
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
        +coeff[ 11]*x11    *x32    
        +coeff[ 12]        *x33    
        +coeff[ 13]*x11*x21    *x41
        +coeff[ 14]*x11        *x42
        +coeff[ 15]*x12            
        +coeff[ 16]*x13            
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x22*x31    
        +coeff[ 18]    *x21*x32    
        +coeff[ 19]    *x21*x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41+coeff[  8]*x31+coeff[ 11]*x32+coeff[ 13]*x21*x41+coeff[ 14]*x42+2*coeff[ 15]*x11+3*coeff[ 16]*x12,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  9]*x31+coeff[ 13]*x11*x41+2*coeff[ 17]*x21*x31+coeff[ 18]*x32+coeff[ 19]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x41+coeff[  8]*x11+coeff[  9]*x21+2*coeff[ 10]*x31+2*coeff[ 11]*x31*x11+3*coeff[ 12]*x32+coeff[ 17]*x22+2*coeff[ 18]*x31*x21+coeff[ 19]*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  7]*x31+coeff[ 13]*x11*x21+2*coeff[ 14]*x41*x11+coeff[ 19]*x21*x31,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.1120875E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.68707927E-05,-0.53315889E-01, 0.22664212E-01, 0.62649786E-01,
        -0.34839489E-01, 0.18685372E-03,-0.16131722E-03,-0.78612618E-04,
        -0.14436597E-04, 0.91251168E-05, 0.40244053E-04,-0.53106280E-04,
         0.13086965E-05,-0.48568381E-06,-0.22450038E-05,-0.22883896E-04,
        -0.54063770E-04, 0.50871720E-04, 0.79825895E-05, 0.17157827E-04,
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
    float x34 = x33*x3;
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
        +coeff[  8]            *x42
        +coeff[  9]    *x21*x32    
        +coeff[ 10]    *x21    *x42
        +coeff[ 11]            *x43
        +coeff[ 12]*x13*x21        
        +coeff[ 13]*x11*x23        
        +coeff[ 14]*x13        *x41
        +coeff[ 15]        *x34*x41
        +coeff[ 16]*x12            
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x11*x21        
        +coeff[ 18]    *x22        
        +coeff[ 19]    *x21*x31    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x31+coeff[  7]*x41+3*coeff[ 12]*x12*x21+coeff[ 13]*x23+3*coeff[ 14]*x12*x41+2*coeff[ 16]*x11+coeff[ 17]*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  9]*x32+coeff[ 10]*x42+coeff[ 12]*x13+3*coeff[ 13]*x22*x11+coeff[ 17]*x11+2*coeff[ 18]*x21+coeff[ 19]*x31,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x11+2*coeff[  6]*x31+2*coeff[  9]*x31*x21+4*coeff[ 15]*x33*x41+coeff[ 19]*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+2*coeff[  8]*x41+2*coeff[ 10]*x41*x21+3*coeff[ 11]*x42+coeff[ 14]*x13+coeff[ 15]*x34,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
