#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 12;
    float avdat=  0.9218322E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 13]={
         0.75835569E-05,-0.38779324E+02, 0.30936418E+01, 0.64461708E+02,
        -0.50582080E+01,-0.28448820E-01, 0.42489488E-01, 0.51747024E-01,
        -0.78203090E-01,-0.75933547E-03, 0.47631240E-02,-0.63708583E-02,
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
    float x32 = x31*x3;
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
        +coeff[  9]        *x32    
        +coeff[ 10]    *x21    *x41
        +coeff[ 11]            *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  7]*x31+coeff[ 10]*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x21+coeff[  8]*x41+2*coeff[  9]*x31,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  8]*x31+coeff[ 10]*x21+2*coeff[ 11]*x41,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 10;
    float avdat=  0.1076995E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 11]={
        -0.42522582E-02,-0.30777843E+01,-0.38781406E+02, 0.50321293E+01,
         0.64465561E+02, 0.84879726E-01, 0.25110701E-01, 0.12346289E-01,
        -0.82747024E-02,-0.93430981E-01,
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
        +coeff[  5]        *x32    
        +coeff[  6]*x12            
        +coeff[  7]    *x21    *x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]    *x22        
        +coeff[  9]*x11    *x31    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  6]*x11+coeff[  9]*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  7]*x41+2*coeff[  8]*x21,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  5]*x31+coeff[  9]*x11,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1834182E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.11935770E-05,-0.17095175E-01,-0.54636925E-01, 0.26253697E-01,
         0.64810492E-01,-0.40204144E-04, 0.62261184E-04, 0.73376141E-04,
         0.59741680E-04,-0.94394985E-04, 0.39012308E-03, 0.18443618E-04,
         0.45297392E-06, 0.27927925E-04,-0.93489696E-04,-0.14008139E-03,
        -0.79707224E-05, 0.58449100E-03,-0.82114100E-03,-0.17559512E-05,
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
        +coeff[  6]    *x21*x31    
        +coeff[  7]        *x32    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x11        *x41
        +coeff[  9]        *x31*x41
        +coeff[ 10]        *x33    
        +coeff[ 11]        *x31*x42
        +coeff[ 12]*x14            
        +coeff[ 13]*x12            
        +coeff[ 14]*x11    *x31    
        +coeff[ 15]*x13            
        +coeff[ 16]*x11*x22        
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12    *x31    
        +coeff[ 18]*x11    *x32    
        +coeff[ 19]    *x21*x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  8]*x41+4*coeff[ 12]*x13+2*coeff[ 13]*x11+coeff[ 14]*x31+3*coeff[ 15]*x12+coeff[ 16]*x22+2*coeff[ 17]*x11*x31+coeff[ 18]*x32,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  6]*x31+2*coeff[ 16]*x21*x11+coeff[ 19]*x32,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  6]*x21+2*coeff[  7]*x31+coeff[  9]*x41+3*coeff[ 10]*x32+coeff[ 11]*x42+coeff[ 14]*x11+coeff[ 17]*x12+2*coeff[ 18]*x31*x11+2*coeff[ 19]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  8]*x11+coeff[  9]*x31+2*coeff[ 11]*x41*x31,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.1661241E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.45297834E-05,-0.54645501E-01, 0.17302604E-01, 0.64803869E-01,
        -0.26602386E-01,-0.27588960E-04, 0.10025123E-03, 0.11744020E-03,
         0.20751288E-04,-0.27920953E-04,-0.18337256E-04,-0.81735852E-05,
         0.10851027E-05,-0.39030073E-06,-0.18578897E-05,-0.10789393E-03,
         0.19009109E-03,-0.89750945E-04,-0.21479420E-03, 0.78985058E-05,
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
        +coeff[  5]*x12            
        +coeff[  6]*x11    *x31    
        +coeff[  7]*x11        *x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]    *x21    *x41
        +coeff[  9]            *x42
        +coeff[ 10]        *x32*x41
        +coeff[ 11]            *x43
        +coeff[ 12]*x13*x21        
        +coeff[ 13]*x11*x23        
        +coeff[ 14]*x13        *x41
        +coeff[ 15]*x11*x21        
        +coeff[ 16]    *x21*x31    
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]        *x32    
        +coeff[ 18]        *x31*x41
        +coeff[ 19]*x12*x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  5]*x11+coeff[  6]*x31+coeff[  7]*x41+3*coeff[ 12]*x12*x21+coeff[ 13]*x23+3*coeff[ 14]*x12*x41+coeff[ 15]*x21+2*coeff[ 19]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  8]*x41+coeff[ 12]*x13+3*coeff[ 13]*x22*x11+coeff[ 15]*x11+coeff[ 16]*x31+coeff[ 19]*x12,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  6]*x11+2*coeff[ 10]*x31*x41+coeff[ 16]*x21+2*coeff[ 17]*x31+coeff[ 18]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+coeff[  8]*x21+2*coeff[  9]*x41+coeff[ 10]*x32+3*coeff[ 11]*x42+coeff[ 14]*x13+coeff[ 18]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
