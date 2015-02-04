#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 10;
    float avdat=  0.1191767E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 11]={
         0.33222139E-04,-0.39416534E+02, 0.65456841E+02, 0.86727107E+00,
        -0.14435229E+01,-0.15901225E+01, 0.26438158E+01,-0.29742345E-01,
         0.49378525E-01, 0.56314259E-02,
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
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;

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
        +coeff[  9]        *x31*x44
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  4]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]*x31+coeff[  7],2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  5]*x21+coeff[  6]*x41+coeff[  9]*x44,2)
        +dx4*dx4*pow(0+coeff[  4]*x11+coeff[  6]*x31+coeff[  8]+4*coeff[  9]*x43*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4782038E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.34589302E-01,-0.39092152E+02, 0.64887482E+02, 0.30870292E+00,
        -0.11154264E+01,-0.88635027E+00, 0.10022482E+01,-0.93715941E-03,
         0.25946356E-03,-0.34537494E-01,-0.40457165E+00, 0.12258588E+01,
         0.24603691E-01,-0.36185205E-01, 0.48397776E-01,-0.11011421E+00,
        -0.24892015E-01,-0.79524554E-02, 0.11463103E+00, 0.10580366E-01,
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
        +coeff[  9]        *x31    
        +coeff[ 10]*x12            
        +coeff[ 11]*x11    *x31    
        +coeff[ 12]*x11            
        +coeff[ 13]*x11*x21*x31    
        +coeff[ 14]    *x21*x32    
        +coeff[ 15]        *x32*x41
        +coeff[ 16]*x12        *x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]    *x22    *x41
        +coeff[ 18]*x11    *x31*x41
        +coeff[ 19]    *x21    *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  7]*x13+2*coeff[ 10]*x11+coeff[ 11]*x31+coeff[ 12]+coeff[ 13]*x21*x31+2*coeff[ 16]*x11*x41+coeff[ 18]*x31*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  3]*x21+coeff[  4]*x41+coeff[ 13]*x11*x31+coeff[ 14]*x32+2*coeff[ 17]*x21*x41+coeff[ 19]*x42,2)
        +dx3*dx3*pow(0+2*coeff[  5]*x31+4*coeff[  8]*x33+coeff[  9]+coeff[ 11]*x11+coeff[ 13]*x11*x21+2*coeff[ 14]*x31*x21+2*coeff[ 15]*x31*x41+coeff[ 18]*x11*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  4]*x21+2*coeff[  6]*x41+coeff[ 15]*x32+coeff[ 16]*x12+coeff[ 17]*x22+coeff[ 18]*x11*x31+2*coeff[ 19]*x41*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1004010E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.17238733E-03, 0.64394772E-02,-0.49919263E-02,-0.27214340E-02,
         0.16223796E-01,-0.13074216E-01, 0.39000700E-02,-0.75487830E-01,
        -0.75713021E+00, 0.84095903E-01, 0.25366268E+01,-0.29726195E+01,
         0.11964098E+01, 0.45413155E-01, 0.64331770E-01,-0.52402634E-01,
        -0.10982844E+00, 0.68441957E-01, 0.43932419E-01,-0.67155078E-01,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;

//                 function

    float v_target_theta                            =avdat
        +coeff[  0]                
        +coeff[  1]            *x41
        +coeff[  2]*x12            
        +coeff[  3]    *x22        
        +coeff[  4]*x11    *x31    
        +coeff[  5]        *x32    
        +coeff[  6]    *x21    *x41
        +coeff[  7]*x12*x21        
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]    *x23        
        +coeff[  9]*x12        *x41
        +coeff[ 10]    *x22    *x41
        +coeff[ 11]    *x21    *x42
        +coeff[ 12]            *x43
        +coeff[ 13]*x14*x21        
        +coeff[ 14]*x12*x23        
        +coeff[ 15]*x14        *x41
        +coeff[ 16]*x12*x22    *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x24    *x41
        +coeff[ 18]*x12*x21    *x42
        +coeff[ 19]    *x23    *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  2]*x11+coeff[  4]*x31+2*coeff[  7]*x11*x21+2*coeff[  9]*x11*x41+4*coeff[ 13]*x13*x21+2*coeff[ 14]*x11*x23+4*coeff[ 15]*x13*x41+2*coeff[ 16]*x11*x22*x41+2*coeff[ 18]*x11*x21*x42,2)
        +dx2*dx2*pow(0+2*coeff[  3]*x21+coeff[  6]*x41+coeff[  7]*x12+3*coeff[  8]*x22+2*coeff[ 10]*x21*x41+coeff[ 11]*x42+coeff[ 13]*x14+3*coeff[ 14]*x22*x12+2*coeff[ 16]*x21*x12*x41+4*coeff[ 17]*x23*x41+coeff[ 18]*x12*x42+3*coeff[ 19]*x22*x42,2)
        +dx3*dx3*pow(0+coeff[  4]*x11+2*coeff[  5]*x31,2)
        +dx4*dx4*pow(0+coeff[  1]+coeff[  6]*x21+coeff[  9]*x12+coeff[ 10]*x22+2*coeff[ 11]*x41*x21+3*coeff[ 12]*x42+coeff[ 15]*x14+coeff[ 16]*x12*x22+coeff[ 17]*x24+2*coeff[ 18]*x41*x12*x21+2*coeff[ 19]*x41*x23,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1182557E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.30169936E-07,-0.58040299E-01, 0.70596866E-01, 0.56541124E-02,
        -0.10317817E-01,-0.95822429E-02, 0.17370474E-01,-0.19255868E-03,
         0.32192009E-03, 0.12300612E-04,-0.94753799E-04,-0.57728703E-05,
         0.79150207E-03, 0.34391916E-04,-0.14829087E-03, 0.19474704E-03,
         0.15244866E-03, 0.24156134E-03,-0.76907064E-03,-0.33765318E-03,
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
    float x43 = x42*x4;
    float x44 = x43*x4;

//                 function

    float v_target_phi                              =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]        *x31    
        +coeff[  3]*x11*x21        
        +coeff[  4]    *x21*x31    
        +coeff[  5]*x11        *x41
        +coeff[  6]        *x31*x41
        +coeff[  7]    *x21        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]            *x41
        +coeff[  9]        *x31*x44
        +coeff[ 10]        *x33*x41
        +coeff[ 11]*x11*x23    *x41
        +coeff[ 12]        *x31*x42
        +coeff[ 13]*x13*x21        
        +coeff[ 14]    *x22*x31*x41
        +coeff[ 15]    *x21*x31*x42
        +coeff[ 16]    *x22*x31    
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x11*x21    *x41
        +coeff[ 18]    *x21*x31*x41
        +coeff[ 19]*x11        *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[ 11]*x23*x41+3*coeff[ 13]*x12*x21+coeff[ 17]*x21*x41+coeff[ 19]*x42,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  7]+3*coeff[ 11]*x22*x11*x41+coeff[ 13]*x13+2*coeff[ 14]*x21*x31*x41+coeff[ 15]*x31*x42+2*coeff[ 16]*x21*x31+coeff[ 17]*x11*x41+coeff[ 18]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x44+3*coeff[ 10]*x32*x41+coeff[ 12]*x42+coeff[ 14]*x22*x41+coeff[ 15]*x21*x42+coeff[ 16]*x22+coeff[ 18]*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  8]+4*coeff[  9]*x43*x31+coeff[ 10]*x33+coeff[ 11]*x11*x23+2*coeff[ 12]*x41*x31+coeff[ 14]*x22*x31+2*coeff[ 15]*x41*x21*x31+coeff[ 17]*x11*x21+coeff[ 18]*x21*x31+2*coeff[ 19]*x41*x11,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
