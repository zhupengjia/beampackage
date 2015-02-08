#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 11;
    float avdat=  0.2201240E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
         0.67881010E-04,-0.39964260E+02, 0.66041382E+02, 0.63049191E+00,
         0.19264846E+01,-0.38750049E-01, 0.64423032E-01,-0.11569226E+01,
        -0.10513455E+01, 0.62664226E-02,-0.26887422E-02,
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
        +coeff[  4]        *x31*x41
        +coeff[  5]    *x21        
        +coeff[  6]            *x41
        +coeff[  7]    *x21*x31    
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]*x11        *x41
        +coeff[  9]        *x31*x42
        +coeff[ 10]*x11*x22        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  8]*x41+coeff[ 10]*x22,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]+coeff[  7]*x31+2*coeff[ 10]*x21*x11,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x41+coeff[  7]*x21+coeff[  9]*x42,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  6]+coeff[  8]*x11+2*coeff[  9]*x41*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 16;
    float avdat=  0.3932708E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 17]={
        -0.23244059E-01,-0.39801453E+02, 0.65756973E+02,-0.80524737E+00,
         0.22054447E+00,-0.70372802E+00, 0.72883314E+00,-0.53798401E-03,
         0.32382429E-01,-0.47571059E-01,-0.31098703E+00, 0.95505810E+00,
        -0.14992652E-01, 0.21124313E-01, 0.24860302E-01,-0.34966238E-01,
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
    float x41 = x4;
    float x42 = x41*x4;

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]    *x21    *x41
        +coeff[  4]    *x22        
        +coeff[  5]        *x32    
        +coeff[  6]            *x42
        +coeff[  7]*x14            
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x11            
        +coeff[  9]        *x31    
        +coeff[ 10]*x12            
        +coeff[ 11]*x11    *x31    
        +coeff[ 12]*x11*x21*x31    
        +coeff[ 13]    *x21*x32    
        +coeff[ 14]*x11    *x31*x41
        +coeff[ 15]        *x32*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  7]*x13+coeff[  8]+2*coeff[ 10]*x11+coeff[ 11]*x31+coeff[ 12]*x21*x31+coeff[ 14]*x31*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  4]*x21+coeff[ 12]*x11*x31+coeff[ 13]*x32,2)
        +dx3*dx3*pow(0+2*coeff[  5]*x31+coeff[  9]+coeff[ 11]*x11+coeff[ 12]*x11*x21+2*coeff[ 13]*x31*x21+coeff[ 14]*x11*x41+2*coeff[ 15]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  6]*x41+coeff[ 14]*x11*x31+coeff[ 15]*x32,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.4882914E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.49076270E-05, 0.43599964E-02, 0.27224857E-02,-0.18256271E-02,
        -0.37071882E-02, 0.26175713E-02,-0.93309544E-01,-0.11176981E+01,
         0.10767224E+00, 0.38983159E+01,-0.46311903E+01, 0.18576928E+01,
         0.30259892E-01, 0.14241798E+01,-0.34917053E-01,-0.50472488E+01,
         0.45640618E-01, 0.59997973E+01,-0.44785324E-01,-0.23888776E+01,
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
        +coeff[  5]    *x21    *x41
        +coeff[  6]*x12*x21        
        +coeff[  7]    *x23        
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x12        *x41
        +coeff[  9]    *x22    *x41
        +coeff[ 10]    *x21    *x42
        +coeff[ 11]            *x43
        +coeff[ 12]*x14*x21        
        +coeff[ 13]*x12*x23        
        +coeff[ 14]*x14        *x41
        +coeff[ 15]*x12*x22    *x41
        +coeff[ 16]    *x24    *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12*x21    *x42
        +coeff[ 18]    *x23    *x42
        +coeff[ 19]*x12        *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  2]*x11+coeff[  4]*x31+2*coeff[  6]*x11*x21+2*coeff[  8]*x11*x41+4*coeff[ 12]*x13*x21+2*coeff[ 13]*x11*x23+4*coeff[ 14]*x13*x41+2*coeff[ 15]*x11*x22*x41+2*coeff[ 17]*x11*x21*x42+2*coeff[ 19]*x11*x43,2)
        +dx2*dx2*pow(0+2*coeff[  3]*x21+coeff[  5]*x41+coeff[  6]*x12+3*coeff[  7]*x22+2*coeff[  9]*x21*x41+coeff[ 10]*x42+coeff[ 12]*x14+3*coeff[ 13]*x22*x12+2*coeff[ 15]*x21*x12*x41+4*coeff[ 16]*x23*x41+coeff[ 17]*x12*x42+3*coeff[ 18]*x22*x42,2)
        +dx3*dx3*pow(0+coeff[  4]*x11,2)
        +dx4*dx4*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  8]*x12+coeff[  9]*x22+2*coeff[ 10]*x41*x21+3*coeff[ 11]*x42+coeff[ 14]*x14+coeff[ 15]*x12*x22+coeff[ 16]*x24+2*coeff[ 17]*x41*x12*x21+2*coeff[ 18]*x41*x23+3*coeff[ 19]*x42*x12,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.2307011E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.57662498E-06,-0.57270821E-01, 0.69199115E-01, 0.38127187E-02,
        -0.69424794E-02,-0.64060888E-02, 0.11626129E-01, 0.38729355E-03,
        -0.22974081E-03, 0.60425013E-04,-0.25433019E-04, 0.11185611E-03,
        -0.11575720E-03, 0.21594441E-03,-0.64206620E-05,-0.38888684E-04,
         0.29305363E-03,-0.41932668E-03,-0.55507771E-04, 0.17474817E-05,
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
        +coeff[  9]        *x31*x42
        +coeff[ 10]*x11*x22        
        +coeff[ 11]    *x21*x33    
        +coeff[ 12]    *x21*x31*x42
        +coeff[ 13]        *x31*x43
        +coeff[ 14]        *x32*x41
        +coeff[ 15]*x13*x21        
        +coeff[ 16]*x11    *x32*x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]        *x33*x41
        +coeff[ 18]*x11        *x43
        +coeff[ 19]    *x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[ 10]*x22+3*coeff[ 15]*x12*x21+coeff[ 16]*x32*x41+coeff[ 18]*x43,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  8]+2*coeff[ 10]*x21*x11+coeff[ 11]*x33+coeff[ 12]*x31*x42+coeff[ 15]*x13+coeff[ 19]*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x42+3*coeff[ 11]*x32*x21+coeff[ 12]*x21*x42+coeff[ 13]*x43+2*coeff[ 14]*x31*x41+2*coeff[ 16]*x31*x11*x41+3*coeff[ 17]*x32*x41,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  7]+2*coeff[  9]*x41*x31+2*coeff[ 12]*x41*x21*x31+3*coeff[ 13]*x42*x31+coeff[ 14]*x32+coeff[ 16]*x11*x32+coeff[ 17]*x33+3*coeff[ 18]*x42*x11+coeff[ 19]*x21,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
