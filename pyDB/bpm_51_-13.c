#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  7;
    float avdat=  0.3007177E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[  8]={
        -0.46120376E-05,-0.38581341E+02, 0.64340477E+02, 0.39774841E+00,
         0.12053976E+01,-0.72694546E+00,-0.65956688E+00,
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
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  6]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]*x31,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x41+coeff[  5]*x21,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  6]*x11,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 10;
    float avdat=  0.2411373E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 11]={
        -0.18069027E-01,-0.38414059E+02, 0.64050293E+02,-0.50189286E+00,
         0.45033842E+00,-0.17961334E+00, 0.13912536E+00, 0.52676702E+00,
        -0.36718246E+00,-0.11600584E-02,
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
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]    *x21    *x41
        +coeff[  4]            *x42
        +coeff[  5]*x12            
        +coeff[  6]    *x22        
        +coeff[  7]*x11    *x31    
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]        *x32    
        +coeff[  9]    *x21*x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  5]*x11+coeff[  7]*x31,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  6]*x21+coeff[  9]*x32,2)
        +dx3*dx3*pow(0+coeff[  7]*x11+2*coeff[  8]*x31+2*coeff[  9]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.6651161E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.44088309E-04,-0.55288166E-01, 0.65972432E-01, 0.80440976E-02,
        -0.64574825E-02,-0.32066917E-02, 0.78165968E-03,-0.24847121E-02,
         0.31041615E-02,-0.55268680E-04, 0.35207657E-04,-0.92519047E-04,
         0.13428133E-03, 0.19687651E-03,-0.28562266E-03,-0.48116766E-04,
         0.27202480E-03,-0.29600170E-03,-0.95545271E-04, 0.13568929E-03,
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
    float avdat=  0.4717440E-05;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.22136508E-08,-0.57051171E-01, 0.68739526E-01, 0.29036242E-02,
        -0.52484390E-02,-0.48699905E-02, 0.87702405E-02, 0.48985377E-04,
        -0.29563525E-04,-0.52388987E-04,-0.13547015E-03, 0.25870850E-04,
         0.94975148E-04,-0.12325608E-03,-0.46532750E-03, 0.32356987E-03,
        -0.47052356E-04, 0.16514871E-03, 0.16880434E-03, 0.54728098E-04,
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
