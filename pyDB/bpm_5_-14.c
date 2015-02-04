#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.1773937E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
         0.11646981E-04,-0.38575089E+02, 0.64334366E+02, 0.39419940E+00,
         0.11983566E+01,-0.21602634E-01, 0.35868648E-01,-0.72199035E+00,
        -0.65469795E+00,
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
        +coeff[  5]    *x21        
        +coeff[  6]            *x41
        +coeff[  7]    *x21*x31    
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]*x11        *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  8]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]+coeff[  7]*x31,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x41+coeff[  7]*x21,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  6]+coeff[  8]*x11,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 12;
    float avdat=  0.2615463E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 13]={
        -0.17202731E-01,-0.38405807E+02, 0.64043571E+02,-0.51218194E+00,
         0.45607549E+00,-0.18224409E+00, 0.14367910E+00, 0.53213865E+00,
        -0.36968824E+00,-0.16505003E-01, 0.12754231E-01,-0.11608150E-02,
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
        +coeff[  9]        *x31    
        +coeff[ 10]*x11            
        +coeff[ 11]    *x21*x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  5]*x11+coeff[  7]*x31+coeff[ 10],2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  6]*x21+coeff[ 11]*x32,2)
        +dx3*dx3*pow(0+coeff[  7]*x11+2*coeff[  8]*x31+coeff[  9]+2*coeff[ 11]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.6671510E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.44779401E-04,-0.55326767E-01, 0.66030882E-01, 0.80072172E-02,
        -0.64353873E-02,-0.32284202E-02, 0.79085230E-03, 0.21381651E-03,
        -0.34172161E-03,-0.24692789E-02, 0.31165371E-02,-0.10573350E-03,
         0.15218141E-03, 0.20170306E-03,-0.29046493E-03,-0.51404983E-04,
         0.29568883E-03,-0.32255700E-03,-0.99929821E-04, 0.14051344E-03,
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
        +coeff[  7]*x11            
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]        *x31    
        +coeff[  9]*x12            
        +coeff[ 10]            *x42
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
        +dx1*dx1*pow(0+coeff[  3]*x31+coeff[  7]+2*coeff[  9]*x11+coeff[ 11]*x21*x31+2*coeff[ 15]*x11*x41+coeff[ 16]*x31*x41+2*coeff[ 18]*x11*x22+2*coeff[ 19]*x11*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  5]*x41+2*coeff[  6]*x21+coeff[ 11]*x11*x31+coeff[ 12]*x32+coeff[ 13]*x32*x41+2*coeff[ 18]*x21*x12+coeff[ 19]*x12*x41,2)
        +dx3*dx3*pow(0+coeff[  3]*x11+2*coeff[  4]*x31+coeff[  8]+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+2*coeff[ 13]*x31*x21*x41+2*coeff[ 14]*x31*x42+coeff[ 16]*x11*x41+2*coeff[ 17]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  5]*x21+2*coeff[ 10]*x41+coeff[ 13]*x21*x32+2*coeff[ 14]*x41*x32+coeff[ 15]*x12+coeff[ 16]*x11*x31+coeff[ 17]*x32+coeff[ 19]*x12*x21,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1158103E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.98307371E-08,-0.57091441E-01, 0.68797648E-01, 0.28572998E-02,
        -0.52101314E-02,-0.48145559E-02, 0.87269414E-02,-0.15381153E-03,
         0.25757897E-03,-0.57707704E-04,-0.91345719E-04, 0.29423754E-04,
        -0.12371347E-03, 0.17624837E-03,-0.37085010E-04,-0.68630681E-04,
         0.89999914E-04, 0.86580940E-05,-0.12805433E-04, 0.90709582E-04,
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
        +coeff[  9]*x11        *x42
        +coeff[ 10]        *x33*x41
        +coeff[ 11]*x11*x22        
        +coeff[ 12]    *x21*x31*x41
        +coeff[ 13]        *x31*x42
        +coeff[ 14]*x13*x21        
        +coeff[ 15]    *x22*x31*x41
        +coeff[ 16]    *x21*x31*x42
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x11*x21*x31    
        +coeff[ 18]    *x21*x32    
        +coeff[ 19]*x12*x21*x31    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[  9]*x42+coeff[ 11]*x22+3*coeff[ 14]*x12*x21+coeff[ 17]*x21*x31+2*coeff[ 19]*x11*x21*x31,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  7]+2*coeff[ 11]*x21*x11+coeff[ 12]*x31*x41+coeff[ 14]*x13+2*coeff[ 15]*x21*x31*x41+coeff[ 16]*x31*x42+coeff[ 17]*x11*x31+coeff[ 18]*x32+coeff[ 19]*x12*x31,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+3*coeff[ 10]*x32*x41+coeff[ 12]*x21*x41+coeff[ 13]*x42+coeff[ 15]*x22*x41+coeff[ 16]*x21*x42+coeff[ 17]*x11*x21+2*coeff[ 18]*x31*x21+coeff[ 19]*x12*x21,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  8]+2*coeff[  9]*x41*x11+coeff[ 10]*x33+coeff[ 12]*x21*x31+2*coeff[ 13]*x41*x31+coeff[ 15]*x22*x31+2*coeff[ 16]*x41*x21*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
