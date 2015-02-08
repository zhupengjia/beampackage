#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.2168621E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
         0.68857335E-04,-0.39155697E+02, 0.65064911E+02,-0.10592775E+01,
         0.17628694E+01, 0.57712805E+00,-0.96148050E+00,-0.35512444E-01,
         0.58989104E-01,
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
        +coeff[  3]    *x21*x31    
        +coeff[  4]        *x31*x41
        +coeff[  5]*x11*x21        
        +coeff[  6]*x11        *x41
        +coeff[  7]    *x21        
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]            *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x31+coeff[  5]*x11+coeff[  7],2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  3]*x21+coeff[  4]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  6]*x11+coeff[  8],2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 16;
    float avdat=  0.3957247E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 17]={
        -0.24037272E-01,-0.39005726E+02, 0.64802704E+02, 0.66658729E+00,
         0.20311919E+00, 0.27329957E-01,-0.39420936E-01,-0.58310479E+00,
        -0.73849553E+00, 0.14566478E-02,-0.45085372E-02, 0.25216362E-02,
        -0.26587474E+00, 0.80660069E+00, 0.46335957E-02,-0.64587002E-02,
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
        +coeff[  3]            *x42
        +coeff[  4]    *x22        
        +coeff[  5]*x11            
        +coeff[  6]        *x31    
        +coeff[  7]        *x32    
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]    *x21    *x41
        +coeff[  9]*x14            
        +coeff[ 10]*x13    *x31    
        +coeff[ 11]*x12    *x32    
        +coeff[ 12]*x12            
        +coeff[ 13]*x11    *x31    
        +coeff[ 14]*x11*x21*x31    
        +coeff[ 15]    *x21*x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  5]+4*coeff[  9]*x13+3*coeff[ 10]*x12*x31+2*coeff[ 11]*x11*x32+2*coeff[ 12]*x11+coeff[ 13]*x31+coeff[ 14]*x21*x31,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  8]*x41+coeff[ 14]*x11*x31+coeff[ 15]*x32,2)
        +dx3*dx3*pow(0+coeff[  6]+2*coeff[  7]*x31+coeff[ 10]*x13+2*coeff[ 11]*x31*x12+coeff[ 13]*x11+coeff[ 14]*x11*x21+2*coeff[ 15]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+2*coeff[  3]*x41+coeff[  8]*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1411881E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.47373142E-05, 0.43588234E-02, 0.27191939E-02,-0.18280419E-02,
        -0.37026971E-02, 0.26203997E-02,-0.93293555E-01,-0.11177913E+01,
         0.10765440E+00, 0.38986669E+01,-0.46316066E+01, 0.18578516E+01,
         0.30255878E-01, 0.14239165E+01,-0.34911972E-01,-0.50463147E+01,
         0.45635439E-01, 0.59986792E+01,-0.44780765E-01,-0.23884270E+01,
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
    float avdat=  0.2306944E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.68881825E-07,-0.57274848E-01, 0.69202870E-01, 0.38690062E-02,
        -0.70175147E-02,-0.65140794E-02, 0.11760153E-01, 0.38778083E-03,
        -0.23018305E-03, 0.26180225E-04, 0.11273946E-03, 0.14055864E-03,
        -0.17463984E-04,-0.68420736E-05, 0.27533701E-04,-0.40552914E-04,
         0.10284865E-03,-0.21352051E-03, 0.30051859E-03,-0.42662484E-03,
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
    float x24 = x23*x2;
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
        +coeff[  7]            *x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]    *x21        
        +coeff[  9]        *x31*x44
        +coeff[ 10]    *x21*x33    
        +coeff[ 11]*x11        *x43
        +coeff[ 12]*x11*x24        
        +coeff[ 13]        *x32*x41
        +coeff[ 14]        *x31*x42
        +coeff[ 15]*x13*x21        
        +coeff[ 16]*x11*x23        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x11*x22    *x41
        +coeff[ 18]*x11    *x32*x41
        +coeff[ 19]        *x33*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[ 11]*x43+coeff[ 12]*x24+3*coeff[ 15]*x12*x21+coeff[ 16]*x23+coeff[ 17]*x22*x41+coeff[ 18]*x32*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  8]+coeff[ 10]*x33+4*coeff[ 12]*x23*x11+coeff[ 15]*x13+3*coeff[ 16]*x22*x11+2*coeff[ 17]*x21*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x44+3*coeff[ 10]*x32*x21+2*coeff[ 13]*x31*x41+coeff[ 14]*x42+2*coeff[ 18]*x31*x11*x41+3*coeff[ 19]*x32*x41,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  7]+4*coeff[  9]*x43*x31+3*coeff[ 11]*x42*x11+coeff[ 13]*x32+2*coeff[ 14]*x41*x31+coeff[ 17]*x11*x22+coeff[ 18]*x11*x32+coeff[ 19]*x33,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
