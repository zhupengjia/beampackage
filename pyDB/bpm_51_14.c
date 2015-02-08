#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.3008425E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
         0.16295412E-05,-0.40164864E+02, 0.66248512E+02, 0.47588089E+00,
        -0.79104871E+00, 0.14446131E+01,-0.86985356E+00,-0.50191237E-02,
         0.83220787E-02,
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
        +coeff[  4]*x11        *x41
        +coeff[  5]        *x31*x41
        +coeff[  6]    *x21*x31    
        +coeff[  7]    *x21        
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]            *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  4]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  6]*x31+coeff[  7],2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  5]*x41+coeff[  6]*x21,2)
        +dx4*dx4*pow(0+coeff[  4]*x11+coeff[  5]*x31+coeff[  8],2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 16;
    float avdat=  0.4132241E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 17]={
        -0.16794538E-01,-0.39952530E+02, 0.65885239E+02, 0.16044335E+00,
        -0.58896416E+00, 0.89415599E-03,-0.70411747E-03, 0.53485501E+00,
        -0.25202361E+00, 0.76012760E+00,-0.55377364E+00, 0.53636413E-02,
        -0.19039864E-02,-0.66133770E-02,-0.69080773E-02, 0.90705641E-02,
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
    float x43 = x42*x4;

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]    *x22        
        +coeff[  4]    *x21    *x41
        +coeff[  5]        *x34    
        +coeff[  6]*x14    *x32    
        +coeff[  7]            *x42
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x12            
        +coeff[  9]*x11    *x31    
        +coeff[ 10]        *x32    
        +coeff[ 11]*x11*x21*x31    
        +coeff[ 12]    *x21*x32    
        +coeff[ 13]        *x32*x41
        +coeff[ 14]    *x21    *x42
        +coeff[ 15]            *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  6]*x13*x32+2*coeff[  8]*x11+coeff[  9]*x31+coeff[ 11]*x21*x31,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  3]*x21+coeff[  4]*x41+coeff[ 11]*x11*x31+coeff[ 12]*x32+coeff[ 14]*x42,2)
        +dx3*dx3*pow(0+4*coeff[  5]*x33+2*coeff[  6]*x31*x14+coeff[  9]*x11+2*coeff[ 10]*x31+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+2*coeff[ 13]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  4]*x21+2*coeff[  7]*x41+coeff[ 13]*x32+2*coeff[ 14]*x41*x21+3*coeff[ 15]*x42,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.5730039E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.44548062E-04,-0.55265389E-01, 0.65944292E-01, 0.80328537E-02,
        -0.64526871E-02,-0.31829779E-02, 0.76938095E-03,-0.24788182E-02,
         0.30933616E-02,-0.55291173E-04, 0.35207260E-04,-0.66684792E-04,
         0.11408583E-03,-0.13537292E-03, 0.29317543E-03,-0.40793908E-03,
        -0.16512869E-03, 0.56658148E-04,-0.26311542E-03, 0.48199052E-03,
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
        +coeff[ 13]        *x32*x41
        +coeff[ 14]    *x22*x32    
        +coeff[ 15]    *x21*x32*x41
        +coeff[ 16]*x12        *x42
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12        *x41
        +coeff[ 18]*x12*x22        
        +coeff[ 19]*x12*x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  3]*x31+2*coeff[  7]*x11+coeff[ 10]+coeff[ 11]*x21*x31+2*coeff[ 16]*x11*x42+2*coeff[ 17]*x11*x41+2*coeff[ 18]*x11*x22+2*coeff[ 19]*x11*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  5]*x41+2*coeff[  6]*x21+coeff[ 11]*x11*x31+coeff[ 12]*x32+2*coeff[ 14]*x21*x32+coeff[ 15]*x32*x41+2*coeff[ 18]*x21*x12+coeff[ 19]*x12*x41,2)
        +dx3*dx3*pow(0+coeff[  3]*x11+2*coeff[  4]*x31+coeff[  9]+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+2*coeff[ 13]*x31*x41+2*coeff[ 14]*x31*x22+2*coeff[ 15]*x31*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  5]*x21+2*coeff[  8]*x41+coeff[ 13]*x32+coeff[ 15]*x21*x32+2*coeff[ 16]*x41*x12+coeff[ 17]*x12+coeff[ 19]*x12*x21,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4721577E-05;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.22899476E-08,-0.57018682E-01, 0.68700247E-01, 0.28747006E-02,
        -0.52170772E-02,-0.48084971E-02, 0.87039536E-02, 0.49076083E-04,
        -0.29616836E-04,-0.11834480E-03, 0.17026623E-03,-0.12077423E-03,
         0.17133245E-03,-0.93920971E-05, 0.21334254E-03,-0.30437051E-03,
        -0.36757323E-04, 0.78788369E-04, 0.64007283E-04,-0.92326591E-04,
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
        +coeff[  9]    *x21*x31*x41
        +coeff[ 10]        *x31*x42
        +coeff[ 11]*x12*x21*x31    
        +coeff[ 12]*x11*x21*x32    
        +coeff[ 13]*x11*x22    *x41
        +coeff[ 14]*x11    *x32*x41
        +coeff[ 15]        *x33*x41
        +coeff[ 16]*x11*x21    *x42
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]        *x31*x43
        +coeff[ 18]*x11*x21    *x41
        +coeff[ 19]*x11        *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+2*coeff[ 11]*x11*x21*x31+coeff[ 12]*x21*x32+coeff[ 13]*x22*x41+coeff[ 14]*x32*x41+coeff[ 16]*x21*x42+coeff[ 18]*x21*x41+coeff[ 19]*x42,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  8]+coeff[  9]*x31*x41+coeff[ 11]*x12*x31+coeff[ 12]*x11*x32+2*coeff[ 13]*x21*x11*x41+coeff[ 16]*x11*x42+coeff[ 18]*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x21*x41+coeff[ 10]*x42+coeff[ 11]*x12*x21+2*coeff[ 12]*x31*x11*x21+2*coeff[ 14]*x31*x11*x41+3*coeff[ 15]*x32*x41+coeff[ 17]*x43,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  7]+coeff[  9]*x21*x31+2*coeff[ 10]*x41*x31+coeff[ 13]*x11*x22+coeff[ 14]*x11*x32+coeff[ 15]*x33+2*coeff[ 16]*x41*x11*x21+3*coeff[ 17]*x42*x31+coeff[ 18]*x11*x21+2*coeff[ 19]*x41*x11,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
