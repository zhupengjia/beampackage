#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.1806671E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
         0.93715616E-05,-0.40188801E+02, 0.66279083E+02, 0.47418493E+00,
        -0.78925103E+00, 0.14428089E+01, 0.43076396E-01,-0.86806732E+00,
        -0.25906313E-01,
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
        +coeff[  6]            *x41
        +coeff[  7]    *x21*x31    
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]    *x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  4]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  7]*x31+coeff[  8],2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  5]*x41+coeff[  7]*x21,2)
        +dx4*dx4*pow(0+coeff[  4]*x11+coeff[  5]*x31+coeff[  6],2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 17;
    float avdat=  0.4371242E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 18]={
        -0.15930772E-01,-0.39976284E+02, 0.65916893E+02, 0.16754286E+00,
        -0.60648340E+00, 0.89390529E-03,-0.75810181E-03, 0.54550362E+00,
         0.18907074E-01,-0.26300058E-01,-0.25050047E+00, 0.75648826E+00,
        -0.55160445E+00, 0.53009270E-02,-0.15048308E-02,-0.70095873E-02,
         0.15583074E-02,
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
        +coeff[  8]*x11            
        +coeff[  9]        *x31    
        +coeff[ 10]*x12            
        +coeff[ 11]*x11    *x31    
        +coeff[ 12]        *x32    
        +coeff[ 13]*x11*x21*x31    
        +coeff[ 14]    *x21*x32    
        +coeff[ 15]        *x32*x41
        +coeff[ 16]            *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  6]*x13*x32+coeff[  8]+2*coeff[ 10]*x11+coeff[ 11]*x31+coeff[ 13]*x21*x31,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  3]*x21+coeff[  4]*x41+coeff[ 13]*x11*x31+coeff[ 14]*x32,2)
        +dx3*dx3*pow(0+4*coeff[  5]*x33+2*coeff[  6]*x31*x14+coeff[  9]+coeff[ 11]*x11+2*coeff[ 12]*x31+coeff[ 13]*x11*x21+2*coeff[ 14]*x31*x21+2*coeff[ 15]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  4]*x21+2*coeff[  7]*x41+coeff[ 15]*x32+3*coeff[ 16]*x42,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.5733608E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.46951602E-04,-0.55295590E-01, 0.65994442E-01, 0.80680642E-02,
        -0.64741704E-02,-0.29697241E-02, 0.68090885E-03, 0.21418027E-03,
        -0.34228229E-03,-0.24964300E-02, 0.29658927E-02,-0.12186050E-03,
         0.17129991E-03,-0.16980096E-03, 0.27021661E-03,-0.37957408E-03,
         0.11117624E-03,-0.31189770E-05, 0.83330830E-04,-0.54035314E-04,
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
        +coeff[ 11]*x12*x21        
        +coeff[ 12]*x11*x21*x31    
        +coeff[ 13]        *x32*x41
        +coeff[ 14]    *x22*x32    
        +coeff[ 15]    *x21*x32*x41
        +coeff[ 16]*x12        *x42
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]        *x31*x41
        +coeff[ 18]*x12        *x41
        +coeff[ 19]*x12*x22        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  3]*x31+coeff[  7]+2*coeff[  9]*x11+2*coeff[ 11]*x11*x21+coeff[ 12]*x21*x31+2*coeff[ 16]*x11*x42+2*coeff[ 18]*x11*x41+2*coeff[ 19]*x11*x22,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  5]*x41+2*coeff[  6]*x21+coeff[ 11]*x12+coeff[ 12]*x11*x31+2*coeff[ 14]*x21*x32+coeff[ 15]*x32*x41+2*coeff[ 19]*x21*x12,2)
        +dx3*dx3*pow(0+coeff[  3]*x11+2*coeff[  4]*x31+coeff[  8]+coeff[ 12]*x11*x21+2*coeff[ 13]*x31*x41+2*coeff[ 14]*x31*x22+2*coeff[ 15]*x31*x21*x41+coeff[ 17]*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  5]*x21+2*coeff[ 10]*x41+coeff[ 13]*x32+coeff[ 15]*x21*x32+2*coeff[ 16]*x41*x12+coeff[ 17]*x31+coeff[ 18]*x12,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1157476E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.94611536E-08,-0.57057861E-01, 0.68758003E-01, 0.28738573E-02,
        -0.52170106E-02,-0.48073600E-02, 0.87039396E-02,-0.15353673E-03,
         0.25746989E-03,-0.25942376E-04, 0.66342647E-04,-0.11897594E-03,
         0.16908266E-03,-0.93717226E-05, 0.20723784E-03,-0.29718652E-03,
        -0.36698006E-04, 0.78622630E-04,-0.29279843E-05,-0.17670296E-04,
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
        +coeff[  7]    *x21        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]            *x41
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
        +coeff[ 18]    *x21*x32    
        +coeff[ 19]*x11*x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+2*coeff[ 11]*x11*x21*x31+coeff[ 12]*x21*x32+coeff[ 13]*x22*x41+coeff[ 14]*x32*x41+coeff[ 16]*x21*x42+coeff[ 19]*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  7]+coeff[  9]*x31*x41+coeff[ 11]*x12*x31+coeff[ 12]*x11*x32+2*coeff[ 13]*x21*x11*x41+coeff[ 16]*x11*x42+coeff[ 18]*x32+coeff[ 19]*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x21*x41+coeff[ 10]*x42+coeff[ 11]*x12*x21+2*coeff[ 12]*x31*x11*x21+2*coeff[ 14]*x31*x11*x41+3*coeff[ 15]*x32*x41+coeff[ 17]*x43+2*coeff[ 18]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  8]+coeff[  9]*x21*x31+2*coeff[ 10]*x41*x31+coeff[ 13]*x11*x22+coeff[ 14]*x11*x32+coeff[ 15]*x33+2*coeff[ 16]*x41*x11*x21+3*coeff[ 17]*x42*x31+coeff[ 19]*x11*x21,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
