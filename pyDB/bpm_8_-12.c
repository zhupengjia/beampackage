#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  7;
    float avdat=  0.9417130E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[  8]={
        -0.50860294E-03,-0.38092712E+02, 0.28792775E+01, 0.63647144E+02,
        -0.47275791E+01, 0.11293771E-01,-0.22315288E-01,
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
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11*x21        
        +coeff[  6]        *x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  6]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 11;
    float avdat=  0.8412593E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
        -0.44745891E-02,-0.28658254E+01,-0.38094673E+02, 0.47056007E+01,
         0.63650642E+02, 0.84199153E-01, 0.24980927E-01,-0.20951441E-01,
         0.55171638E-02,-0.92755951E-01, 0.20084566E-01,
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
        +coeff[  5]        *x32    
        +coeff[  6]*x12            
        +coeff[  7]    *x21    *x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]    *x22        
        +coeff[  9]*x11    *x31    
        +coeff[ 10]            *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  6]*x11+coeff[  9]*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  7]*x41+2*coeff[  8]*x21,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  5]*x31+coeff[  9]*x11,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x21+2*coeff[ 10]*x41,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1922371E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.12011561E-05,-0.16649608E-01,-0.54781962E-01, 0.25727727E-01,
         0.65033726E-01,-0.22788810E-04, 0.11704432E-04, 0.75501426E-04,
        -0.11052881E-04,-0.14394645E-03, 0.11892091E-03, 0.38174760E-04,
        -0.10804214E-03, 0.44361994E-04,-0.17600076E-04,-0.15798520E-05,
         0.21850346E-04, 0.82153056E-05,-0.10936083E-04, 0.14312148E-04,
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
    float x43 = x42*x4;

//                 function

    float v_target_theta                            =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]        *x31*x41
        +coeff[  6]*x11*x21        
        +coeff[  7]        *x32    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x14            
        +coeff[  9]*x11*x21*x32*x41
        +coeff[ 10]        *x33*x42
        +coeff[ 11]*x12            
        +coeff[ 12]*x11    *x31    
        +coeff[ 13]*x12*x22*x31    
        +coeff[ 14]*x11    *x32    
        +coeff[ 15]    *x21*x32    
        +coeff[ 16]        *x33    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x21    *x42
        +coeff[ 18]            *x43
        +coeff[ 19]*x13    *x31    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  6]*x21+4*coeff[  8]*x13+coeff[  9]*x21*x32*x41+2*coeff[ 11]*x11+coeff[ 12]*x31+2*coeff[ 13]*x11*x22*x31+coeff[ 14]*x32+3*coeff[ 19]*x12*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  6]*x11+coeff[  9]*x11*x32*x41+2*coeff[ 13]*x21*x12*x31+coeff[ 15]*x32+coeff[ 17]*x42,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x41+2*coeff[  7]*x31+2*coeff[  9]*x31*x11*x21*x41+3*coeff[ 10]*x32*x42+coeff[ 12]*x11+coeff[ 13]*x12*x22+2*coeff[ 14]*x31*x11+2*coeff[ 15]*x31*x21+3*coeff[ 16]*x32+coeff[ 19]*x13,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  5]*x31+coeff[  9]*x11*x21*x32+2*coeff[ 10]*x41*x33+2*coeff[ 17]*x41*x21+3*coeff[ 18]*x42,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.1505907E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.44532935E-05,-0.54786190E-01, 0.16848592E-01, 0.65020286E-01,
        -0.26067747E-01,-0.26980359E-04, 0.98178825E-04, 0.11658196E-03,
         0.20702686E-04,-0.27733373E-04, 0.14769178E-05,-0.53532062E-06,
        -0.24610181E-05, 0.89546356E-05,-0.35132998E-04,-0.10654630E-03,
         0.18860119E-03,-0.88087501E-04,-0.21380052E-03, 0.83109444E-05,
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
        +coeff[ 10]*x13*x21        
        +coeff[ 11]*x11*x23        
        +coeff[ 12]*x13        *x41
        +coeff[ 13]    *x22*x32*x41
        +coeff[ 14]        *x32*x43
        +coeff[ 15]*x11*x21        
        +coeff[ 16]    *x21*x31    
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]        *x32    
        +coeff[ 18]        *x31*x41
        +coeff[ 19]*x12*x22    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  5]*x11+coeff[  6]*x31+coeff[  7]*x41+3*coeff[ 10]*x12*x21+coeff[ 11]*x23+3*coeff[ 12]*x12*x41+coeff[ 15]*x21+2*coeff[ 19]*x11*x22*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  8]*x41+coeff[ 10]*x13+3*coeff[ 11]*x22*x11+2*coeff[ 13]*x21*x32*x41+coeff[ 15]*x11+coeff[ 16]*x31+2*coeff[ 19]*x21*x12*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  6]*x11+2*coeff[ 13]*x31*x22*x41+2*coeff[ 14]*x31*x43+coeff[ 16]*x21+2*coeff[ 17]*x31+coeff[ 18]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+coeff[  8]*x21+2*coeff[  9]*x41+coeff[ 12]*x13+coeff[ 13]*x22*x32+3*coeff[ 14]*x42*x32+coeff[ 18]*x31+coeff[ 19]*x12*x22,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
