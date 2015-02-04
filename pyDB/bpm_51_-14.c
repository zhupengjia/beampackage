#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  7;
    float avdat=  0.3007129E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[  8]={
        -0.22757356E-05,-0.38552334E+02, 0.64305519E+02, 0.39603245E+00,
         0.12004460E+01,-0.72397202E+00,-0.65670151E+00,
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
    int ncoeff= 11;
    float avdat=  0.2377510E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
        -0.18015277E-01,-0.38384773E+02, 0.64014427E+02,-0.49983913E+00,
         0.44827333E+00,-0.18193060E+00, 0.13868949E+00, 0.53132701E+00,
        -0.36910897E+00,-0.11038817E-02, 0.14630270E-02,
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
    float x43 = x42*x4;

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
        +coeff[ 10]            *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  5]*x11+coeff[  7]*x31,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  6]*x21+coeff[  9]*x32,2)
        +dx3*dx3*pow(0+coeff[  7]*x11+2*coeff[  8]*x31+2*coeff[  9]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41+3*coeff[ 10]*x42,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.6668012E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.42813514E-04,-0.55292066E-01, 0.65974712E-01, 0.78757759E-02,
        -0.63568014E-02,-0.31988875E-02, 0.77923335E-03,-0.24145676E-02,
         0.31012301E-02,-0.55237026E-04, 0.35174260E-04,-0.92705028E-04,
         0.13489898E-03,-0.28140348E-03,-0.17379814E-04, 0.19525512E-03,
        -0.24933997E-03,-0.95761061E-04, 0.32486999E-03, 0.33393319E-05,
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
        +coeff[ 13]        *x32*x42
        +coeff[ 14]*x12        *x41
        +coeff[ 15]*x11    *x31*x41
        +coeff[ 16]        *x32*x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12*x22        
        +coeff[ 18]*x11*x21*x31*x41
        +coeff[ 19]    *x21    *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  3]*x31+2*coeff[  7]*x11+coeff[ 10]+coeff[ 11]*x21*x31+2*coeff[ 14]*x11*x41+coeff[ 15]*x31*x41+2*coeff[ 17]*x11*x22+coeff[ 18]*x21*x31*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  5]*x41+2*coeff[  6]*x21+coeff[ 11]*x11*x31+coeff[ 12]*x32+2*coeff[ 17]*x21*x12+coeff[ 18]*x11*x31*x41+coeff[ 19]*x42,2)
        +dx3*dx3*pow(0+coeff[  3]*x11+2*coeff[  4]*x31+coeff[  9]+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+2*coeff[ 13]*x31*x42+coeff[ 15]*x11*x41+2*coeff[ 16]*x31*x41+coeff[ 18]*x11*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  5]*x21+2*coeff[  8]*x41+2*coeff[ 13]*x41*x32+coeff[ 14]*x12+coeff[ 15]*x11*x31+coeff[ 16]*x32+coeff[ 18]*x11*x21*x31+2*coeff[ 19]*x41*x21,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4717260E-05;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17988E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.20264230E-08,-0.57051808E-01, 0.68740286E-01, 0.29039837E-02,
        -0.52488404E-02,-0.48703328E-02, 0.87706130E-02, 0.48983551E-04,
        -0.29562560E-04,-0.52334068E-04,-0.13526103E-03, 0.25844050E-04,
         0.94865289E-04,-0.12306486E-03,-0.46388691E-03, 0.32258284E-03,
        -0.46998826E-04, 0.16489938E-03, 0.16827168E-03, 0.54627326E-04,
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
