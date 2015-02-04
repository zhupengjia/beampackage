#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 17;
    float avdat=  0.1732959E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 18]={
        -0.17467410E-03,-0.39805229E+02, 0.65795189E+02, 0.80619186E+00,
         0.23993089E+01,-0.43502867E-01, 0.71488850E-01,-0.14582486E+01,
        -0.13272983E+01, 0.33826433E-01,-0.47716647E-01,-0.14799167E-01,
        -0.27837517E-01, 0.70728928E-01,-0.69723693E-02, 0.94177807E-02,
         0.13596930E-02,
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
        +coeff[  9]*x11*x22        
        +coeff[ 10]*x11*x21    *x41
        +coeff[ 11]    *x22*x31    
        +coeff[ 12]    *x21*x31*x41
        +coeff[ 13]        *x31*x42
        +coeff[ 14]*x12    *x31    
        +coeff[ 15]*x11    *x32    
        +coeff[ 16]        *x32*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  8]*x41+coeff[  9]*x22+coeff[ 10]*x21*x41+2*coeff[ 14]*x11*x31+coeff[ 15]*x32,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]+coeff[  7]*x31+2*coeff[  9]*x21*x11+coeff[ 10]*x11*x41+2*coeff[ 11]*x21*x31+coeff[ 12]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x41+coeff[  7]*x21+coeff[ 11]*x22+coeff[ 12]*x21*x41+coeff[ 13]*x42+coeff[ 14]*x12+2*coeff[ 15]*x31*x11+2*coeff[ 16]*x31*x42,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  6]+coeff[  8]*x11+coeff[ 10]*x11*x21+coeff[ 12]*x21*x31+2*coeff[ 13]*x41*x31+2*coeff[ 16]*x41*x32,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.5615679E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.33428401E-01,-0.39186367E+02, 0.64752258E+02,-0.90559894E+00,
         0.81947458E+00,-0.29058534E+00, 0.24660850E+00, 0.87939525E+00,
        -0.62416047E+00, 0.34063622E-01, 0.28426916E-03,-0.25121702E-03,
         0.14097884E-01,-0.17696610E-01,-0.18527970E-01,-0.12718515E-01,
         0.72349273E-01,-0.99398255E-01, 0.17617563E-01,-0.51359745E-03,
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
        +coeff[ 10]        *x33    
        +coeff[ 11]*x14    *x31    
        +coeff[ 12]*x11            
        +coeff[ 13]        *x31    
        +coeff[ 14]*x12*x21        
        +coeff[ 15]    *x22    *x41
        +coeff[ 16]*x11    *x31*x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]        *x32*x41
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]*x11*x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  5]*x11+coeff[  7]*x31+4*coeff[ 11]*x13*x31+coeff[ 12]+2*coeff[ 14]*x11*x21+coeff[ 16]*x31*x41+coeff[ 19]*x21,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  6]*x21+coeff[  9]*x32+coeff[ 14]*x12+2*coeff[ 15]*x21*x41+coeff[ 18]*x42+coeff[ 19]*x11,2)
        +dx3*dx3*pow(0+coeff[  7]*x11+2*coeff[  8]*x31+2*coeff[  9]*x31*x21+3*coeff[ 10]*x32+coeff[ 11]*x14+coeff[ 13]+coeff[ 16]*x11*x41+2*coeff[ 17]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41+coeff[ 15]*x22+coeff[ 16]*x11*x31+coeff[ 17]*x32+2*coeff[ 18]*x41*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1206173E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.98631273E-04,-0.50751101E-01, 0.58757514E-01,-0.50309966E-02,
         0.16722104E-02, 0.16218727E-01,-0.12964579E-01,-0.66042193E-02,
         0.40482040E-03,-0.64044894E-03, 0.62300740E-02, 0.39293809E-03,
        -0.56033832E-03,-0.84412168E-04, 0.12615279E-03,-0.21514601E-03,
        -0.19717780E-03, 0.27428154E-03, 0.14514776E-03,-0.56012409E-05,
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
        +coeff[  3]*x12            
        +coeff[  4]    *x22        
        +coeff[  5]*x11    *x31    
        +coeff[  6]        *x32    
        +coeff[  7]    *x21    *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x11            
        +coeff[  9]        *x31    
        +coeff[ 10]            *x42
        +coeff[ 11]    *x21*x32*x41
        +coeff[ 12]        *x32*x42
        +coeff[ 13]*x12*x21        
        +coeff[ 14]*x11*x21*x31    
        +coeff[ 15]*x11    *x31*x41
        +coeff[ 16]*x12*x22        
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12*x21    *x41
        +coeff[ 18]*x12        *x41
        +coeff[ 19]    *x21*x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  3]*x11+coeff[  5]*x31+coeff[  8]+2*coeff[ 13]*x11*x21+coeff[ 14]*x21*x31+coeff[ 15]*x31*x41+2*coeff[ 16]*x11*x22+2*coeff[ 17]*x11*x21*x41+2*coeff[ 18]*x11*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  7]*x41+coeff[ 11]*x32*x41+coeff[ 13]*x12+coeff[ 14]*x11*x31+2*coeff[ 16]*x21*x12+coeff[ 17]*x12*x41+coeff[ 19]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  5]*x11+2*coeff[  6]*x31+coeff[  9]+2*coeff[ 11]*x31*x21*x41+2*coeff[ 12]*x31*x42+coeff[ 14]*x11*x21+coeff[ 15]*x11*x41+coeff[ 19]*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  7]*x21+2*coeff[ 10]*x41+coeff[ 11]*x21*x32+2*coeff[ 12]*x41*x32+coeff[ 15]*x11*x31+coeff[ 17]*x12*x21+coeff[ 18]*x12+coeff[ 19]*x21*x31,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.7720074E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.54136845E-09, 0.52938256E-02, 0.58978773E-02,-0.10591876E-01,
        -0.97948126E-02, 0.17525565E-01,-0.11334425E+01,-0.50522786E-01,
         0.38591793E+01, 0.56711297E-01,-0.45441222E+01, 0.18336517E+01,
         0.22788195E-01, 0.30587792E-01, 0.18230096E+01,-0.24661165E-01,
        -0.35293087E-01,-0.63278222E+01, 0.74122758E+01,-0.29173980E+01,
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
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;

//                 function

    float v_target_phi                              =avdat
        +coeff[  0]                
        +coeff[  1]        *x31    
        +coeff[  2]*x11*x21        
        +coeff[  3]    *x21*x31    
        +coeff[  4]*x11        *x41
        +coeff[  5]        *x31*x41
        +coeff[  6]*x13            
        +coeff[  7]*x11*x22        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]*x12    *x31    
        +coeff[  9]    *x22*x31    
        +coeff[ 10]*x11    *x32    
        +coeff[ 11]        *x33    
        +coeff[ 12]*x13*x22        
        +coeff[ 13]*x11*x24        
        +coeff[ 14]*x14    *x31    
        +coeff[ 15]*x12*x22*x31    
        +coeff[ 16]    *x24*x31    
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x13    *x32    
        +coeff[ 18]*x12    *x33    
        +coeff[ 19]*x11    *x34    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  2]*x21+coeff[  4]*x41+3*coeff[  6]*x12+coeff[  7]*x22+2*coeff[  8]*x11*x31+coeff[ 10]*x32+3*coeff[ 12]*x12*x22+coeff[ 13]*x24+4*coeff[ 14]*x13*x31+2*coeff[ 15]*x11*x22*x31+3*coeff[ 17]*x12*x32+2*coeff[ 18]*x11*x33+coeff[ 19]*x34,2)
        +dx2*dx2*pow(0+coeff[  2]*x11+coeff[  3]*x31+2*coeff[  7]*x21*x11+2*coeff[  9]*x21*x31+2*coeff[ 12]*x21*x13+4*coeff[ 13]*x23*x11+2*coeff[ 15]*x21*x12*x31+4*coeff[ 16]*x23*x31,2)
        +dx3*dx3*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[  8]*x12+coeff[  9]*x22+2*coeff[ 10]*x31*x11+3*coeff[ 11]*x32+coeff[ 14]*x14+coeff[ 15]*x12*x22+coeff[ 16]*x24+2*coeff[ 17]*x31*x13+3*coeff[ 18]*x32*x12+4*coeff[ 19]*x33*x11,2)
        +dx4*dx4*pow(0+coeff[  4]*x11+coeff[  5]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
