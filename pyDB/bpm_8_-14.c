#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  7;
    float avdat=  0.9440846E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[  8]={
        -0.50532544E-03,-0.38005749E+02, 0.28525734E+01, 0.63543930E+02,
        -0.46862478E+01, 0.11253360E-01,-0.22252690E-01,
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
    float avdat=  0.8106563E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
        -0.44665970E-02,-0.28393984E+01,-0.38007645E+02, 0.46647635E+01,
         0.63547329E+02, 0.83910994E-01, 0.24801873E-01,-0.20669114E-01,
         0.53995606E-02,-0.92305668E-01, 0.19914169E-01,
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
    float avdat=  0.1933514E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.11816728E-05,-0.16592231E-01,-0.54799728E-01, 0.25661219E-01,
         0.65061189E-01,-0.22804239E-04, 0.11664303E-04, 0.75199670E-04,
         0.83474890E-06, 0.15572034E-04, 0.29174469E-04,-0.97006399E-04,
        -0.29890537E-05, 0.76131851E-05,-0.85346792E-05, 0.86172686E-05,
        -0.50242475E-05,-0.16742346E-05, 0.80589216E-05,-0.10740502E-04,
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
        +coeff[  9]    *x21*x33*x41
        +coeff[ 10]*x12            
        +coeff[ 11]*x11    *x31    
        +coeff[ 12]*x11*x22        
        +coeff[ 13]        *x31*x42
        +coeff[ 14]*x12*x22*x31    
        +coeff[ 15]*x11    *x34    
        +coeff[ 16]*x12    *x31    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x21*x32    
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]            *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  6]*x21+4*coeff[  8]*x13+2*coeff[ 10]*x11+coeff[ 11]*x31+coeff[ 12]*x22+2*coeff[ 14]*x11*x22*x31+coeff[ 15]*x34+2*coeff[ 16]*x11*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  6]*x11+coeff[  9]*x33*x41+2*coeff[ 12]*x21*x11+2*coeff[ 14]*x21*x12*x31+coeff[ 17]*x32+coeff[ 18]*x42,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x41+2*coeff[  7]*x31+3*coeff[  9]*x32*x21*x41+coeff[ 11]*x11+coeff[ 13]*x42+coeff[ 14]*x12*x22+4*coeff[ 15]*x33*x11+coeff[ 16]*x12+2*coeff[ 17]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  5]*x31+coeff[  9]*x21*x33+2*coeff[ 13]*x41*x31+2*coeff[ 18]*x41*x21+3*coeff[ 19]*x42,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.1485738E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.44302815E-05,-0.54803751E-01, 0.16792927E-01, 0.65047376E-01,
        -0.26001647E-01,-0.26703245E-04, 0.97386677E-04, 0.11546986E-03,
         0.20625052E-04,-0.27614296E-04, 0.14686339E-05,-0.54227354E-06,
        -0.24800022E-05,-0.22513494E-04,-0.10570916E-03, 0.18744243E-03,
        -0.87514934E-04,-0.21219048E-03,-0.11979147E-07, 0.83260593E-05,
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
        +coeff[ 13]    *x21*x32*x42
        +coeff[ 14]*x11*x21        
        +coeff[ 15]    *x21*x31    
        +coeff[ 16]        *x32    
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]        *x31*x41
        +coeff[ 18]*x12*x21        
        +coeff[ 19]*x12*x23        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  5]*x11+coeff[  6]*x31+coeff[  7]*x41+3*coeff[ 10]*x12*x21+coeff[ 11]*x23+3*coeff[ 12]*x12*x41+coeff[ 14]*x21+2*coeff[ 18]*x11*x21+2*coeff[ 19]*x11*x23,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  8]*x41+coeff[ 10]*x13+3*coeff[ 11]*x22*x11+coeff[ 13]*x32*x42+coeff[ 14]*x11+coeff[ 15]*x31+coeff[ 18]*x12+3*coeff[ 19]*x22*x12,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  6]*x11+2*coeff[ 13]*x31*x21*x42+coeff[ 15]*x21+2*coeff[ 16]*x31+coeff[ 17]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+coeff[  8]*x21+2*coeff[  9]*x41+coeff[ 12]*x13+2*coeff[ 13]*x41*x21*x32+coeff[ 17]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
