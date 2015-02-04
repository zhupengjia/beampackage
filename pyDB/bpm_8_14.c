#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 12;
    float avdat=  0.8971587E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 13]={
         0.52621639E-04,-0.39550629E+02, 0.33418801E+01, 0.65375999E+02,
        -0.54386673E+01,-0.29866450E-01, 0.43991320E-01, 0.54301340E-01,
        -0.81049487E-01,-0.74321742E-03, 0.50709872E-02,-0.67824512E-02,
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
    float x42 = x41*x4;

//                 function

    float v_target_x                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11*x21        
        +coeff[  6]*x11        *x41
        +coeff[  7]    *x21*x31    
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]        *x31*x41
        +coeff[  9]*x11    *x31    
        +coeff[ 10]    *x21    *x41
        +coeff[ 11]            *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41+coeff[  9]*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  7]*x31+coeff[ 10]*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x21+coeff[  8]*x41+coeff[  9]*x11,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  8]*x31+coeff[ 10]*x21+2*coeff[ 11]*x41,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 10;
    float avdat=  0.1329268E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 11]={
        -0.42762267E-02,-0.33233316E+01,-0.39552856E+02, 0.54079723E+01,
         0.65380249E+02, 0.86148746E-01, 0.25654672E-01, 0.12124118E-01,
        -0.80780024E-02,-0.95122755E-01,
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
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  6]*x11+coeff[  9]*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  7]*x41+2*coeff[  8]*x21,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  5]*x31+coeff[  9]*x11,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1734684E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.10383814E-05,-0.17609218E-01,-0.54472655E-01, 0.26865251E-01,
         0.64559095E-01,-0.99555655E-04,-0.42767766E-04, 0.65233973E-04,
         0.53671687E-06, 0.17385042E-04, 0.23171931E-04,-0.22439483E-05,
         0.64791425E-05,-0.99823328E-05,-0.81296690E-04, 0.66637731E-04,
         0.62696308E-04,-0.10206275E-04, 0.14699649E-04,-0.17294996E-05,
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
        +coeff[ 11]*x11*x22        
        +coeff[ 12]        *x31*x42
        +coeff[ 13]*x12*x22*x31    
        +coeff[ 14]*x11    *x31    
        +coeff[ 15]    *x21*x31    
        +coeff[ 16]*x11        *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12    *x31    
        +coeff[ 18]*x11    *x32    
        +coeff[ 19]    *x21*x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  6]*x21+4*coeff[  8]*x13+2*coeff[ 10]*x11+coeff[ 11]*x22+2*coeff[ 13]*x11*x22*x31+coeff[ 14]*x31+coeff[ 16]*x41+2*coeff[ 17]*x11*x31+coeff[ 18]*x32,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  6]*x11+coeff[  9]*x33*x41+2*coeff[ 11]*x21*x11+2*coeff[ 13]*x21*x12*x31+coeff[ 15]*x31+coeff[ 19]*x32,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x41+2*coeff[  7]*x31+3*coeff[  9]*x32*x21*x41+coeff[ 12]*x42+coeff[ 13]*x12*x22+coeff[ 14]*x11+coeff[ 15]*x21+coeff[ 17]*x12+2*coeff[ 18]*x31*x11+2*coeff[ 19]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  5]*x31+coeff[  9]*x21*x33+2*coeff[ 12]*x41*x31+coeff[ 16]*x11,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.1827296E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.49725122E-05,-0.54483049E-01, 0.17810287E-01, 0.64555556E-01,
        -0.27208183E-01,-0.28090870E-04, 0.10164021E-03, 0.11772666E-03,
        -0.69177600E-04, 0.76219487E-06,-0.24877909E-06,-0.12770261E-05,
        -0.21605796E-04,-0.10915773E-03,-0.31168238E-04, 0.19107375E-03,
        -0.90612390E-04, 0.92673385E-04,-0.21445785E-03, 0.79848214E-05,
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
        +coeff[  8]            *x42
        +coeff[  9]*x13*x21        
        +coeff[ 10]*x11*x23        
        +coeff[ 11]*x13        *x41
        +coeff[ 12]    *x21*x32*x42
        +coeff[ 13]*x11*x21        
        +coeff[ 14]    *x22        
        +coeff[ 15]    *x21*x31    
        +coeff[ 16]        *x32    
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]    *x21    *x41
        +coeff[ 18]        *x31*x41
        +coeff[ 19]*x12*x23        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  5]*x11+coeff[  6]*x31+coeff[  7]*x41+3*coeff[  9]*x12*x21+coeff[ 10]*x23+3*coeff[ 11]*x12*x41+coeff[ 13]*x21+2*coeff[ 19]*x11*x23,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  9]*x13+3*coeff[ 10]*x22*x11+coeff[ 12]*x32*x42+coeff[ 13]*x11+2*coeff[ 14]*x21+coeff[ 15]*x31+coeff[ 17]*x41+3*coeff[ 19]*x22*x12,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  6]*x11+2*coeff[ 12]*x31*x21*x42+coeff[ 15]*x21+2*coeff[ 16]*x31+coeff[ 18]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+2*coeff[  8]*x41+coeff[ 11]*x13+2*coeff[ 12]*x41*x21*x32+coeff[ 17]*x21+coeff[ 18]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
