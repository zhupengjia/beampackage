#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 12;
    float avdat=  0.6622339E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 13]={
        -0.26988884E-04,-0.39077824E+02, 0.64968704E+02, 0.16204052E+01,
         0.53668213E+00,-0.97893369E+00,-0.88897568E+00,-0.11097169E-01,
         0.18341077E-01, 0.13730576E-01,-0.44999760E-01, 0.37094247E-01,
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
    float x22 = x21*x2;
    float x31 = x3;
    float x41 = x4;
    float x42 = x41*x4;

//                 function

    float v_target_x                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]        *x31    
        +coeff[  3]        *x31*x41
        +coeff[  4]*x11*x21        
        +coeff[  5]    *x21*x31    
        +coeff[  6]*x11        *x41
        +coeff[  7]    *x21        
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]            *x41
        +coeff[  9]    *x22*x31    
        +coeff[ 10]    *x21*x31*x41
        +coeff[ 11]        *x31*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  4]*x21+coeff[  6]*x41,2)
        +dx2*dx2*pow(0+coeff[  4]*x11+coeff[  5]*x31+coeff[  7]+2*coeff[  9]*x21*x31+coeff[ 10]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  3]*x41+coeff[  5]*x21+coeff[  9]*x22+coeff[ 10]*x21*x41+coeff[ 11]*x42,2)
        +dx4*dx4*pow(0+coeff[  3]*x31+coeff[  6]*x11+coeff[  8]+coeff[ 10]*x21*x31+2*coeff[ 11]*x41*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 18;
    float avdat=  0.3408799E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 19]={
        -0.23684982E-01,-0.38773140E+02, 0.64447250E+02, 0.95516807E-04,
         0.20226325E+00,-0.70920932E+00, 0.62342417E+00,-0.26164120E+00,
         0.75595945E+00,-0.52237791E+00,-0.14951042E-01,-0.70953160E-02,
         0.46765227E-02, 0.32572094E-02,-0.66269371E-02, 0.89619458E-02,
         0.20708866E-01,-0.11561035E-01,
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
    float x44 = x43*x4;

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]            *x44
        +coeff[  4]    *x22        
        +coeff[  5]    *x21    *x41
        +coeff[  6]            *x42
        +coeff[  7]*x12            
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x11    *x31    
        +coeff[  9]        *x32    
        +coeff[ 10]        *x33    
        +coeff[ 11]        *x32*x41
        +coeff[ 12]*x13            
        +coeff[ 13]*x12*x21        
        +coeff[ 14]    *x22    *x41
        +coeff[ 15]    *x21    *x42
        +coeff[ 16]*x11    *x34    
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]*x14    *x33    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  7]*x11+coeff[  8]*x31+3*coeff[ 12]*x12+2*coeff[ 13]*x11*x21+coeff[ 16]*x34+4*coeff[ 17]*x13*x33,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  5]*x41+coeff[ 13]*x12+2*coeff[ 14]*x21*x41+coeff[ 15]*x42,2)
        +dx3*dx3*pow(0+coeff[  8]*x11+2*coeff[  9]*x31+3*coeff[ 10]*x32+2*coeff[ 11]*x31*x41+4*coeff[ 16]*x33*x11+3*coeff[ 17]*x32*x14,2)
        +dx4*dx4*pow(0+coeff[  2]+4*coeff[  3]*x43+coeff[  5]*x21+2*coeff[  6]*x41+coeff[ 11]*x32+coeff[ 14]*x22+2*coeff[ 15]*x41*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.8808760E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.63606152E-04,-0.54414712E-01, 0.64651877E-01,-0.33755519E-02,
         0.98321971E-03, 0.10880431E-01,-0.41380376E-02,-0.87120356E-02,
         0.40503307E-02, 0.10605837E-03,-0.16909219E-03,-0.18023585E-03,
         0.25975477E-03,-0.74016680E-04, 0.47396490E-03,-0.53468643E-03,
         0.15455901E-03,-0.25349576E-03, 0.46741367E-04, 0.72734474E-05,
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
        +coeff[  6]    *x21    *x41
        +coeff[  7]        *x32    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]            *x42
        +coeff[  9]*x11            
        +coeff[ 10]        *x31    
        +coeff[ 11]*x11*x21*x31    
        +coeff[ 12]    *x21*x32    
        +coeff[ 13]*x12        *x41
        +coeff[ 14]*x11    *x31*x41
        +coeff[ 15]        *x32*x41
        +coeff[ 16]    *x22*x32    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x21*x32*x41
        +coeff[ 18]*x12        *x42
        +coeff[ 19]    *x21    *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  3]*x11+coeff[  5]*x31+coeff[  9]+coeff[ 11]*x21*x31+2*coeff[ 13]*x11*x41+coeff[ 14]*x31*x41+2*coeff[ 18]*x11*x42,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  6]*x41+coeff[ 11]*x11*x31+coeff[ 12]*x32+2*coeff[ 16]*x21*x32+coeff[ 17]*x32*x41+coeff[ 19]*x42,2)
        +dx3*dx3*pow(0+coeff[  5]*x11+2*coeff[  7]*x31+coeff[ 10]+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+coeff[ 14]*x11*x41+2*coeff[ 15]*x31*x41+2*coeff[ 16]*x31*x22+2*coeff[ 17]*x31*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  6]*x21+2*coeff[  8]*x41+coeff[ 13]*x12+coeff[ 14]*x11*x31+coeff[ 15]*x32+coeff[ 17]*x21*x32+2*coeff[ 18]*x41*x12+2*coeff[ 19]*x41*x21,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.3949403E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.22433214E-07,-0.57536207E-01, 0.69544472E-01, 0.39079185E-02,
        -0.70607434E-02,-0.65358677E-02, 0.11760601E-01,-0.78453369E-04,
         0.13009892E-03, 0.60636394E-04,-0.25678826E-04,-0.24945661E-03,
         0.33126355E-03, 0.18597531E-03, 0.12217718E-03,-0.16630210E-03,
        -0.13511562E-03,-0.97161428E-04, 0.71908662E-03,-0.72206103E-03,
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
    float x23 = x22*x2;
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
        +coeff[  9]        *x31*x42
        +coeff[ 10]*x11*x22        
        +coeff[ 11]*x11*x21*x32    
        +coeff[ 12]    *x21*x33    
        +coeff[ 13]        *x31*x43
        +coeff[ 14]*x11*x23        
        +coeff[ 15]*x11*x22    *x41
        +coeff[ 16]*x12    *x31*x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]    *x22*x31*x41
        +coeff[ 18]*x11    *x32*x41
        +coeff[ 19]        *x33*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[ 10]*x22+coeff[ 11]*x21*x32+coeff[ 14]*x23+coeff[ 15]*x22*x41+2*coeff[ 16]*x11*x31*x41+coeff[ 18]*x32*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  7]+2*coeff[ 10]*x21*x11+coeff[ 11]*x11*x32+coeff[ 12]*x33+3*coeff[ 14]*x22*x11+2*coeff[ 15]*x21*x11*x41+2*coeff[ 17]*x21*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x42+2*coeff[ 11]*x31*x11*x21+3*coeff[ 12]*x32*x21+coeff[ 13]*x43+coeff[ 16]*x12*x41+coeff[ 17]*x22*x41+2*coeff[ 18]*x31*x11*x41+3*coeff[ 19]*x32*x41,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  8]+2*coeff[  9]*x41*x31+3*coeff[ 13]*x42*x31+coeff[ 15]*x11*x22+coeff[ 16]*x12*x31+coeff[ 17]*x22*x31+coeff[ 18]*x11*x32+coeff[ 19]*x33,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
