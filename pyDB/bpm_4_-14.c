#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 12;
    float avdat=  0.6620286E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 13]={
        -0.24476411E-04,-0.39048553E+02, 0.64933327E+02, 0.16138657E+01,
         0.53434956E+00,-0.97489798E+00,-0.88521814E+00,-0.11000644E-01,
         0.18203694E-01, 0.13614306E-01,-0.44694360E-01, 0.36879323E-01,
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
    int ncoeff= 16;
    float avdat=  0.3363785E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 17]={
        -0.23711629E-01,-0.38743790E+02, 0.64412399E+02, 0.19906969E+00,
        -0.69993073E+00, 0.61671156E+00,-0.25775033E+00, 0.74523944E+00,
        -0.51482654E+00,-0.31444186E-04,-0.70121721E-02, 0.46217460E-02,
        -0.59979833E-02, 0.31709950E-02,-0.82673514E-02, 0.10751640E-01,
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

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]    *x22        
        +coeff[  4]    *x21    *x41
        +coeff[  5]            *x42
        +coeff[  6]*x12            
        +coeff[  7]*x11    *x31    
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]        *x32    
        +coeff[  9]        *x33    
        +coeff[ 10]        *x32*x41
        +coeff[ 11]*x11            
        +coeff[ 12]        *x31    
        +coeff[ 13]*x12*x21        
        +coeff[ 14]    *x22    *x41
        +coeff[ 15]    *x21    *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  6]*x11+coeff[  7]*x31+coeff[ 11]+2*coeff[ 13]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  3]*x21+coeff[  4]*x41+coeff[ 13]*x12+2*coeff[ 14]*x21*x41+coeff[ 15]*x42,2)
        +dx3*dx3*pow(0+coeff[  7]*x11+2*coeff[  8]*x31+3*coeff[  9]*x32+2*coeff[ 10]*x31*x41+coeff[ 12],2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  4]*x21+2*coeff[  5]*x41+coeff[ 10]*x32+coeff[ 14]*x22+2*coeff[ 15]*x41*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.8831159E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.63517000E-04,-0.54420207E-01, 0.64658977E-01,-0.33842751E-02,
         0.99590537E-03, 0.10901488E-01,-0.41688480E-02,-0.87245563E-02,
         0.40693353E-02, 0.10524398E-03,-0.16820524E-03,-0.17282559E-03,
         0.25211199E-03,-0.85905252E-04, 0.49350725E-03,-0.54222369E-03,
         0.14565203E-03,-0.24240933E-03, 0.45940724E-04, 0.71404893E-05,
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
    float avdat=  0.3949477E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.22365155E-07,-0.57537340E-01, 0.69545835E-01, 0.39084232E-02,
        -0.70613115E-02,-0.65363743E-02, 0.11761174E-01,-0.78458492E-04,
         0.13010419E-03, 0.60552680E-04,-0.25645471E-04,-0.24869209E-03,
         0.33031113E-03, 0.18545261E-03, 0.12184097E-03,-0.16586507E-03,
        -0.13531296E-03,-0.96866192E-04, 0.71840012E-03,-0.72091015E-03,
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
