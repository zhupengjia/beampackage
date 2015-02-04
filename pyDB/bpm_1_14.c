#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 11;
    float avdat=  0.1754370E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
        -0.35898942E-04,-0.41396221E+02, 0.67710220E+02, 0.96926558E+00,
        -0.15984436E+01,-0.17515659E+01, 0.28849988E+01, 0.85750714E-01,
        -0.52165996E-01,-0.19577140E-01, 0.27551534E-01,
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

//                 function

    float v_target_x                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]        *x31    
        +coeff[  3]*x11*x21        
        +coeff[  4]*x11        *x41
        +coeff[  5]    *x21*x31    
        +coeff[  6]        *x31*x41
        +coeff[  7]            *x41
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]    *x21        
        +coeff[  9]*x11*x22        
        +coeff[ 10]*x11*x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  4]*x41+coeff[  9]*x22+coeff[ 10]*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]*x31+coeff[  8]+2*coeff[  9]*x21*x11+coeff[ 10]*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  5]*x21+coeff[  6]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]*x11+coeff[  6]*x31+coeff[  7]+coeff[ 10]*x11*x21,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.2504251E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.30604592E-01,-0.40611912E+02, 0.66402824E+02, 0.29850423E+00,
        -0.11002679E+01, 0.25276293E-02,-0.19836328E-02, 0.99828529E+00,
        -0.35657451E-01,-0.43530267E+00, 0.13466977E+01,-0.99793547E+00,
         0.25377380E-01, 0.34270197E-01,-0.17630886E-01,-0.63824854E-02,
        -0.13336268E-01, 0.85981928E-01,-0.10870022E+00, 0.18289633E-01,
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
        +coeff[  8]        *x31    
        +coeff[  9]*x12            
        +coeff[ 10]*x11    *x31    
        +coeff[ 11]        *x32    
        +coeff[ 12]*x11            
        +coeff[ 13]    *x21*x32    
        +coeff[ 14]*x12*x21        
        +coeff[ 15]*x12        *x41
        +coeff[ 16]    *x22    *x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]*x11    *x31*x41
        +coeff[ 18]        *x32*x41
        +coeff[ 19]    *x21    *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  6]*x13*x32+2*coeff[  9]*x11+coeff[ 10]*x31+coeff[ 12]+2*coeff[ 14]*x11*x21+2*coeff[ 15]*x11*x41+coeff[ 17]*x31*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  3]*x21+coeff[  4]*x41+coeff[ 13]*x32+coeff[ 14]*x12+2*coeff[ 16]*x21*x41+coeff[ 19]*x42,2)
        +dx3*dx3*pow(0+4*coeff[  5]*x33+2*coeff[  6]*x31*x14+coeff[  8]+coeff[ 10]*x11+2*coeff[ 11]*x31+2*coeff[ 13]*x31*x21+coeff[ 17]*x11*x41+2*coeff[ 18]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  4]*x21+2*coeff[  7]*x41+coeff[ 15]*x12+coeff[ 16]*x22+coeff[ 17]*x11*x31+coeff[ 18]*x32+2*coeff[ 19]*x41*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1025925E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.99108620E-04,-0.50634447E-01, 0.58621772E-01,-0.50448515E-02,
         0.16193050E-02, 0.16216347E-01,-0.12940251E-01,-0.64892112E-02,
         0.40956857E-03,-0.64579095E-03, 0.61683203E-02,-0.65604155E-03,
         0.21154194E-02, 0.88584970E-03,-0.94714114E-05, 0.34428402E-03,
        -0.10989191E-02,-0.17177727E-02, 0.90081139E-05,-0.14146050E-04,
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
        +coeff[ 11]    *x22*x32    
        +coeff[ 12]    *x21*x32*x41
        +coeff[ 13]*x12        *x42
        +coeff[ 14]    *x21*x32    
        +coeff[ 15]*x12*x22        
        +coeff[ 16]*x12*x21    *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]        *x32*x42
        +coeff[ 18]    *x22*x31    
        +coeff[ 19]*x11*x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  3]*x11+coeff[  5]*x31+coeff[  8]+2*coeff[ 13]*x11*x42+2*coeff[ 15]*x11*x22+2*coeff[ 16]*x11*x21*x41+coeff[ 19]*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  7]*x41+2*coeff[ 11]*x21*x32+coeff[ 12]*x32*x41+coeff[ 14]*x32+2*coeff[ 15]*x21*x12+coeff[ 16]*x12*x41+2*coeff[ 18]*x21*x31+coeff[ 19]*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  5]*x11+2*coeff[  6]*x31+coeff[  9]+2*coeff[ 11]*x31*x22+2*coeff[ 12]*x31*x21*x41+2*coeff[ 14]*x31*x21+2*coeff[ 17]*x31*x42+coeff[ 18]*x22,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  7]*x21+2*coeff[ 10]*x41+coeff[ 12]*x21*x32+2*coeff[ 13]*x41*x12+coeff[ 16]*x12*x21+2*coeff[ 17]*x41*x32+coeff[ 19]*x11*x21,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.7703677E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.12564089E-09, 0.52821618E-02, 0.58486769E-02,-0.10534083E-01,
        -0.97396988E-02, 0.17461266E-01,-0.11311375E+01,-0.50413437E-01,
         0.38513763E+01, 0.56586403E-01,-0.45349979E+01, 0.18299983E+01,
         0.22728330E-01, 0.30526116E-01, 0.18193710E+01,-0.24593156E-01,
        -0.35222426E-01,-0.63152509E+01, 0.73976283E+01,-0.29116678E+01,
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
