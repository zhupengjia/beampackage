#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 17;
    float avdat=  0.1743476E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 18]={
        -0.17636981E-03,-0.40583023E+02, 0.66730476E+02, 0.88647437E+00,
        -0.14606659E+01,-0.16024895E+01, 0.26380181E+01,-0.47752146E-01,
         0.78487180E-01, 0.33229034E-01,-0.46933763E-01,-0.14068650E-01,
        -0.28934056E-01, 0.70947282E-01,-0.78048818E-02, 0.10551029E-01,
         0.11877671E-02,
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
        +coeff[  4]*x11        *x41
        +coeff[  5]    *x21*x31    
        +coeff[  6]        *x31*x41
        +coeff[  7]    *x21        
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]            *x41
        +coeff[  9]*x11*x22        
        +coeff[ 10]*x11*x21    *x41
        +coeff[ 11]    *x22*x31    
        +coeff[ 12]    *x21*x31*x41
        +coeff[ 13]        *x31*x42
        +coeff[ 14]*x12    *x31    
        +coeff[ 15]*x11    *x32    
        +coeff[ 16]*x11    *x31*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  4]*x41+coeff[  9]*x22+coeff[ 10]*x21*x41+2*coeff[ 14]*x11*x31+coeff[ 15]*x32+coeff[ 16]*x31*x42,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]*x31+coeff[  7]+2*coeff[  9]*x21*x11+coeff[ 10]*x11*x41+2*coeff[ 11]*x21*x31+coeff[ 12]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  5]*x21+coeff[  6]*x41+coeff[ 11]*x22+coeff[ 12]*x21*x41+coeff[ 13]*x42+coeff[ 14]*x12+2*coeff[ 15]*x31*x11+coeff[ 16]*x11*x42,2)
        +dx4*dx4*pow(0+coeff[  4]*x11+coeff[  6]*x31+coeff[  8]+coeff[ 10]*x11*x21+coeff[ 12]*x21*x31+2*coeff[ 13]*x41*x31+2*coeff[ 16]*x41*x11*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.4025083E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.32111365E-01,-0.39886692E+02, 0.65563179E+02,-0.10011243E+01,
         0.90763301E+00,-0.80648327E+00,-0.76016056E-03, 0.11819611E-03,
        -0.36110654E+00, 0.27183130E+00, 0.11079873E+01,-0.26308103E-01,
         0.19471901E-01, 0.35428759E-01,-0.19442575E-01,-0.13427667E-01,
         0.75223878E-01,-0.10322651E+00, 0.18335944E-01,-0.55389776E-03,
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
        +coeff[  3]    *x21    *x41
        +coeff[  4]            *x42
        +coeff[  5]        *x32    
        +coeff[  6]*x14            
        +coeff[  7]        *x34    
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x12            
        +coeff[  9]    *x22        
        +coeff[ 10]*x11    *x31    
        +coeff[ 11]        *x31    
        +coeff[ 12]*x11            
        +coeff[ 13]    *x21*x32    
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
        +dx1*dx1*pow(0+4*coeff[  6]*x13+2*coeff[  8]*x11+coeff[ 10]*x31+coeff[ 12]+2*coeff[ 14]*x11*x21+coeff[ 16]*x31*x41+coeff[ 19]*x21,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  9]*x21+coeff[ 13]*x32+coeff[ 14]*x12+2*coeff[ 15]*x21*x41+coeff[ 18]*x42+coeff[ 19]*x11,2)
        +dx3*dx3*pow(0+2*coeff[  5]*x31+4*coeff[  7]*x33+coeff[ 10]*x11+coeff[ 11]+2*coeff[ 13]*x31*x21+coeff[ 16]*x11*x41+2*coeff[ 17]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41+coeff[ 15]*x22+coeff[ 16]*x11*x31+coeff[ 17]*x32+2*coeff[ 18]*x41*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1117651E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.99764933E-04,-0.50699815E-01, 0.58698010E-01,-0.50125537E-02,
         0.15764602E-02, 0.16160807E-01,-0.12924738E-01,-0.63843061E-02,
         0.39825638E-03,-0.63266262E-03, 0.61046197E-02,-0.63638022E-03,
         0.89668925E-03, 0.89368224E-03,-0.12596308E-02, 0.27908969E-04,
        -0.50039322E-04, 0.23908231E-04,-0.61413166E-05,-0.32098735E-04,
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
        +coeff[ 11]*x11*x21*x31*x41
        +coeff[ 12]    *x21*x32*x41
        +coeff[ 13]*x11    *x31*x42
        +coeff[ 14]        *x32*x42
        +coeff[ 15]*x11*x21*x31    
        +coeff[ 16]        *x32*x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x11*x21    *x41
        +coeff[ 18]*x11        *x42
        +coeff[ 19]        *x31*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  3]*x11+coeff[  5]*x31+coeff[  8]+coeff[ 11]*x21*x31*x41+coeff[ 13]*x31*x42+coeff[ 15]*x21*x31+coeff[ 17]*x21*x41+coeff[ 18]*x42,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  7]*x41+coeff[ 11]*x11*x31*x41+coeff[ 12]*x32*x41+coeff[ 15]*x11*x31+coeff[ 17]*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  5]*x11+2*coeff[  6]*x31+coeff[  9]+coeff[ 11]*x11*x21*x41+2*coeff[ 12]*x31*x21*x41+coeff[ 13]*x11*x42+2*coeff[ 14]*x31*x42+coeff[ 15]*x11*x21+2*coeff[ 16]*x31*x41+coeff[ 19]*x42,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  7]*x21+2*coeff[ 10]*x41+coeff[ 11]*x11*x21*x31+coeff[ 12]*x21*x32+2*coeff[ 13]*x41*x11*x31+2*coeff[ 14]*x41*x32+coeff[ 16]*x32+coeff[ 17]*x11*x21+2*coeff[ 18]*x41*x11+2*coeff[ 19]*x41*x31,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.7711709E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.35720810E-09, 0.52877953E-02, 0.58736377E-02,-0.10563425E-01,
        -0.97676162E-02, 0.17493853E-01,-0.11322616E+01,-0.50466456E-01,
         0.38551836E+01, 0.56646805E-01,-0.45394521E+01, 0.18317829E+01,
         0.22757124E-01, 0.30556185E-01, 0.18211429E+01,-0.24625771E-01,
        -0.35256881E-01,-0.63213725E+01, 0.74047604E+01,-0.29144576E+01,
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
