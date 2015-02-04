#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.6676099E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
        -0.23387236E-04,-0.39861843E+02, 0.65917770E+02,-0.10740857E+01,
         0.17790698E+01, 0.58905059E+00,-0.97671878E+00, 0.20039359E-01,
        -0.12104545E-01,
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
        +coeff[  3]    *x21*x31    
        +coeff[  4]        *x31*x41
        +coeff[  5]*x11*x21        
        +coeff[  6]*x11        *x41
        +coeff[  7]            *x41
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]    *x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x31+coeff[  5]*x11+coeff[  8],2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  3]*x21+coeff[  4]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  6]*x11+coeff[  7],2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4571130E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.22821445E-01,-0.39518459E+02, 0.65333115E+02, 0.67580760E+00,
         0.21418124E+00,-0.64283490E+00,-0.76141554E+00, 0.25919394E-02,
        -0.75578876E-02, 0.46052765E-02,-0.30748791E+00, 0.90598822E+00,
         0.60408171E-02,-0.83330376E-02, 0.18341703E-01,-0.98282220E-02,
         0.37695020E-01,-0.52101631E-01,-0.89309392E-02, 0.11589214E-01,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]            *x42
        +coeff[  4]    *x22        
        +coeff[  5]        *x32    
        +coeff[  6]    *x21    *x41
        +coeff[  7]*x14            
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x13    *x31    
        +coeff[  9]*x12    *x32    
        +coeff[ 10]*x12            
        +coeff[ 11]*x11    *x31    
        +coeff[ 12]*x11            
        +coeff[ 13]        *x31    
        +coeff[ 14]    *x21*x32    
        +coeff[ 15]*x12*x21        
        +coeff[ 16]*x11    *x31*x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]        *x32*x41
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]            *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  7]*x13+3*coeff[  8]*x12*x31+2*coeff[  9]*x11*x32+2*coeff[ 10]*x11+coeff[ 11]*x31+coeff[ 12]+2*coeff[ 15]*x11*x21+coeff[ 16]*x31*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  6]*x41+coeff[ 14]*x32+coeff[ 15]*x12+coeff[ 18]*x42,2)
        +dx3*dx3*pow(0+2*coeff[  5]*x31+coeff[  8]*x13+2*coeff[  9]*x31*x12+coeff[ 11]*x11+coeff[ 13]+2*coeff[ 14]*x31*x21+coeff[ 16]*x11*x41+2*coeff[ 17]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+2*coeff[  3]*x41+coeff[  6]*x21+coeff[ 16]*x11*x31+coeff[ 17]*x32+2*coeff[ 18]*x41*x21+3*coeff[ 19]*x42,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.8208974E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.59646551E-04,-0.54393612E-01, 0.64627349E-01,-0.33002095E-02,
         0.11229175E-02, 0.10709785E-01,-0.44625797E-02,-0.86077312E-02,
         0.42447313E-02, 0.10586443E-03,-0.16888343E-03,-0.19640998E-03,
         0.28225780E-03,-0.46756095E-03,-0.95556796E-04, 0.55147818E-03,
        -0.60149119E-03,-0.16241729E-03, 0.54406968E-03, 0.66117777E-05,
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
        +coeff[ 13]        *x32*x42
        +coeff[ 14]*x12        *x41
        +coeff[ 15]*x11    *x31*x41
        +coeff[ 16]        *x32*x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12*x22        
        +coeff[ 18]*x11*x21*x31*x41
        +coeff[ 19]    *x22    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  3]*x11+coeff[  5]*x31+coeff[  9]+coeff[ 11]*x21*x31+2*coeff[ 14]*x11*x41+coeff[ 15]*x31*x41+2*coeff[ 17]*x11*x22+coeff[ 18]*x21*x31*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  6]*x41+coeff[ 11]*x11*x31+coeff[ 12]*x32+2*coeff[ 17]*x21*x12+coeff[ 18]*x11*x31*x41+2*coeff[ 19]*x21*x41,2)
        +dx3*dx3*pow(0+coeff[  5]*x11+2*coeff[  7]*x31+coeff[ 10]+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+2*coeff[ 13]*x31*x42+coeff[ 15]*x11*x41+2*coeff[ 16]*x31*x41+coeff[ 18]*x11*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  6]*x21+2*coeff[  8]*x41+2*coeff[ 13]*x41*x32+coeff[ 14]*x12+coeff[ 15]*x11*x31+coeff[ 16]*x32+coeff[ 18]*x11*x21*x31+coeff[ 19]*x22,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.3947516E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.22348891E-07,-0.57507142E-01, 0.69509491E-01, 0.39617335E-02,
        -0.71257828E-02,-0.65887882E-02, 0.11824809E-01,-0.78381629E-04,
         0.13004107E-03, 0.62197192E-04,-0.26274030E-04,-0.25839644E-03,
         0.34174198E-03, 0.11750825E-03, 0.11776925E-03,-0.26552597E-03,
        -0.11908164E-03, 0.69350947E-03,-0.71351830E-03, 0.72032897E-04,
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
        +coeff[ 13]*x11        *x43
        +coeff[ 14]*x11*x23        
        +coeff[ 15]*x11*x22    *x41
        +coeff[ 16]*x12    *x31*x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x11    *x32*x41
        +coeff[ 18]        *x33*x41
        +coeff[ 19]    *x21*x31*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[ 10]*x22+coeff[ 11]*x21*x32+coeff[ 13]*x43+coeff[ 14]*x23+coeff[ 15]*x22*x41+2*coeff[ 16]*x11*x31*x41+coeff[ 17]*x32*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  7]+2*coeff[ 10]*x21*x11+coeff[ 11]*x11*x32+coeff[ 12]*x33+3*coeff[ 14]*x22*x11+2*coeff[ 15]*x21*x11*x41+coeff[ 19]*x31*x42,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x42+2*coeff[ 11]*x31*x11*x21+3*coeff[ 12]*x32*x21+coeff[ 16]*x12*x41+2*coeff[ 17]*x31*x11*x41+3*coeff[ 18]*x32*x41+coeff[ 19]*x21*x42,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  8]+2*coeff[  9]*x41*x31+3*coeff[ 13]*x42*x11+coeff[ 15]*x11*x22+coeff[ 16]*x12*x31+coeff[ 17]*x11*x32+coeff[ 18]*x33+2*coeff[ 19]*x41*x21*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
