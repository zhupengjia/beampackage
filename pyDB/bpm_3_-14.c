#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 12;
    float avdat=  0.1237758E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 13]={
        -0.65833668E-03,-0.37899681E+02, 0.37500329E+01, 0.63367073E+02,
        -0.61604438E+01,-0.36709167E-01, 0.55649545E-01, 0.66012412E-01,
        -0.10113548E+00,-0.39087194E-02, 0.58258213E-02,-0.55220462E-02,
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
        +coeff[  9]    *x21*x32*x41
        +coeff[ 10]*x12*x22        
        +coeff[ 11]*x12*x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41+2*coeff[ 10]*x11*x22+2*coeff[ 11]*x11*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  7]*x31+coeff[  9]*x32*x41+2*coeff[ 10]*x21*x12+coeff[ 11]*x12*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x21+coeff[  8]*x41+2*coeff[  9]*x31*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  8]*x31+coeff[  9]*x21*x32+coeff[ 11]*x12*x21,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 11;
    float avdat= -0.1397799E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
        -0.68252617E-02,-0.37336247E+01,-0.37898949E+02, 0.61336670E+01,
         0.63365860E+02,-0.14324167E+00, 0.12400502E+00, 0.54383161E-02,
        -0.22973416E-01, 0.41376181E-01, 0.23528269E-01,
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
        +coeff[  5]*x11    *x31    
        +coeff[  6]        *x32    
        +coeff[  7]    *x22        
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]    *x21    *x41
        +coeff[  9]*x12            
        +coeff[ 10]            *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x31+2*coeff[  9]*x11,2)
        +dx2*dx2*pow(0+coeff[  2]+2*coeff[  7]*x21+coeff[  8]*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x11+2*coeff[  6]*x31,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  8]*x21+2*coeff[ 10]*x41,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1794413E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.33807939E-05,-0.21758484E-01,-0.53595312E-01, 0.33636115E-01,
         0.63103825E-01,-0.73540403E-04, 0.11159950E-03,-0.17881418E-03,
        -0.14938906E-03, 0.11590434E-03, 0.11188471E-03, 0.66711131E-04,
         0.19361749E-04, 0.48538499E-04, 0.56002220E-04,-0.84363337E-05,
        -0.11035437E-03,-0.34471000E-05,-0.42447150E-05, 0.13748117E-05,
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
        +coeff[  5]*x11*x21        
        +coeff[  6]*x11        *x41
        +coeff[  7]        *x31*x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x11    *x31    
        +coeff[  9]    *x21*x31    
        +coeff[ 10]        *x32    
        +coeff[ 11]        *x33    
        +coeff[ 12]        *x31*x42
        +coeff[ 13]*x12            
        +coeff[ 14]*x13            
        +coeff[ 15]*x11*x22        
        +coeff[ 16]*x12    *x31    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x21*x32    
        +coeff[ 18]            *x43
        +coeff[ 19]    *x22        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41+coeff[  8]*x31+2*coeff[ 13]*x11+3*coeff[ 14]*x12+coeff[ 15]*x22+2*coeff[ 16]*x11*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  9]*x31+2*coeff[ 15]*x21*x11+coeff[ 17]*x32+2*coeff[ 19]*x21,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x41+coeff[  8]*x11+coeff[  9]*x21+2*coeff[ 10]*x31+3*coeff[ 11]*x32+coeff[ 12]*x42+coeff[ 16]*x12+2*coeff[ 17]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  7]*x31+2*coeff[ 12]*x41*x31+3*coeff[ 18]*x42,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.9089716E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.71265426E-05,-0.53588562E-01, 0.22016559E-01, 0.63069411E-01,
        -0.34079585E-01, 0.18348916E-03,-0.15902171E-03, 0.14918389E-03,
         0.21905709E-04,-0.29508383E-04,-0.75731323E-05, 0.19679997E-05,
        -0.66016338E-06,-0.32636347E-05,-0.18550261E-04,-0.52896754E-04,
        -0.13791193E-03, 0.24384132E-03,-0.27339277E-03, 0.76947781E-05,
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
        +coeff[  5]*x11    *x31    
        +coeff[  6]        *x32    
        +coeff[  7]*x11        *x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]    *x21    *x41
        +coeff[  9]            *x42
        +coeff[ 10]        *x32*x41
        +coeff[ 11]*x13*x21        
        +coeff[ 12]*x11*x23        
        +coeff[ 13]*x13        *x41
        +coeff[ 14]        *x32*x43
        +coeff[ 15]*x12            
        +coeff[ 16]*x11*x21        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]    *x21*x31    
        +coeff[ 18]        *x31*x41
        +coeff[ 19]*x12*x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x31+coeff[  7]*x41+3*coeff[ 11]*x12*x21+coeff[ 12]*x23+3*coeff[ 13]*x12*x41+2*coeff[ 15]*x11+coeff[ 16]*x21+2*coeff[ 19]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  8]*x41+coeff[ 11]*x13+3*coeff[ 12]*x22*x11+coeff[ 16]*x11+coeff[ 17]*x31+coeff[ 19]*x12,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x11+2*coeff[  6]*x31+2*coeff[ 10]*x31*x41+2*coeff[ 14]*x31*x43+coeff[ 17]*x21+coeff[ 18]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+coeff[  8]*x21+2*coeff[  9]*x41+coeff[ 10]*x32+coeff[ 13]*x13+3*coeff[ 14]*x42*x32+coeff[ 18]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
