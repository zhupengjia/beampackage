#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 13;
    float avdat=  0.1234657E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 14]={
        -0.30487232E-03,-0.38077770E+02, 0.38235295E+01, 0.63576630E+02,
        -0.62741103E+01,-0.37088349E-01, 0.56043938E-01, 0.66768467E-01,
        -0.10197682E+00,-0.15610778E-02,-0.21038583E-03, 0.21905401E-02,
        -0.38106705E-02,
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
    float x32 = x31*x3;
    float x33 = x32*x3;
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;

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
        +coeff[  9]        *x34    
        +coeff[ 10]            *x44
        +coeff[ 11]    *x22        
        +coeff[ 12]            *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  7]*x31+2*coeff[ 11]*x21,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x21+coeff[  8]*x41+4*coeff[  9]*x33,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  8]*x31+4*coeff[ 10]*x43+2*coeff[ 12]*x41,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 11;
    float avdat= -0.1338627E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 12]={
        -0.68365140E-02,-0.38063147E+01,-0.38077145E+02, 0.62459488E+01,
         0.63575623E+02,-0.14382723E+00, 0.12445895E+00, 0.52609686E-02,
        -0.22571528E-01, 0.41571561E-01, 0.23285935E-01,
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
    float avdat=  0.1763839E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.28496361E-05,-0.21911189E-01,-0.53529706E-01, 0.33815194E-01,
         0.63000545E-01,-0.74273390E-04, 0.11769223E-03, 0.11259916E-03,
        -0.15200586E-03, 0.11418912E-03,-0.18126934E-03, 0.73147305E-04,
        -0.84953543E-04, 0.12501205E-03, 0.49223923E-04, 0.60929411E-04,
        -0.12033424E-03,-0.35253711E-05, 0.51842162E-04,-0.76856595E-04,
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
        +coeff[  5]*x11*x21        
        +coeff[  6]    *x21*x31    
        +coeff[  7]*x11        *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x11    *x31    
        +coeff[  9]        *x32    
        +coeff[ 10]        *x31*x41
        +coeff[ 11]        *x33    
        +coeff[ 12]    *x21*x31*x41
        +coeff[ 13]        *x31*x42
        +coeff[ 14]*x12            
        +coeff[ 15]*x13            
        +coeff[ 16]*x12    *x31    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x21*x32    
        +coeff[ 18]*x11*x21    *x41
        +coeff[ 19]*x11        *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  7]*x41+coeff[  8]*x31+2*coeff[ 14]*x11+3*coeff[ 15]*x12+2*coeff[ 16]*x11*x31+coeff[ 18]*x21*x41+coeff[ 19]*x42,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  6]*x31+coeff[ 12]*x31*x41+coeff[ 17]*x32+coeff[ 18]*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  6]*x21+coeff[  8]*x11+2*coeff[  9]*x31+coeff[ 10]*x41+3*coeff[ 11]*x32+coeff[ 12]*x21*x41+coeff[ 13]*x42+coeff[ 16]*x12+2*coeff[ 17]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+coeff[ 10]*x31+coeff[ 12]*x21*x31+2*coeff[ 13]*x41*x31+coeff[ 18]*x11*x21+2*coeff[ 19]*x41*x11,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.9602526E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.69351017E-05,-0.53525135E-01, 0.22176098E-01, 0.62971555E-01,
        -0.34262396E-01, 0.18447482E-03,-0.15972098E-03, 0.15108255E-03,
        -0.15551228E-04,-0.15870062E-04,-0.90797239E-05, 0.17016279E-05,
        -0.63876314E-06,-0.29104119E-05,-0.53239022E-04,-0.13930067E-03,
         0.86352420E-05, 0.24615086E-03,-0.27665950E-03, 0.78345838E-05,
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
        +coeff[  8]            *x42
        +coeff[  9]    *x21*x32    
        +coeff[ 10]            *x43
        +coeff[ 11]*x13*x21        
        +coeff[ 12]*x11*x23        
        +coeff[ 13]*x13        *x41
        +coeff[ 14]*x12            
        +coeff[ 15]*x11*x21        
        +coeff[ 16]    *x22        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]    *x21*x31    
        +coeff[ 18]        *x31*x41
        +coeff[ 19]*x12*x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x31+coeff[  7]*x41+3*coeff[ 11]*x12*x21+coeff[ 12]*x23+3*coeff[ 13]*x12*x41+2*coeff[ 14]*x11+coeff[ 15]*x21+2*coeff[ 19]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  9]*x32+coeff[ 11]*x13+3*coeff[ 12]*x22*x11+coeff[ 15]*x11+2*coeff[ 16]*x21+coeff[ 17]*x31+coeff[ 19]*x12,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x11+2*coeff[  6]*x31+2*coeff[  9]*x31*x21+coeff[ 17]*x21+coeff[ 18]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+2*coeff[  8]*x41+3*coeff[ 10]*x42+coeff[ 13]*x13+coeff[ 18]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
