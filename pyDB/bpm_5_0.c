#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.1790319E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
         0.63137586E-05,-0.39382095E+02, 0.65306908E+02, 0.43441033E+00,
        -0.72221911E+00,-0.79527116E+00, 0.13208510E+01, 0.39553825E-01,
        -0.23821346E-01,
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
        +coeff[  3]*x11*x21        
        +coeff[  4]*x11        *x41
        +coeff[  5]    *x21*x31    
        +coeff[  6]        *x31*x41
        +coeff[  7]            *x41
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]    *x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  4]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]*x31+coeff[  8],2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  5]*x21+coeff[  6]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]*x11+coeff[  6]*x31+coeff[  7],2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 16;
    float avdat=  0.3526662E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 17]={
        -0.16621113E-01,-0.39190487E+02, 0.64979134E+02, 0.15603992E+00,
        -0.56043983E+00,-0.46052605E+00, 0.50152624E+00,-0.43107473E-03,
         0.53009666E-04,-0.20796459E-01,-0.21615823E+00, 0.64443403E+00,
         0.15324904E-01, 0.45976071E-02,-0.63788076E-02, 0.14582445E-02,
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

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]    *x22        
        +coeff[  4]    *x21    *x41
        +coeff[  5]        *x32    
        +coeff[  6]            *x42
        +coeff[  7]*x14            
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]        *x34    
        +coeff[  9]        *x31    
        +coeff[ 10]*x12            
        +coeff[ 11]*x11    *x31    
        +coeff[ 12]*x11            
        +coeff[ 13]*x11*x21*x31    
        +coeff[ 14]    *x21*x32    
        +coeff[ 15]            *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  7]*x13+2*coeff[ 10]*x11+coeff[ 11]*x31+coeff[ 12]+coeff[ 13]*x21*x31,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  3]*x21+coeff[  4]*x41+coeff[ 13]*x11*x31+coeff[ 14]*x32,2)
        +dx3*dx3*pow(0+2*coeff[  5]*x31+4*coeff[  8]*x33+coeff[  9]+coeff[ 11]*x11+coeff[ 13]*x11*x21+2*coeff[ 14]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  4]*x21+2*coeff[  6]*x41+3*coeff[ 15]*x42,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.6202480E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.45493704E-04,-0.55311322E-01, 0.66010721E-01, 0.80056107E-02,
        -0.64322841E-02,-0.31208822E-02, 0.74379134E-03, 0.21637636E-03,
        -0.34485047E-03,-0.24715587E-02, 0.30562545E-02,-0.64829459E-04,
         0.10978914E-03,-0.12935656E-03,-0.31915744E-03, 0.45420384E-03,
         0.45561741E-03,-0.64743758E-03, 0.54771735E-04,-0.16959303E-05,
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
        +coeff[  3]*x11    *x31    
        +coeff[  4]        *x32    
        +coeff[  5]    *x21    *x41
        +coeff[  6]    *x22        
        +coeff[  7]*x11            
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]        *x31    
        +coeff[  9]*x12            
        +coeff[ 10]            *x42
        +coeff[ 11]*x11*x21*x31    
        +coeff[ 12]    *x21*x32    
        +coeff[ 13]        *x32*x41
        +coeff[ 14]*x11*x21*x31*x41
        +coeff[ 15]    *x21*x32*x41
        +coeff[ 16]*x11    *x31*x42
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]        *x32*x42
        +coeff[ 18]*x12        *x41
        +coeff[ 19]*x11*x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  3]*x31+coeff[  7]+2*coeff[  9]*x11+coeff[ 11]*x21*x31+coeff[ 14]*x21*x31*x41+coeff[ 16]*x31*x42+2*coeff[ 18]*x11*x41+coeff[ 19]*x21,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  5]*x41+2*coeff[  6]*x21+coeff[ 11]*x11*x31+coeff[ 12]*x32+coeff[ 14]*x11*x31*x41+coeff[ 15]*x32*x41+coeff[ 19]*x11,2)
        +dx3*dx3*pow(0+coeff[  3]*x11+2*coeff[  4]*x31+coeff[  8]+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+2*coeff[ 13]*x31*x41+coeff[ 14]*x11*x21*x41+2*coeff[ 15]*x31*x21*x41+coeff[ 16]*x11*x42+2*coeff[ 17]*x31*x42,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  5]*x21+2*coeff[ 10]*x41+coeff[ 13]*x32+coeff[ 14]*x11*x21*x31+coeff[ 15]*x21*x32+2*coeff[ 16]*x41*x11*x31+2*coeff[ 17]*x41*x32+coeff[ 18]*x12,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1157775E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18013E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.10655735E-07,-0.57073589E-01, 0.68776883E-01, 0.28866159E-02,
        -0.52297642E-02,-0.48527080E-02, 0.87516392E-02,-0.15369957E-03,
         0.25751023E-03,-0.24750456E-04, 0.63995205E-04,-0.12249475E-03,
         0.17248289E-03, 0.15560185E-03, 0.21245069E-03,-0.30277797E-03,
        -0.43908594E-03, 0.31329971E-03,-0.23957589E-05,-0.17305043E-04,
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
        +coeff[  9]    *x21*x31*x41
        +coeff[ 10]        *x31*x42
        +coeff[ 11]*x12*x21*x31    
        +coeff[ 12]*x11*x21*x32    
        +coeff[ 13]    *x22*x31*x41
        +coeff[ 14]*x11    *x32*x41
        +coeff[ 15]        *x33*x41
        +coeff[ 16]    *x21*x31*x42
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]        *x31*x43
        +coeff[ 18]*x11*x21*x31    
        +coeff[ 19]*x11*x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+2*coeff[ 11]*x11*x21*x31+coeff[ 12]*x21*x32+coeff[ 14]*x32*x41+coeff[ 18]*x21*x31+coeff[ 19]*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  7]+coeff[  9]*x31*x41+coeff[ 11]*x12*x31+coeff[ 12]*x11*x32+2*coeff[ 13]*x21*x31*x41+coeff[ 16]*x31*x42+coeff[ 18]*x11*x31+coeff[ 19]*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x21*x41+coeff[ 10]*x42+coeff[ 11]*x12*x21+2*coeff[ 12]*x31*x11*x21+coeff[ 13]*x22*x41+2*coeff[ 14]*x31*x11*x41+3*coeff[ 15]*x32*x41+coeff[ 16]*x21*x42+coeff[ 17]*x43+coeff[ 18]*x11*x21,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  8]+coeff[  9]*x21*x31+2*coeff[ 10]*x41*x31+coeff[ 13]*x22*x31+coeff[ 14]*x11*x32+coeff[ 15]*x33+2*coeff[ 16]*x41*x21*x31+3*coeff[ 17]*x42*x31+coeff[ 19]*x11*x21,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
