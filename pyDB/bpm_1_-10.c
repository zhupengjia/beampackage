#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 17;
    float avdat=  0.1735133E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 18]={
        -0.17353994E-03,-0.39966061E+02, 0.65988586E+02, 0.82277334E+00,
         0.24487102E+01,-0.44373125E-01, 0.72930992E-01,-0.14880599E+01,
        -0.13549029E+01, 0.33331286E-01,-0.47123618E-01,-0.14235701E-01,
        -0.28480794E-01, 0.70704229E-01,-0.71445196E-02, 0.96489033E-02,
         0.13692006E-02,
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
        +coeff[  4]        *x31*x41
        +coeff[  5]    *x21        
        +coeff[  6]            *x41
        +coeff[  7]    *x21*x31    
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]*x11        *x41
        +coeff[  9]*x11*x22        
        +coeff[ 10]*x11*x21    *x41
        +coeff[ 11]    *x22*x31    
        +coeff[ 12]    *x21*x31*x41
        +coeff[ 13]        *x31*x42
        +coeff[ 14]*x12    *x31    
        +coeff[ 15]*x11    *x32    
        +coeff[ 16]        *x32*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  8]*x41+coeff[  9]*x22+coeff[ 10]*x21*x41+2*coeff[ 14]*x11*x31+coeff[ 15]*x32,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]+coeff[  7]*x31+2*coeff[  9]*x21*x11+coeff[ 10]*x11*x41+2*coeff[ 11]*x21*x31+coeff[ 12]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x41+coeff[  7]*x21+coeff[ 11]*x22+coeff[ 12]*x21*x41+coeff[ 13]*x42+coeff[ 14]*x12+2*coeff[ 15]*x31*x11+2*coeff[ 16]*x31*x42,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  6]+coeff[  8]*x11+coeff[ 10]*x11*x21+coeff[ 12]*x21*x31+2*coeff[ 13]*x41*x31+2*coeff[ 16]*x41*x32,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.5276875E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.33035144E-01,-0.39331837E+02, 0.64922562E+02,-0.92426658E+00,
         0.83690047E+00, 0.25150564E+00,-0.81284180E-01, 0.69105968E-01,
        -0.68296021E+00,-0.73341921E-01, 0.17557897E+00,-0.10123075E+00,
        -0.19194622E-01, 0.96198004E+00,-0.16017871E-03, 0.11098673E-03,
        -0.50299703E-02, 0.10815259E-01, 0.15034854E-01,-0.31930739E+00,
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
        +coeff[  5]    *x22        
        +coeff[  6]*x11    *x33    
        +coeff[  7]        *x34    
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]        *x32    
        +coeff[  9]*x14    *x32    
        +coeff[ 10]*x13    *x33    
        +coeff[ 11]*x12    *x34    
        +coeff[ 12]        *x31    
        +coeff[ 13]*x11    *x31    
        +coeff[ 14]*x13            
        +coeff[ 15]*x12    *x31    
        +coeff[ 16]    *x21*x32    
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]*x14            
        +coeff[ 18]*x11            
        +coeff[ 19]*x12            
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  6]*x33+4*coeff[  9]*x13*x32+3*coeff[ 10]*x12*x33+2*coeff[ 11]*x11*x34+coeff[ 13]*x31+3*coeff[ 14]*x12+2*coeff[ 15]*x11*x31+4*coeff[ 17]*x13+coeff[ 18]+2*coeff[ 19]*x11,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  5]*x21+coeff[ 16]*x32,2)
        +dx3*dx3*pow(0+3*coeff[  6]*x32*x11+4*coeff[  7]*x33+2*coeff[  8]*x31+2*coeff[  9]*x31*x14+3*coeff[ 10]*x32*x13+4*coeff[ 11]*x33*x12+coeff[ 12]+coeff[ 13]*x11+coeff[ 15]*x12+2*coeff[ 16]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1187868E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.98607008E-04,-0.50741822E-01, 0.58746818E-01,-0.50295992E-02,
         0.16727285E-02, 0.16214082E-01,-0.12960944E-01,-0.66102063E-02,
         0.40581363E-03,-0.64184266E-03, 0.62360736E-02, 0.39695684E-03,
        -0.56960213E-03, 0.78185076E-04,-0.53216878E-04,-0.49920214E-04,
        -0.19526058E-03, 0.27454420E-03,-0.29715757E-05,-0.55305841E-05,
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
        +coeff[ 11]    *x21*x32*x41
        +coeff[ 12]        *x32*x42
        +coeff[ 13]*x11*x21*x31    
        +coeff[ 14]    *x21*x32    
        +coeff[ 15]*x11    *x31*x41
        +coeff[ 16]*x12*x22        
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12*x21    *x41
        +coeff[ 18]        *x31*x41
        +coeff[ 19]    *x21*x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  3]*x11+coeff[  5]*x31+coeff[  8]+coeff[ 13]*x21*x31+coeff[ 15]*x31*x41+2*coeff[ 16]*x11*x22+2*coeff[ 17]*x11*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  7]*x41+coeff[ 11]*x32*x41+coeff[ 13]*x11*x31+coeff[ 14]*x32+2*coeff[ 16]*x21*x12+coeff[ 17]*x12*x41+coeff[ 19]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  5]*x11+2*coeff[  6]*x31+coeff[  9]+2*coeff[ 11]*x31*x21*x41+2*coeff[ 12]*x31*x42+coeff[ 13]*x11*x21+2*coeff[ 14]*x31*x21+coeff[ 15]*x11*x41+coeff[ 18]*x41+coeff[ 19]*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  7]*x21+2*coeff[ 10]*x41+coeff[ 11]*x21*x32+2*coeff[ 12]*x41*x32+coeff[ 15]*x11*x31+coeff[ 17]*x12*x21+coeff[ 18]*x31+coeff[ 19]*x21*x31,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.7718283E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18023E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.90319530E-09, 0.52924980E-02, 0.58929841E-02,-0.10586272E-01,
        -0.97893178E-02, 0.17519360E-01,-0.11331893E+01,-0.50510630E-01,
         0.38583229E+01, 0.56697324E-01,-0.45431218E+01, 0.18332515E+01,
         0.22781638E-01, 0.30581018E-01, 0.18226070E+01,-0.24653697E-01,
        -0.35285328E-01,-0.63264303E+01, 0.74106526E+01,-0.29167624E+01,
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
