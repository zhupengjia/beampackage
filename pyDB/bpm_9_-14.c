#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.2136013E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
         0.69862159E-04,-0.38346039E+02, 0.64086540E+02, 0.15986624E+01,
         0.52331877E+00,-0.32166649E-01, 0.53429574E-01,-0.96115023E+00,
        -0.87106371E+00,
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
        +coeff[  3]        *x31*x41
        +coeff[  4]*x11*x21        
        +coeff[  5]    *x21        
        +coeff[  6]            *x41
        +coeff[  7]    *x21*x31    
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]*x11        *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  4]*x21+coeff[  8]*x41,2)
        +dx2*dx2*pow(0+coeff[  4]*x11+coeff[  5]+coeff[  7]*x31,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  3]*x41+coeff[  7]*x21,2)
        +dx4*dx4*pow(0+coeff[  3]*x31+coeff[  6]+coeff[  8]*x11,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 13;
    float avdat=  0.3892798E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 14]={
        -0.24757134E-01,-0.38211613E+02, 0.63850468E+02, 0.39737643E-04,
         0.18470141E+00,-0.66929352E+00, 0.60278416E+00, 0.22382032E-01,
        -0.31391360E-01,-0.21915385E+00, 0.65348411E+00,-0.45998740E+00,
        -0.10626253E-02,
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
        +coeff[  7]*x11            
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]        *x31    
        +coeff[  9]*x12            
        +coeff[ 10]*x11    *x31    
        +coeff[ 11]        *x32    
        +coeff[ 12]    *x21*x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  7]+2*coeff[  9]*x11+coeff[ 10]*x31,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  5]*x41+coeff[ 12]*x32,2)
        +dx3*dx3*pow(0+coeff[  8]+coeff[ 10]*x11+2*coeff[ 11]*x31+2*coeff[ 12]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+4*coeff[  3]*x43+coeff[  5]*x21+2*coeff[  6]*x41,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.7706756E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.43979262E-05, 0.43578837E-02, 0.27158663E-02,-0.18280480E-02,
        -0.36980701E-02, 0.26206921E-02,-0.93288720E-01,-0.11178124E+01,
         0.10764994E+00, 0.38987522E+01,-0.46317053E+01, 0.18578880E+01,
         0.30254643E-01, 0.14237638E+01,-0.34909949E-01,-0.50457749E+01,
         0.45628585E-01, 0.59980397E+01,-0.44775221E-01,-0.23881733E+01,
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
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;

//                 function

    float v_target_theta                            =avdat
        +coeff[  0]                
        +coeff[  1]            *x41
        +coeff[  2]*x12            
        +coeff[  3]    *x22        
        +coeff[  4]*x11    *x31    
        +coeff[  5]    *x21    *x41
        +coeff[  6]*x12*x21        
        +coeff[  7]    *x23        
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x12        *x41
        +coeff[  9]    *x22    *x41
        +coeff[ 10]    *x21    *x42
        +coeff[ 11]            *x43
        +coeff[ 12]*x14*x21        
        +coeff[ 13]*x12*x23        
        +coeff[ 14]*x14        *x41
        +coeff[ 15]*x12*x22    *x41
        +coeff[ 16]    *x24    *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12*x21    *x42
        +coeff[ 18]    *x23    *x42
        +coeff[ 19]*x12        *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  2]*x11+coeff[  4]*x31+2*coeff[  6]*x11*x21+2*coeff[  8]*x11*x41+4*coeff[ 12]*x13*x21+2*coeff[ 13]*x11*x23+4*coeff[ 14]*x13*x41+2*coeff[ 15]*x11*x22*x41+2*coeff[ 17]*x11*x21*x42+2*coeff[ 19]*x11*x43,2)
        +dx2*dx2*pow(0+2*coeff[  3]*x21+coeff[  5]*x41+coeff[  6]*x12+3*coeff[  7]*x22+2*coeff[  9]*x21*x41+coeff[ 10]*x42+coeff[ 12]*x14+3*coeff[ 13]*x22*x12+2*coeff[ 15]*x21*x12*x41+4*coeff[ 16]*x23*x41+coeff[ 17]*x12*x42+3*coeff[ 18]*x22*x42,2)
        +dx3*dx3*pow(0+coeff[  4]*x11,2)
        +dx4*dx4*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  8]*x12+coeff[  9]*x22+2*coeff[ 10]*x41*x21+3*coeff[ 11]*x42+coeff[ 14]*x14+coeff[ 15]*x12*x22+coeff[ 16]*x24+2*coeff[ 17]*x41*x12*x21+2*coeff[ 18]*x41*x23+3*coeff[ 19]*x42*x12,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.2306969E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18018E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.66938888E-07,-0.57274122E-01, 0.69201767E-01, 0.38210528E-02,
        -0.69664959E-02,-0.64595118E-02, 0.11696957E-01, 0.38595410E-03,
        -0.22888249E-03, 0.36780028E-04, 0.51412826E-04, 0.20718587E-05,
         0.13245939E-03,-0.27215610E-04,-0.54600673E-05, 0.27717808E-04,
         0.96146425E-04,-0.20031213E-03, 0.10810144E-03,-0.20455763E-03,
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
    float x44 = x43*x4;

//                 function

    float v_target_phi                              =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]        *x31    
        +coeff[  3]*x11*x21        
        +coeff[  4]    *x21*x31    
        +coeff[  5]*x11        *x41
        +coeff[  6]        *x31*x41
        +coeff[  7]            *x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]    *x21        
        +coeff[  9]        *x31*x44
        +coeff[ 10]*x11*x21*x32    
        +coeff[ 11]    *x21*x33    
        +coeff[ 12]*x11        *x43
        +coeff[ 13]*x11*x23    *x41
        +coeff[ 14]    *x21*x32    
        +coeff[ 15]        *x31*x42
        +coeff[ 16]*x11*x23        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x11*x22    *x41
        +coeff[ 18]*x12    *x31*x41
        +coeff[ 19]*x11    *x32*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[ 10]*x21*x32+coeff[ 12]*x43+coeff[ 13]*x23*x41+coeff[ 16]*x23+coeff[ 17]*x22*x41+2*coeff[ 18]*x11*x31*x41+coeff[ 19]*x32*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  8]+coeff[ 10]*x11*x32+coeff[ 11]*x33+3*coeff[ 13]*x22*x11*x41+coeff[ 14]*x32+3*coeff[ 16]*x22*x11+2*coeff[ 17]*x21*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x44+2*coeff[ 10]*x31*x11*x21+3*coeff[ 11]*x32*x21+2*coeff[ 14]*x31*x21+coeff[ 15]*x42+coeff[ 18]*x12*x41+2*coeff[ 19]*x31*x11*x41,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  7]+4*coeff[  9]*x43*x31+3*coeff[ 12]*x42*x11+coeff[ 13]*x11*x23+2*coeff[ 15]*x41*x31+coeff[ 17]*x11*x22+coeff[ 18]*x12*x31+coeff[ 19]*x11*x32,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
