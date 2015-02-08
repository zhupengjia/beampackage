#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 12;
    float avdat=  0.6633394E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 13]={
        -0.22446113E-04,-0.39239883E+02, 0.65164581E+02, 0.16527215E+01,
         0.54724091E+00,-0.99826092E+00,-0.90673774E+00,-0.11272131E-01,
         0.18647004E-01, 0.13955253E-01,-0.45747057E-01, 0.37718121E-01,
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
    int ncoeff= 20;
    float avdat=  0.3655740E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.23521552E-01,-0.38928326E+02, 0.64632324E+02, 0.23908507E-03,
         0.20355858E+00,-0.71635467E+00, 0.63139242E+00, 0.41022335E-03,
        -0.13457376E-02, 0.26551276E-02,-0.95455401E-03,-0.63943135E-03,
        -0.26902753E+00, 0.78250575E+00,-0.54479527E+00,-0.65980749E-02,
        -0.79425871E-02, 0.30801559E-02, 0.49849884E-02, 0.22052047E-02,
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
        +coeff[  7]        *x34    
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x12*x21    *x41
        +coeff[  9]*x11*x21*x31*x41
        +coeff[ 10]    *x21*x32*x41
        +coeff[ 11]*x14    *x32    
        +coeff[ 12]*x12            
        +coeff[ 13]*x11    *x31    
        +coeff[ 14]        *x32    
        +coeff[ 15]        *x31    
        +coeff[ 16]        *x32*x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]*x12        *x43
        +coeff[ 18]*x11            
        +coeff[ 19]*x12*x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  8]*x11*x21*x41+coeff[  9]*x21*x31*x41+4*coeff[ 11]*x13*x32+2*coeff[ 12]*x11+coeff[ 13]*x31+2*coeff[ 17]*x11*x43+coeff[ 18]+2*coeff[ 19]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  5]*x41+coeff[  8]*x12*x41+coeff[  9]*x11*x31*x41+coeff[ 10]*x32*x41+coeff[ 19]*x12,2)
        +dx3*dx3*pow(0+4*coeff[  7]*x33+coeff[  9]*x11*x21*x41+2*coeff[ 10]*x31*x21*x41+2*coeff[ 11]*x31*x14+coeff[ 13]*x11+2*coeff[ 14]*x31+coeff[ 15]+2*coeff[ 16]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+4*coeff[  3]*x43+coeff[  5]*x21+2*coeff[  6]*x41+coeff[  8]*x12*x21+coeff[  9]*x11*x21*x31+coeff[ 10]*x21*x32+coeff[ 16]*x32+3*coeff[ 17]*x42*x12,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.8684783E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.63728905E-04,-0.54413211E-01, 0.64651951E-01,-0.33863313E-02,
         0.99002256E-03, 0.10906124E-01,-0.41570030E-02,-0.87274052E-02,
         0.40634037E-02, 0.10525627E-03,-0.16821490E-03,-0.17666815E-03,
         0.25822711E-03,-0.97097676E-04, 0.52534137E-03,-0.56655589E-03,
         0.15465649E-03,-0.25718217E-03, 0.48893424E-04, 0.61955830E-05,
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
    float avdat=  0.3949004E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.22601707E-07,-0.57530001E-01, 0.69536991E-01, 0.39723776E-02,
        -0.71384152E-02,-0.66005946E-02, 0.11839512E-01,-0.78439523E-04,
         0.13008750E-03, 0.61065599E-04,-0.25850466E-04,-0.25288039E-03,
         0.33551134E-03, 0.11199110E-03, 0.11248745E-03,-0.25475572E-03,
        -0.13306399E-03, 0.71935286E-03,-0.72546105E-03, 0.70682967E-04,
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
