#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  7;
    float avdat=  0.9390720E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[  8]={
        -0.50346670E-03,-0.38187920E+02, 0.29086261E+01, 0.63760128E+02,
        -0.47729683E+01, 0.11348484E-01,-0.22395801E-01,
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
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11*x21        
        +coeff[  6]        *x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  6]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff=  7;
    float avdat=  0.8745710E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[  8]={
        -0.26567599E-02,-0.28948627E+01,-0.38189888E+02, 0.47504463E+01,
         0.63763664E+02, 0.28973518E-01,-0.13956552E-01,
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
    float x31 = x3;
    float x32 = x31*x3;
    float x41 = x4;

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]        *x32    
        +coeff[  6]*x12            
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  6]*x11,2)
        +dx2*dx2*pow(0+coeff[  2],2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  5]*x31,2)
        +dx4*dx4*pow(0+coeff[  4],2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1910166E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.12474261E-05,-0.16716288E-01,-0.54762572E-01, 0.25805537E-01,
         0.65003783E-01,-0.22767312E-04, 0.11699444E-04, 0.77126715E-04,
         0.85477723E-06, 0.57372927E-04, 0.71828217E-04, 0.29555753E-04,
        -0.98982222E-04, 0.62735321E-05, 0.11059568E-04,-0.56454719E-04,
        -0.65981294E-04,-0.16918485E-05, 0.85108904E-05,-0.11231980E-04,
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

    float v_target_theta                            =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]        *x31*x41
        +coeff[  6]*x11*x21        
        +coeff[  7]        *x32    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x14            
        +coeff[  9]*x11    *x34    
        +coeff[ 10]        *x33*x42
        +coeff[ 11]*x12            
        +coeff[ 12]*x11    *x31    
        +coeff[ 13]*x13            
        +coeff[ 14]*x13*x22        
        +coeff[ 15]*x12    *x33    
        +coeff[ 16]*x11*x21*x32*x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x21*x32    
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]            *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  6]*x21+4*coeff[  8]*x13+coeff[  9]*x34+2*coeff[ 11]*x11+coeff[ 12]*x31+3*coeff[ 13]*x12+3*coeff[ 14]*x12*x22+2*coeff[ 15]*x11*x33+coeff[ 16]*x21*x32*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  6]*x11+2*coeff[ 14]*x21*x13+coeff[ 16]*x11*x32*x41+coeff[ 17]*x32+coeff[ 18]*x42,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x41+2*coeff[  7]*x31+4*coeff[  9]*x33*x11+3*coeff[ 10]*x32*x42+coeff[ 12]*x11+3*coeff[ 15]*x32*x12+2*coeff[ 16]*x31*x11*x21*x41+2*coeff[ 17]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  5]*x31+2*coeff[ 10]*x41*x33+coeff[ 16]*x11*x21*x32+2*coeff[ 18]*x41*x21+3*coeff[ 19]*x42,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.1527863E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18004E+02,-0.18006E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.44750327E-05,-0.54766893E-01, 0.16912552E-01, 0.64990528E-01,
        -0.26143156E-01,-0.27084376E-04, 0.98556076E-04, 0.11745966E-03,
         0.20737911E-04,-0.27803528E-04, 0.14350483E-05,-0.50386745E-06,
        -0.24007461E-05, 0.34955420E-03,-0.24096087E-03,-0.10721003E-03,
         0.18957442E-03,-0.88414898E-04,-0.21518274E-03,-0.12701652E-03,
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

//                 function

    float v_target_phi                              =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x12            
        +coeff[  6]*x11    *x31    
        +coeff[  7]*x11        *x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]    *x21    *x41
        +coeff[  9]            *x42
        +coeff[ 10]*x13*x21        
        +coeff[ 11]*x11*x23        
        +coeff[ 12]*x13        *x41
        +coeff[ 13]    *x22*x32*x41
        +coeff[ 14]    *x21*x32*x42
        +coeff[ 15]*x11*x21        
        +coeff[ 16]    *x21*x31    
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]        *x32    
        +coeff[ 18]        *x31*x41
        +coeff[ 19]    *x23*x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  5]*x11+coeff[  6]*x31+coeff[  7]*x41+3*coeff[ 10]*x12*x21+coeff[ 11]*x23+3*coeff[ 12]*x12*x41+coeff[ 15]*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  8]*x41+coeff[ 10]*x13+3*coeff[ 11]*x22*x11+2*coeff[ 13]*x21*x32*x41+coeff[ 14]*x32*x42+coeff[ 15]*x11+coeff[ 16]*x31+3*coeff[ 19]*x22*x32,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  6]*x11+2*coeff[ 13]*x31*x22*x41+2*coeff[ 14]*x31*x21*x42+coeff[ 16]*x21+2*coeff[ 17]*x31+coeff[ 18]*x41+2*coeff[ 19]*x31*x23,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+coeff[  8]*x21+2*coeff[  9]*x41+coeff[ 12]*x13+coeff[ 13]*x22*x32+2*coeff[ 14]*x41*x21*x32+coeff[ 18]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
