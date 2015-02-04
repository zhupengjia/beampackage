#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 19;
    float avdat= -0.4738893E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 20]={
         0.32039545E-02,-0.37667664E+02, 0.55687985E+01, 0.62923306E+02,
        -0.91373014E+01,-0.59005085E-01, 0.84228024E-01, 0.11180373E+00,
         0.14285327E+00,-0.16336235E+00,-0.10769519E+00,-0.44085946E-01,
        -0.21345116E-01,-0.47225036E-01, 0.60630105E-01,-0.37672182E-02,
         0.26158907E-02, 0.82978699E-02,-0.10464624E-01,
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
        +coeff[  8]    *x21    *x41
        +coeff[  9]        *x31*x41
        +coeff[ 10]            *x42
        +coeff[ 11]        *x32    
        +coeff[ 12]*x12            
        +coeff[ 13]    *x22        
        +coeff[ 14]*x11    *x31    
        +coeff[ 15]    *x21*x32    
        +coeff[ 16]*x12        *x41
    ;
    v_target_x                                =v_target_x                                
        +coeff[ 17]    *x21    *x42
        +coeff[ 18]            *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41+2*coeff[ 12]*x11+coeff[ 14]*x31+2*coeff[ 16]*x11*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  7]*x31+coeff[  8]*x41+2*coeff[ 13]*x21+coeff[ 15]*x32+coeff[ 17]*x42,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x21+coeff[  9]*x41+2*coeff[ 11]*x31+coeff[ 14]*x11+2*coeff[ 15]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  8]*x21+coeff[  9]*x31+2*coeff[ 10]*x41+coeff[ 16]*x12+2*coeff[ 17]*x41*x21+3*coeff[ 18]*x42,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 13;
    float avdat=  0.2298697E+02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 14]={
        -0.92038466E-02,-0.55143533E+01,-0.37707043E+02, 0.90549812E+01,
         0.63002014E+02, 0.17375329E+00, 0.54285366E-01, 0.19979041E-01,
        -0.12949491E-01,-0.19480489E+00, 0.24914071E-01,-0.19344132E-01,
         0.17603704E-02,
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

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]        *x32    
        +coeff[  6]*x12            
        +coeff[  7]    *x21    *x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]    *x22        
        +coeff[  9]*x11    *x31    
        +coeff[ 10]*x11        *x41
        +coeff[ 11]*x11*x21        
        +coeff[ 12]        *x33*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  6]*x11+coeff[  9]*x31+coeff[ 10]*x41+coeff[ 11]*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  7]*x41+2*coeff[  8]*x21+coeff[ 11]*x11,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  5]*x31+coeff[  9]*x11+3*coeff[ 12]*x32*x42,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x21+coeff[ 10]*x11+2*coeff[ 12]*x41*x33,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4107630E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.16897910E-04,-0.31854764E-01,-0.50046127E-01, 0.49166918E-01,
         0.57370473E-01,-0.34874006E-05,-0.38419239E-03, 0.39652895E-03,
         0.19584215E-04, 0.10462097E-02,-0.72132098E-03, 0.12686379E-03,
         0.49204882E-04,-0.14483975E-04,-0.44881690E-05, 0.55236287E-04,
        -0.72881456E-04,-0.45550289E-03, 0.16512313E-04,-0.13135867E-04,
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
    float x33 = x32*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;

//                 function

    float v_target_theta                            =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]            *x44
        +coeff[  6]    *x22        
        +coeff[  7]*x11    *x32*x42
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]    *x21*x31    
        +coeff[  9]    *x21    *x41
        +coeff[ 10]            *x42
        +coeff[ 11]*x13*x22        
        +coeff[ 12]        *x33*x43
        +coeff[ 13]*x11*x21        
        +coeff[ 14]*x11    *x31    
        +coeff[ 15]    *x22    *x41
        +coeff[ 16]    *x21    *x42
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x12*x21*x31*x41
        +coeff[ 18]        *x31*x44
        +coeff[ 19]*x13*x23        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  7]*x32*x42+3*coeff[ 11]*x12*x22+coeff[ 13]*x21+coeff[ 14]*x31+2*coeff[ 17]*x11*x21*x31*x41+3*coeff[ 19]*x12*x23,2)
        +dx2*dx2*pow(0+coeff[  2]+2*coeff[  6]*x21+coeff[  8]*x31+coeff[  9]*x41+2*coeff[ 11]*x21*x13+coeff[ 13]*x11+2*coeff[ 15]*x21*x41+coeff[ 16]*x42+coeff[ 17]*x12*x31*x41+3*coeff[ 19]*x22*x13,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  7]*x31*x11*x42+coeff[  8]*x21+3*coeff[ 12]*x32*x43+coeff[ 14]*x11+coeff[ 17]*x12*x21*x41+coeff[ 18]*x44,2)
        +dx4*dx4*pow(0+coeff[  4]+4*coeff[  5]*x43+2*coeff[  7]*x41*x11*x32+coeff[  9]*x21+2*coeff[ 10]*x41+3*coeff[ 12]*x42*x33+coeff[ 15]*x22+2*coeff[ 16]*x41*x21+coeff[ 17]*x12*x21*x31+4*coeff[ 18]*x43*x31,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.2393771E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.54827582E-04,-0.50168227E-01, 0.32434233E-01, 0.57446141E-01,
        -0.50132450E-01, 0.12975635E-02,-0.10211422E-02,-0.19223765E-03,
         0.65754261E-03,-0.12802718E-03,-0.41736261E-03,-0.56413759E-03,
        -0.89418878E-04, 0.67656278E-03,-0.64300984E-03, 0.80768135E-04,
        -0.36181962E-04, 0.38343238E-04,-0.11483623E-03,-0.14989791E-04,
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
        +coeff[  5]    *x21    *x41
        +coeff[  6]            *x42
        +coeff[  7]*x12            
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]*x11    *x31    
        +coeff[  9]        *x31*x43
        +coeff[ 10]    *x22        
        +coeff[ 11]        *x32    
        +coeff[ 12]        *x32*x41
        +coeff[ 13]    *x21    *x42
        +coeff[ 14]            *x43
        +coeff[ 15]*x11*x23        
        +coeff[ 16]*x11        *x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x12*x21        
        +coeff[ 18]    *x23        
        +coeff[ 19]*x11        *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  7]*x11+coeff[  8]*x31+coeff[ 15]*x23+coeff[ 16]*x41+2*coeff[ 17]*x11*x21+coeff[ 19]*x42,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x41+2*coeff[ 10]*x21+coeff[ 13]*x42+3*coeff[ 15]*x22*x11+coeff[ 17]*x12+3*coeff[ 18]*x22,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  8]*x11+coeff[  9]*x43+2*coeff[ 11]*x31+2*coeff[ 12]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  5]*x21+2*coeff[  6]*x41+3*coeff[  9]*x42*x31+coeff[ 12]*x32+2*coeff[ 13]*x41*x21+3*coeff[ 14]*x42+coeff[ 16]*x11+2*coeff[ 19]*x41*x11,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
