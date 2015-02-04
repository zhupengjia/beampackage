#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 16;
    float avdat=  0.1208478E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 17]={
         0.43641485E-04,-0.40235130E+02, 0.66453529E+02, 0.94636196E+00,
        -0.15779893E+01, 0.28885868E+01,-0.17340155E+01,-0.32455273E-01,
         0.53916696E-01, 0.13452040E-01,-0.57596741E-02, 0.16814124E-01,
        -0.32554887E-01,-0.22769528E-02, 0.35528064E-01,-0.19191816E-01,
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

    float v_target_x                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]        *x31    
        +coeff[  3]*x11*x21        
        +coeff[  4]*x11        *x41
        +coeff[  5]        *x31*x41
        +coeff[  6]    *x21*x31    
        +coeff[  7]    *x21        
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]            *x41
        +coeff[  9]        *x31*x42
        +coeff[ 10]*x11*x22        
        +coeff[ 11]*x11    *x32    
        +coeff[ 12]        *x33    
        +coeff[ 13]    *x21*x33    
        +coeff[ 14]*x11    *x34    
        +coeff[ 15]*x14    *x33    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  4]*x41+coeff[ 10]*x22+coeff[ 11]*x32+coeff[ 14]*x34+4*coeff[ 15]*x13*x33,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  6]*x31+coeff[  7]+2*coeff[ 10]*x21*x11+coeff[ 13]*x33,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  5]*x41+coeff[  6]*x21+coeff[  9]*x42+2*coeff[ 11]*x31*x11+3*coeff[ 12]*x32+3*coeff[ 13]*x32*x21+4*coeff[ 14]*x33*x11+3*coeff[ 15]*x32*x14,2)
        +dx4*dx4*pow(0+coeff[  4]*x11+coeff[  5]*x31+coeff[  8]+2*coeff[  9]*x41*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 17;
    float avdat=  0.4730060E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 18]={
        -0.33230975E-01,-0.39878941E+02, 0.65828423E+02, 0.33496416E+00,
        -0.12156650E+01, 0.22799745E-02,-0.18007526E-02, 0.10954813E+01,
         0.28831398E-01,-0.41306622E-01,-0.47678977E+00, 0.14590634E+01,
        -0.10745407E+01,-0.12966401E-01, 0.57562059E-02, 0.51306998E-02,
        -0.21242858E-02,
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
        +coeff[  5]        *x34    
        +coeff[  6]*x14    *x32    
        +coeff[  7]            *x42
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x11            
        +coeff[  9]        *x31    
        +coeff[ 10]*x12            
        +coeff[ 11]*x11    *x31    
        +coeff[ 12]        *x32    
        +coeff[ 13]        *x32*x41
        +coeff[ 14]*x12*x21        
        +coeff[ 15]            *x43
        +coeff[ 16]    *x23        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  6]*x13*x32+coeff[  8]+2*coeff[ 10]*x11+coeff[ 11]*x31+2*coeff[ 14]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  3]*x21+coeff[  4]*x41+coeff[ 14]*x12+3*coeff[ 16]*x22,2)
        +dx3*dx3*pow(0+4*coeff[  5]*x33+2*coeff[  6]*x31*x14+coeff[  9]+coeff[ 11]*x11+2*coeff[ 12]*x31+2*coeff[ 13]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  4]*x21+2*coeff[  7]*x41+coeff[ 13]*x32+3*coeff[ 15]*x42,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.8357406E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.17293135E-03, 0.64404523E-02,-0.49785962E-02,-0.27202333E-02,
         0.16203512E-01,-0.13070500E-01, 0.38982583E-02,-0.75483769E-01,
        -0.75719357E+00, 0.84087446E-01, 0.25368462E+01,-0.29729149E+01,
         0.11965495E+01, 0.45415934E-01, 0.64324751E-01,-0.52405261E-01,
        -0.10980803E+00, 0.68466730E-01, 0.43921210E-01,-0.67179501E-01,
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
        +coeff[  5]        *x32    
        +coeff[  6]    *x21    *x41
        +coeff[  7]*x12*x21        
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]    *x23        
        +coeff[  9]*x12        *x41
        +coeff[ 10]    *x22    *x41
        +coeff[ 11]    *x21    *x42
        +coeff[ 12]            *x43
        +coeff[ 13]*x14*x21        
        +coeff[ 14]*x12*x23        
        +coeff[ 15]*x14        *x41
        +coeff[ 16]*x12*x22    *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x24    *x41
        +coeff[ 18]*x12*x21    *x42
        +coeff[ 19]    *x23    *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  2]*x11+coeff[  4]*x31+2*coeff[  7]*x11*x21+2*coeff[  9]*x11*x41+4*coeff[ 13]*x13*x21+2*coeff[ 14]*x11*x23+4*coeff[ 15]*x13*x41+2*coeff[ 16]*x11*x22*x41+2*coeff[ 18]*x11*x21*x42,2)
        +dx2*dx2*pow(0+2*coeff[  3]*x21+coeff[  6]*x41+coeff[  7]*x12+3*coeff[  8]*x22+2*coeff[ 10]*x21*x41+coeff[ 11]*x42+coeff[ 13]*x14+3*coeff[ 14]*x22*x12+2*coeff[ 16]*x21*x12*x41+4*coeff[ 17]*x23*x41+coeff[ 18]*x12*x42+3*coeff[ 19]*x22*x42,2)
        +dx3*dx3*pow(0+coeff[  4]*x11+2*coeff[  5]*x31,2)
        +dx4*dx4*pow(0+coeff[  1]+coeff[  6]*x21+coeff[  9]*x12+coeff[ 10]*x22+2*coeff[ 11]*x41*x21+3*coeff[ 12]*x42+coeff[ 15]*x14+coeff[ 16]*x12*x22+coeff[ 17]*x24+2*coeff[ 18]*x41*x12*x21+2*coeff[ 19]*x41*x23,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1182643E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.28554062E-07,-0.58046218E-01, 0.70603020E-01, 0.56321002E-02,
        -0.10298776E-01,-0.94945095E-02, 0.17273918E-01,-0.19226034E-03,
         0.32171339E-03, 0.12722681E-04,-0.95335039E-04,-0.43583527E-05,
         0.49968826E-03, 0.35025103E-04,-0.19939958E-04,-0.76008597E-04,
         0.16389017E-03, 0.19305696E-03,-0.27545996E-03,-0.33832504E-03,
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
    float x24 = x23*x2;
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
        +coeff[  7]    *x21        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]            *x41
        +coeff[  9]        *x31*x44
        +coeff[ 10]        *x33*x41
        +coeff[ 11]*x11*x24        
        +coeff[ 12]        *x31*x42
        +coeff[ 13]*x13*x21        
        +coeff[ 14]*x11*x22    *x41
        +coeff[ 15]*x11*x21    *x42
        +coeff[ 16]        *x31*x43
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x11*x22        
        +coeff[ 18]*x11*x21    *x41
        +coeff[ 19]    *x21*x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[ 11]*x24+3*coeff[ 13]*x12*x21+coeff[ 14]*x22*x41+coeff[ 15]*x21*x42+coeff[ 17]*x22+coeff[ 18]*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  7]+4*coeff[ 11]*x23*x11+coeff[ 13]*x13+2*coeff[ 14]*x21*x11*x41+coeff[ 15]*x11*x42+2*coeff[ 17]*x21*x11+coeff[ 18]*x11*x41+coeff[ 19]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x44+3*coeff[ 10]*x32*x41+coeff[ 12]*x42+coeff[ 16]*x43+coeff[ 19]*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  8]+4*coeff[  9]*x43*x31+coeff[ 10]*x33+2*coeff[ 12]*x41*x31+coeff[ 14]*x11*x22+2*coeff[ 15]*x41*x11*x21+3*coeff[ 16]*x42*x31+coeff[ 18]*x11*x21+coeff[ 19]*x21*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
