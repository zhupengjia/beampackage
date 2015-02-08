#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.1175047E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
         0.30512228E-04,-0.38595741E+02, 0.64459351E+02, 0.78696525E+00,
         0.23984890E+01,-0.27052755E-01,-0.14435251E+01,-0.13084998E+01,
         0.44864040E-01,
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
        +coeff[  4]        *x31*x41
        +coeff[  5]    *x21        
        +coeff[  6]    *x21*x31    
        +coeff[  7]*x11        *x41
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]            *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  7]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]+coeff[  6]*x31,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x41+coeff[  6]*x21,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  7]*x11+coeff[  8],2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4701689E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.35758674E-01,-0.38308296E+02, 0.63948692E+02,-0.10145059E+01,
         0.90846586E+00, 0.28206238E+00, 0.55594830E-03, 0.48947628E-02,
        -0.13178140E-01, 0.91903545E-02,-0.97970793E-03, 0.20437224E-01,
        -0.27849616E-01,-0.33731869E+00, 0.10034401E+01,-0.70581752E+00,
        -0.27928010E-02, 0.59767622E-02, 0.12322493E-02,-0.73418193E-02,
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
        +coeff[  6]        *x34    
        +coeff[  7]*x12*x21    *x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x11*x21*x31*x41
        +coeff[  9]    *x21*x32*x41
        +coeff[ 10]*x14    *x32    
        +coeff[ 11]*x11            
        +coeff[ 12]        *x31    
        +coeff[ 13]*x12            
        +coeff[ 14]*x11    *x31    
        +coeff[ 15]        *x32    
        +coeff[ 16]    *x21*x32    
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]*x12*x21        
        +coeff[ 18]    *x22    *x41
        +coeff[ 19]*x11    *x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  7]*x11*x21*x41+coeff[  8]*x21*x31*x41+4*coeff[ 10]*x13*x32+coeff[ 11]+2*coeff[ 13]*x11+coeff[ 14]*x31+2*coeff[ 17]*x11*x21+coeff[ 19]*x31*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  5]*x21+coeff[  7]*x12*x41+coeff[  8]*x11*x31*x41+coeff[  9]*x32*x41+coeff[ 16]*x32+coeff[ 17]*x12+2*coeff[ 18]*x21*x41,2)
        +dx3*dx3*pow(0+4*coeff[  6]*x33+coeff[  8]*x11*x21*x41+2*coeff[  9]*x31*x21*x41+2*coeff[ 10]*x31*x14+coeff[ 12]+coeff[ 14]*x11+2*coeff[ 15]*x31+2*coeff[ 16]*x31*x21+coeff[ 19]*x11*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41+coeff[  7]*x12*x21+coeff[  8]*x11*x21*x31+coeff[  9]*x21*x32+coeff[ 18]*x22+coeff[ 19]*x11*x31,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1036501E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.17127546E-03, 0.64421971E-02,-0.49357903E-02,-0.27216307E-02,
         0.16081879E-01,-0.12984480E-01, 0.39000839E-02,-0.75513422E-01,
        -0.75698102E+00, 0.84122397E-01, 0.25360603E+01,-0.29719479E+01,
         0.11961522E+01, 0.45426164E-01, 0.64358599E-01,-0.52416906E-01,
        -0.10987433E+00, 0.68462946E-01, 0.43955121E-01,-0.67176826E-01,
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
    float avdat=  0.1182578E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.29666761E-07,-0.58045961E-01, 0.70603915E-01, 0.58026146E-02,
        -0.10503924E-01,-0.97431354E-02, 0.17571192E-01,-0.19288399E-03,
         0.32217248E-03, 0.23795436E-03,-0.33142214E-03, 0.34524556E-03,
        -0.27130128E-03,-0.17241924E-03,-0.24564812E-03,-0.94499334E-03,
         0.65790198E-03, 0.33197796E-03, 0.34236294E-03, 0.10584341E-03,
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
        +coeff[  9]*x11*x22        
        +coeff[ 10]*x11*x21    *x41
        +coeff[ 11]        *x31*x42
        +coeff[ 12]        *x33*x41
        +coeff[ 13]    *x22*x31    
        +coeff[ 14]*x13*x21        
        +coeff[ 15]    *x22*x31*x41
        +coeff[ 16]    *x21*x31*x42
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x12*x21*x31    
        +coeff[ 18]    *x23*x31    
        +coeff[ 19]*x13        *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[  9]*x22+coeff[ 10]*x21*x41+3*coeff[ 14]*x12*x21+2*coeff[ 17]*x11*x21*x31+3*coeff[ 19]*x12*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  7]+2*coeff[  9]*x21*x11+coeff[ 10]*x11*x41+2*coeff[ 13]*x21*x31+coeff[ 14]*x13+2*coeff[ 15]*x21*x31*x41+coeff[ 16]*x31*x42+coeff[ 17]*x12*x31+3*coeff[ 18]*x22*x31,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[ 11]*x42+3*coeff[ 12]*x32*x41+coeff[ 13]*x22+coeff[ 15]*x22*x41+coeff[ 16]*x21*x42+coeff[ 17]*x12*x21+coeff[ 18]*x23,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  8]+coeff[ 10]*x11*x21+2*coeff[ 11]*x41*x31+coeff[ 12]*x33+coeff[ 15]*x22*x31+2*coeff[ 16]*x41*x21*x31+coeff[ 19]*x13,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
