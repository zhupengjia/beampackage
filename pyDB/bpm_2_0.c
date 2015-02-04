#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 13;
    float avdat= -0.5082599E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 14]={
         0.26555189E-02,-0.38373188E+02, 0.60369120E+01, 0.63729485E+02,
        -0.98587456E+01,-0.63315675E-01, 0.89078210E-01, 0.11973865E+00,
         0.36383349E-01,-0.17280684E+00,-0.48869036E-01, 0.13130589E-01,
        -0.18522762E-01,
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
    float x32 = x31*x3;
    float x41 = x4;
    float x42 = x41*x4;

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
        +coeff[ 11]*x11    *x31    
        +coeff[ 12]        *x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41+coeff[ 11]*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  7]*x31+coeff[  8]*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x21+coeff[  9]*x41+coeff[ 11]*x11+2*coeff[ 12]*x31,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  8]*x21+coeff[  9]*x31+2*coeff[ 10]*x41,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 14;
    float avdat=  0.2355136E+02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 15]={
        -0.92371814E-02,-0.59728336E+01,-0.38411163E+02, 0.97595301E+01,
         0.63807404E+02, 0.17383678E+00, 0.54168776E-01, 0.10452395E-01,
        -0.19476195E+00, 0.25994845E-01,-0.20171273E-01,-0.37897723E-02,
         0.38344611E-02,-0.19413864E-02,
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
    float x24 = x23*x2;
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
        +coeff[  5]        *x32    
        +coeff[  6]*x12            
        +coeff[  7]            *x42
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x11    *x31    
        +coeff[  9]*x11        *x41
        +coeff[ 10]*x11*x21        
        +coeff[ 11]    *x22        
        +coeff[ 12]*x11        *x42
        +coeff[ 13]*x11*x24        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  6]*x11+coeff[  8]*x31+coeff[  9]*x41+coeff[ 10]*x21+coeff[ 12]*x42+coeff[ 13]*x24,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[ 10]*x11+2*coeff[ 11]*x21+4*coeff[ 13]*x23*x11,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  5]*x31+coeff[  8]*x11,2)
        +dx4*dx4*pow(0+coeff[  4]+2*coeff[  7]*x41+coeff[  9]*x11+2*coeff[ 12]*x41*x11,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.3872641E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.93905282E-05,-0.32732274E-01,-0.49419153E-01, 0.50161030E-01,
         0.56410898E-01,-0.18060621E-03, 0.12047253E-04, 0.82629180E-04,
         0.13486644E-03, 0.11311271E-03,-0.22278407E-04,-0.33870394E-05,
        -0.30787793E-03, 0.28843433E-03, 0.15105748E-04,-0.24280040E-04,
        -0.29861827E-04, 0.54322383E-04,-0.76951719E-05,-0.72916351E-04,
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
        +coeff[  7]        *x31*x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]        *x31*x42
        +coeff[  9]    *x23    *x43
        +coeff[ 10]*x11*x21        
        +coeff[ 11]*x11*x22        
        +coeff[ 12]*x11    *x32    
        +coeff[ 13]        *x33    
        +coeff[ 14]*x12            
        +coeff[ 15]*x11    *x31    
        +coeff[ 16]*x11        *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x13            
        +coeff[ 18]    *x21*x32    
        +coeff[ 19]*x11*x21    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[ 10]*x21+coeff[ 11]*x22+coeff[ 12]*x32+2*coeff[ 14]*x11+coeff[ 15]*x31+coeff[ 16]*x41+3*coeff[ 17]*x12+coeff[ 19]*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+2*coeff[  6]*x21+3*coeff[  9]*x22*x43+coeff[ 10]*x11+2*coeff[ 11]*x21*x11+coeff[ 18]*x32+coeff[ 19]*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x41+coeff[  8]*x42+2*coeff[ 12]*x31*x11+3*coeff[ 13]*x32+coeff[ 15]*x11+2*coeff[ 18]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+4*coeff[  5]*x43+coeff[  7]*x31+2*coeff[  8]*x41*x31+3*coeff[  9]*x42*x23+coeff[ 16]*x11+coeff[ 19]*x11*x21,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.2466951E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.55398301E-04,-0.49585205E-01, 0.33369187E-01, 0.56544334E-01,
        -0.51211026E-01, 0.14203591E-02,-0.11258411E-02,-0.29953205E-03,
         0.19615279E-03, 0.15759794E-03,-0.45337967E-03,-0.16610850E-03,
         0.20096102E-03, 0.14332185E-02,-0.99841726E-03,-0.47480029E-04,
        -0.99242265E-04,-0.51962264E-03, 0.38797117E-03,-0.56032010E-03,
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
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;

//                 function

    float v_target_phi                              =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]    *x21    *x41
        +coeff[  6]            *x42
        +coeff[  7]        *x34    
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]*x13    *x33    
        +coeff[  9]*x11*x21        
        +coeff[ 10]    *x22        
        +coeff[ 11]*x11        *x41
        +coeff[ 12]    *x21*x32    
        +coeff[ 13]    *x22    *x41
        +coeff[ 14]    *x21    *x42
        +coeff[ 15]        *x31*x41
        +coeff[ 16]*x12*x21        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]    *x23        
        +coeff[ 18]*x11    *x31*x41
        +coeff[ 19]        *x32*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+3*coeff[  8]*x12*x33+coeff[  9]*x21+coeff[ 11]*x41+2*coeff[ 16]*x11*x21+coeff[ 18]*x31*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x41+coeff[  9]*x11+2*coeff[ 10]*x21+coeff[ 12]*x32+2*coeff[ 13]*x21*x41+coeff[ 14]*x42+coeff[ 16]*x12+3*coeff[ 17]*x22,2)
        +dx3*dx3*pow(0+coeff[  3]+4*coeff[  7]*x33+3*coeff[  8]*x32*x13+2*coeff[ 12]*x31*x21+coeff[ 15]*x41+coeff[ 18]*x11*x41+2*coeff[ 19]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  5]*x21+2*coeff[  6]*x41+coeff[ 11]*x11+coeff[ 13]*x22+2*coeff[ 14]*x41*x21+coeff[ 15]*x31+coeff[ 18]*x11*x31+coeff[ 19]*x32,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
