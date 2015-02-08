#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 19;
    float avdat= -0.4818825E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 20]={
         0.34093824E-02,-0.37834362E+02, 0.56771369E+01, 0.63114079E+02,
        -0.93046389E+01,-0.60064539E-01, 0.85521728E-01, 0.11377884E+00,
         0.14729758E+00,-0.16582236E+00,-0.11128739E+00,-0.46117913E-01,
        -0.22114443E-01,-0.48591252E-01, 0.63099302E-01,-0.38528959E-02,
         0.26315276E-02, 0.92558302E-02,-0.11654451E-01,
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
    int ncoeff= 14;
    float avdat=  0.2312269E+02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 15]={
        -0.91260346E-02,-0.56214800E+01,-0.37873535E+02, 0.92197695E+01,
         0.63192749E+02, 0.17314428E+00, 0.53827368E-01, 0.19516537E-01,
        -0.12626247E-01,-0.19374527E+00, 0.25211716E-01,-0.19609187E-01,
         0.32987485E-02,-0.17337612E-02,
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
        +coeff[  7]    *x21    *x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]    *x22        
        +coeff[  9]*x11    *x31    
        +coeff[ 10]*x11        *x41
        +coeff[ 11]*x11*x21        
        +coeff[ 12]*x11        *x42
        +coeff[ 13]*x11*x24        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  6]*x11+coeff[  9]*x31+coeff[ 10]*x41+coeff[ 11]*x21+coeff[ 12]*x42+coeff[ 13]*x24,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  7]*x41+2*coeff[  8]*x21+coeff[ 11]*x11+4*coeff[ 13]*x23*x11,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  5]*x31+coeff[  9]*x11,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x21+coeff[ 10]*x11+2*coeff[ 12]*x41*x11,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4052441E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.86391356E-05,-0.32046527E-01,-0.49881358E-01, 0.49372900E-01,
         0.57125710E-01,-0.17761676E-03, 0.11963940E-04, 0.39931090E-03,
         0.11129084E-03,-0.29352261E-04, 0.33920012E-04, 0.15636248E-03,
         0.14828316E-03, 0.35929155E-04,-0.38374928E-05, 0.11315462E-03,
        -0.23098396E-03,-0.23351469E-03,-0.90619296E-05,-0.25940986E-03,
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
        +coeff[  7]        *x31*x42
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]    *x23    *x43
        +coeff[  9]*x11*x21        
        +coeff[ 10]    *x21*x31    
        +coeff[ 11]*x11*x22        
        +coeff[ 12]        *x33    
        +coeff[ 13]        *x33*x43
        +coeff[ 14]*x12            
        +coeff[ 15]*x13            
        +coeff[ 16]*x12    *x31    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x11*x21    *x41
        +coeff[ 18]    *x22    *x41
        +coeff[ 19]    *x21*x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  9]*x21+coeff[ 11]*x22+2*coeff[ 14]*x11+3*coeff[ 15]*x12+2*coeff[ 16]*x11*x31+coeff[ 17]*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+2*coeff[  6]*x21+3*coeff[  8]*x22*x43+coeff[  9]*x11+coeff[ 10]*x31+2*coeff[ 11]*x21*x11+coeff[ 17]*x11*x41+2*coeff[ 18]*x21*x41+coeff[ 19]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x42+coeff[ 10]*x21+3*coeff[ 12]*x32+3*coeff[ 13]*x32*x43+coeff[ 16]*x12+coeff[ 19]*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+4*coeff[  5]*x43+2*coeff[  7]*x41*x31+3*coeff[  8]*x42*x23+3*coeff[ 13]*x42*x33+coeff[ 17]*x11*x21+coeff[ 18]*x22+coeff[ 19]*x21*x31,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.2411477E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.57151268E-04,-0.50035413E-01, 0.32656122E-01, 0.57240423E-01,
        -0.50390825E-01, 0.13513581E-02,-0.10669561E-02,-0.19811901E-03,
         0.67710271E-03, 0.16201747E-03,-0.43304617E-03,-0.58035133E-03,
        -0.16474398E-03,-0.92789087E-04, 0.70333824E-03,-0.66810625E-03,
        -0.58238391E-04, 0.39804821E-04,-0.11951922E-03,-0.18894834E-04,
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
        +coeff[  9]*x11*x21        
        +coeff[ 10]    *x22        
        +coeff[ 11]        *x32    
        +coeff[ 12]*x11        *x41
        +coeff[ 13]        *x32*x41
        +coeff[ 14]    *x21    *x42
        +coeff[ 15]            *x43
        +coeff[ 16]        *x31*x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x12*x21        
        +coeff[ 18]    *x23        
        +coeff[ 19]        *x31*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  7]*x11+coeff[  8]*x31+coeff[  9]*x21+coeff[ 12]*x41+2*coeff[ 17]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x41+coeff[  9]*x11+2*coeff[ 10]*x21+coeff[ 14]*x42+coeff[ 17]*x12+3*coeff[ 18]*x22,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  8]*x11+2*coeff[ 11]*x31+2*coeff[ 13]*x31*x41+coeff[ 16]*x41+coeff[ 19]*x42,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  5]*x21+2*coeff[  6]*x41+coeff[ 12]*x11+coeff[ 13]*x32+2*coeff[ 14]*x41*x21+3*coeff[ 15]*x42+coeff[ 16]*x31+2*coeff[ 19]*x41*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
