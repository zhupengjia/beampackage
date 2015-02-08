#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 19;
    float avdat= -0.4751081E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 20]={
         0.32439802E-02,-0.37693146E+02, 0.55852671E+01, 0.62952484E+02,
        -0.91627703E+01,-0.59051778E-01, 0.84410697E-01, 0.11200695E+00,
         0.14309627E+00,-0.16375269E+00,-0.10800450E+00,-0.44419769E-01,
        -0.21519341E-01,-0.47243133E-01, 0.61100990E-01,-0.37480250E-02,
         0.26012689E-02, 0.85785547E-02,-0.10795593E-01,
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
    float avdat=  0.2300784E+02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 15]={
        -0.91933832E-02,-0.55303373E+01,-0.37732525E+02, 0.90797615E+01,
         0.63031212E+02, 0.17328939E+00, 0.54065462E-01, 0.20011730E-01,
        -0.12977852E-01,-0.19416219E+00, 0.25011599E-01,-0.19447805E-01,
         0.48788274E-02,-0.24000464E-02,
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
        +coeff[ 12]*x11    *x32*x42
        +coeff[ 13]*x14*x22*x31    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  6]*x11+coeff[  9]*x31+coeff[ 10]*x41+coeff[ 11]*x21+coeff[ 12]*x32*x42+4*coeff[ 13]*x13*x22*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  7]*x41+2*coeff[  8]*x21+coeff[ 11]*x11+2*coeff[ 13]*x21*x14*x31,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  5]*x31+coeff[  9]*x11+2*coeff[ 12]*x31*x11*x42+coeff[ 13]*x14*x22,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x21+coeff[ 10]*x11+2*coeff[ 12]*x41*x11*x32,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4099201E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.16954731E-04,-0.31888224E-01,-0.50024614E-01, 0.49205549E-01,
         0.57337236E-01,-0.31614036E-05,-0.38236938E-03, 0.39646396E-03,
         0.19927998E-04, 0.10428403E-02,-0.72032789E-03, 0.12647112E-03,
         0.49811730E-04,-0.14388170E-04,-0.45363495E-05, 0.55903478E-04,
        -0.73667143E-04,-0.45471510E-03, 0.16662350E-04,-0.13354766E-04,
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
    float avdat= -0.2396496E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.55214903E-04,-0.50147232E-01, 0.32468237E-01, 0.57413798E-01,
        -0.50172087E-01, 0.13040095E-02,-0.10274021E-02,-0.19304297E-03,
         0.66023716E-03,-0.12813251E-03,-0.41885441E-03,-0.56640035E-03,
        -0.89954046E-04, 0.68215438E-03,-0.64816378E-03, 0.80801990E-04,
        -0.36091315E-04, 0.38591665E-04,-0.11583963E-03,-0.14938358E-04,
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
