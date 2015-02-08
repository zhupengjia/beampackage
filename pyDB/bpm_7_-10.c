#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 15;
    float avdat=  0.1178972E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 16]={
         0.37057216E-04,-0.38784473E+02, 0.64687134E+02, 0.80597597E+00,
         0.24562647E+01,-0.27648419E-01,-0.14781024E+01,-0.13403569E+01,
         0.45877095E-01, 0.14152665E-01,-0.40826220E-01,-0.69843838E-02,
         0.90779634E-02,-0.52698860E-02, 0.38859911E-01,
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
    float x32 = x31*x3;
    float x33 = x32*x3;
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
        +coeff[  6]    *x21*x31    
        +coeff[  7]*x11        *x41
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]            *x41
        +coeff[  9]    *x22*x31    
        +coeff[ 10]    *x21*x31*x41
        +coeff[ 11]*x11    *x32    
        +coeff[ 12]        *x33    
        +coeff[ 13]*x11*x21    *x41
        +coeff[ 14]        *x31*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  7]*x41+coeff[ 11]*x32+coeff[ 13]*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]+coeff[  6]*x31+2*coeff[  9]*x21*x31+coeff[ 10]*x31*x41+coeff[ 13]*x11*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x41+coeff[  6]*x21+coeff[  9]*x22+coeff[ 10]*x21*x41+2*coeff[ 11]*x31*x11+3*coeff[ 12]*x32+coeff[ 14]*x42,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  7]*x11+coeff[  8]+coeff[ 10]*x21*x31+coeff[ 13]*x11*x21+2*coeff[ 14]*x41*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4732498E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.33969257E-01,-0.38492970E+02, 0.64169716E+02,-0.10382326E+01,
         0.93062860E+00, 0.28840563E+00,-0.12207831E+01, 0.78137803E+01,
        -0.37371315E+02, 0.67029175E+02, 0.15978609E+02, 0.21446101E-01,
        -0.29450102E-01,-0.68514895E+00, 0.18483069E+01,-0.53439453E+02,
        -0.67614685E-02, 0.78594666E-02,-0.51588370E-02, 0.12783625E-02,
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
        +coeff[  6]        *x32    
        +coeff[  7]*x14            
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x13    *x31    
        +coeff[  9]*x12    *x32    
        +coeff[ 10]        *x34    
        +coeff[ 11]*x11            
        +coeff[ 12]        *x31    
        +coeff[ 13]*x12            
        +coeff[ 14]*x11    *x31    
        +coeff[ 15]*x11    *x33    
        +coeff[ 16]    *x21*x32    
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]*x12*x21        
        +coeff[ 18]*x12        *x41
        +coeff[ 19]    *x22    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  7]*x13+3*coeff[  8]*x12*x31+2*coeff[  9]*x11*x32+coeff[ 11]+2*coeff[ 13]*x11+coeff[ 14]*x31+coeff[ 15]*x33+2*coeff[ 17]*x11*x21+2*coeff[ 18]*x11*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  5]*x21+coeff[ 16]*x32+coeff[ 17]*x12+2*coeff[ 19]*x21*x41,2)
        +dx3*dx3*pow(0+2*coeff[  6]*x31+coeff[  8]*x13+2*coeff[  9]*x31*x12+4*coeff[ 10]*x33+coeff[ 12]+coeff[ 14]*x11+3*coeff[ 15]*x32*x11+2*coeff[ 16]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41+coeff[ 18]*x12+coeff[ 19]*x22,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.8162675E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.17099638E-03, 0.64412681E-02,-0.49014878E-02,-0.27223020E-02,
         0.16007632E-01,-0.12945268E-01, 0.39007787E-02,-0.75507231E-01,
        -0.75702214E+00, 0.84116250E-01, 0.25362561E+01,-0.29722130E+01,
         0.11962638E+01, 0.45424718E-01, 0.64290375E-01,-0.52414879E-01,
        -0.10972859E+00, 0.68452939E-01, 0.43876544E-01,-0.67167021E-01,
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
    float avdat=  0.1182564E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.29500598E-07,-0.58044918E-01, 0.70602678E-01, 0.57981173E-02,
        -0.10499128E-01,-0.97391028E-02, 0.17567212E-01,-0.19280807E-03,
         0.32211203E-03, 0.23950190E-03,-0.33356113E-03, 0.34776059E-03,
        -0.27367100E-03,-0.17366522E-03,-0.24763870E-03,-0.96372236E-03,
         0.67080028E-03, 0.33465566E-03, 0.34925333E-03, 0.10691698E-03,
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
