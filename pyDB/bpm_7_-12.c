#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.1176918E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
         0.39342427E-04,-0.38687878E+02, 0.64571449E+02, 0.79579461E+00,
         0.24256828E+01,-0.27302086E-01,-0.14597168E+01,-0.13234000E+01,
         0.45306094E-01,
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
    float avdat=  0.4717299E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.35573743E-01,-0.38395794E+02, 0.64054413E+02,-0.10250223E+01,
         0.91864556E+00, 0.28480151E+00,-0.85188881E-01, 0.70048667E-01,
         0.20946451E-01,-0.28640117E-01,-0.73696208E+00,-0.11151000E+00,
         0.24953368E+00,-0.13792077E+00, 0.10491796E+01, 0.14799912E-01,
        -0.35485005E+00,-0.65669906E-02, 0.75881872E-02,-0.49504815E-02,
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
        +coeff[  6]*x11    *x33    
        +coeff[  7]        *x34    
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x11            
        +coeff[  9]        *x31    
        +coeff[ 10]        *x32    
        +coeff[ 11]*x14    *x32    
        +coeff[ 12]*x13    *x33    
        +coeff[ 13]*x12    *x34    
        +coeff[ 14]*x11    *x31    
        +coeff[ 15]*x14            
        +coeff[ 16]*x12            
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]    *x21*x32    
        +coeff[ 18]*x12*x21        
        +coeff[ 19]*x12        *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  6]*x33+coeff[  8]+4*coeff[ 11]*x13*x32+3*coeff[ 12]*x12*x33+2*coeff[ 13]*x11*x34+coeff[ 14]*x31+4*coeff[ 15]*x13+2*coeff[ 16]*x11+2*coeff[ 18]*x11*x21+2*coeff[ 19]*x11*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+coeff[  3]*x41+2*coeff[  5]*x21+coeff[ 17]*x32+coeff[ 18]*x12,2)
        +dx3*dx3*pow(0+3*coeff[  6]*x32*x11+4*coeff[  7]*x33+coeff[  9]+2*coeff[ 10]*x31+2*coeff[ 11]*x31*x14+3*coeff[ 12]*x32*x13+4*coeff[ 13]*x33*x12+coeff[ 14]*x11+2*coeff[ 17]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  3]*x21+2*coeff[  4]*x41+coeff[ 19]*x12,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.9313676E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.17097403E-03, 0.64415508E-02,-0.49076839E-02,-0.27224668E-02,
         0.16022112E-01,-0.12953234E-01, 0.39008611E-02,-0.75503819E-01,
        -0.75710100E+00, 0.84112041E-01, 0.25365117E+01,-0.29724948E+01,
         0.11963677E+01, 0.45419347E-01, 0.64323284E-01,-0.52408017E-01,
        -0.10979564E+00, 0.68452410E-01, 0.43909468E-01,-0.67164510E-01,
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
    float avdat=  0.1182570E-02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18027E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.29930181E-07,-0.58045428E-01, 0.70603281E-01, 0.58005964E-02,
        -0.10501818E-01,-0.97414143E-02, 0.17569594E-01,-0.19284803E-03,
         0.32214398E-03, 0.23875826E-03,-0.33253626E-03, 0.34653503E-03,
        -0.27269832E-03,-0.17305676E-03,-0.24685063E-03,-0.95431949E-03,
         0.66433102E-03, 0.33359166E-03, 0.34578887E-03, 0.10646420E-03,
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
