#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.6731912E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
        -0.22483178E-04,-0.40674686E+02, 0.66900352E+02, 0.64354694E+00,
         0.19439750E+01,-0.11730441E+01,-0.10679458E+01, 0.21938613E-01,
        -0.13257444E-01,
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
        +coeff[  5]    *x21*x31    
        +coeff[  6]*x11        *x41
        +coeff[  7]            *x41
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]    *x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  6]*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  5]*x31+coeff[  8],2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x41+coeff[  5]*x21,2)
        +dx4*dx4*pow(0+coeff[  4]*x31+coeff[  6]*x11+coeff[  7],2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 19;
    float avdat=  0.5689921E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 20]={
        -0.21988930E-01,-0.40292496E+02, 0.66252861E+02, 0.74494427E+00,
         0.23691931E+00,-0.76903200E+00,-0.44932740E-03,-0.35665786E+00,
         0.10639725E+01,-0.84048021E+00, 0.72799828E-02,-0.10406538E-01,
        -0.27789472E-01, 0.37874371E-01,-0.19520964E-01,-0.79827541E-02,
         0.90814710E-01,-0.88432491E-01, 0.10716435E-01,
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
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]            *x42
        +coeff[  4]    *x22        
        +coeff[  5]        *x32    
        +coeff[  6]*x14            
        +coeff[  7]*x12            
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x11    *x31    
        +coeff[  9]    *x21    *x41
        +coeff[ 10]*x11            
        +coeff[ 11]        *x31    
        +coeff[ 12]*x11*x21*x31    
        +coeff[ 13]    *x21*x32    
        +coeff[ 14]*x12        *x41
        +coeff[ 15]    *x22    *x41
        +coeff[ 16]*x11    *x31*x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]        *x32*x41
        +coeff[ 18]    *x21    *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+4*coeff[  6]*x13+2*coeff[  7]*x11+coeff[  8]*x31+coeff[ 10]+coeff[ 12]*x21*x31+2*coeff[ 14]*x11*x41+coeff[ 16]*x31*x41,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  9]*x41+coeff[ 12]*x11*x31+coeff[ 13]*x32+2*coeff[ 15]*x21*x41+coeff[ 18]*x42,2)
        +dx3*dx3*pow(0+2*coeff[  5]*x31+coeff[  8]*x11+coeff[ 11]+coeff[ 12]*x11*x21+2*coeff[ 13]*x31*x21+coeff[ 16]*x11*x41+2*coeff[ 17]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+2*coeff[  3]*x41+coeff[  9]*x21+coeff[ 14]*x12+coeff[ 15]*x22+coeff[ 16]*x11*x31+coeff[ 17]*x32+2*coeff[ 18]*x41*x21,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.7587097E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.59738351E-04,-0.54363795E-01, 0.64592645E-01,-0.32853298E-02,
         0.10996829E-02, 0.10672694E-01,-0.44115721E-02,-0.85832439E-02,
         0.42187045E-02, 0.10729975E-03,-0.17068036E-03, 0.21523732E-03,
         0.52299572E-04,-0.29892343E-03,-0.53925382E-03, 0.14669778E-03,
        -0.18515075E-03, 0.62380353E-03,-0.17857879E-03, 0.56214103E-05,
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
        +coeff[ 13]        *x32*x41
        +coeff[ 14]        *x32*x42
        +coeff[ 15]*x12        *x41
        +coeff[ 16]*x12*x22        
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]*x11*x21*x31*x41
        +coeff[ 18]*x12*x21        
        +coeff[ 19]    *x22    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  3]*x11+coeff[  5]*x31+coeff[  9]+coeff[ 11]*x21*x31+2*coeff[ 15]*x11*x41+2*coeff[ 16]*x11*x22+coeff[ 17]*x21*x31*x41+2*coeff[ 18]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  6]*x41+coeff[ 11]*x11*x31+coeff[ 12]*x32+2*coeff[ 16]*x21*x12+coeff[ 17]*x11*x31*x41+coeff[ 18]*x12+2*coeff[ 19]*x21*x41,2)
        +dx3*dx3*pow(0+coeff[  5]*x11+2*coeff[  7]*x31+coeff[ 10]+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+2*coeff[ 13]*x31*x41+2*coeff[ 14]*x31*x42+coeff[ 17]*x11*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  6]*x21+2*coeff[  8]*x41+coeff[ 13]*x32+2*coeff[ 14]*x41*x32+coeff[ 15]*x12+coeff[ 17]*x11*x21*x31+coeff[ 19]*x22,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.3945712E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.22313438E-07,-0.57482023E-01, 0.69476612E-01, 0.38894997E-02,
        -0.70322398E-02,-0.64990618E-02, 0.11714112E-01,-0.78317462E-04,
         0.12999844E-03, 0.27847206E-03, 0.10082289E-03, 0.11770023E-03,
        -0.11770229E-03, 0.21960362E-03,-0.41015108E-04, 0.31066695E-03,
        -0.44311525E-03,-0.56738987E-04,-0.14584018E-03,-0.18763977E-03,
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
        +coeff[ 11]    *x21*x33    
        +coeff[ 12]    *x21*x31*x42
        +coeff[ 13]        *x31*x43
        +coeff[ 14]*x13*x21        
        +coeff[ 15]*x11    *x32*x41
        +coeff[ 16]        *x33*x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x11        *x43
        +coeff[ 18]*x11*x21    *x41
        +coeff[ 19]    *x21*x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[ 10]*x22+3*coeff[ 14]*x12*x21+coeff[ 15]*x32*x41+coeff[ 17]*x43+coeff[ 18]*x21*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  7]+2*coeff[ 10]*x21*x11+coeff[ 11]*x33+coeff[ 12]*x31*x42+coeff[ 14]*x13+coeff[ 18]*x11*x41+coeff[ 19]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x42+3*coeff[ 11]*x32*x21+coeff[ 12]*x21*x42+coeff[ 13]*x43+2*coeff[ 15]*x31*x11*x41+3*coeff[ 16]*x32*x41+coeff[ 19]*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  8]+2*coeff[  9]*x41*x31+2*coeff[ 12]*x41*x21*x31+3*coeff[ 13]*x42*x31+coeff[ 15]*x11*x32+coeff[ 16]*x33+3*coeff[ 17]*x42*x11+coeff[ 18]*x11*x21+coeff[ 19]*x21*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
