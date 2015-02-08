#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 12;
    float avdat=  0.6626595E+00;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 13]={
        -0.31075251E-04,-0.39139900E+02, 0.65043732E+02, 0.16327300E+01,
         0.54067707E+00,-0.98628187E+00,-0.89573318E+00,-0.11147065E-01,
         0.18435029E-01, 0.13803155E-01,-0.45273315E-01, 0.37333041E-01,
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
    float x41 = x4;
    float x42 = x41*x4;

//                 function

    float v_target_x                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]        *x31    
        +coeff[  3]        *x31*x41
        +coeff[  4]*x11*x21        
        +coeff[  5]    *x21*x31    
        +coeff[  6]*x11        *x41
        +coeff[  7]    *x21        
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]            *x41
        +coeff[  9]    *x22*x31    
        +coeff[ 10]    *x21*x31*x41
        +coeff[ 11]        *x31*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  4]*x21+coeff[  6]*x41,2)
        +dx2*dx2*pow(0+coeff[  4]*x11+coeff[  5]*x31+coeff[  7]+2*coeff[  9]*x21*x31+coeff[ 10]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  3]*x41+coeff[  5]*x21+coeff[  9]*x22+coeff[ 10]*x21*x41+coeff[ 11]*x42,2)
        +dx4*dx4*pow(0+coeff[  3]*x31+coeff[  6]*x11+coeff[  8]+coeff[ 10]*x21*x31+2*coeff[ 11]*x41*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.3503791E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.23606796E-01,-0.38834766E+02, 0.64520309E+02,-0.11162079E-03,
         0.20286299E+00,-0.71306205E+00, 0.62767470E+00, 0.59531006E-03,
         0.44718729E-02,-0.11564376E-01, 0.75565781E-02,-0.70197828E-03,
        -0.26565644E+00, 0.76972777E+00,-0.53349334E+00,-0.65519037E-02,
        -0.73752012E-02, 0.49814908E-02, 0.33525629E-02, 0.12585418E-02,
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
    float x44 = x43*x4;

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]    *x21        
        +coeff[  2]            *x41
        +coeff[  3]            *x44
        +coeff[  4]    *x22        
        +coeff[  5]    *x21    *x41
        +coeff[  6]            *x42
        +coeff[  7]        *x34    
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x12*x21    *x41
        +coeff[  9]*x11*x21*x31*x41
        +coeff[ 10]    *x21*x32*x41
        +coeff[ 11]*x14    *x32    
        +coeff[ 12]*x12            
        +coeff[ 13]*x11    *x31    
        +coeff[ 14]        *x32    
        +coeff[ 15]        *x31    
        +coeff[ 16]        *x32*x41
    ;
    v_target_y                                =v_target_y                                
        +coeff[ 17]*x11            
        +coeff[ 18]*x12*x21        
        +coeff[ 19]    *x22    *x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  8]*x11*x21*x41+coeff[  9]*x21*x31*x41+4*coeff[ 11]*x13*x32+2*coeff[ 12]*x11+coeff[ 13]*x31+coeff[ 17]+2*coeff[ 18]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  5]*x41+coeff[  8]*x12*x41+coeff[  9]*x11*x31*x41+coeff[ 10]*x32*x41+coeff[ 18]*x12+2*coeff[ 19]*x21*x41,2)
        +dx3*dx3*pow(0+4*coeff[  7]*x33+coeff[  9]*x11*x21*x41+2*coeff[ 10]*x31*x21*x41+2*coeff[ 11]*x31*x14+coeff[ 13]*x11+2*coeff[ 14]*x31+coeff[ 15]+2*coeff[ 16]*x31*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+4*coeff[  3]*x43+coeff[  5]*x21+2*coeff[  6]*x41+coeff[  8]*x12*x21+coeff[  9]*x11*x21*x31+coeff[ 10]*x21*x32+coeff[ 16]*x32+coeff[ 19]*x22,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.8761276E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.63698084E-04,-0.54413874E-01, 0.64651005E-01,-0.33838956E-02,
         0.98515861E-03, 0.10900175E-01,-0.41407836E-02,-0.87237982E-02,
         0.40509258E-02, 0.10427052E-03,-0.16686678E-03,-0.18207880E-03,
         0.26160257E-03,-0.84314910E-04, 0.50021510E-03,-0.55061298E-03,
         0.15406412E-03,-0.25380895E-03, 0.47863869E-04, 0.74822738E-05,
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
        +coeff[ 13]*x12        *x41
        +coeff[ 14]*x11    *x31*x41
        +coeff[ 15]        *x32*x41
        +coeff[ 16]    *x22*x32    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x21*x32*x41
        +coeff[ 18]*x12        *x42
        +coeff[ 19]    *x21    *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+2*coeff[  3]*x11+coeff[  5]*x31+coeff[  9]+coeff[ 11]*x21*x31+2*coeff[ 13]*x11*x41+coeff[ 14]*x31*x41+2*coeff[ 18]*x11*x42,2)
        +dx2*dx2*pow(0+coeff[  1]+2*coeff[  4]*x21+coeff[  6]*x41+coeff[ 11]*x11*x31+coeff[ 12]*x32+2*coeff[ 16]*x21*x32+coeff[ 17]*x32*x41+coeff[ 19]*x42,2)
        +dx3*dx3*pow(0+coeff[  5]*x11+2*coeff[  7]*x31+coeff[ 10]+coeff[ 11]*x11*x21+2*coeff[ 12]*x31*x21+coeff[ 14]*x11*x41+2*coeff[ 15]*x31*x41+2*coeff[ 16]*x31*x22+2*coeff[ 17]*x31*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  2]+coeff[  6]*x21+2*coeff[  8]*x41+coeff[ 13]*x12+coeff[ 14]*x11*x31+coeff[ 15]*x32+coeff[ 17]*x21*x32+2*coeff[ 18]*x41*x12+2*coeff[ 19]*x41*x21,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.3949251E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18000E+02,-0.18017E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.22604956E-07,-0.57533816E-01, 0.69541596E-01, 0.39736410E-02,
        -0.71398546E-02,-0.66019176E-02, 0.11841082E-01,-0.78448902E-04,
         0.13009558E-03, 0.60807200E-04,-0.25747526E-04,-0.25091408E-03,
         0.33308318E-03, 0.11082672E-03, 0.11138435E-03,-0.25238469E-03,
        -0.13445887E-03, 0.71972772E-03,-0.72385429E-03, 0.70236463E-04,
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
        +coeff[ 11]*x11*x21*x32    
        +coeff[ 12]    *x21*x33    
        +coeff[ 13]*x11        *x43
        +coeff[ 14]*x11*x23        
        +coeff[ 15]*x11*x22    *x41
        +coeff[ 16]*x12    *x31*x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]*x11    *x32*x41
        +coeff[ 18]        *x33*x41
        +coeff[ 19]    *x21*x31*x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  3]*x21+coeff[  5]*x41+coeff[ 10]*x22+coeff[ 11]*x21*x32+coeff[ 13]*x43+coeff[ 14]*x23+coeff[ 15]*x22*x41+2*coeff[ 16]*x11*x31*x41+coeff[ 17]*x32*x41,2)
        +dx2*dx2*pow(0+coeff[  3]*x11+coeff[  4]*x31+coeff[  7]+2*coeff[ 10]*x21*x11+coeff[ 11]*x11*x32+coeff[ 12]*x33+3*coeff[ 14]*x22*x11+2*coeff[ 15]*x21*x11*x41+coeff[ 19]*x31*x42,2)
        +dx3*dx3*pow(0+coeff[  2]+coeff[  4]*x21+coeff[  6]*x41+coeff[  9]*x42+2*coeff[ 11]*x31*x11*x21+3*coeff[ 12]*x32*x21+coeff[ 16]*x12*x41+2*coeff[ 17]*x31*x11*x41+3*coeff[ 18]*x32*x41+coeff[ 19]*x21*x42,2)
        +dx4*dx4*pow(0+coeff[  5]*x11+coeff[  6]*x31+coeff[  8]+2*coeff[  9]*x41*x31+3*coeff[ 13]*x42*x11+coeff[ 15]*x11*x22+coeff[ 16]*x12*x31+coeff[ 17]*x11*x32+coeff[ 18]*x33+2*coeff[ 19]*x41*x21*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
