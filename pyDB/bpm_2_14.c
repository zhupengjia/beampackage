#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.5436230E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.48330333E-02,-0.39070366E+02, 0.65111666E+01, 0.64522881E+02,
        -0.10584745E+02,-0.66846877E-01, 0.92868358E-01,-0.12514404E-02,
        -0.59959650E-01, 0.12669688E+00,-0.18089432E+00,-0.28673466E-01,
         0.83851524E-01, 0.18238951E+00,-0.13787946E+00,-0.62927589E-01,
         0.15296741E-01,-0.19613069E-01,-0.45022960E-02, 0.17920663E-02,
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
    float x34 = x33*x3;
    float x41 = x4;
    float x42 = x41*x4;
    float x43 = x42*x4;
    float x44 = x43*x4;

//                 function

    float v_target_x                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11*x21        
        +coeff[  6]*x11        *x41
        +coeff[  7]            *x44
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]    *x22        
        +coeff[  9]    *x21*x31    
        +coeff[ 10]        *x31*x41
        +coeff[ 11]*x12            
        +coeff[ 12]*x11    *x31    
        +coeff[ 13]    *x21    *x41
        +coeff[ 14]            *x42
        +coeff[ 15]        *x32    
        +coeff[ 16]    *x21    *x42
    ;
    v_target_x                                =v_target_x                                
        +coeff[ 17]            *x43
        +coeff[ 18]    *x21*x34    
        +coeff[ 19]*x12*x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41+2*coeff[ 11]*x11+coeff[ 12]*x31+2*coeff[ 19]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+2*coeff[  8]*x21+coeff[  9]*x31+coeff[ 13]*x41+coeff[ 16]*x42+coeff[ 18]*x34+coeff[ 19]*x12,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  9]*x21+coeff[ 10]*x41+coeff[ 12]*x11+2*coeff[ 15]*x31+4*coeff[ 18]*x33*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+4*coeff[  7]*x43+coeff[ 10]*x31+coeff[ 13]*x21+2*coeff[ 14]*x41+2*coeff[ 16]*x41*x21+3*coeff[ 17]*x42,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff= 14;
    float avdat=  0.2408240E+02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 15]={
        -0.89786351E-02,-0.64418197E+01,-0.39106441E+02, 0.10476825E+02,
         0.64599350E+02, 0.17293526E+00, 0.53395424E-01,-0.92865136E-02,
        -0.19310041E+00, 0.26956381E-01, 0.15047486E-01,-0.20857375E-01,
         0.64898413E-02,-0.14244443E-02,
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
    float x44 = x43*x4;

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]        *x32    
        +coeff[  6]*x12            
        +coeff[  7]    *x22        
    ;
    v_target_y                                =v_target_y                                
        +coeff[  8]*x11    *x31    
        +coeff[  9]*x11        *x41
        +coeff[ 10]    *x21    *x41
        +coeff[ 11]*x11*x21        
        +coeff[ 12]        *x33*x44
        +coeff[ 13]*x13*x22        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  6]*x11+coeff[  8]*x31+coeff[  9]*x41+coeff[ 11]*x21+3*coeff[ 13]*x12*x22,2)
        +dx2*dx2*pow(0+coeff[  2]+2*coeff[  7]*x21+coeff[ 10]*x41+coeff[ 11]*x11+2*coeff[ 13]*x21*x13,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  5]*x31+coeff[  8]*x11+3*coeff[ 12]*x32*x44,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  9]*x11+coeff[ 10]*x21+4*coeff[ 12]*x43*x33,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.3636806E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.18171037E-04,-0.33632044E-01,-0.48830267E-01, 0.51181894E-01,
         0.55511553E-01, 0.10376864E-02,-0.71643846E-03, 0.84303014E-04,
        -0.38171947E-03, 0.32458785E-04,-0.97753464E-05,-0.35960995E-04,
         0.13937657E-03, 0.23666424E-04,-0.30256197E-04, 0.10619048E-03,
        -0.21665166E-03,-0.86750870E-05, 0.54478038E-04,-0.71678005E-04,
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

    float v_target_theta                            =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]    *x21    *x41
        +coeff[  6]            *x42
        +coeff[  7]        *x31*x42
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]    *x22        
        +coeff[  9]    *x21*x31    
        +coeff[ 10]        *x32    
        +coeff[ 11]*x11*x22        
        +coeff[ 12]        *x33    
        +coeff[ 13]        *x33*x41
        +coeff[ 14]*x11*x21        
        +coeff[ 15]*x13            
        +coeff[ 16]*x12    *x31    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x21*x32    
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]            *x43
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[ 11]*x22+coeff[ 14]*x21+3*coeff[ 15]*x12+2*coeff[ 16]*x11*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x41+2*coeff[  8]*x21+coeff[  9]*x31+2*coeff[ 11]*x21*x11+coeff[ 14]*x11+coeff[ 17]*x32+coeff[ 18]*x42,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x42+coeff[  9]*x21+2*coeff[ 10]*x31+3*coeff[ 12]*x32+3*coeff[ 13]*x32*x41+coeff[ 16]*x12+2*coeff[ 17]*x31*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  5]*x21+2*coeff[  6]*x41+2*coeff[  7]*x41*x31+coeff[ 13]*x33+2*coeff[ 18]*x41*x21+3*coeff[ 19]*x42,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.2533993E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.54472788E-04,-0.48984755E-01, 0.34211226E-01, 0.55629451E-01,
        -0.52157454E-01, 0.12500786E-02,-0.98803639E-03,-0.19616443E-03,
         0.66991278E-03,-0.40076350E-03,-0.57313283E-03,-0.16901022E-03,
        -0.12833595E-03, 0.67123526E-03,-0.63895324E-03, 0.15008760E-03,
         0.25113031E-05,-0.11358857E-03,-0.28569078E-04, 0.70400987E-04,
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
        +coeff[  9]    *x22        
        +coeff[ 10]        *x32    
        +coeff[ 11]*x11        *x41
        +coeff[ 12]        *x32*x41
        +coeff[ 13]    *x21    *x42
        +coeff[ 14]            *x43
        +coeff[ 15]*x11*x21        
        +coeff[ 16]*x12*x21        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]    *x23        
        +coeff[ 18]        *x31*x41
        +coeff[ 19]*x11*x21*x31    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  7]*x11+coeff[  8]*x31+coeff[ 11]*x41+coeff[ 15]*x21+2*coeff[ 16]*x11*x21+coeff[ 19]*x21*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x41+2*coeff[  9]*x21+coeff[ 13]*x42+coeff[ 15]*x11+coeff[ 16]*x12+3*coeff[ 17]*x22+coeff[ 19]*x11*x31,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  8]*x11+2*coeff[ 10]*x31+2*coeff[ 12]*x31*x41+coeff[ 18]*x41+coeff[ 19]*x11*x21,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  5]*x21+2*coeff[  6]*x41+coeff[ 11]*x11+coeff[ 12]*x32+2*coeff[ 13]*x41*x21+3*coeff[ 14]*x42+coeff[ 18]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
