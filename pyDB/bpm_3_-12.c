#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff=  9;
    float avdat=  0.1236288E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 10]={
        -0.87424723E-03,-0.37984703E+02, 0.37850506E+01, 0.63467129E+02,
        -0.62146215E+01,-0.36824390E-01, 0.55777323E-01, 0.66302627E-01,
        -0.10146715E+00,
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
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11*x21        
        +coeff[  6]*x11        *x41
        +coeff[  7]    *x21*x31    
    ;
    v_target_x                                =v_target_x                                
        +coeff[  8]        *x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_x                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  6]*x41,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  7]*x31,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  7]*x21+coeff[  8]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  6]*x11+coeff[  8]*x31,2)
        );
    dx[0]=dv_target_x                                ;
    return v_target_x                                ;
}
#include <math.h>
float target_y                                (float *x,float *dx,int m){
    int ncoeff=  7;
    float avdat= -0.1369434E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[  8]={
        -0.46303957E-02,-0.37682531E+01,-0.37984043E+02, 0.61871805E+01,
         0.63466038E+02,-0.45517422E-01, 0.66252351E-01,
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

//                 function

    float v_target_y                                =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11    *x31    
        +coeff[  6]        *x32    
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x31,2)
        +dx2*dx2*pow(0+coeff[  2],2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x11+2*coeff[  6]*x31,2)
        +dx4*dx4*pow(0+coeff[  4],2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.1779821E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
        -0.27864760E-05,-0.21831026E-01,-0.53563274E-01, 0.33721723E-01,
         0.63052408E-01,-0.74007359E-04, 0.11695299E-03, 0.11237301E-03,
        -0.15084111E-03, 0.11312551E-03,-0.18040874E-03, 0.70809867E-04,
        -0.23755081E-04, 0.48896942E-04, 0.59155391E-04, 0.63917028E-05,
        -0.11672908E-03,-0.49579521E-04,-0.35025319E-05, 0.78409248E-04,
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

//                 function

    float v_target_theta                            =avdat
        +coeff[  0]                
        +coeff[  1]*x11            
        +coeff[  2]    *x21        
        +coeff[  3]        *x31    
        +coeff[  4]            *x41
        +coeff[  5]*x11*x21        
        +coeff[  6]    *x21*x31    
        +coeff[  7]*x11        *x41
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[  8]*x11    *x31    
        +coeff[  9]        *x32    
        +coeff[ 10]        *x31*x41
        +coeff[ 11]        *x33    
        +coeff[ 12]*x11        *x42
        +coeff[ 13]*x12            
        +coeff[ 14]*x13            
        +coeff[ 15]*x11*x22        
        +coeff[ 16]*x12    *x31    
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x22*x31    
        +coeff[ 18]    *x21*x32    
        +coeff[ 19]    *x21*x31*x41
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x21+coeff[  7]*x41+coeff[  8]*x31+coeff[ 12]*x42+2*coeff[ 13]*x11+3*coeff[ 14]*x12+coeff[ 15]*x22+2*coeff[ 16]*x11*x31,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x11+coeff[  6]*x31+2*coeff[ 15]*x21*x11+2*coeff[ 17]*x21*x31+coeff[ 18]*x32+coeff[ 19]*x31*x41,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  6]*x21+coeff[  8]*x11+2*coeff[  9]*x31+coeff[ 10]*x41+3*coeff[ 11]*x32+coeff[ 16]*x12+coeff[ 17]*x22+2*coeff[ 18]*x31*x21+coeff[ 19]*x21*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+coeff[ 10]*x31+2*coeff[ 12]*x41*x11+coeff[ 19]*x21*x31,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.9335596E-03;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18005E+02,-0.18008E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18003E+02, 0.18001E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.71016507E-05,-0.53558350E-01, 0.22094907E-01, 0.63022785E-01,
        -0.34167029E-01, 0.18387387E-03,-0.15927899E-03, 0.15010464E-03,
         0.21519936E-04,-0.28960487E-04,-0.15163333E-04,-0.85638758E-05,
         0.17771067E-05,-0.64242579E-06,-0.30360932E-05,-0.53047861E-04,
        -0.13856114E-03, 0.24500999E-03,-0.27505876E-03, 0.74903905E-05,
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
        +coeff[  5]*x11    *x31    
        +coeff[  6]        *x32    
        +coeff[  7]*x11        *x41
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[  8]    *x21    *x41
        +coeff[  9]            *x42
        +coeff[ 10]    *x21*x32    
        +coeff[ 11]            *x43
        +coeff[ 12]*x13*x21        
        +coeff[ 13]*x11*x23        
        +coeff[ 14]*x13        *x41
        +coeff[ 15]*x12            
        +coeff[ 16]*x11*x21        
    ;
    v_target_phi                              =v_target_phi                              
        +coeff[ 17]    *x21*x31    
        +coeff[ 18]        *x31*x41
        +coeff[ 19]*x12*x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  5]*x31+coeff[  7]*x41+3*coeff[ 12]*x12*x21+coeff[ 13]*x23+3*coeff[ 14]*x12*x41+2*coeff[ 15]*x11+coeff[ 16]*x21+2*coeff[ 19]*x11*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  8]*x41+coeff[ 10]*x32+coeff[ 12]*x13+3*coeff[ 13]*x22*x11+coeff[ 16]*x11+coeff[ 17]*x31+coeff[ 19]*x12,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  5]*x11+2*coeff[  6]*x31+2*coeff[ 10]*x31*x21+coeff[ 17]*x21+coeff[ 18]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x11+coeff[  8]*x21+2*coeff[  9]*x41+3*coeff[ 11]*x42+coeff[ 14]*x13+coeff[ 18]*x31,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
