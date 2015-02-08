#include <math.h>
float target_x                                (float *x,float *dx,int m){
    int ncoeff= 19;
    float avdat= -0.4776969E+01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 20]={
         0.32949860E-02,-0.37747269E+02, 0.56203742E+01, 0.63014431E+02,
        -0.92170153E+01,-0.59460890E-01, 0.84875040E-01, 0.11269446E+00,
         0.14468814E+00,-0.16456492E+00,-0.10924708E+00,-0.45048472E-01,
        -0.21761464E-01,-0.47755767E-01, 0.61873309E-01,-0.37752711E-02,
         0.26132148E-02, 0.88863885E-02,-0.11188229E-01,
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
    int ncoeff= 12;
    float avdat=  0.2305200E+02;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 13]={
        -0.91789924E-02,-0.55652819E+01,-0.37786579E+02, 0.91338167E+01,
         0.63093155E+02, 0.17337595E+00, 0.54116964E-01, 0.19807277E-01,
        -0.12827802E-01,-0.19429608E+00, 0.25058651E-01,-0.19494768E-01,
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
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_y                                =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  6]*x11+coeff[  9]*x31+coeff[ 10]*x41+coeff[ 11]*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  7]*x41+2*coeff[  8]*x21+coeff[ 11]*x11,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  5]*x31+coeff[  9]*x11,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  7]*x21+coeff[ 10]*x11,2)
        );
    dx[0]=dv_target_y                                ;
    return v_target_y                                ;
}
#include <math.h>
float target_theta                            (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat=  0.4081295E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.16888196E-04,-0.31955551E-01,-0.49978644E-01, 0.49285531E-01,
         0.57266228E-01,-0.82093269E-04,-0.37598956E-03, 0.10438695E-03,
         0.17833844E-03, 0.11237516E-04, 0.43667180E-04,-0.37514030E-04,
        -0.10914430E-03,-0.73155893E-05,-0.47155727E-05, 0.10177912E-02,
        -0.69231266E-03, 0.56362042E-04,-0.74226868E-04,-0.17008388E-04,
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
    float x24 = x23*x2;
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
        +coeff[  8]    *x22    *x44
        +coeff[  9]    *x21*x31    
        +coeff[ 10]        *x33*x41
        +coeff[ 11]*x13*x22        
        +coeff[ 12]    *x24    *x44
        +coeff[ 13]*x11*x21        
        +coeff[ 14]*x11    *x31    
        +coeff[ 15]    *x21    *x41
        +coeff[ 16]            *x42
    ;
    v_target_theta                            =v_target_theta                            
        +coeff[ 17]    *x22    *x41
        +coeff[ 18]    *x21    *x42
        +coeff[ 19]*x13*x21        
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_theta                            =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+coeff[  7]*x32*x42+3*coeff[ 11]*x12*x22+coeff[ 13]*x21+coeff[ 14]*x31+3*coeff[ 19]*x12*x21,2)
        +dx2*dx2*pow(0+coeff[  2]+2*coeff[  6]*x21+2*coeff[  8]*x21*x44+coeff[  9]*x31+2*coeff[ 11]*x21*x13+4*coeff[ 12]*x23*x44+coeff[ 13]*x11+coeff[ 15]*x41+2*coeff[ 17]*x21*x41+coeff[ 18]*x42+coeff[ 19]*x13,2)
        +dx3*dx3*pow(0+coeff[  3]+2*coeff[  7]*x31*x11*x42+coeff[  9]*x21+3*coeff[ 10]*x32*x41+coeff[ 14]*x11,2)
        +dx4*dx4*pow(0+coeff[  4]+4*coeff[  5]*x43+2*coeff[  7]*x41*x11*x32+4*coeff[  8]*x43*x22+coeff[ 10]*x33+4*coeff[ 12]*x43*x24+coeff[ 15]*x21+2*coeff[ 16]*x41+coeff[ 17]*x22+2*coeff[ 18]*x41*x21,2)
        );
    dx[0]=dv_target_theta                            ;
    return v_target_theta                            ;
}
#include <math.h>
float target_phi                              (float *x,float *dx,int m){
    int ncoeff= 20;
    float avdat= -0.2402258E-01;
    float xmin[10]={
        -0.15000E+02,-0.15000E+02,-0.18003E+02,-0.18010E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10]={
         0.15000E+02, 0.15000E+02, 0.18010E+02, 0.18002E+02, 0.00000E+00,
         0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10]={0};
    float coeff[ 21]={
         0.56029367E-04,-0.50102632E-01, 0.32540374E-01, 0.57345070E-01,
        -0.50256159E-01, 0.13251291E-02,-0.10446544E-02,-0.19545112E-03,
         0.66784420E-03, 0.16239920E-03,-0.42537740E-03,-0.57248981E-03,
        -0.16451985E-03,-0.91015499E-04, 0.69025345E-03,-0.65593584E-03,
        -0.59342521E-04, 0.39044033E-04,-0.11718085E-03,-0.14780680E-04,
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
        +coeff[ 19]*x11        *x42
        ;
//                 error


    float dxall=0+dx1+dx2+dx3+dx4;
    float dv_target_phi                              =dxall==0?0:sqrt(
        +dx1*dx1*pow(0+coeff[  1]+2*coeff[  7]*x11+coeff[  8]*x31+coeff[  9]*x21+coeff[ 12]*x41+2*coeff[ 17]*x11*x21+coeff[ 19]*x42,2)
        +dx2*dx2*pow(0+coeff[  2]+coeff[  5]*x41+coeff[  9]*x11+2*coeff[ 10]*x21+coeff[ 14]*x42+coeff[ 17]*x12+3*coeff[ 18]*x22,2)
        +dx3*dx3*pow(0+coeff[  3]+coeff[  8]*x11+2*coeff[ 11]*x31+2*coeff[ 13]*x31*x41+coeff[ 16]*x41,2)
        +dx4*dx4*pow(0+coeff[  4]+coeff[  5]*x21+2*coeff[  6]*x41+coeff[ 12]*x11+coeff[ 13]*x32+2*coeff[ 14]*x41*x21+3*coeff[ 15]*x42+coeff[ 16]*x31+2*coeff[ 19]*x41*x11,2)
        );
    dx[0]=dv_target_phi                              ;
    return v_target_phi                              ;
}
