//this file contains several functions to boost the code

#include <math.h>

//for slow raster function
extern "C" float slrasterfunc(float x,float *par)
{
    //par:am phase,phase,am period,period,max am,center am,clock rate
    float pi=3.141592653589793;
    float am,dp;
    float time=x/par[6];
    float tmp1=time+par[0];
    float tmp2=par[2]*4;
    dp=tmp1-tmp2*(int)(tmp1/tmp2);

    if(dp<=par[2]) am=sqrt(dp);
    else if(dp>par[2] && dp<=par[2]*2) am=sqrt(2*par[2]-dp);
    else if(dp>par[2]*2 && dp<=par[2]*3) am=-sqrt(-2*par[2]+dp);
    else am=-sqrt(4*par[2]-dp);
    return par[4]/sqrt(par[2])*am*sin(2*pi*time/par[3]+par[1])+par[5];
}
