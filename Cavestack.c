#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

class avestack{
public:
    avestack(int sizeinput);
    ~avestack();
    bool isfull();
    bool isempty();
    int getsize();
    void empty();
    virtual void push(float x,float y);
    virtual float getave(int xy);
    virtual float getrms(int xy);
    float isqrt(float x);
protected:
    int start_pointer;
    int end_pointer;
    float *astack[2];
    float tot[2];
    int size;
    int i,pointer,buffersize;
    float totrms,ave;
};

///////////////////////extern///////////////////////////////////////////////////////////////////
extern "C"{
    avestack* avestack_new(int sizeinput){return new avestack(sizeinput);}
    void avestack_del(avestack* astack){delete astack;}
    bool avestack_isfull(avestack* astack){return astack->isfull();}
    bool avestack_isempty(avestack* astack){return astack->isempty();}
    int avestack_getsize(avestack* astack){return astack->getsize();}
    void avestack_empty(avestack* astack){astack->empty();}
    void avestack_push(avestack* astack,float x,float y){astack->push(x,y);}
    float avestack_ave(avestack* astack,int xy){return astack->getave(xy);}
    float avestack_rms(avestack* astack,int xy){return astack->getrms(xy);}
}

/////////////////////avestack///////////////////////////////////////////////////////////////////
avestack::avestack(int sizeinput){
    size=sizeinput;
    start_pointer=0;
    end_pointer=0;
    buffersize=0;
    for(i=0;i<2;i++){
	tot[i]=0;
	astack[i]=new float[size+1];
    }
}

avestack::~avestack(){
    for(i=0;i<2;i++){
	delete[] astack[i];
    }
}

int avestack::getsize(){
    return buffersize;
}

bool avestack::isempty(){
    return buffersize==0;
}

bool avestack::isfull(){
    return buffersize>=size;
}

void avestack::empty(){
	for(i=0;i<2;i++) tot[i]=0;
	start_pointer=0;
	end_pointer=0;
	buffersize=0;
}

void avestack::push(float x,float y=0){
    if(buffersize>=size){
	for(i=0;i<2;i++) tot[i]-=astack[i][end_pointer];
	end_pointer=end_pointer<size?end_pointer+1:0;
    }
    else buffersize++;
    astack[0][start_pointer]=x;
    astack[1][start_pointer]=y;
    tot[0]+=x;
    tot[1]+=y;
    start_pointer=start_pointer<size?start_pointer+1:0;
}

float avestack::getave(int xy=0){
    return xy>1?isqrt(tot[0]*tot[0]+tot[1]*tot[1])/buffersize:tot[xy]/buffersize;
}

float avestack::getrms(int xy=0){
    float purev;
    ave=getave(xy);
    totrms=0;
    for(i=0;i<buffersize;i++){
	pointer=(end_pointer+i>size)?end_pointer+i-size-1:end_pointer+i;
	purev=xy>1?isqrt(astack[0][pointer]*astack[0][pointer]+astack[1][pointer]*astack[1][pointer])-ave:astack[xy][pointer]-ave;
	totrms+=purev*purev;
    }
    return isqrt(totrms/buffersize);
}

float avestack::isqrt(float x){
    float xhalf = 0.5f*x;
    int i = *(int*)&x;
    i = 0x5f375a86- (i>>1);
    x = *(float*)&i;
    x = x*(1.5f-xhalf*x*x);
    return 1/x;
}