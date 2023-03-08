#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

using namespace std;

double g=1000;                       //garivtational accleration
double H=3;                          //hight
double lf=1;                         //half shoulder length
double lh=1;                         //half hip length
double h=1;                          //half width
double l=5;                          //half body length
double L1=1*l, L2=0*l;               //front, hind part body length
double L=L1+L2;                      //total body length
double D=sqrt(L*L+(lh+lf)*(lh+lf));  //diagonal length
double freq=sqrt(g/H);               //inverted pendulum frequency
double Tswc;                         //typical swing time

// double m   = 30;     // mouse mass in g
double I   = L*L/3;  	// moment of inertia over mass
double tau = .25;

double DM[5]     = {0,1*L,1*L,1*L,1*L};  //max leg length

// double xcom(double x0,double y0,double vx,double vy,double theta,double delta)  //distance of extrapolated center of mass from the diagonal, delta is the body bending angle
// {
// 	double xh=x0-L2*cos(theta), yh=y0-L2*sin(theta);	//pelvis
// 	double x2=xh+h*sin(theta), y2=yh-h*cos(theta);
	
// 	theta+=delta;	//delta is the body bending angle

// 	double xf=x0+L1*cos(theta), yf=y0+L1*sin(theta);	//shoulder
// 	double x1=xf-h*sin(theta), y1=yf+h*cos(theta);

// 	double a=y1-y2,b=x2-x1,ab=sqrt(a*a+b*b);
// 	a/=ab; b/=ab;

// 	if(h<0) { a=-a; b=-b; }
	
// 	double vz=vx*a+vy*b;
	
// 	double cor=0;   //-L1/2*sin(delta/2);
// 	double xc=x0+cor*cos(theta+(M_PI-delta)/2), yc=y0+cor*sin(theta+(M_PI-delta)/2);
	
// 	return (xc-x1)*a+(yc-y1)*b+vz/freq;
// }

double rnd1() {return ((double) rand() / (RAND_MAX)) * 2 - 1;}
double rnd2() {srand (time(NULL)); return double(rand())/RAND_MAX;}

void update_F(double xfh[], double yfh[], double xc, double yc, double F[][5]) {
    double u1=xfh[1]-xc,u2=xfh[2]-xc,u3=xfh[3]-xc,u4=xfh[4]-xc;
    double v1=yfh[1]-yc,v2=yfh[2]-yc,v3=yfh[3]-yc,v4=yfh[4]-yc;

    double xi=u1+u2+u3+u4, eta=v1+v2+v3+v4;
    double xi2=u1*u1+u2*u2+u3*u3+u4*u4, eta2=v1*v1+v2*v2+v3*v3+v4*v4;
    double xieta=u1*v1+u2*v2+u3*v3+u4*v4;
    double W=4*(xi2*eta2-xieta*xieta)-xi*(xi*eta2-xieta*eta)+eta*(xi*xieta-xi2*eta);
    double lam1=-(xi2*eta2-xieta*xieta)/W;
    double lam2=(xi*eta2-xieta*eta)/W;
    double lam3=-(xi*xieta-xi2*eta)/W;

    F[1][0]=-lam1-lam2*u1-lam3*v1;
    F[2][0]=-lam1-lam2*u2-lam3*v2;
    F[3][0]=-lam1-lam2*u3-lam3*v3;
    F[4][0]=-lam1-lam2*u4-lam3*v4;

    double W4=(u2*v3-u3*v2-u1*v3+u3*v1+u1*v2-u2*v1);
    F[1][4]=(u2*v3-u3*v2)/W4;
    F[2][4]=(-u1*v3+u3*v1)/W4;
    F[3][4]=(u1*v2-u2*v1)/W4;

    double W3=(u2*v4-u4*v2-u1*v4+u4*v1+u1*v2-u2*v1);
    F[1][3]=(u2*v4-u4*v2)/W3;
    F[2][3]=(-u1*v4+u4*v1)/W3;
    F[4][3]=(u1*v2-u2*v1)/W3;

    double W2=(u3*v4-u4*v3-u1*v4+u4*v1+u1*v3-u3*v1);
    F[1][2]=(u3*v4-u4*v3)/W2;
    F[3][2]=(-u1*v4+u4*v1)/W2;
    F[4][2]=(u1*v3-u3*v1)/W2;

    double W1=(u3*v4-u4*v3-u2*v4+u4*v2+u2*v3-u3*v2);
    F[2][1]=(u3*v4-u4*v3)/W1;
    F[3][1]=(-u2*v4+u4*v2)/W1;
    F[4][1]=(u2*v3-u3*v2)/W1;
}

int main()
{
	// int Fl[5]={0,1,1,1,0};
	// int count=0;
	// for(int i=1;i<=4;i++) {count+=Fl[i];}

	// cout<<count<<xcom(0,0,0,0,0,0)<<endl;

	// double x[5]={0,2,4,6,8}, y[5]={0,1,3,5,7};
	// for(int i=1;i<=4;i++){cout<<"set object "<< i <<" circle front at "<<x[i]<<','<<y[i]<<" size "<<.2<<" fc \"brown\""<<endl;}

	// int time = 9;
	// if (time < 10) {
	// cout << "Good morning.";
	// } else if (time < 20) {
	// cout << "Good day.";
	// } else {
	// cout << "Good evening.";
	// }
	// srand (time(NULL)); 
	// int i=0;
	// while(i<=5){
	// 	cout<<rnd1()<<'\t'<<rnd1()<<endl;
	// 	i++;
	// }
	double vv, vexp;//speed in magnitude, expected speed
	double omega;	//body angular velocity
	double theta;	//body angle

	double xc, yc;	//position of COM
    double vx, vy;	//velocity of COM
    double xf, yf, xh, yh;	//middile of shoulder and hip positions
    double xfh[5]={},yfh[5]={},xfh0[5]={},yfh0[5]={},xfh1[5]={},yfh1[5]={};
    double &xfl=xfh[1], &xfr=xfh[2], &xhl=xfh[3], &xhr=xfh[4];
    double &yfl=yfh[1], &yfr=yfh[2], &yhl=yfh[3], &yhr=yfh[4];			//limb positions
    double &xfl0=xfh0[1], &xfr0=xfh0[2], &xhl0=xfh0[3], &xhr0=xfh0[4];
    double &yfl0=yfh0[1], &yfr0=yfh0[2], &yhl0=yfh0[3], &yhr0=yfh0[4];	//joint positions
    double &xfl1=xfh1[1], &xfr1=xfh1[2], &xhl1=xfh1[3], &xhr1=xfh1[4];
    double &yfl1=yfh1[1], &yfr1=yfh1[2], &yhl1=yfh1[3], &yhr1=yfh1[4];	//corrected landing positions
	xc=0; yc=0;
	theta=M_PI/4;

	omega=0;

	xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
	xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);
	xfl0=xf-lf*sin(theta), xfr0=xf+lf*sin(theta);
	xhl0=xh-lh*sin(theta), xhr0=xh+lh*sin(theta);
	yfl0=yf+lf*cos(theta), yfr0=yf-lf*cos(theta);
	yhl0=yh+lh*cos(theta), yhr0=yh-lh*cos(theta);

	
	return 0;
}