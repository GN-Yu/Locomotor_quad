#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

using namespace std;

double g    = 1000;             // garivtational accleration in cm/s/s
double H    = 1.5;              // mouse hight in cm
double L    = 5;                // mouse half-length in cm
double h    = 1.5;              // mouse half-width in cm
double D    = sqrt(L*L+4*h*h);  // diagonal length
double freq = sqrt(g/H);        // frequency
double m    = 25;               // mouse mass in g
double I    = L*L/3;            // moment of inertia over mass
double F0   = 10;               // locomotor speed
double tau  = 0.2;              // lambda = 1/tau
double L1   = 1*L, L2 = .2*L;   //hind legs are close to the COM

double xcom(double x0,double y0,double vx,double vy,double theta,double delta)  //distance of extrapolated center of mass from the diagonal, delta is the body bending angle
{
	double xh=x0-L2*cos(theta), yh=y0-L2*sin(theta);	//pelvis
	double x2=xh+h*sin(theta), y2=yh-h*cos(theta);
	
	theta+=delta;	//delta is the body bending angle

	double xf=x0+L1*cos(theta), yf=y0+L1*sin(theta);	//shoulder
	double x1=xf-h*sin(theta), y1=yf+h*cos(theta);

	double a=y1-y2,b=x2-x1,ab=sqrt(a*a+b*b);
	a/=ab; b/=ab;

	if(h<0) { a=-a; b=-b; }
	
	double vz=vx*a+vy*b;
	
	double cor=0;   //-L1/2*sin(delta/2);
	double xc=x0+cor*cos(theta+(M_PI-delta)/2), yc=y0+cor*sin(theta+(M_PI-delta)/2);
	
	return (xc-x1)*a+(yc-y1)*b+vz/freq;
}

double rnd1() {return ((double) rand() / (RAND_MAX)) * 2 - 1;}
double rnd2() {srand (time(NULL)); return double(rand())/RAND_MAX;}

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
	srand (time(NULL)); 
	int i=0;
	while(i<=5){
		cout<<rnd1()<<'\t'<<rnd1()<<endl;
		i++;
	}
	
	
	return 0;
}