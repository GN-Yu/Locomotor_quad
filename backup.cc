#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

double g=1000; // garivtational accleration in cm/s/s
double H=1.5; // mouse hight in cm
double L=5; // mouse half-length in cm
double h=1.5; // mouse half-width in cm
double D=sqrt(L*L+4*h*h); // diagonal length
double freq=sqrt(g/H); // frequency

// double m=25; // mouse mass in g
double I=L*L/3; // moment of inertia over mass

double F0=10;	// locomotor speed
double tau=0.2;	// lambda = 1/tau

double L1=1*L, L2=.2*L; //hind legs are close to the COM

double excom(double x0,double y0,double vx,double vy,double h,double theta,double delta)  //distance of extrapolated center of mass from the diagonal, delta is the body bending angle
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

int main()
{
	double xc=0, yc=0, theta=0;
	double vx=3*F0*cos(theta), vy=3*F0*sin(theta), omega=0; //why 3? because there are 3 legs on the ground (means the speed is actually 3*F0)
	
	double xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
	double xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);

	double xfl0=xf-h*sin(theta), xfr0=xf+h*sin(theta);
	double xhl0=xh-h*sin(theta), xhr0=xh+h*sin(theta);
	
	double yfl0=yf+h*cos(theta), yfr0=yf-h*cos(theta);
	double yhl0=yh+h*cos(theta), yhr0=yh-h*cos(theta);
	
	double xfl=xfl0, yfl=yfl0;
	double xfr=xfr0, yfr=yfr0;
	double xhl=xhl0, yhl=yhl0;
	double xhr=xhr0, yhr=yhr0;  //initial positions

	int Fl[5]={0,0,0,0,0};	//means which foot is lifted
	int fl1,fl2,fl3,fl4;	//means which foot could be lifted
	
	double T=2, dt=.001;
	
	for(double t=0;t<=T;t+=dt)
	{
		if(Fl[1]+Fl[2]+Fl[3]+Fl[4]>2) break;

		double u1=xfl-xc,u2=xfr-xc,u3=xhl-xc,u4=xhr-xc;
		double v1=yfl-yc,v2=yfr-yc,v3=yhl-yc,v4=yhr-yc;

		double W4=(u2*v3-u3*v2-u1*v3+u3*v1+u1*v2-u2*v1);
		double F14=(u2*v3-u3*v2)/W4;
		double F24=(-u1*v3+u3*v1)/W4;
		double F34=(u1*v2-u2*v1)/W4;
		double F44=0;
			
		double W3=(u2*v4-u4*v2-u1*v4+u4*v1+u1*v2-u2*v1);
		double F13=(u2*v4-u4*v2)/W3;
		double F23=(-u1*v4+u4*v1)/W3;
		double F43=(u1*v2-u2*v1)/W3;
		double F33=0;
			
		double W2=(u3*v4-u4*v3-u1*v4+u4*v1+u1*v3-u3*v1);
		double F12=(u3*v4-u4*v3)/W2;
		double F32=(-u1*v4+u4*v1)/W2;
		double F42=(u1*v3-u3*v1)/W2;
		double F22=0;
			
		double W1=(u3*v4-u4*v3-u2*v4+u4*v2+u2*v3-u3*v2);
		double F21=(u3*v4-u4*v3)/W1;
		double F31=(-u2*v4+u4*v2)/W1;
		double F41=(u2*v3-u3*v2)/W1;
		double F11=0;	
        
        // double W14;
        // double F114=0;
        // double F414=0;

        // double W23;
        // double F223=0;
        // double F323=0;
        // solutions of equilibrium when one or two of the limbs are lifted (cross products of each force and its arm)
		
		static int isw=0;
		
		fl1=F21>=0 && F31>=0 && F41>=0; //change to (original plus "or 'diagonal on the groud and excom>=0")
		fl2=F12>=0 && F32>=0 && F42>=0;
		fl3=F13>=0 && F23>=0 && F43>=0;
		fl4=F14>=0 && F24>=0 && F34>=0;
		
		double d[5]={0,
			hypot(xfl0-xfl,yfl0-yfl),hypot(xfr0-xfr,yfr0-yfr),hypot(xhl0-xhl,yhl0-yhl),hypot(xhr0-xhr,yhr0-yhr)};
		double dmax=0,dsec=0;
		int imax=0, isec=0;
		for(int i=1;i<5;i++) if(d[i]>=dmax) { dsec=dmax; dmax=d[i]; isec=imax; imax=i; }
		
		static double tpre=0,P=0.07;
		
		if(t-tpre<P && isw==1 && fl1) { xfl=xfl0; yfl=yfl0; }
		else if(t-tpre<P && isw==2 && fl2) { xfr=xfr0; yfr=yfr0; }
		else if(t-tpre<P && isw==3 && fl3) { xhl=xhl0; yhl=yhl0; }
		else if(t-tpre<P && isw==4 && fl4) { xhr=xhr0; yhr=yhr0; }
		else if(d[1]==dmax && fl1) { isw=1; xfl=xfl0; yfl=yfl0; tpre=t; }
		else if(d[2]==dmax && fl2) { isw=2; xfr=xfr0; yfr=yfr0; tpre=t; }
		else if(d[3]==dmax && fl3) { isw=3; xhl=xhl0; yhl=yhl0; tpre=t; }
		else if(d[4]==dmax && fl4) { isw=4; xhr=xhr0; yhr=yhr0; tpre=t; }
		else if(d[1]==dsec && fl1) { isw=1; xfl=xfl0; yfl=yfl0; tpre=t; }
		else if(d[2]==dsec && fl2) { isw=2; xfr=xfr0; yfr=yfr0; tpre=t; }
		else if(d[3]==dsec && fl3) { isw=3; xhl=xhl0; yhl=yhl0; tpre=t; }
		else if(d[4]==dsec && fl4) { isw=4; xhr=xhr0; yhr=yhr0; tpre=t; }


		double Fx[5]={},Fy[5]={};   //forces, F[0] not used
		double delta=0;
		double FL=(1+delta)*F0, FR=(1-delta)*F0;    //biased forces for different sides
		
		// if(isw!=1) if(d[1]==0) { Fx[1]=FL*cos(theta); Fy[1]=FL*sin(theta); } 
		// 			else { Fx[1]=FL*(xfl0-xfl)/d[1]; Fy[1]=FL*(yfl0-yfl)/d[1]; }
		// if(isw!=2) if(d[2]==0) { Fx[2]=FR*cos(theta); Fy[2]=FR*sin(theta); } 
		// 			else { Fx[2]=FR*(xfr0-xfr)/d[2]; Fy[2]=FR*(yfr0-yfr)/d[2]; }
		// if(isw!=3) if(d[3]==0) { Fx[3]=FL*cos(theta); Fy[3]=FL*sin(theta); } 
		// 			else { Fx[3]=FL*(xhl0-xhl)/d[3]; Fy[3]=FL*(yhl0-yhl)/d[3]; }
		// if(isw!=4) if(d[4]==0) { Fx[4]=FR*cos(theta); Fy[4]=FR*sin(theta); } 
		// 			else { Fx[4]=FR*(xhr0-xhr)/d[4]; Fy[4]=FR*(yhr0-yhr)/d[4]; }

		if(isw!=1) { Fx[1]=FL*cos(theta); Fy[1]=FL*sin(theta); }
		if(isw!=2) { Fx[2]=FR*cos(theta); Fy[2]=FR*sin(theta); }
		if(isw!=3) { Fx[3]=FL*cos(theta); Fy[3]=FL*sin(theta); }
		if(isw!=4) { Fx[4]=FR*cos(theta); Fy[4]=FR*sin(theta); }

		xc+=vx*dt;
		yc+=vy*dt;
		theta+=omega*dt;

		vx+=(Fx[1]+Fx[2]+Fx[3]+Fx[4]-vx)*dt/tau;
		vy+=(Fy[1]+Fy[2]+Fy[3]+Fy[4]-vy)*dt/tau;
		
		double M[5]={0,
			(xfl-xc)*Fy[1]-(yfl-yc)*Fx[1],
			(xfr-xc)*Fy[2]-(yfr-yc)*Fx[2],
			(xhl-xc)*Fy[3]-(yhl-yc)*Fx[3],
			(xhr-xc)*Fy[4]-(yhr-yc)*Fx[4]};
		
		omega+=((M[1]+M[2]+M[3]+M[4])/I-omega)*dt/tau;
		
		
		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);

		xfl0=xf-h*sin(theta), xfr0=xf+h*sin(theta);
		xhl0=xh-h*sin(theta), xhr0=xh+h*sin(theta);
	
		yfl0=yf+h*cos(theta), yfr0=yf-h*cos(theta);
		yhl0=yh+h*cos(theta), yhr0=yh-h*cos(theta);
		

		cout<<t<<'\t'<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<theta<<'\t';
		// switch(isw)
		// {
		// 	case 1: cout<<F11<<'\t'<<F21<<'\t'<<F31<<'\t'<<F41<<endl; break;
		// 	case 2: cout<<F12<<'\t'<<F22<<'\t'<<F32<<'\t'<<F42<<endl; break;
		// 	case 3: cout<<F13<<'\t'<<F23<<'\t'<<F33<<'\t'<<F43<<endl; break;
		// 	case 4: cout<<F14<<'\t'<<F24<<'\t'<<F34<<'\t'<<F44<<endl; break;
		// }
		
		cerr<<"set xrange[-2:80]"<<endl;
		cerr<<"set yrange[-10:10]"<<endl;
		cerr<<"unset object"<<endl;
		cerr<<"set object 1 circle front at "<<xfl<<','<<yfl<<" size "<<.2<<" fc \"brown\""<<endl;
		cerr<<"set object 2 circle front at "<<xfr<<','<<yfr<<" size "<<.2<<" fc \"blue\""<<endl;
		cerr<<"set object 3 circle front at "<<xhl<<','<<yhl<<" size "<<.2<<" fc \"red\""<<endl;
		cerr<<"set object 4 circle front at "<<xhr<<','<<yhr<<" size "<<.2<<" fc \"green\""<<endl;
		cerr<<"set object 5 circle front at "<<xc<<','<<yc<<" size "<<.2<<" fc \"black\""<<endl;

		cerr<<"unset arrow"<<endl;
		if(isw!=1 && isw!=4) cerr<<"set arrow 1 from "<<xfl<<','<<yfl<<" to "<<xhr<<','<<yhr<<" nohead lt 0"<<endl;
		if(isw!=2 && isw!=3) cerr<<"set arrow 2 from "<<xfr<<','<<yfr<<" to "<<xhl<<','<<yhl<<" nohead lt 0"<<endl;
		if(isw!=2 && isw!=1) cerr<<"set arrow 3 from "<<xfr<<','<<yfr<<" to "<<xfl<<','<<yfl<<" nohead lt 0"<<endl;
		if(isw!=3 && isw!=4) cerr<<"set arrow 4 from "<<xhr<<','<<yhr<<" to "<<xhl<<','<<yhl<<" nohead lt 0"<<endl;
		cerr<<"set arrow 5 from "<<xf<<','<<yf<<" to "<<xfl<<','<<yfl<<" nohead"<<endl;
		cerr<<"set arrow 6 from "<<xf<<','<<yf<<" to "<<xfr<<','<<yfr<<" nohead"<<endl;
		cerr<<"set arrow 7 from "<<xh<<','<<yh<<" to "<<xhl<<','<<yhl<<" nohead"<<endl;
		cerr<<"set arrow 8 from "<<xh<<','<<yh<<" to "<<xhr<<','<<yhr<<" nohead"<<endl;
		cerr<<"set arrow 9 from "<<xh<<','<<yh<<" to "<<xf<<','<<yf<<" lw 2 nohead"<<endl;
		
		cerr<<"plot \"dat\" u 2:3 w d"<<endl;
	}
	return 0;
}