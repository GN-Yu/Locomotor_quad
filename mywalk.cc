#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

double g=1000; // garivtational accleration in cm/s/s

/* define biomechanics of a mouse */
double H=1.5; // mouse hight in cm
double L=5; // mouse half-length in cm
double h=1.5; // mouse half-width in cm
// double m=25; // mouse mass in g

double L1=1*L, L2=.2*L; //hind legs are closer to the COM
double fm=0.8*L;
double hm=1.1*L;	//maximum extension of the front and hind limbs
double P=0.5;	//swing period of a leg
double vs=30;	//typical swing speed

double D=sqrt(L*L+4*h*h); // diagonal length
double freq=sqrt(g/H); // frequency
double I=L*L/3; // moment of inertia over mass
double tll = 0.9;	//threshold of lifting a leg

/* distance of extrapolated center of mass from the supporting diagonal */
double dds(double x0,double y0,double s1x,double s1y,double s2x,double s2y)	//postitive: xcom has already crossed the diagonal; negative otherwise
{
	return ((x0-s2x) * (s2y-s1y) + (s1x-s2x) * (y0-s2y)) / sqrt((s1x-s2x)*(s1x-s2x) + (s1y-s2y)*(s1y-s2y));
}

int main()
{
	/* initial values */
	double delta=0;	//to be used as body bending angle
	double xc=0, yc=0, theta=0;
	double F0=10;	// locomotor speed
	double tau=0.1;	// lambda = 1/tau

	double vx=3*F0*cos(theta), vy=3*F0*sin(theta), omega=0; //why 3? because there are 3 legs on the ground (means the speed is actually 3*F0)

	double xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
	double xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);

	double x0[5],y0[5],x[5],y[5];
	x0[1]=xf-h*sin(theta), x0[2]=xf+h*sin(theta);
	x0[3]=xh-h*sin(theta), x0[4]=xh+h*sin(theta);

	y0[1]=yf+h*cos(theta), y0[2]=yf-h*cos(theta);
	y0[3]=yh+h*cos(theta), y0[4]=yh-h*cos(theta);

	x[1]=x0[1], y[1]=y0[1];
	x[2]=x0[2], y[2]=y0[2];
	x[3]=x0[3], y[3]=y0[3];
	x[4]=x0[4], y[4]=y0[4];

	double xcomx = xc+vx/freq;
	double xcomy = yc+vy/freq;

	/* lifting decisions */
	int Fl[5]={0,1,0,0,0};	//means which foot is lifted, true for lifted
	int fl[5]={0,1,0,0,0};	//means which foot could be lifted

	double T=2, dt=.0005;

	for(double t=0;t<=T;t+=dt)
	{
		int count=0;	//number of legs that are swinging
		for(int i=1;i<=4;i++) {count+=Fl[i];}
		//if(count>2) break;

		fl[1]=0, fl[2]=0, fl[3]=0, fl[4]=0;
		double u1=x[1]-xc,u2=x[2]-xc,u3=x[3]-xc,u4=x[4]-xc;
		double v1=y[1]-yc,v2=y[2]-yc,v3=y[3]-yc,v4=y[4]-yc;
		double F41pre;
		double F32pre;
		// if(count==0)
		// {
		// 	double xi=u1+u2+u3+u4, eta=v1+v2+v3+v4;
		// 	double xi2=u1*u1+u2*u2+u3*u3+u4*u4, eta2=v1*v1+v2*v2+v3*v3+v4*v4;
		// 	double xieta=u1*v1+u2*v2+u3*v3+u4*v4;
		// 	double WW=4*(xi2*eta2-xieta*xieta)-xi*(xi*eta2-xieta*eta)+eta*(xi*xieta-xi2*eta);
		// 	double lam1=-(xi2*eta2-xieta*xieta)/WW;
		// 	double lam2=(xi*eta2-xieta*eta)/WW;
		// 	double lam3=-(xi*xieta-xi2*eta)/WW;
		// 	double F10=-lam1-lam2*u1-lam3*v1;
		// 	double F20=-lam1-lam2*u2-lam3*v2;
		// 	double F30=-lam1-lam2*u3-lam3*v3;
		// 	double F40=-lam1-lam2*u4-lam3*v4;
		// 	if(F10<=0) fl[1]=1;
		// 	if(F20<=0) fl[2]=1;
		// 	if(F30<=0) fl[3]=1;
		// 	if(F40<=0) fl[4]=1;
		// }
		// else if(!Fl[1] && Fl[2] && !Fl[4] && dds(xcomx,xcomy,x[1],y[1],x[4],y[4])>=0) {fl[1]=0; fl[2]=1; fl[3]=1; fl[4]=0;}
		// else if(Fl[1] && !Fl[2] && !Fl[3] && dds(xcomx,xcomy,x[2],y[2],x[3],y[3])>=0) {fl[1]=1; fl[2]=0; fl[3]=0; fl[4]=1;}
		if(count<=1)
		{
			double W4=(u2*v3-u3*v2-u1*v3+u3*v1+u1*v2-u2*v1);
			double F14=(u2*v3-u3*v2)/W4;
			double F24=(-u1*v3+u3*v1)/W4;
			double F34=(u1*v2-u2*v1)/W4;
			// double F44=0;

			double W3=(u2*v4-u4*v2-u1*v4+u4*v1+u1*v2-u2*v1);
			double F13=(u2*v4-u4*v2)/W3;
			double F23=(-u1*v4+u4*v1)/W3;
			double F43=(u1*v2-u2*v1)/W3;
			// double F33=0;

			double W2=(u3*v4-u4*v3-u1*v4+u4*v1+u1*v3-u3*v1);
			double F12=(u3*v4-u4*v3)/W2;
			double F32=(-u1*v4+u4*v1)/W2;
			double F42=(u1*v3-u3*v1)/W2;
			// double F22=0;

			double W1=(u3*v4-u4*v3-u2*v4+u4*v2+u2*v3-u3*v2);
			double F21=(u3*v4-u4*v3)/W1;
			double F31=(-u2*v4+u4*v2)/W1;
			double F41=(u2*v3-u3*v2)/W1;
			// double F11=0;

			if(!Fl[2] && !Fl[3] && !Fl[4] && F21>=0 && F31>=0 && F41>=0) 
			{
				fl[1]=1;
				if(Fl[1] && F41<F41pre && F41<tll) fl[4] = 1;
			}
			if(!Fl[1] && !Fl[3] && !Fl[4] && F12>=0 && F32>=0 && F42>=0) 
			{
				fl[2]=1;
				if(Fl[2] && F32<F32pre && F32<tll) fl[3] = 1;
			}
			if(!Fl[1] && !Fl[2] && !Fl[4] && F13>=0 && F23>=0 && F43>=0) fl[3]=1;
			if(!Fl[1] && !Fl[2] && !Fl[3] && F14>=0 && F24>=0 && F34>=0) fl[4]=1;

			F41pre = F41;
			F32pre = F32;
		}

		double d[5]={0,hypot(x0[1]-x[1],y0[1]-y[1]),hypot(x0[2]-x[2],y0[2]-y[2]),hypot(x0[3]-x[3],y0[3]-y[3]),hypot(x0[4]-x[4],y0[4]-y[4])};
		for(int i=1;i<=2;i++) if(d[i]>=fm) fl[i]=1;
		for(int i=3;i<=4;i++) if(d[i]>=hm) fl[i]=1;	//a leg should be lifted when fully extended

		static double tpre[5]={0,0,0,0,0};	//record the time when each leg is lifted

		for(int i=1;i<=4;i++)	//leg lifting
		{
			if(!Fl[i] && fl[i]) {Fl[i]=1; tpre[i]=t; count++; break;} //lift a leg that can be lifted
		}


		/* update variables (move a little) */
		double Fx[5]={},Fy[5]={};   //forces
		double beta=0;
		double FL=(1+beta)*F0, FR=(1-beta)*F0;    //biased forces for different sides

		for(int i=1;i<=4;i++)
		{
			if(Fl[i])
			{
				if(d[i]!=.0) {x[i]+=dt*vs*(x0[i]-x[i]); y[i]+=dt*vs*(y0[i]-y[i]);}	//movement of swinging legs
				Fx[i]=0;
				Fy[i]=0;
			}
			else { Fx[i]=FL*cos(theta); Fy[i]=FL*sin(theta); }	//forces of standing legs pointing along the body orientation
			//else if(d[i]) {Fx[i]=FR*(x0[i]-x[i])/d[i]; Fy[i]=FR*(y0[i]-y[i])/d[i]; }	//forces of standing legs pointing to the shoulder/pelvis joints
		}

		xc+=vx*dt;
		yc+=vy*dt;
		theta+=omega*dt;

		xcomx = xc+vx/freq;
		xcomy = yc+vy/freq;

		double hgx=0, hgy=0;	//"horizontal projection" of gravity pointing to the direction perpendicular to the diagonal
		double s1x, s1y, s2x, s2y, dd;
		int temps;
		if(count<=1 || count==4) {hgx=0; hgy=0;}
		else if(count==2)
		{
			s1x=0, s1y=0 ,s2x=0 ,s2y=0;
			temps=0;
			for(int i=1;i<=4;i++)	//recognize the supporting legs
			{
				if(!Fl[i] && temps<1) {s1x=x[i]; s1y=y[i]; temps++; continue;}
				if(!Fl[i] && temps<2) {s2x=x[i]; s2y=y[i]; temps++;}
			}
			hgx=(s2y-s1y) * (g/H) * dds(xc,yc,s1x,s1y,s2x,s2y) / sqrt((s1y-s2y)*(s1y-s2y)+(s1x-s2x)*(s1x-s2x));
			hgy=(s1x-s2x) * (g/H) * dds(xc,yc,s1x,s1y,s2x,s2y) / sqrt((s1y-s2y)*(s1y-s2y)+(s1x-s2x)*(s1x-s2x));	//wtong with direction
		}
		// else if(count==3)


		vx+=(Fx[1]+Fx[2]+Fx[3]+Fx[4]-vx)*dt/tau;
		vy+=(Fy[1]+Fy[2]+Fy[3]+Fy[4]-vy)*dt/tau;
		// vx+=(Fx[1]+Fx[2]+Fx[3]+Fx[4]+tau*hgx-vx)*dt/tau;
		// vy+=(Fy[1]+Fy[2]+Fy[3]+Fy[4]+tau*hgy-vy)*dt/tau;

		double M[5]={0,
			(x[1]-xc)*Fy[1]-(y[1]-yc)*Fx[1],
			(x[2]-xc)*Fy[2]-(y[2]-yc)*Fx[2],
			(x[3]-xc)*Fy[3]-(y[3]-yc)*Fx[3],
			(x[4]-xc)*Fy[4]-(y[4]-yc)*Fx[4]};

		omega+=((M[1]+M[2]+M[3]+M[4])/I-omega)*dt/tau;

		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);

		x0[1]=xf-h*sin(theta), x0[2]=xf+h*sin(theta);
		x0[3]=xh-h*sin(theta), x0[4]=xh+h*sin(theta);

		y0[1]=yf+h*cos(theta), y0[2]=yf-h*cos(theta);
		y0[3]=yh+h*cos(theta), y0[4]=yh-h*cos(theta);


		/* data output */
		cout<<t<<'\t'<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<theta<<'\t';
		// switch(isw)
		// {
		// 	case 1: cout<<F11<<'\t'<<F21<<'\t'<<F31<<'\t'<<F41<<endl; break;
		// 	case 2: cout<<F12<<'\t'<<F22<<'\t'<<F32<<'\t'<<F42<<endl; break;
		// 	case 3: cout<<F13<<'\t'<<F23<<'\t'<<F33<<'\t'<<F43<<endl; break;
		// 	case 4: cout<<F14<<'\t'<<F24<<'\t'<<F34<<'\t'<<F44<<endl; break;
		// }
		cout<<" after lifting: "<<Fl[1]<<'\t'<<Fl[2]<<'\t'<<Fl[3]<<'\t'<<Fl[4]<<'\t';

		cerr<<"set xrange[-10:40]"<<endl;
		cerr<<"set yrange[-5:5]"<<endl;
		cerr<<"unset object"<<endl;
		cerr<<"set object 1 circle front at "<<x[1]<<','<<y[1]<<" size "<<.2<<" fc \"brown\""<<endl;
		cerr<<"set object 2 circle front at "<<x[2]<<','<<y[2]<<" size "<<.2<<" fc \"blue\""<<endl;
		cerr<<"set object 3 circle front at "<<x[3]<<','<<y[3]<<" size "<<.2<<" fc \"red\""<<endl;
		cerr<<"set object 4 circle front at "<<x[4]<<','<<y[4]<<" size "<<.2<<" fc \"green\""<<endl;
		cerr<<"set object 5 circle front at "<<xc<<','<<yc<<" size "<<.2<<" fc \"black\""<<endl;

		cerr<<"unset arrow"<<endl;
		if(!Fl[1] && !Fl[4]) cerr<<"set arrow 1 from "<<x[1]<<','<<y[1]<<" to "<<x[4]<<','<<y[4]<<" nohead lt 0"<<endl;
		if(!Fl[2] && !Fl[3]) cerr<<"set arrow 2 from "<<x[2]<<','<<y[2]<<" to "<<x[3]<<','<<y[3]<<" nohead lt 0"<<endl;
		if(!Fl[2] && !Fl[1]) cerr<<"set arrow 3 from "<<x[2]<<','<<y[2]<<" to "<<x[1]<<','<<y[1]<<" nohead lt 0"<<endl;
		if(!Fl[3] && !Fl[4]) cerr<<"set arrow 4 from "<<x[4]<<','<<y[4]<<" to "<<x[3]<<','<<y[3]<<" nohead lt 0"<<endl;
		cerr<<"set arrow 5 from "<<xf<<','<<yf<<" to "<<x[1]<<','<<y[1]<<" nohead"<<endl;
		cerr<<"set arrow 6 from "<<xf<<','<<yf<<" to "<<x[2]<<','<<y[2]<<" nohead"<<endl;
		cerr<<"set arrow 7 from "<<xh<<','<<yh<<" to "<<x[3]<<','<<y[3]<<" nohead"<<endl;
		cerr<<"set arrow 8 from "<<xh<<','<<yh<<" to "<<x[4]<<','<<y[4]<<" nohead"<<endl;
		cerr<<"set arrow 9 from "<<xh<<','<<yh<<" to "<<xf<<','<<yf<<" lw 2 nohead"<<endl;

		cerr<<"plot \"dat\" u 2:3 w d"<<endl;


		/* landing decisions */
		for(int i=1;i<=4;i++)	//leg landing
		{
			if(Fl[i] && t-tpre[i]>=P) {Fl[i]=0; count--;}	//a leg is not going to swing for too long
		}
		/** sense of supporting? still got problems **/
		double vr, Fr, Frpre;
		if(count==2)
		{
			s1x=0, s1y=0 ,s2x=0 ,s2y=0;
			temps=0;
			for(int i=1;i<=4;i++)	//recognize the supporting legs
			{
				if(!Fl[i] && temps<1) {s1x=x[i]; s1y=y[i]; temps++; continue;}
				if(!Fl[i] && temps<2) {s2x=x[i]; s2y=y[i]; temps++;}
			}
			dd=dds(xc,yc,s1x,s1y,s2x,s2y);	//distance from COM to the support
			vr=(-(s1y-s2y)*vx+(s1x-s2x)*vy) / hypot((s1y-s2y),(s1x-s2x));

			Frpre=Fr;
			Fr=(g/H) * sqrt(H*H-dd*dd) - vr*vr/H;
			if((Fr-Frpre)/dt<0)
			{
				int maxswing = 0;
				double tmax = 0;
				for(int i=1;i<=4;i++)
				{if(Fl[i] && t-tpre[i]>tmax) maxswing = i;}
				Fl[maxswing] = 0;
			}
		}
		cout<<"after landing: "<<Fl[1]<<'\t'<<Fl[2]<<'\t'<<Fl[3]<<'\t'<<Fl[4]<<endl;
	}
	return 0;
}