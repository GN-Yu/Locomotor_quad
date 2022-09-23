#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>

using namespace std;

double g=1000; // garivtational accleration in cm/s/s
double H=3; // mouse hight in cm
double L=5; // mouse half-length in cm
double h=1.5; // mouse half-width in cm
//double D=sqrt(L*L+4*h*h); // diagonal length
//double omega=sqrt(g/H); // frequency

// double m=25; // mouse mass in g
double I=L*L/3; // moment of inertia over mass

double F0=10; // locomotor speed
double tau=.25;

double L1=1*L, L2=0*L;

double EL = -58, C = .01, gL = 2.5, v12mp = -40, kmp = -6, v12hp = -55, khp = 10, tauhp = .1, gNaP = 5, ENa = 50;
double EK=-90, Ee=0, Ei=-75, v12n=-30, kn=-4, gKdr=5;

inline double ninf(double v) { return 1./(1.+exp((v-v12n)/kn)); }
inline double mpinf(double v) { return 1./(1.+exp((v-v12mp)/kmp)); }
inline double hpinf(double v) { return 1./(1.+exp((v-v12hp)/khp)); }
inline double tauinf(double v) { return tauhp/cosh((v-v12hp)/khp); }

inline double f(double v) { return v<-50.?0.:((v+50.)/50.); }



int main(int argc,char** argv)
{
	double T=2, dt=.00001, DT=.01;
	double Gu=.15, Gv=.9;
	double taus=.01;
	double F1=10, F2=10;
	int ini=0;
	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i],"-T")==0) T=atof(argv[++i]);
		else if(strcmp(argv[i],"-dt")==0) DT=atof(argv[++i]);
		else if(strcmp(argv[i],"-v")==0) { F1=atof(argv[++i]); F2=atof(argv[++i]); }
		else if(strcmp(argv[i],"-u")==0) Gu=atof(argv[++i]);
		else if(strcmp(argv[i],"-w")==0) Gv=atof(argv[++i]);
		else if(strcmp(argv[i],"-ts")==0) taus=atof(argv[++i]);
		else if(strcmp(argv[i],"-el")==0) EL=atof(argv[++i]);
		else if(strcmp(argv[i],"-i")==0) ini=1;
		else return 1;
	}
	
	double xc=0, yc=0, theta=M_PI/4;
	double vx=3*F0*cos(theta), vy=3*F0*sin(theta), omega=0;
	static double v[5]={0,-30,-60,-60,-60}, hp[5]={0,.7,.4,.3,.6};
	double xfl, xfr, xhl, xhr;
	double yfl, yfr, yhl, yhr;
	if(ini) ifstream("ini")>>xc>>yc>>vx>>vy>>theta>>omega
			>>xfl>>xfr>>xhl>>xhr
			>>yfl>>yfr>>yhl>>yhr
			>>v[1]>>v[2]>>v[3]>>v[4]>>hp[1]>>hp[2]>>hp[3]>>hp[4];
//	xc=yc=0;
	
	double xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
	double xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);

	double xfl0=xf-h*sin(theta), xfr0=xf+h*sin(theta);
	double xhl0=xh-h*sin(theta), xhr0=xh+h*sin(theta);
	
	double yfl0=yf+h*cos(theta), yfr0=yf-h*cos(theta);
	double yhl0=yh+h*cos(theta), yhr0=yh-h*cos(theta);

	if(!ini)
	{
		xfl=xfl0, yfl=yfl0;
		xfr=xfr0, yfr=yfr0;
		xhl=xhl0-.1, yhl=yhl0-.1;
		xhr=xhr0, yhr=yhr0;
	}
	
	
	for(double t=0;t<=T;t+=dt)
	{
		F0=F1+(F2-F1)*t/T;
		
		double u1=xfl-xc,u2=xfr-xc,u3=xhl-xc,u4=xhr-xc;
		double v1=yfl-yc,v2=yfr-yc,v3=yhl-yc,v4=yhr-yc;
		double F[5][5]={};

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

		int fl[5]={1,1,1,1,1};
		for(int i=1;i<=4;i++) for(int k=1;k<=4;k++) if(i!=k) fl[i]=fl[i]&&(F[k][i]>=0);
		
		static int sw[5]={0,1,0,0,0};
		for(int k=1;k<=4;k++) sw[k]=f(v[k])>.01;
		
		int nsw=0;
		for(int i=1;i<=4;i++) if(sw[i]) nsw++;
		
		double G[5]={},Gtot=0;
		double Zx=0,Zy=0;
		if(nsw==0) for(int k=1;k<=4;k++) G[k]=F[k][0];
		else if(nsw==1) { for(int k=1;k<=4;k++) if(sw[k]) for(int i=1;i<=4;i++) G[i]=F[i][k]; }
		else if(nsw==2)
		{ 
					int k1=1; while(sw[k1]) k1++;
					double x1,y1;
					switch(k1) 
					{
						case 1: x1=xfl; y1=yfl; break;
						case 2: x1=xfr; y1=yfr; break;
						case 3: x1=xhl; y1=yhl; break;
						case 4: x1=xhr; y1=yhr; break;
					}
					int k2=k1+1; while(sw[k2]) k2++;
					double x2,y2;
					switch(k2) 
					{
						case 1: x2=xfl; y2=yfl; break;
						case 2: x2=xfr; y2=yfr; break;
						case 3: x2=xhl; y2=yhl; break;
						case 4: x2=xhr; y2=yhr; break;
					}
					double a=y1-y2,b=x2-x1,ab=hypot(a,b); a/=ab; b/=ab;
					double z=(xc-x2)*a+(yc-y2)*b;
					double vz=vx*a+vy*b;
					Zx=z*a; Zy=z*b;
					Gtot=1-z*z/H/H;
					if(Gtot<0) { cerr<<"game over"; break; }
					Gtot=sqrt(Gtot)-vz*vz/H/g;
					double x0=xc-z*a, y0=yc-z*b;
					double d12=hypot(x2-x1,y2-y1);
					double r2=((xc-x2)*(x1-x2)+(yc-y2)*(y1-y2))/d12;
					double r1=((xc-x1)*(x1-x2)+(yc-y1)*(y1-y2))/d12;
					G[k1]=r2/(r2-r1)*Gtot;
					G[k2]=r1/(r1-r2)*Gtot; 
		}
		
		double d[5]={0,
			hypot(xfl0-xfl,yfl0-yfl),hypot(xfr0-xfr,yfr0-yfr),
			hypot(xhl0-xhl,yhl0-yhl),hypot(xhr0-xhr,yhr0-yhr)};
		double dmax=0,dsec=0;
		int imax=0, isec=0;
		for(int i=1;i<5;i++) if(d[i]>=dmax) { dsec=dmax; dmax=d[i]; isec=imax; imax=i; }

		if(sw[1]) { xfl+=(xfl0-xfl)*dt/taus; yfl+=(yfl0-yfl)*dt/taus; }
		if(sw[2]) { xfr+=(xfr0-xfr)*dt/taus; yfr+=(yfr0-yfr)*dt/taus; }
		if(sw[3]) { xhl+=(xhl0-xhl)*dt/taus; yhl+=(yhl0-yhl)*dt/taus; }
		if(sw[4]) { xhr+=(xhr0-xhr)*dt/taus; yhr+=(yhr0-yhr)*dt/taus; }

		double Fx[5]={},Fy[5]={};
		double delta=0;
		double FL=(1+delta)*F0, FR=(1-delta)*F0;
		
		if(!sw[1] && d[1]!=0) { Fx[1]=FL*(xfl0-xfl)/d[1]; Fy[1]=FL*(yfl0-yfl)/d[1]; }
		if(!sw[2] && d[2]!=0) { Fx[2]=FR*(xfr0-xfr)/d[2]; Fy[2]=FR*(yfr0-yfr)/d[2]; }
		if(!sw[3] && d[3]!=0) { Fx[3]=FL*(xhl0-xhl)/d[3]; Fy[3]=FL*(yhl0-yhl)/d[3]; }
		if(!sw[4] && d[4]!=0) { Fx[4]=FR*(xhr0-xhr)/d[4]; Fy[4]=FR*(yhr0-yhr)/d[4]; }

		xc+=vx*dt;
		yc+=vy*dt;
		theta+=omega*dt;

		vx+=(Fx[1]+Fx[2]+Fx[3]+Fx[4]-vx+tau*g*Zx/H)*dt/tau;
		vy+=(Fy[1]+Fy[2]+Fy[3]+Fy[4]-vy+tau*g*Zy/H)*dt/tau;
		
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
		
		double gLoad=2, gLength=1, gBal=2;
		double P=1, DM[5]={0,L,L,2*L,2*L};
		static double GS=0; 
		double Gpre=GS;
		GS=0; for(int i=1;i<=4;i++) GS+=G[i];

		double ul=2,lr=2;
		double b[4][4]={0,lr,ul,0,
						lr,0,0,ul,
						ul,0,0,lr,
						0,ul,lr,0};

		double diag=0; //1.5;
		double a[4][4]={0,0,0,diag,
						0,0,diag,0,
						0,diag,0,0,
						diag,0,0,0};
						
		double bal=(Gpre-GS)/dt;
		for(int k=1;k<=4;k++)
		{
			double n=ninf(v[k]); n*=n;
			double vp=-gKdr*n*n*(v[k]-EK)-gL*(v[k]-EL)-gNaP*mpinf(v[k])*hp[k]*(v[k]-ENa);
			
			for(int i=1;i<=4;i++) vp-=b[i-1][k-1]*f(v[i])*(v[k]-Ei)+a[i-1][k-1]*f(v[i])*(v[k]-Ee);
			
			vp-=gLoad*(G[k]>0?pow(G[k],.5):0)*(v[k]-Ei);
			
			double lfb=d[k]/DM[k]; lfb*=lfb;
			vp-=gLength*(lfb*lfb)*(v[k]-Ee);
			
			vp-=gBal*((Gpre<Gv && bal>0)?(bal<10?bal:10):0)*(v[k]-Ei);
			
			v[k]+=vp/C*dt;
			double hpi=hpinf(v[k]);
			hp[k]+=(hpi-hp[k])/tauhp*dt;
		}
		
	if(int(t/dt)%int(DT/dt)==0)
	{
		cout<<t<<'\t'<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<Gtot;
		for(int k=1;k<=4;k++) cout<<'\t'<<v[k];
		
		for(int k=1;k<=4;k++) cout<<'\t'<<G[k];  

		for(int i=1;i<=4;i++) cout<<'\t'<<(sw[i]?i:0);
		cout<<endl;
		
		cerr<<"unset object"<<endl;
		cerr<<"set object 1 circle front at "<<xfl<<','<<yfl<<" size "<<.2<<" fc \"brown\""<<endl;
		cerr<<"set object 2 circle front at "<<xfr<<','<<yfr<<" size "<<.2<<" fc \"blue\""<<endl;
		cerr<<"set object 3 circle front at "<<xhl<<','<<yhl<<" size "<<.2<<" fc \"red\""<<endl;
		cerr<<"set object 4 circle front at "<<xhr<<','<<yhr<<" size "<<.2<<" fc \"green\""<<endl;
		cerr<<"set object 5 circle front at "<<xc<<','<<yc<<" size "<<.2<<" fc \"black\""<<endl;


		cerr<<"unset arrow"<<endl;
		if(!sw[1] && !sw[4]) cerr<<"set arrow 1 from "<<xfl<<','<<yfl<<" to "<<xhr<<','<<yhr<<" nohead dt 2"<<endl;
		if(!sw[2] && !sw[3]) cerr<<"set arrow 2 from "<<xfr<<','<<yfr<<" to "<<xhl<<','<<yhl<<" nohead dt 2"<<endl;
		if(!sw[2] && !sw[1]) cerr<<"set arrow 3 from "<<xfr<<','<<yfr<<" to "<<xfl<<','<<yfl<<" nohead dt 2"<<endl;
		if(!sw[3] && !sw[4]) cerr<<"set arrow 4 from "<<xhr<<','<<yhr<<" to "<<xhl<<','<<yhl<<" nohead dt 2"<<endl;
		cerr<<"set arrow 5 from "<<xf<<','<<yf<<" to "<<xfl<<','<<yfl<<" nohead"<<endl;
		cerr<<"set arrow 6 from "<<xf<<','<<yf<<" to "<<xfr<<','<<yfr<<" nohead"<<endl;
		cerr<<"set arrow 7 from "<<xh<<','<<yh<<" to "<<xhl<<','<<yhl<<" nohead"<<endl;
		cerr<<"set arrow 8 from "<<xh<<','<<yh<<" to "<<xhr<<','<<yhr<<" nohead"<<endl;
		cerr<<"set arrow 9 from "<<xh<<','<<yh<<" to "<<xf<<','<<yf<<" lt -1 lw 5 nohead"<<endl;
		
		cerr<<"plot \"dat\" u 2:3 w d"<<endl;
	}
		
/*		static double tpre[5]={};
		
		int flcont=0;
		

		static double GP[5]={};
		for(int k=1;k<=4;k++) if(!sw[k]) if((GP[k]>=Gu && G[k]<Gu) || d[k]>DM[k]) 
		{ 
			sw[k]=1; 
			tpre[k]=t; flcont=1; 
		}
//		if(flcont) continue;
		
//		double tmax=.2;
//		for(int k=1;k<=4;k++) if(sw[k] && t-tpre[k]>=tmax) { sw[k]=0; flcont=1; }
//		if(flcont) continue;
		
		for(int k=1;k<=4;k++) GP[k]=G[k];

		double tsw[5]={};
		for(int k=1;k<=4;k++) if(sw[k] && ((Gpre<Gv && GS<Gpre) || GS<=0.)) tsw[k]=t-tpre[k];
		int kmax=0;
		for(int k=1;k<=4;k++) if(tsw[k]>tsw[kmax]) kmax=k;
		sw[kmax]=0;

//		if(GS<=0.) for(int k=1;k<=4;k++) sw[k]=0;

*/
	}
	ofstream("ini")<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<theta<<'\t'<<omega<<endl
		<<xfl<<'\t'<<xfr<<'\t'<<xhl<<'\t'<<xhr<<endl
		<<yfl<<'\t'<<yfr<<'\t'<<yhl<<'\t'<<yhr<<endl
		<<v[1]<<'\t'<<v[2]<<'\t'<<v[3]<<'\t'<<v[4]<<endl
		<<hp[1]<<'\t'<<hp[2]<<'\t'<<hp[3]<<'\t'<<hp[4]<<endl;
	return 0;
}

