#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>

using namespace std;

double g=1000; // garivtational accleration in cm/s/s
double H=3; // mouse COM hight in cm
double L=5; // mouse half-length in cm
double h=1; // mouse half-width in cm
//double D=sqrt(L*L+4*h*h); // diagonal length
double omg=sqrt(g/H); // frequency

// double m=25; // mouse mass in g
double I=L*L/3; // moment of inertia over mass

double F0=10; // propulsive force in cm/s
double tau=.25;  // mass/lambda in s

double L1=1*L, L2=0*L;

int main(int argc,char** argv)
{
	double T=2, dt=.0001, DT=.01;
	double Gu1=.15, Gu2=Gu1, Gv1=.9, Gv2=Gv1;
	double taus=.01;
	double F1=10, F2=10;
	double T1=.1, T2=.1;
	int ini=0, ui=0, iu=0, lri=0;
	double kr=1, kkk=1;
	double tmin=0;
	char fdur[100]="dur";
	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i],"-T")==0) T=atof(argv[++i]);
		else if(strcmp(argv[i],"-dt")==0) DT=atof(argv[++i]);
		else if(strcmp(argv[i],"-v")==0) { F1=atof(argv[++i]); F2=atof(argv[++i]); }
		else if(strcmp(argv[i],"-u")==0) { Gu1=atof(argv[++i]); Gu2=atof(argv[++i]); }
		else if(strcmp(argv[i],"-w")==0) { Gv1=atof(argv[++i]); Gv2=atof(argv[++i]); }
		else if(strcmp(argv[i],"-s")==0) { T1=atof(argv[++i]); T2=atof(argv[++i]); }
		else if(strcmp(argv[i],"-ts")==0) taus=atof(argv[++i]);
		else if(strcmp(argv[i],"-r")==0) kr=atof(argv[++i]);
		else if(strcmp(argv[i],"-k")==0) kkk=atof(argv[++i]);
		else if(strcmp(argv[i],"-m")==0) tmin=atof(argv[++i]);
		else if(strcmp(argv[i],"-i")==0) ini=1;
		else if(strcmp(argv[i],"-t")==0) ui=1;
		else if(strcmp(argv[i],"-q")==0) iu=1;
		else if(strcmp(argv[i],"-lr")==0) lri=1;
		else if(strcmp(argv[i],"-o")==0) strcpy(fdur,argv[++i]);
		else return 1;
	}
	
	double xc=0, yc=0, theta=M_PI/4;
	double vx=F1*cos(theta), vy=F1*sin(theta), omega=0;
	double xfh[5]={},yfh[5]={},xfh0[5]={},yfh0[5]={};
	double &xfl=xfh[1], &xfr=xfh[2], &xhl=xfh[3], &xhr=xfh[4];
	double &yfl=yfh[1], &yfr=yfh[2], &yhl=yfh[3], &yhr=yfh[4];
	double &xfl0=xfh0[1], &xfr0=xfh0[2], &xhl0=xfh0[3], &xhr0=xfh0[4];
	double &yfl0=yfh0[1], &yfr0=yfh0[2], &yhl0=yfh0[3], &yhr0=yfh0[4];
	int sw[5]={0,1,0,0,0};
	double tpre[5]={};
	if(ini) 
	{
		ifstream("ini")>>xc>>yc>>vx>>vy>>theta>>omega
			>>xfl>>xfr>>xhl>>xhr
			>>yfl>>yfr>>yhl>>yhr
			>>sw[1]>>sw[2]>>sw[3]>>sw[4]
			>>tpre[1]>>tpre[2]>>tpre[3]>>tpre[4]>>F0;
		xfl-=xc; xfr-=xc; xhl-=xc; xhr-=xc;
		yfl-=yc; yfr-=yc; yhl-=yc; yhr-=yc;
		xc=0; yc=0;
	}
			
	
	double xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
	double xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);

	xfl0=xf-h*sin(theta); xfr0=xf+h*sin(theta);
	xhl0=xh-h*sin(theta); xhr0=xh+h*sin(theta);
	
	yfl0=yf+h*cos(theta); yfr0=yf-h*cos(theta);
	yhl0=yh+h*cos(theta); yhr0=yh-h*cos(theta);

	if(!ini)
	{
		xfl=xfl0, yfl=yfl0;
		xfr=xfr0, yfr=yfr0;
		xhl=xhl0-.1, yhl=yhl0-.1;
		xhr=xhr0, yhr=yhr0;
	}
	
	ofstream out(fdur);
	
	for(double t=0;t<=T;t+=dt)
	{
		int gameover=0;
		
		F0=F1+(F2-F1)*t/T;
		double Gu=Gu1+(Gu2-Gu1)*t/T;
		double Gv=Gv1+(Gv2-Gv1)*t/T;
		
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
		
		int nsw=0;
		for(int i=1;i<=4;i++) if(sw[i]) nsw++;
		
		double G[5]={},Gtot=0;
		double Zx=0,Zy=0,cosphi=1;
		if(nsw==0) for(int k=1;k<=4;k++) G[k]=F[k][0];
		else if(nsw==1) { for(int k=1;k<=4;k++) if(sw[k]) for(int i=1;i<=4;i++) G[i]=F[i][k]; }
		else if(nsw==2)
		{ 
					int k1=1; while(sw[k1]) k1++;
					double x1=xfh[k1],y1=yfh[k1];

					int k2=k1+1; while(sw[k2]) k2++;
					double x2=xfh[k2],y2=yfh[k2];

					double a=y1-y2,b=x2-x1,ab=hypot(a,b); a/=ab; b/=ab;
					double z=(xc-x2)*a+(yc-y2)*b;
					double vz=vx*a+vy*b;
					Zx=z*a; Zy=z*b;
					cosphi=1-z*z/H/H;
					if(cosphi<0) { gameover=1; cosphi=0; }
					cosphi=sqrt(cosphi);
					Gtot=cosphi-vz*vz/H/g;
					double d12=hypot(x2-x1,y2-y1);
					double r2=((xc-x2)*(x1-x2)+(yc-y2)*(y1-y2))/d12;
					double r1=((xc-x1)*(x1-x2)+(yc-y1)*(y1-y2))/d12;
					G[k1]=r2/(r2-r1)*Gtot;
					G[k2]=r1/(r1-r2)*Gtot; 
		}
		else if(nsw==3)
		{ 
					int k1=1; while(sw[k1]) k1++;
					double x1=xfh[k1],y1=yfh[k1];
					double a=xc-x1,b=yc-y1,ab=hypot(a,b); a/=ab; b/=ab;
					double z=(xc-x1)*a+(yc-y1)*b;
					double vz=vx*a+vy*b;
					Zx=-z*a; Zy=-z*b;
					cosphi=1-z*z/H/H;
					if(cosphi<0) { gameover=1; cosphi=0; }
					Gtot=sqrt(cosphi)-vz*vz/H/g;
					G[k1]=Gtot;
		}

		double d[5]={0,
			hypot(xfl0-xfl,yfl0-yfl),hypot(xfr0-xfr,yfr0-yfr),
			hypot(xhl0-xhl,yhl0-yhl),hypot(xhr0-xhr,yhr0-yhr)};

		for(int k=1;k<=4;k++) 
			if(sw[k] && d[k]>0.) { xfh[k]+=(xfh0[k]-xfh[k])*dt/d[k]/taus; 
									yfh[k]+=(yfh0[k]-yfh[k])*dt/d[k]/taus; }
		
	
		double Fx[5]={},Fy[5]={};
		double delta=0; //kr*omega;
		double FL=(1+delta)*F0/2, FR=(1-delta)*F0/2;
		
		for(int k=1;k<=4;k++) if(!sw[k] && d[k]>0.) 
		{ 
			Fx[k]=FL*(xfh0[k]-xfh[k])/d[k]; Fy[k]=FL*(yfh0[k]-yfh[k])/d[k]; 
//			Fx[k]=FL*cos(theta); Fy[k]=FL*sin(theta); 
		}

		xc+=vx*dt;
		yc+=vy*dt;
		theta+=omega*dt;
		
		vx+=(Fx[1]+Fx[2]+Fx[3]+Fx[4]-vx+tau*g*Zx/H)*dt/tau;
		vy+=(Fy[1]+Fy[2]+Fy[3]+Fy[4]-vy+tau*g*Zy/H)*dt/tau;
		
		double M=0; for(int k=1;k<=4;k++) M+=(xfh[k]-xc)*Fy[k]-(yfh[k]-yc)*Fx[k];

		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);

		
		xfl0=xf-h*sin(theta), xfr0=xf+h*sin(theta);
		xhl0=xh-h*sin(theta), xhr0=xh+h*sin(theta);
	
		yfl0=yf+h*cos(theta), yfr0=yf-h*cos(theta);
		yhl0=yh+h*cos(theta), yhr0=yh-h*cos(theta);
		
		static double phpre[5]={};
		double phi[5]={};
		double MA=0;
		for(int k=1;k<=4;k++)
		{
			double hh=hypot(yfh0[k]-yfh[k],xfh0[k]-xfh[k]);
			phi[k]=hh*sin(atan2(yfh0[k]-yfh[k],xfh0[k]-xfh[k])-theta);
			MA+=kr*(phi[k]-phpre[k])/dt;
			phpre[k]=phi[k];
		}
		
		omega+=(M/I-omega+MA)*dt/tau;
		
		
	if(int(t/dt)%int(DT/dt)==0)
	{
		cout<<t<<'\t'<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<Gtot;

		for(int k=1;k<=4;k++) cout<<'\t'<<G[k];  

		for(int i=1;i<=4;i++) cout<<'\t'<<(sw[i]?i:0);
		cout<<endl;
		
		cerr<<"unset object"<<endl;
		for(int i=1;i<=4;i++) cerr<<"set object circle front at "<<xfh[i]<<','<<yfh[i]<<",0 size "<<.2<<endl;

		cerr<<"unset arrow"<<endl;
		for(int i=1;i<=4;i++)
		{
			for(int k=1;k<=4;k++) if(!sw[i] && !sw[k])
				cerr<<"set arrow from "<<xfh[i]<<','<<yfh[i]<<" to "<<xfh[k]<<','<<yfh[k]<<" nohead"<<endl;
			
			double x0=xfh0[i]-L/2*cos(theta), y0=yfh0[i]-L/2*sin(theta);
			double a=x0-xfh[i], b=y0-yfh[i];
			double d=hypot(a,b), R=hypot(H,d), S=hypot(H,L/2);
			double q=L/2>d?(sqrt(L*L/4-d*d)/2):0;
			if(a*cos(theta)+b*sin(theta)<0) d=-d;
			if(i<=2) d=-d;
			double r=-d/2+q*H/R, z=H/2-q*d/R;
			double x1=x0, y1=y0;
			x1+=fabs(d)>.01?r*a/d:r*cos(theta); y1+=fabs(d)>.01?r*b/d:r*sin(theta);
			cerr<<"set arrow from "<<x0<<','<<y0<<","<<H<<" to "<<x1<<','<<y1<<','<<z<<" nohead"<<endl;
			cerr<<"set arrow from "<<x1<<','<<y1<<','<<z<<" to "<<xfh[i]<<','<<yfh[i]<<",0 nohead"<<endl;
		}
		
		double xt=xh-(xf-xh), yt=yh-(yf-yh);
		cerr<<"set arrow from "<<xt<<','<<yt<<","<<H<<" to "<<xf<<','<<yf<<","<<H<<" lt -1 lw 5 nohead"<<endl;

		double x1=xhl0-L/2*cos(theta), y1=yhl0-L/2*sin(theta);
		double x2=xhr0-L/2*cos(theta), y2=yhr0-L/2*sin(theta);
		double x3=xfl0-L/2*cos(theta), y3=yfl0-L/2*sin(theta);
		double x4=xfr0-L/2*cos(theta), y4=yfr0-L/2*sin(theta);
		cerr<<"set arrow from "<<x1<<','<<y1<<","<<H<<" to "<<x2<<','<<y2<<","<<H<<" lt -1 lw 3 nohead"<<endl;
		cerr<<"set arrow from "<<x3<<','<<y3<<","<<H<<" to "<<x4<<','<<y4<<","<<H<<" lt -1 lw 3 nohead"<<endl;
		
		cerr<<"set xrange ["<<(xc-20)<<':'<<(xc+20)<<']'<<endl;
		cerr<<"set yrange ["<<(yc-20)<<':'<<(yc+20)<<']'<<endl;
		cerr<<"set title \"V="<<int(hypot(vx,vy))<<'\"'<<endl;
		cerr<<"splot \"dat\" u 2:3:(0) w d"<<endl;
	}
		if(gameover) return 1;

		double tmax=T1+(T2-T1)*t/T;
		for(int k=1;k<=4;k++) if(sw[k] && (t-tpre[k]>=tmax))
		{
			sw[k]=0; out<<t<<'\t'<<F0<<'\t'<<k<<'\t'<<sw[k]<<'\t'<<(t-tpre[k])<<'\t'<<5<<endl; tpre[k]=t; nsw--;
		}
		
		double DM[5]={0,L,L,kkk*L,kkk*L};

		static double GP[5]={};
		int lr[5]={0,2,1,4,3};
		for(int k=1;k<=4;k++) if(!sw[k]) if(G[k]>0 && GP[k]>G[k] && G[k]<Gu) 
		{ 
			if((lri && sw[lr[k]]) || t-tpre[k]<tmin) continue;
			sw[k]=1; out<<t<<'\t'<<F0<<'\t'<<k<<'\t'<<sw[k]<<'\t'<<(t-tpre[k])<<'\t'<<4<<endl; tpre[k]=t; nsw++;
//			if(ui && k<=2 && sw[k+2]) 
//			{ 
//				sw[k+2]=0; 
//				out<<t<<'\t'<<F0<<'\t'<<(k+2)<<'\t'<<sw[k+2]<<'\t'<<(t-tpre[k+2])<<'\t'<<3<<endl; tpre[k+2]=t; 
//			}
		}
		for(int k=1;k<=4;k++) if(!sw[k] && t-tpre[k]>tmin) if(G[k]<0 || d[k]>DM[k]) 
		{ 
			sw[k]=1; out<<t<<'\t'<<F0<<'\t'<<k<<'\t'<<sw[k]<<'\t'<<(t-tpre[k])<<'\t'<<(G[k]<0?0:1)<<endl; tpre[k]=t; nsw++;
			int ka=lr[k];
			if(lri && sw[ka]) { sw[ka]=0; out<<t<<'\t'<<F0<<'\t'<<ka<<'\t'<<sw[ka]<<'\t'<<(t-tpre[ka])<<'\t'<<2<<endl; tpre[ka]=t; nsw--; }
			if(ui && k<=2) if(sw[k+2]) 
			{ 
				sw[k+2]=0; nsw--;
				out<<t<<'\t'<<F0<<'\t'<<(k+2)<<'\t'<<sw[k+2]<<'\t'<<(t-tpre[k+2])<<'\t'<<3<<endl; tpre[k+2]=t; 
			}
			if(iu && k>2) if(sw[k-2]) 
			{ 
				sw[k-2]=0; nsw--;
				out<<t<<'\t'<<F0<<'\t'<<(k-2)<<'\t'<<sw[k-2]<<'\t'<<(t-tpre[k-2])<<'\t'<<3<<endl; tpre[k-2]=t; 
			}
		}
		
		for(int k=1;k<=4;k++) GP[k]=G[k];
		
		static double Gpre=0; 
		double bal=Gtot>0?(Gpre-Gtot)/dt:0;
		Gpre=Gtot;
		
		double tsw[5]={};
		for(int k=1;k<=4;k++) if(sw[k] && (nsw>2 || bal>Gv)) tsw[k]=t-tpre[k];
		int kmax=0;
		for(int k=1;k<=4;k++) if(tsw[k]>tsw[kmax]) kmax=k;
		sw[kmax]=0; if(kmax && t-tpre[kmax]>tmin) 
		{
			out<<t<<'\t'<<F0<<'\t'<<kmax<<'\t'<<sw[kmax]<<'\t'<<(t-tpre[kmax])<<'\t'<<6<<endl; tpre[kmax]=t;
			nsw--;
		}
	}
	for(int k=1;k<=4;k++) tpre[k]-=T;
	ofstream("ini")<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<theta<<'\t'<<omega<<endl
		<<xfl<<'\t'<<xfr<<'\t'<<xhl<<'\t'<<xhr<<endl
		<<yfl<<'\t'<<yfr<<'\t'<<yhl<<'\t'<<yhr<<endl
		<<sw[1]<<'\t'<<sw[2]<<'\t'<<sw[3]<<'\t'<<sw[4]<<endl
		<<tpre[1]<<'\t'<<tpre[2]<<'\t'<<tpre[3]<<'\t'<<tpre[4]<<'\t'<<F0<<endl;
	return 0;
}

