#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>

using namespace std;

//units: length - cm, time - s, mass - g

double g=1000;	//garivtational accleration
double H=3;		//hight
double lf=1.5;	//half shoulder length
double lh=1.5;	//half hip length
double h=1.5;	//half width
double l=4.1;	//half body length
double L1=1*l, L2=0*l;		//front, hind part body length
double L=L1+L2;				//total body length
double D=sqrt(L*L+(lh+lf)*(lh+lf)); 	//diagonal length
double freq=sqrt(g/H);			//inverted pendulum frequency
double Tswc;				//swing time

// double m=30; // mouse mass in g
double I=L*L/3; // moment of inertia over mass
double tau=.25;

double DM[5]={0,L,L,1.3*L,1.3*L};	//max leg length
static double phi[5]={0,0.1,-M_PI,-M_PI,-M_PI}; //swing phase indicator
static double Omega[5]={}; //swing phase indicator angular velocity

int contraside[5]={0,2,1,4,3};
int homoside[5]={0,3,4,1,2};

string colors[5]={"black","brown","blue","red","green"};
string inifiles[4]={"","CPGini","CPGslowini","CPGfastini"};


int main(int argc,char** argv)
{
	//output settings
	ofstream out_timers("strideinfo");
	ofstream out_polar("ani_isolators.plot");

	cerr<<"set xtics 0,5,1000"<<endl;
	cerr<<"set ytics 0,5,2000"<<endl;
	cerr<<"set size ratio -1"<<endl;
	cerr<<"set grid"<<endl;
	cerr<<"set xlabel 'x'"<<endl;
	cerr<<"set ylabel 'y'"<<endl;


	//variable declarations
	double T, dt, DT;
	dt=.0001;
	int ini=0;
	int newini=0;
	int cor=0;

	double Guini, Gufin; 
	double Gvini, Gvfin;	
	double Gu, Gv;	//thresholds
	double Fini, Ffin;
	double F0;		//forces in the unit cm/s, original force F/lambda
	double xc, yc;	//position of COM
	double vx, vy;	//velocity of COM
	double vv, vexp;//speed in magnitude, expected speed
	double omega;	//body angular velocity
	double theta;	//body angle

	double xf, yf, xh, yh;	//middile of shoulder and hip positions
	double xfh[5]={},yfh[5]={},xfh0[5]={},yfh0[5]={},xfh1[5]={},yfh1[5]={};
	double &xfl=xfh[1], &xfr=xfh[2], &xhl=xfh[3], &xhr=xfh[4];
	double &yfl=yfh[1], &yfr=yfh[2], &yhl=yfh[3], &yhr=yfh[4];			//limb positions
	double &xfl0=xfh0[1], &xfr0=xfh0[2], &xhl0=xfh0[3], &xhr0=xfh0[4];
	double &yfl0=yfh0[1], &yfr0=yfh0[2], &yhl0=yfh0[3], &yhr0=yfh0[4];	//joint positions
	double &xfl1=xfh1[1], &xfr1=xfh1[2], &xhl1=xfh1[3], &xhr1=xfh1[4];
	double &yfl1=yfh1[1], &yfr1=yfh1[2], &yhl1=yfh1[3], &yhr1=yfh1[4];	//corrected landing positions

	double d[5]={};	//the distance from a limb to its shoulder or hip

	int nsw;		//number of swinging limbs
	double G[5]={};	//loads of each foot
	double Gtot;	//total load
	double Zx, Zy;	//distance to diagonal when two limbs supporting
	int flcont=0;
	int flcontpre=0;
	static double Gpre=0; //total supporting force of the previous time step
	static double GP[5]={};
	static double tswpre[5]={0,0,0,0,0};
	double tsw[5]={};

	int sw[5]={};

	double kv;
	double kvini;
	double kvfin;
    double Tswcini;
	double Tswcfin;

	double stdOmega=M_PI/0.05;
	for(int i=1;i<=4;i++) Omega[i]=stdOmega;
	

	//initialize parameters
	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i],"-T")==0) T=atof(argv[++i]);
		else if(strcmp(argv[i],"-fps")==0) DT=1/atof(argv[++i]);
		//else if(strcmp(argv[i],"-v")==0) F0=atof(argv[++i]);
		else if(strcmp(argv[i],"-v0")==0) Fini=atof(argv[++i]);
		else if(strcmp(argv[i],"-v1")==0) Ffin=atof(argv[++i]);
		else if(strcmp(argv[i],"-Gu0")==0) Guini=atof(argv[++i]);
		else if(strcmp(argv[i],"-Gu1")==0) Gufin=atof(argv[++i]);
		else if(strcmp(argv[i],"-Gv0")==0) Gvini=atof(argv[++i]);
		else if(strcmp(argv[i],"-Gv1")==0) Gvfin=atof(argv[++i]);
		else if(strcmp(argv[i],"-kv0")==0) kvini=atof(argv[++i]);
		else if(strcmp(argv[i],"-kv1")==0) kvfin=atof(argv[++i]);
        else if(strcmp(argv[i],"-tsw0")==0) Tswcini=atof(argv[++i]);
		else if(strcmp(argv[i],"-tsw1")==0) Tswcfin=atof(argv[++i]);
		else if(strcmp(argv[i],"-i")==0) {ini=1;}
		else if(strcmp(argv[i],"-slowi")==0) {ini=2;}
		else if(strcmp(argv[i],"-fasti")==0) {ini=3;}
		else if(strcmp(argv[i],"-newi")==0) {newini=1;}
		else if(strcmp(argv[i],"-newslowi")==0) {newini=2;}
		else if(strcmp(argv[i],"-newfasti")==0) {newini=3;}
		else if(strcmp(argv[i],"-cor")==0) {cor=1;}
		else return 1;
	}
	F0=Fini;
	Gu=Guini;
	Gv=Gvini;
	kv=kvini;

	if(ini==0)
	{
		xc=0; yc=0;
		theta=M_PI/4;
		vx=3*F0*cos(theta); vy=3*F0*sin(theta);
		omega=0;

		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);
		xfl0=xf-lf*sin(theta), xfr0=xf+lf*sin(theta);
		xhl0=xh-lh*sin(theta), xhr0=xh+lh*sin(theta);
		yfl0=yf+lf*cos(theta), yfr0=yf-lf*cos(theta);
		yhl0=yh+lh*cos(theta), yhr0=yh-lh*cos(theta);

		xfl=xfl0; yfl=yfl0;
		xfr=xfr0; yfr=yfr0;
		xhl=xhl0-.1; yhl=yhl0-.1;
		xhr=xhr0; yhr=yhr0;
	}
	else 
	{
		ifstream(inifiles[ini])>>xc>>yc>>vx>>vy>>theta>>omega>>Gpre
		>>xfl>>xfr>>xhl>>xhr>>yfl>>yfr>>yhl>>yhr
		>>phi[1]>>phi[2]>>phi[3]>>phi[4];

		for(int i=1;i<=4;i++) {xfh[i]-=xc; yfh[i]-=yc;}
		xc=0; yc=0;	//reset to (0,0) position
		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);
		xfl0=xf-lf*sin(theta), xfr0=xf+lf*sin(theta);
		xhl0=xh-lh*sin(theta), xhr0=xh+lh*sin(theta);
		yfl0=yf+lf*cos(theta), yfr0=yf-lf*cos(theta);
		yhl0=yh+lh*cos(theta), yhr0=yh-lh*cos(theta);
	}
	vv=sqrt(vx*vx+vy*vy);
	vexp=vv;


	//step cycles
	for(double t=0;t<=T;t+=dt)
	{
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
		
		nsw=0;
		for(int i=1;i<=4;i++) G[i]=0;
		Gtot=0;
		Zx=0; Zy=0;

		for(int i=1;i<=4;i++) {sw[i]=(sin(phi[i])>0?1:0);}
		for(int i=1;i<=4;i++) if(sw[i]) nsw++;	//swing legs count

		if(nsw==0) {for(int k=1;k<=4;k++) G[k]=F[k][0];}
		else if(nsw==1) {for(int k=1;k<=4;k++) if(sw[k]) for(int i=1;i<=4;i++) G[i]=F[i][k];}
		else if(nsw==2)
		{ 
            int k1=1; while(sw[k1]) k1++;
            double x1=xfh[k1],y1=yfh[k1];
            int k2=k1+1; while(sw[k2]) k2++;
            double x2=xfh[k2],y2=yfh[k2];
//			cerr<<k1<<'\t'<<k2<<endl;
            double a=y1-y2,b=x2-x1,ab=hypot(a,b); a/=ab; b/=ab;
            double z=(xc-x2)*a+(yc-y2)*b;
            Zx=z*a; Zy=z*b;
			double vz=vx*a+vy*b;
            Gtot=1-z*z/H/H;
			if(Gtot<0) { cerr<<"pisdez"; break; }
			Gtot=sqrt(Gtot)-vz*vz/H/g;
            double x0=xc-z*a, y0=yc-z*b;
            double d12=hypot(x2-x1,y2-y1);
            double r2=((xc-x2)*(x1-x2)+(yc-y2)*(y1-y2))/d12;
            double r1=((xc-x1)*(x1-x2)+(yc-y1)*(y1-y2))/d12;
            G[k1]=r2/(r2-r1)*Gtot;
            G[k2]=r1/(r1-r2)*Gtot;
//			cerr<<G[k1]<<'\t'<<G[k2]<<endl; return 1; 
		}
		else if(nsw==3){for(int k=1;k<=4;k++) if(!sw[k]) G[k]=1;}
		// {
		// 	int k1=1; while(sw[k1]) k1++;
		// 	double x1=xfh[k1],y1=yfh[k1];
		// 	double a=xc-x1,b=yc-y1,ab=hypot(a,b); a/=ab; b/=ab;
		// 	double z=(xc-x1)*a+(yc-y1)*b;
		// 	double vz=vx*a+vy*b;
		// 	Zx=-z*a; Zy=-z*b;
		// 	Gtot=1-z*z/H/H;
		// 	// if(Gtot<0) { cerr<<"pisdez"; break; }
		// 	Gtot=sqrt(Gtot)-vz*vz/H/g;
		// 	G[k1]=Gtot;
		// }
		// {for(int k=1;k<=4;k++) if(!sw[k]) G[k]=1;}
		// else if(nsw==4){
		// 	cerr<<"flying"<<'\t'<<flcont<<'\t'<<Gu<<endl; 
		// 	cerr<<G[1]<<'\t'<<G[2]<<'\t'<<G[3]<<'\t'<<G[4]<<endl; 
		// 	cerr<<GP[1]<<'\t'<<GP[2]<<'\t'<<GP[3]<<'\t'<<GP[4]<<endl;
		// 	cerr<<tsw[1]<<'\t'<<tsw[2]<<'\t'<<tsw[3]<<'\t'<<tsw[4]<<endl;
		// 	break;}
		
		
		//calculate corrections of limb positionings when swing
		for(int i=1;i<=4;i++) {d[i]=hypot(xfh0[i]-xfh[i],yfh0[i]-yfh[i]);}
		double vtor=vx*(yfl0-yhr0)/D - vy*(xfl0-xhr0)/D;
		double vtol=-vx*(yfr0-yhl0)/D + vy*(xfr0-xhl0)/D;
		double deltar=(L/2-L2)-vtor*D/2/h/freq+2*vexp/freq/(1+exp(freq*L/vexp));
		double deltal=(L/2-L2)-vtol*D/2/h/freq+2*vexp/freq/(1+exp(freq*L/vexp));

		xfl1=xfl0-deltar*cos(theta); yfl1=yfl0-deltar*sin(theta);
		xhr1=xhr0-deltar*cos(theta); yhr1=yhr0-deltar*sin(theta);
		xfr1=xfr0-deltal*cos(theta); yfr1=yfr0-deltal*sin(theta);
		xhl1=xhl0-deltal*cos(theta); yhl1=yhl0-deltal*sin(theta);

		double taus=.001;
		double Vswc=(L+vv*Tswc)/Tswc;

		//swing limb movement
		for(int i=1;i<=4;i++) if(sw[i])
		{
			if(cor) {xfh[i]+=(xfh1[i]-xfh[i])*dt/taus; yfh[i]+=(yfh1[i]-yfh[i])*dt/taus;}
			else {xfh[i]+=(xfh0[i]-xfh[i])*dt/taus; yfh[i]+=(yfh0[i]-yfh[i])*dt/taus;}
		}

		//horizontal force calculations
		F0=Fini+(Ffin-Fini)*t/T;	//F0 update
		double Fx[5]={},Fy[5]={};
		// double delta=0*omega;
		// double FL=(1+delta)*F0, FR=(1-delta)*F0;
		double Fbase[5];
		for(int i=1;i<=4;i++) Fbase[i]=F0;
		for(int i=1;i<=4;i++) if(!sw[i] && d[i]>0.) {Fx[i]=Fbase[i]*(xfh0[i]-xfh[i])/d[i]; Fy[i]=Fbase[i]*(yfh0[i]-yfh[i])/d[i];}

		//update variables
		xc+=vx*dt;
		yc+=vy*dt;
		theta+=omega*dt;
		vx+=((Fx[1]+Fx[2]+Fx[3]+Fx[4])/tau-vx/tau+g*Zx/H)*dt;
		vy+=((Fy[1]+Fy[2]+Fy[3]+Fy[4])/tau-vy/tau+g*Zy/H)*dt;
		vv=sqrt(vx*vx+vy*vy);

		double tauvv=1;
		vexp+=(vv-vexp)*dt/tauvv;

		xf=xc+L1*cos(theta); yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta); yh=yc-L2*sin(theta);
		xfl0=xf-lf*sin(theta), xfr0=xf+lf*sin(theta);
		xhl0=xh-lh*sin(theta), xhr0=xh+lh*sin(theta);
		yfl0=yf+lf*cos(theta), yfr0=yf-lf*cos(theta);
		yhl0=yh+lh*cos(theta), yhr0=yh-lh*cos(theta);

		double M=0;
		for(int k=1;k<=4;k++) {M+=(xfh[k]-xc)*Fy[k]-(yfh[k]-yc)*Fx[k];}
		omega+=(M/I-omega)*dt/tau;

		Gu=Guini+(Gufin-Guini)*t/T;
		Gv=Gvini+(Gvfin-Gvini)*t/T;
		kv=kvini+(kvfin-kvini)*t/T;
        Tswc=Tswcini+(Tswcfin-Tswcini)*t/T;
		Tswc=0.5;

		//phase isolator updates
		stdOmega=M_PI/Tswc;
		double Omega_signal[5]={};
		for(int i=1;i<=4;i++)
		{
			Omega[i]=stdOmega;
			Omega_signal[i]=0;
		}

		double small_number=10;
		double normal_number=30;
        double large_number=100;
		double very_large_number=1000;

		// // for(int i=1;i<=4;i++) if(!sw[i]) if(GP[i]>G[i] && G[i]<Gu)
		// // {
		// // 	// if(sw[contraside[i]]==1) {Omega[i]+=A*sin(phi[i]-phi[contraside[i]]);}
		// // 	// if(sw[homoside[i]]==1) {Omega[i]+=B*sin(phi[i]-phi[homoside[i]]);}
		// // 	Omega[i]+=stdOmega;
		// // }	//weak lifting conditions

		// for(int i=1;i<=4;i++) if(!sw[i])
		// {
		// 	if(G[i]<=0) 
		// 	{
		// 		Omega_signal[i]+=very_large_number;
		// 	}
		// 	else if(G[i]>Gu){Omega_signal[i]+= -1*stdOmega;}
		// 	else {Omega_signal[i]+=2*stdOmega;}
		// 	// else {omega_signal[i]=stdOmega-stdOmega*sqrt(G[i]);}
		// 	// Omega_signal[i]+=stdOmega*(1-G[i]/Gu);
		// 	Omega_signal[i]+=pow(d[i]/DM[i],8);
		// }

		// // for(int i=1;i<=4;i++) if(!sw[i]) if(d[i]>DM[i])
		// // {
		// // 	Omega_signal[i]+=3*stdOmega;
		// // }

		// if((Gtot!=0. && -(Gtot-Gpre)/dt>kv)) for(int i=1;i<=4;i++) if(sw[i])
		// {
		// 	Omega_signal[i]+=4*stdOmega;
		// }

		// // if(Gtot!=0. && (Gtot-Gpre)/dt<0) for(int i=1;i<=4;i++) if(sw[i])
		// // {
		// // 	Omega_signal[i]+= -kv*(Gtot-Gpre)/dt;
		// // }

		for(int i=1;i<=4;i++) if(!sw[i]) if(GP[i]>G[i] && G[i]<Gu)
		{
			if(sw[contraside[i]]==1 || sw[homoside[i]]==1) {Omega_signal[i]+=small_number;}
			else{Omega_signal[i]+=very_large_number;}
		}	//weak lifting conditions --> add a bit to omega

		for(int i=1;i<=4;i++) if(!sw[i]) if(G[i]<0 || d[i]>DM[i])
		{
			Omega_signal[i]+=very_large_number;
			Omega_signal[contraside[i]]+=very_large_number*sin(phi[contraside[i]]-phi[i]);
			Omega_signal[homoside[i]]+=very_large_number*sin(phi[homoside[i]]-phi[i]);
		}	//strong lifting conditions --> add very large number to omega
		
		if(flcontpre==0 && ((Gtot!=0 && -(Gtot-Gpre)/dt>kv) || nsw>2))
		{
			for(int i=1;i<=4;i++) if(sw[i]) {Omega_signal[i]+=normal_number;}
		}	//stop swing a leg for lossing balance --> add not so large number to omega

		for(int i=0;i<=4;i++) 
		{
			// Omega_signal[i]+=large_number*sin(phi[i]-phi[contraside[i]]);
			// Omega_signal[i]+=large_number*sin(phi[i]-phi[homoside[i]]);
			if(Omega[i]+Omega_signal[i]>0)	{phi[i]+=(Omega[i]+Omega_signal[i])*dt;}
			if(phi[i]>M_PI) phi[i]-=2*M_PI;
			if(phi[i]<-M_PI) phi[i]+=2*M_PI;
		}

		for(int i=1;i<=4;i++) {GP[i]=G[i];}
		Gpre=Gtot;

		
		//timers
		static int swpre[5];
		static double ttswpre[5]={};
		double ttsw[5]={};
		static double ttstpre[5]={};
		double ttst[5]={};
		static double stridepre=0;
		static double stridetime=0;
		for(int k=1;k<=4;k++) if(swpre[k] && !sw[k])	//stop swing tracker
		{
			ttstpre[k]=t; 
			ttsw[k]=t-ttswpre[k];
			if(k==1) { stridetime=t-stridepre; stridepre=t; }
			if(stridepre!=0) out_timers<<t<<'\t'<<vv<<'\t'<<k<<'\t'<<ttsw[k]<<'\t'<<-1<<'\t'<<1/stridetime<<endl;
		}
		for(int k=1;k<=4;k++) if(!swpre[k] && sw[k])	//start swing tracker
		{
			ttst[k]=t-ttstpre[k];
			ttswpre[k]=t;
			out_timers<<t<<'\t'<<vv<<'\t'<<k<<'\t'<<-1<<'\t'<<ttst[k]<<'\t'<<1/stridetime<<endl;
		}
		for(int k=1;k<=4;k++) swpre[k]=sw[k];


		//output data
		// if(int(t/dt)%int(DT/dt)==0)
		cout<<t<<'\t'<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<Gtot;
		for(int i=1;i<=4;i++) cout<<'\t'<<G[i];
		for(int i=1;i<=4;i++) cout<<'\t'<<(sw[i]?i:0);
		for(int i=1;i<=4;i++) cout<<'\t'<<phi[i];
		for(int i=1;i<=4;i++) cout<<'\t'<<Omega[i]+Omega_signal[i];
		cout<<endl;
		// for(int i=1;i<=4;i++) cout<<'\t'<<Omega[i];
		// cout<<endl;
		
		cerr<<"unset object"<<endl;
		cerr<<"unset arrow"<<endl;
		cerr<<"set xrange ["<<xc-20<<':'<<xc+20<<"]"<<endl;
		cerr<<"set yrange ["<<yc-20<<':'<<yc+20<<"]"<<endl;
		cerr<<"set title \"V="<<int(vv)<<'\"'<<endl;
		for(int i=1;i<=4;i++) {cerr<<"set object "<<i<<" circle front at "<<xfh[i]<<','<<yfh[i]<<" size "<<.2
								<<" fc \""<<colors[i]<<"\" fs "<<(sw[i] ? "empty" : "solid")<<endl;}

		cerr<<"set object 11 circle front at "<<xc<<','<<yc<<" size "<<.2<<" fc \"black\""<<endl;
		
		if(!sw[1] && !sw[2]) cerr<<"set arrow 1 from "<<xfl<<','<<yfl<<" to "<<xfr<<','<<yfr<<" nohead dt 2"<<endl;
		if(!sw[1] && !sw[3]) cerr<<"set arrow 2 from "<<xfl<<','<<yfl<<" to "<<xhl<<','<<yhl<<" nohead dt 2"<<endl;
		if(!sw[1] && !sw[4]) cerr<<"set arrow 3 from "<<xfl<<','<<yfl<<" to "<<xhr<<','<<yhr<<" nohead dt 2"<<endl;
		if(!sw[2] && !sw[3]) cerr<<"set arrow 4 from "<<xfr<<','<<yfr<<" to "<<xhl<<','<<yhl<<" nohead dt 2"<<endl;
		if(!sw[2] && !sw[4]) cerr<<"set arrow 5 from "<<xfr<<','<<yfr<<" to "<<xhr<<','<<yhr<<" nohead dt 2"<<endl;
		if(!sw[3] && !sw[4]) cerr<<"set arrow 6 from "<<xhl<<','<<yhl<<" to "<<xhr<<','<<yhr<<" nohead dt 2"<<endl;
		cerr<<"set arrow 7 from "<<xfl0<<','<<yfl0<<" to "<<xfl<<','<<yfl<<" nohead"<<endl;
		cerr<<"set arrow 8 from "<<xfr0<<','<<yfr0<<" to "<<xfr<<','<<yfr<<" nohead"<<endl;
		cerr<<"set arrow 9 from "<<xhl0<<','<<yhl0<<" to "<<xhl<<','<<yhl<<" nohead"<<endl;
		cerr<<"set arrow 10 from "<<xhr0<<','<<yhr0<<" to "<<xhr<<','<<yhr<<" nohead"<<endl;
		cerr<<"set arrow 11 from "<<xh<<','<<yh<<" to "<<xf<<','<<yf<<" lt -1 lw 5 nohead"<<endl;
		cerr<<"set arrow 12 from "<<xfl0<<','<<yfl0<<" to "<<xfr0<<','<<yfr0<<" nohead dt 3"<<endl;
		cerr<<"set arrow 13 from "<<xhl0<<','<<yhl0<<" to "<<xhr0<<','<<yhr0<<" nohead dt 3"<<endl;
		cerr<<"set arrow 14 from "<<xhl0<<','<<yhl0<<" to "<<xfl0<<','<<yfl0<<" dt 3"<<endl;
		cerr<<"set arrow 15 from "<<xhr0<<','<<yhr0<<" to "<<xfr0<<','<<yfr0<<" dt 3"<<endl;
		cerr<<"plot \"dat\" u 2:3 w d"<<endl;
	}

	if(newini) 
	{
		ofstream(inifiles[newini])<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<theta<<'\t'<<omega<<'\t'<<Gpre<<endl
		<<xfl<<'\t'<<xfr<<'\t'<<xhl<<'\t'<<xhr<<endl
		<<yfl<<'\t'<<yfr<<'\t'<<yhl<<'\t'<<yhr<<endl
		<<phi[1]<<'\t'<<phi[2]<<'\t'<<phi[3]<<'\t'<<phi[4]<<endl;
	}

	return 0;
}