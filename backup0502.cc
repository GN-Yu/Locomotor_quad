#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>

using namespace std;

double g=1000; // garivtational accleration in cm/s/s
double H=3; // mouse hight in cm
double l=4.1; // mouse "half"-length in cm
double h=1.5; // mouse half-width in cm
double L1=1*l, L2=0*l;
double L=L1+L2; //	mouse body length
double D=sqrt(L*L+4*h*h); // diagonal length
double w=sqrt(g/H); // frequency
double Tswc;	//typical swing time

// double m=30; // mouse mass in g
double I=L*L/3; // moment of inertia over mass
double tau=.25;

static int sw[5]={0,1,0,0,0};

int main(int argc,char** argv)
{
	ofstream abc("strideinfo");

	double T, dt;
	int ini=0;
	int newini=0;
	int cor=0;

	double Guini, Gufin; 
	double Gvini, Gvfin;	
	double Gu, Gv;	//thresholds
	double Fini, Ffin;
	double F0;		//forces
	double xc, yc;
	double vx, vy;
	double omega;	//angular velocity
	double theta;	//body angle
	double xf, yf, xh, yh;	
	double xfl0, yfl0, xfr0, yfr0, xhl0, yhl0, xhr0, yhr0;	//shoulder and hip positions
	double xfl, xfr, xhl, xhr;
	double yfl, yfr, yhl, yhr;	//limb positions
	double vv, vexp;			//speed in magnitude, expected speed
	int nsw;		//number of swinging limbs
	double G[5]={};	//loads of each foot
	double Gtot;	//total load
	double Zx, Zy;	//distance to diagonal when two limbs supporting
	int flcont=0;
	int flcontpre=0;
	static double GP[5]={};
	double tsw[5]={};

	double kv;
	double kvini;
	double kvfin;
	
	//initialize parameters
	theta=M_PI/4;

	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i],"-T")==0) T=atof(argv[++i]);
		else if(strcmp(argv[i],"-dt")==0) dt=atof(argv[++i]);
		//else if(strcmp(argv[i],"-v")==0) F0=atof(argv[++i]);
		else if(strcmp(argv[i],"-v0")==0) Fini=atof(argv[++i]);
		else if(strcmp(argv[i],"-v1")==0) Ffin=atof(argv[++i]);
		else if(strcmp(argv[i],"-Gu0")==0) Guini=atof(argv[++i]);
		else if(strcmp(argv[i],"-Gu1")==0) Gufin=atof(argv[++i]);
		else if(strcmp(argv[i],"-Gv0")==0) Gvini=atof(argv[++i]);
		else if(strcmp(argv[i],"-Gv1")==0) Gvfin=atof(argv[++i]);
		else if(strcmp(argv[i],"-kv0")==0) kvini=atof(argv[++i]);
		else if(strcmp(argv[i],"-kv1")==0) kvfin=atof(argv[++i]);
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


	//initialize locomotion parameters
	if(ini==0)
	{
		xc=0; yc=0;
		vx=3*F0*cos(theta); vy=3*F0*sin(theta);
		omega=0;

		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);
		xfl0=xf-h*sin(theta), xfr0=xf+h*sin(theta);
		xhl0=xh-h*sin(theta), xhr0=xh+h*sin(theta);
		yfl0=yf+h*cos(theta), yfr0=yf-h*cos(theta);
		yhl0=yh+h*cos(theta), yhr0=yh-h*cos(theta);

		xfl=xfl0; yfl=yfl0;
		xfr=xfr0; yfr=yfr0;
		xhl=xhl0-.1; yhl=yhl0-.1;
		xhr=xhr0; yhr=yhr0;
	}
	else if(ini==1)
	{
		ifstream("ini")>>xc>>yc>>vx>>vy>>theta>>omega
		>>xfl>>xfr>>xhl>>xhr
		>>yfl>>yfr>>yhl>>yhr
		>>sw[1]>>sw[2]>>sw[3]>>sw[4];

		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);
		xfl0=xf-h*sin(theta), xfr0=xf+h*sin(theta);
		xhl0=xh-h*sin(theta), xhr0=xh+h*sin(theta);
		yfl0=yf+h*cos(theta), yfr0=yf-h*cos(theta);
		yhl0=yh+h*cos(theta), yhr0=yh-h*cos(theta);
	}
	else if(ini==2)
	{
		ifstream("slowini")>>xc>>yc>>vx>>vy>>theta>>omega
		>>xfl>>xfr>>xhl>>xhr
		>>yfl>>yfr>>yhl>>yhr
		>>sw[1]>>sw[2]>>sw[3]>>sw[4];

		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);
		xfl0=xf-h*sin(theta), xfr0=xf+h*sin(theta);
		xhl0=xh-h*sin(theta), xhr0=xh+h*sin(theta);
		yfl0=yf+h*cos(theta), yfr0=yf-h*cos(theta);
		yhl0=yh+h*cos(theta), yhr0=yh-h*cos(theta);
	}
	else if(ini==3)
	{
		ifstream("fastini")>>xc>>yc>>vx>>vy>>theta>>omega
		>>xfl>>xfr>>xhl>>xhr
		>>yfl>>yfr>>yhl>>yhr
		>>sw[1]>>sw[2]>>sw[3]>>sw[4];

		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);
		xfl0=xf-h*sin(theta), xfr0=xf+h*sin(theta);
		xhl0=xh-h*sin(theta), xhr0=xh+h*sin(theta);
		yfl0=yf+h*cos(theta), yfr0=yf-h*cos(theta);
		yhl0=yh+h*cos(theta), yhr0=yh-h*cos(theta);
	}
	vv=sqrt(vx*vx+vy*vy);
	vexp=vv;
	Tswc=0.1;


	cerr<<"set xtics 0,5,1000"<<endl;
	cerr<<"set ytics 0,5,1000"<<endl;
	cerr<<"set size ratio -1"<<endl;
	

	for(double t=0;t<=T;t+=dt)	//step cycles
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

		for(int i=1;i<=4;i++) if(sw[i]) nsw++;	//swing legs count

		if(nsw==0) {for(int k=1;k<=4;k++) G[k]=F[k][0];}
		else if(nsw==1) {for(int k=1;k<=4;k++) if(sw[k]) for(int i=1;i<=4;i++) G[i]=F[i][k];}
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
		else if(nsw==4){
			cerr<<"flying"<<'\t'<<flcont<<'\t'<<Gu<<endl; 
			cerr<<G[1]<<'\t'<<G[2]<<'\t'<<G[3]<<'\t'<<G[4]<<endl; 
			cerr<<GP[1]<<'\t'<<GP[2]<<'\t'<<GP[3]<<'\t'<<GP[4]<<endl;
			cerr<<tsw[1]<<'\t'<<tsw[2]<<'\t'<<tsw[3]<<'\t'<<tsw[4]<<endl;
			break;}
		
		
		double d[5]={0,
			hypot(xfl0-xfl,yfl0-yfl),hypot(xfr0-xfr,yfr0-yfr),
			hypot(xhl0-xhl,yhl0-yhl),hypot(xhr0-xhr,yhr0-yhr)};	//the distance from a limb to its shoulder or hip


		// correction of limb positioning when swing

		double vtor=vx*(yfl0-yhr0)/D - vy*(xfl0-xhr0)/D;
		double vtol=-vx*(yfr0-yhl0)/D + vy*(xfr0-xhl0)/D;
		double deltar=(L/2-L2)-vtor*D/2/h/w+2*vexp/w/(1+exp(w*L/vexp));
		double deltal=(L/2-L2)-vtol*D/2/h/w+2*vexp/w/(1+exp(w*L/vexp));

		double xfl1,yfl1,xfr1,yfr1,xhl1,yhl1,xhr1,yhr1;
		xfl1=xfl0-deltar*cos(theta); yfl1=yfl0-deltar*sin(theta);
		xhr1=xhr0-deltar*cos(theta); yhr1=yhr0-deltar*sin(theta);
		xfr1=xfr0-deltal*cos(theta); yfr1=yfr0-deltal*sin(theta);
		xhl1=xhl0-deltal*cos(theta); yhl1=yhl0-deltal*sin(theta);

		double taus=.001;
		double Vswc=(L+vv*Tswc)/Tswc;

		// if(cor)
		// {
		// 	if(sw[1]) { xfl+=(xfl1-xfl)*dt/taus; yfl+=(yfl1-yfl)*dt/taus; }
		// 	if(sw[2]) { xfr+=(xfr1-xfr)*dt/taus; yfr+=(yfr1-yfr)*dt/taus; }
		// 	if(sw[3]) { xhl+=(xhl1-xhl)*dt/taus; yhl+=(yhl1-yhl)*dt/taus; }
		// 	if(sw[4]) { xhr+=(xhr1-xhr)*dt/taus; yhr+=(yhr1-yhr)*dt/taus; }
		// }
		// else 
		// {
		// 	if(d[1]>0.01 && sw[1]) { xfl+=(xfl0-xfl)*dt*Vswc/d[1]; yfl+=(yfl0-yfl)*dt*Vswc/d[1]; }
		// 	if(d[2]>0.01 && sw[2]) { xfr+=(xfr0-xfr)*dt*Vswc/d[2]; yfr+=(yfr0-yfr)*dt*Vswc/d[2]; }
		// 	if(d[3]>0.01 && sw[3]) { xhl+=(xhl0-xhl)*dt*Vswc/d[3]; yhl+=(yhl0-yhl)*dt*Vswc/d[3]; }
		// 	if(d[4]>0.01 && sw[4]) { xhr+=(xhr0-xhr)*dt*Vswc/d[4]; yhr+=(yhr0-yhr)*dt*Vswc/d[4]; }
		// }

		if(cor)
		{
			if(sw[1]) { xfl+=(xfl1-xfl)*dt/taus; yfl+=(yfl1-yfl)*dt/taus; }
			if(sw[2]) { xfr+=(xfr1-xfr)*dt/taus; yfr+=(yfr1-yfr)*dt/taus; }
			if(sw[3]) { xhl+=(xhl1-xhl)*dt/taus; yhl+=(yhl1-yhl)*dt/taus; }
			if(sw[4]) { xhr+=(xhr1-xhr)*dt/taus; yhr+=(yhr1-yhr)*dt/taus; }
		}
		else 
		{
			if(sw[1]) { xfl+=(xfl0-xfl)*dt/taus; yfl+=(yfl0-yfl)*dt/taus; }
			if(sw[2]) { xfr+=(xfr0-xfr)*dt/taus; yfr+=(yfr0-yfr)*dt/taus; }
			if(sw[3]) { xhl+=(xhl0-xhl)*dt/taus; yhl+=(yhl0-yhl)*dt/taus; }
			if(sw[4]) { xhr+=(xhr0-xhr)*dt/taus; yhr+=(yhr0-yhr)*dt/taus; }
		}
		

		double Fx[5]={},Fy[5]={};
		double delta=0.1*omega;

		F0=Fini+(Ffin-Fini)*t/T;	//F0 change

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
		vv=sqrt(vx*vx+vy*vy);

		double tauvv=1;
		vexp+=(vv-vexp)*dt/tauvv;
		
		double M[5]={0,
			(xfl-xc)*Fy[1]-(yfl-yc)*Fx[1],
			(xfr-xc)*Fy[2]-(yfr-yc)*Fx[2],
			(xhl-xc)*Fy[3]-(yhl-yc)*Fx[3],
			(xhr-xc)*Fy[4]-(yhr-yc)*Fx[4]};
		
		omega+=((M[1]+M[2]+M[3]+M[4])/I-omega)*dt/tau;
		
		xf=xc+L1*cos(theta); yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta); yh=yc-L2*sin(theta);

		xfl0=xf-h*sin(theta); xfr0=xf+h*sin(theta);
		xhl0=xh-h*sin(theta); xhr0=xh+h*sin(theta);
		yfl0=yf+h*cos(theta); yfr0=yf-h*cos(theta);
		yhl0=yh+h*cos(theta); yhr0=yh-h*cos(theta);
		

		cout<<t<<'\t'<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<Gtot;
		for(int k=1;k<=4;k++) cout<<'\t'<<G[k];
		for(int i=1;i<=4;i++) cout<<'\t'<<(sw[i]?i:0);
		cout<<'\t'<<deltal<<'\t'<<deltar;
		cout<<endl;
		
		cerr<<"unset object"<<endl;
		cerr<<"set xrange ["<<xc-25<<':'<<xc+25<<"]"<<endl;
		cerr<<"set yrange ["<<yc-25<<':'<<yc+25<<"]"<<endl;
		cerr<<"set object 1 circle front at "<<xfl<<','<<yfl<<" size "<<.2<<" fc \"brown\""<<endl;
		cerr<<"set object 2 circle front at "<<xfr<<','<<yfr<<" size "<<.2<<" fc \"blue\""<<endl;
		cerr<<"set object 3 circle front at "<<xhl<<','<<yhl<<" size "<<.2<<" fc \"red\""<<endl;
		cerr<<"set object 4 circle front at "<<xhr<<','<<yhr<<" size "<<.2<<" fc \"green\""<<endl;
		cerr<<"set object 5 circle front at "<<xc<<','<<yc<<" size "<<.2<<" fc \"black\""<<endl;
		cerr<<"unset arrow"<<endl;
		if(!sw[1] && !sw[2]) cerr<<"set arrow 1 from "<<xfl<<','<<yfl<<" to "<<xfr<<','<<yfr<<" nohead dt 2"<<endl;
		if(!sw[1] && !sw[3]) cerr<<"set arrow 2 from "<<xfl<<','<<yfl<<" to "<<xhl<<','<<yhl<<" nohead dt 2"<<endl;
		if(!sw[1] && !sw[4]) cerr<<"set arrow 3 from "<<xfl<<','<<yfl<<" to "<<xhr<<','<<yhr<<" nohead dt 2"<<endl;
		if(!sw[2] && !sw[3]) cerr<<"set arrow 4 from "<<xfr<<','<<yfr<<" to "<<xhl<<','<<yhl<<" nohead dt 2"<<endl;
		if(!sw[2] && !sw[4]) cerr<<"set arrow 5 from "<<xfr<<','<<yfr<<" to "<<xhr<<','<<yhr<<" nohead dt 2"<<endl;
		if(!sw[3] && !sw[4]) cerr<<"set arrow 6 from "<<xhl<<','<<yhl<<" to "<<xhr<<','<<yhr<<" nohead dt 2"<<endl;
		cerr<<"set arrow 7 from "<<xf<<','<<yf<<" to "<<xfl<<','<<yfl<<" nohead"<<endl;
		cerr<<"set arrow 8 from "<<xf<<','<<yf<<" to "<<xfr<<','<<yfr<<" nohead"<<endl;
		cerr<<"set arrow 9 from "<<xh<<','<<yh<<" to "<<xhl<<','<<yhl<<" nohead"<<endl;
		cerr<<"set arrow 10 from "<<xh<<','<<yh<<" to "<<xhr<<','<<yhr<<" nohead"<<endl;
		cerr<<"set arrow 11 from "<<xh<<','<<yh<<" to "<<xf<<','<<yf<<" lt -1 lw 5 nohead"<<endl;
		cerr<<"plot \"dat\" u 2:3 w d"<<endl;
		
		
		static double tswpre[5]={};
		for(int k=1;k<=4;k++) tsw[k]=0;
		double P=1, DM[5]={0,L,L,1.1*L,1.1*L};
		flcont=0;
		static double GS=0;	//total supporting force
		double Gpre=GS;	//total supporting force of the previous time step
		
		//GS=0; for(int i=1;i<=4;i++) GS+=G[i];
		GS=Gtot;

		Gu=Guini+(Gufin-Guini)*t/T;
		Gv=Gvini+(Gvfin-Gvini)*t/T;
		kv=kvini+(kvfin-kvini)*t/T;

		if(nsw<2) for(int k=1;k<=4;k++) if(!sw[k]) if((k>2 && GP[k]>G[k] && G[k]<Gu) || d[k]>DM[k])
		{
			sw[k]=1;
			if(k<=2) {sw[k+2]=0;}
			tswpre[k]=t;
			flcont++;
		}	//start swinging for nonsupport or overextension

		
		for(int k=1;k<=4;k++) GP[k]=G[k];

		// for(int k=1;k<=4;k++) if(sw[k] && d[k]<0.01) sw[k]=0;
		
		if(flcontpre==0 && ((Gtot!=0 && -(Gtot-Gpre)/dt>kv) || nsw>2))
		{
			int kmax=0;
			for(int k=1;k<=4;k++) if(sw[k]) {tsw[k]=t-tswpre[k];}
			for(int k=1;k<=4;k++) if(tsw[k]>tsw[kmax]) {kmax=k;}
			sw[kmax]=0;
		}	//stop swing a leg for lossing balance

		flcontpre=flcont;

		// for(int k=1;k<=4;k++) if(sw[k] && (-(Gtot-Gpre)/dt>kv || nsw==3)) tsw[k]=t-tswpre[k];
		// //for(int k=1;k<=4;k++) if(sw[k] && ((Gpre<Gv && GS<Gpre) || nsw==3)) tsw[k]=t-tswpre[k];
		// int kmax=0;
		// for(int k=1;k<=4;k++) if(tsw[k]>tsw[kmax]) kmax=k;
		// sw[kmax]=0;

		//for(int k=1;k<=2;k++) if(sw[k]) sw[k+2]=0;

		
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
			if(stridepre!=0) abc<<t<<'\t'<<vv<<'\t'<<k<<'\t'<<ttsw[k]<<'\t'<<-1<<'\t'<<1/stridetime<<endl;
		}
		for(int k=1;k<=4;k++) if(!swpre[k] && sw[k])	//start swing tracker
		{
			ttst[k]=t-ttstpre[k];
			ttswpre[k]=t;
			abc<<t<<'\t'<<vv<<'\t'<<k<<'\t'<<-1<<'\t'<<ttst[k]<<'\t'<<1/stridetime<<endl;
		}
		for(int k=1;k<=4;k++) swpre[k]=sw[k];
	}

	if(newini==1)
	{	
		ofstream("ini")<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<theta<<'\t'<<omega<<endl
		<<xfl<<'\t'<<xfr<<'\t'<<xhl<<'\t'<<xhr<<endl
		<<yfl<<'\t'<<yfr<<'\t'<<yhl<<'\t'<<yhr<<endl
		<<sw[1]<<'\t'<<sw[2]<<'\t'<<sw[3]<<'\t'<<sw[4]<<endl;
	}
	else if(newini==2)
	{	
		ofstream("slowini")<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<theta<<'\t'<<omega<<endl
		<<xfl<<'\t'<<xfr<<'\t'<<xhl<<'\t'<<xhr<<endl
		<<yfl<<'\t'<<yfr<<'\t'<<yhl<<'\t'<<yhr<<endl
		<<sw[1]<<'\t'<<sw[2]<<'\t'<<sw[3]<<'\t'<<sw[4]<<endl;
	}
	else if(newini==3)
	{	
		ofstream("fastini")<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<theta<<'\t'<<omega<<endl
		<<xfl<<'\t'<<xfr<<'\t'<<xhl<<'\t'<<xhr<<endl
		<<yfl<<'\t'<<yfr<<'\t'<<yhl<<'\t'<<yhr<<endl
		<<sw[1]<<'\t'<<sw[2]<<'\t'<<sw[3]<<'\t'<<sw[4]<<endl;
	}
	
	return 0;
}