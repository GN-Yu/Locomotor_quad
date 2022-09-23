#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>

using namespace std;

//units: length - cm, time - s, mass - g

double g=1000;	//garivtational accleration
double H=3;		//hight
double lf=1;	//half shoulder length
double lh=1;	//half hip length
double h=1;	//half width
double l=5;	//half body length
double L1=1*l, L2=0*l;		//front, hind part body length
double L=L1+L2;				//total body length
double D=sqrt(L*L+(lh+lf)*(lh+lf)); 	//diagonal length
double freq=sqrt(g/H);			//inverted pendulum frequency
double Tswc;				//typical swing time

// double m=30; // mouse mass in g
double I=L*L/3; // moment of inertia over mass
double tau=.25;

double DM[5]={0,1*L,1*L,1*L,1*L};	//max leg length
static int sw[5]={0,1,0,0,1};

int contraside[5]={0,2,1,4,3};
int homoside[5]={0,3,4,1,2};

char colors[5][10]={"black","brown","blue","red","green"};
char inifiles[4][50]={"","ini","slowini","fastini"};


int main(int argc,char** argv)
{
	//variable declarations
	double T, dt, DT;
	dt=.00001;
	DT=.01;
	int ini=0;
	int newini=0;
	int cor=0;
	int ani_out=0;
	int stridedata_out=0;
	int cookini=0;
	char fdur[50]="dur";
	char ini_name[50]="ini";

	double Guini,Gufin; 
	double Gvini,Gvfin;	
	double kvini,kvfin;
    double Tswcini,Tswcfin;
	double Fini,Ffin;
	double Gu, Gv, kv;	//thresholds
	double F0;		//forces in the unit cm/s, original force F/lambda
	
	
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

    double Zx, Zy;	//distance to diagonal when two limbs supporting
    double hc=H;    //height of the COM

	double d[5]={};	//the distance from a limb to its shoulder or hip

	int swing_count;		//number of swinging limbs
	double load[5]={};	//loads of each foot
    double F[5][5]={};
	double total_load;	//total load
	static double Gpre=0; //total supporting force of the previous time step
	static double GP[5]={};
	static double tswpre[5]={0,0,0,0,0};
	double tsw[5]={};


	double kr=0.5;	//damp of rotation
	int inhib=0;	//indicator of contra/same side inhibitions

	//initialize parameters
	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i],"-T")==0) T=atof(argv[++i]);
		else if(strcmp(argv[i],"-fps")==0) {ani_out=1;DT=1/atof(argv[++i]);}
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
		else if(strcmp(argv[i],"-kr")==0) kr=atof(argv[++i]);
		else if(strcmp(argv[i],"-inhib")==0) {inhib=1;}
		else if(strcmp(argv[i],"-cor")==0) {cor=1;}
		else if(strcmp(argv[i],"-i")==0) {ini=1;}
		else if(strcmp(argv[i],"-slowi")==0) {ini=2;}
		else if(strcmp(argv[i],"-fasti")==0) {ini=3;}
		else if(strcmp(argv[i],"-newi")==0) {newini=1;}
		else if(strcmp(argv[i],"-newslowi")==0) {newini=2;}
		else if(strcmp(argv[i],"-newfasti")==0) {newini=3;}
		else if(strcmp(argv[i],"-cookini")==0) {strcpy(inifiles[1],argv[++i]); newini=1;}
		//else if(strcmp(argv[i],"-out_strideinfo")) {strcpy(fdur,"strideinfo"); stridedata_out=1;}
		else if(strcmp(argv[i],"-out_swp")==0) {strcpy(fdur,argv[++i]); stridedata_out=1;}
		else if(strcmp(argv[i],"-swpini")==0) {strcpy(inifiles[1],argv[++i]); ini=1;}
		else return 1;
	}

	//output settings
	ofstream out_timers(fdur);
	// ofstream out_timers("strideinfo");

	if(ani_out)
	{
		cerr<<"set xtics 0,5,1000"<<endl;
		cerr<<"set ytics 0,5,2000"<<endl;
        cerr<<"set ztics 0,1,5"<<endl;
		cerr<<"set size ratio -1"<<endl;
		cerr<<"set grid"<<endl;
		cerr<<"set xlabel 'x'"<<endl;
		cerr<<"set ylabel 'y'"<<endl;
        cerr<<"set zlabel 'z'"<<endl;
        cerr<<"set xyplane 0"<<endl;
        cerr<<"set zrange [0:4]"<<endl;
        cerr<<"set view equal xyz"<<endl;
	}


	//initialize locomotion parameters
	F0=Fini;
	Gu=Guini;
	Gv=Gvini;
	kv=kvini;

	if(ini==0)
	{
		xc=0; yc=0;
		theta=M_PI/4;
		vx=2*F0*cos(theta); vy=2*F0*sin(theta);
		omega=0;

		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);
		xfl0=xf-lf*sin(theta), xfr0=xf+lf*sin(theta);
		xhl0=xh-lh*sin(theta), xhr0=xh+lh*sin(theta);
		yfl0=yf+lf*cos(theta), yfr0=yf-lf*cos(theta);
		yhl0=yh+lh*cos(theta), yhr0=yh-lh*cos(theta);

		double qq=L/sqrt(2)/2;
		xfl=xfl0, yfl=yfl0;
		xfr=xfr0-qq, yfr=yfr0-qq;
		xhl=xhl0-qq, yhl=yhl0-qq+.001;
		xhr=xhr0, yhr=yhr0;

		// double pp=L/sqrt(2)/2+0.2;
		// xfl=xfl0, yfl=yfl0;
		// xfr=xfr0-qq, yfr=yfr0-qq;
		// xhl=xhl0-pp, yhl=yhl0-pp;
		// xhr=xhr0, yhr=yhr0;

		// xfl=xfl0, yfl=yfl0;
		// xfr=xfr0-qq, yfr=yfr0-qq;
		// xhl=xhl0, yhl=yhl0;
		// xhr=xhr0-qq, yhr=yhr0-qq;

		// xfl=xfl0-3*cos(theta); yfl=yfl0-3*sin(theta);
		// xfr=xfr0; yfr=yfr0;
		// xhl=xhl0; yhl=yhl0;
		// xhr=xhr0-3*cos(theta); yhr=yhr0-3*sin(theta);
	}
	else
	{
		ifstream(inifiles[ini])>>xc>>yc>>vx>>vy>>theta>>omega>>Gpre
		>>xfl>>xfr>>xhl>>xhr>>yfl>>yfr>>yhl>>yhr
		>>sw[1]>>sw[2]>>sw[3]>>sw[4]
		>>tswpre[1]>>tswpre[2]>>tswpre[3]>>tswpre[4];
	
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
		Zx=0; Zy=0;
		total_load=0;
		swing_count=0;
		for(int i=1;i<=4;i++) load[i]=0;

		for(int i=1;i<=4;i++) if(sw[i]) swing_count++;	//swing legs count

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

        if(swing_count==0)
        {
            for(int k=1;k<=4;k++) load[k]=F[k][0];
            hc=H;
        }
        else if(swing_count==1)
        {
            for(int k=1;k<=4;k++) if(sw[k]) for(int i=1;i<=4;i++) load[i]=F[i][k];
            hc=H;
        }
        else if(swing_count==2)
        { 
            int k1=1; while(sw[k1]) k1++;
            double x1=xfh[k1],y1=yfh[k1];
            int k2=k1+1; while(sw[k2]) k2++;
            double x2=xfh[k2],y2=yfh[k2];
            double a=y1-y2,b=x2-x1,ab=hypot(a,b); a/=ab; b/=ab;
            double z=(xc-x2)*a+(yc-y2)*b;
            double vz=vx*a+vy*b;
            Zx=z*a; Zy=z*b;
            double R=sqrt(z*z+hc*hc);
            total_load=hc/R-vz*vz/R/g;
            double d12=hypot(x2-x1,y2-y1);
            double r2=((xc-x2)*(x1-x2)+(yc-y2)*(y1-y2))/d12;
            double r1=((xc-x1)*(x1-x2)+(yc-y1)*(y1-y2))/d12;
            load[k1]=r2/(r2-r1)*total_load;
            load[k2]=r1/(r1-r2)*total_load;
            hc-=z/hc*vz*dt;
        }
        else if(swing_count==3)
        {
            int k1=1; while(sw[k1]) k1++;
            double x1=xfh[k1],y1=yfh[k1];
            double a=xc-x1,b=yc-y1,ab=hypot(a,b); a/=ab; b/=ab;
            double z=(xc-x1)*a+(yc-y1)*b;
            double vz=vx*a+vy*b;
            Zx=-z*a; Zy=-z*b;
            double R=sqrt(z*z+hc*hc);
            total_load=hc/R-vz*vz/R/g;
            load[k1]=total_load;
            hc-=z/hc*vz*dt;
        }

		if(hc<0) 
		{
			cerr<<"fall"<<endl;
			break;
		}
		else if(total_load<0)
		{
			cerr<<"fly"<<endl;
			break;
		}

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

		double taus=.0001;
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
		double delta=0*omega;
		double FL=(1+delta)*F0, FR=(1-delta)*F0;
		
		if(!sw[1] && d[1]>0.) { Fx[1]=FL*(xfl0-xfl)/d[1]; Fy[1]=FL*(yfl0-yfl)/d[1]; }
		if(!sw[2] && d[2]>0.) { Fx[2]=FR*(xfr0-xfr)/d[2]; Fy[2]=FR*(yfr0-yfr)/d[2]; }
		if(!sw[3] && d[3]>0.) { Fx[3]=FL*(xhl0-xhl)/d[3]; Fy[3]=FL*(yhl0-yhl)/d[3]; }
		if(!sw[4] && d[4]>0.) { Fx[4]=FR*(xhr0-xhr)/d[4]; Fy[4]=FR*(yhr0-yhr)/d[4]; }

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

		static double psipre[5]={};
		double psi[5]={};
		double MA=0;
		for(int k=1;k<=4;k++)
		{
			double hh=hypot(yfh0[k]-yfh[k],xfh0[k]-xfh[k]);
			psi[k]=hh*sin(atan2(yfh0[k]-yfh[k],xfh0[k]-xfh[k])-theta);
			MA+=kr*(psi[k]-psipre[k])/dt; //+kr*kr/4/I*phi[k];
			psipre[k]=psi[k];
		}
		
		omega+=(M/I-omega+MA)*dt/tau;
		// omega+=(M/I-omega)*dt/tau;

		Gu=Guini+(Gufin-Guini)*t/T;
		Gv=Gvini+(Gvfin-Gvini)*t/T;
		kv=kvini+(kvfin-kvini)*t/T;
        Tswc=Tswcini+(Tswcfin-Tswcini)*t/T;

		
		for(int k=1;k<=4;k++) tsw[k]=0;

		for(int k=1;k<=4;k++) if(sw[k]) if(t-tswpre[k]>=Tswc) 
		{
            if(swing_count==2)
			{
				int s=1;
				while(!sw[s] || s==k) {s++;}
				if(F[k][s]<0.) continue;
			}
			else if(swing_count==1 && F[k][0]<0.) continue;
			sw[k]=0;
			swing_count--;
		}	//stop swing a leg after specific time, and skip if it gets negative load after landing

		for(int k=1;k<=4;k++) if(!sw[k]) if(GP[k]>load[k] && load[k]<Gu)
		{
			if(inhib) if(sw[contraside[k]]==1 || sw[homoside[k]]==1) {continue;}
			sw[k]=1;
			swing_count++;
			tswpre[k]=t;
		}	//weak lifting conditions

		for(int k=1;k<=4;k++) if(!sw[k]) if(load[k]<0 || d[k]>DM[k])
		{
			sw[k]=1;
			swing_count++;
            tswpre[k]=t;
			if(inhib)
			{
				if(sw[contraside[k]]==1) {sw[contraside[k]]=0; swing_count--;}
				if(sw[homoside[k]]==1) {sw[homoside[k]]=0; swing_count--;}
			}
		}	//strong lifting conditions
        
		//if((total_load!=0 && -(total_load-Gpre)/dt>kv) || swing_count==4)
		if(total_load!=0 && -(total_load-Gpre)/dt>kv)
		{
			int kmax=0;
			for(int k=1;k<=4;k++) if(sw[k]) {tsw[k]=t-tswpre[k];}
			for(int k=1;k<=4;k++) if(tsw[k]>tsw[kmax]) {kmax=k;}
            if(swing_count==2)
			{
				int s=1; while(!sw[s] || s==kmax) s++;
				if(F[kmax][s]<0.) continue;
			}
			else if(swing_count==1 && F[kmax][0]<0.) continue;
			sw[kmax]=0;
			swing_count--;
		}	//stop swing a leg when lose balance


		for(int k=1;k<=4;k++) GP[k]=load[k];
		Gpre=total_load;

	
		//output data
		if(ani_out && int(t/dt)%int(DT/dt)==0)
		{
			cout<<t<<'\t'<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<total_load;
			for(int i=1;i<=4;i++) cout<<'\t'<<load[i];
			for(int i=1;i<=4;i++) cout<<'\t'<<(sw[i]?i:0);
			// cout<<'\t'<<deltal<<'\t'<<deltar;
			cout<<endl;
			
			cerr<<"unset object"<<endl;
			cerr<<"unset arrow"<<endl;
			cerr<<"set xrange ["<<xc-20<<':'<<xc+20<<"]"<<endl;
			cerr<<"set yrange ["<<yc-20<<':'<<yc+20<<"]"<<endl;
			cerr<<"set title \"V="<<int(vv)<<'\"'<<endl;
			for(int i=1;i<=4;i++) {cerr<<"set object "<<i<<" circle front at "<<xfh[i]<<','<<yfh[i]<<",0 size "<<.2
									<<" fc \""<<colors[i]<<"\" fs "<<(sw[i] ? "empty" : "solid")<<endl;}

			cerr<<"set object 11 circle front at "<<xc<<','<<yc<<','<<hc<<" size "<<.2<<" fc \"black\""<<endl;
			
			if(!sw[1] && !sw[2]) cerr<<"set arrow 1 from "<<xfl<<','<<yfl<<",0 to "<<xfr<<','<<yfr<<",0 nohead dt 2"<<endl;
			if(!sw[1] && !sw[3]) cerr<<"set arrow 2 from "<<xfl<<','<<yfl<<",0 to "<<xhl<<','<<yhl<<",0 nohead dt 2"<<endl;
			if(!sw[1] && !sw[4]) cerr<<"set arrow 3 from "<<xfl<<','<<yfl<<",0 to "<<xhr<<','<<yhr<<",0 nohead dt 2"<<endl;
			if(!sw[2] && !sw[3]) cerr<<"set arrow 4 from "<<xfr<<','<<yfr<<",0 to "<<xhl<<','<<yhl<<",0 nohead dt 2"<<endl;
			if(!sw[2] && !sw[4]) cerr<<"set arrow 5 from "<<xfr<<','<<yfr<<",0 to "<<xhr<<','<<yhr<<",0 nohead dt 2"<<endl;
			if(!sw[3] && !sw[4]) cerr<<"set arrow 6 from "<<xhl<<','<<yhl<<",0 to "<<xhr<<','<<yhr<<",0 nohead dt 2"<<endl;
			cerr<<"set arrow 7 from "<<xfl0<<','<<yfl0<<','<<hc<<" to "<<xfl<<','<<yfl<<",0 nohead"<<endl;
			cerr<<"set arrow 8 from "<<xfr0<<','<<yfr0<<','<<hc<<" to "<<xfr<<','<<yfr<<",0 nohead"<<endl;
			cerr<<"set arrow 9 from "<<xhl0<<','<<yhl0<<','<<hc<<" to "<<xhl<<','<<yhl<<",0 nohead"<<endl;
			cerr<<"set arrow 10 from "<<xhr0<<','<<yhr0<<','<<hc<<" to "<<xhr<<','<<yhr<<",0 nohead"<<endl;
			cerr<<"set arrow 11 from "<<xh<<','<<yh<<','<<hc<<" to "<<xf<<','<<yf<<','<<hc<<" lt -1 lw 5 nohead"<<endl;
			cerr<<"set arrow 12 from "<<xfl0<<','<<yfl0<<','<<hc<<" to "<<xfr0<<','<<yfr0<<','<<hc<<" nohead dt 3"<<endl;
			cerr<<"set arrow 13 from "<<xhl0<<','<<yhl0<<','<<hc<<" to "<<xhr0<<','<<yhr0<<','<<hc<<" nohead dt 3"<<endl;
			cerr<<"set arrow 14 from "<<xhl0<<','<<yhl0<<','<<hc<<" to "<<xfl0<<','<<yfl0<<','<<hc<<" dt 3"<<endl;
			cerr<<"set arrow 15 from "<<xhr0<<','<<yhr0<<','<<hc<<" to "<<xfr0<<','<<yfr0<<','<<hc<<" dt 3"<<endl;
			cerr<<"splot \"dat\" u 2:3:(0) w d"<<endl;
		}

		//falling indicator: if there is at least one leg having negative load for a long time (say 0.2s), it is falling
		// static double fallpre=0;
		// int iter_nonneg;
		// iter_nonneg=1;
		// while(load[iter_nonneg]>=0 && iter_nonneg<4) iter_nonneg++;
		// if(iter_nonneg<4) {if(t-fallpre>0.2) break;}
		// else fallpre=t;

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
			//if(stridedata_out==1 && stridepre!=0) out_timers<<t<<'\t'<<vv<<'\t'<<k<<'\t'<<ttsw[k]<<'\t'<<-1<<'\t'<<1/stridetime<<endl;
		}
		for(int k=1;k<=4;k++) if(!swpre[k] && sw[k])	//start swing tracker
		{
			ttst[k]=t-ttstpre[k];
			ttswpre[k]=t;
			//if(stridepre!=0) out_timers<<Tswc<<'\t'<<F0<<'\t'<<k<<'\t'<<-1<<'\t'<<ttst[k]<<'\t'<<(100*ttst[k]/stridetime)<<endl;
			if(stridepre!=0) out_timers<<t<<'\t'<<vv<<'\t'<<k<<'\t'<<-1<<'\t'<<ttst[k]<<'\t'<<1/stridetime<<endl;
		}
		for(int k=1;k<=4;k++) swpre[k]=sw[k];
	}
	
	if(newini)
	{
		for(int k=1;k<=4;k++) tswpre[k]-=T;
		ofstream(inifiles[newini])<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<theta<<'\t'<<omega<<'\t'<<Gpre<<'\t'
		<<xfl<<'\t'<<xfr<<'\t'<<xhl<<'\t'<<xhr<<'\t'
		<<yfl<<'\t'<<yfr<<'\t'<<yhl<<'\t'<<yhr<<'\t'
		<<sw[1]<<'\t'<<sw[2]<<'\t'<<sw[3]<<'\t'<<sw[4]<<'\t'
		<<tswpre[1]<<'\t'<<tswpre[2]<<'\t'<<tswpre[3]<<'\t'<<tswpre[4]<<endl;
	}
	
	return 0;
}