#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>

using namespace std;

// plot flags
bool plot_parameter_sweep = false;
bool plot_specified_velocity = false;

// quadruped body geometry
// units: length - cm, time - s, mass - g
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
static int sw[5] = {0,1,0,0,0};

int contralateral[5] = {0,2,1,4,3};
int ipsilateral[5]   = {0,3,4,1,2};
int diagnal[5]       = {0,4,3,2,1};

double balance_sensory[5]={0,1,1,1,1};

char colors[5][10]   = {"black","brown","blue","red","green"};
char inifiles[4][50] = {"","ini","slowini","fastini"};

double rnd() { return ((double) rand() / (RAND_MAX)) * 2 - 1; }


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

int positive_load_check(int k,int swing_count,int sw[],double F[][5]){
	//check if limb k has positive load after landing
	int k_pos=1;
	if(swing_count == 2){
		int s=1;
		while(!sw[s] || s==k) s++;
		if(F[k][s]<0.) k_pos=0;
	}
	else if(swing_count==1 && F[k][0]<0.) k_pos=0;

	return k_pos;
}

void update_corrected_swing_target(double theta,double vx,double vy,double vexp,double xfh0[],double yfh0[],double xfh1[],double yfh1[]){
	//calculate corrections of limb positionings when swing
	double vtor=vx*(yfh0[1]-yfh0[4])/D - vy*(xfh0[1]-xfh0[4])/D;
	double vtol=-vx*(yfh0[2]-yfh0[3])/D + vy*(xfh0[2]-xfh0[3])/D;
	double deltar=(L/2-L2)-vtor*D/2/h/freq+2*vexp/freq/(1+exp(freq*L/vexp));
	double deltal=(L/2-L2)-vtol*D/2/h/freq+2*vexp/freq/(1+exp(freq*L/vexp));

	xfh1[1]=xfh0[1]-deltar*cos(theta); yfh1[1]=yfh0[1]-deltar*sin(theta);
	xfh1[4]=xfh0[4]-deltar*cos(theta); yfh1[4]=yfh0[4]-deltar*sin(theta);
	xfh1[2]=xfh0[2]-deltal*cos(theta); yfh1[2]=yfh0[2]-deltal*sin(theta);
	xfh1[3]=xfh0[3]-deltal*cos(theta); yfh1[3]=yfh0[3]-deltal*sin(theta);
}

void plotcmd_head() {
	cerr<<"set terminal qt size 800,800"<<endl;
	cerr<<"set xlabel ''"<<endl;
	cerr<<"set ylabel ''"<<endl;
	cerr<<"set zlabel ''"<<endl;
	cerr<<"set key noautotitle"<<endl;
	cerr<<"set xtics -50,5,1000"<<endl;
	cerr<<"set ytics -50,5,2000"<<endl;
	cerr<<"set ztics 0,1,5"<<endl;
	cerr<<"set size ratio -1"<<endl;
	cerr<<"set grid"<<endl;
	cerr<<"set xyplane 0"<<endl;
	cerr<<"set zrange [0:4]"<<endl;
	cerr<<"set view 0,0,1,1"<<endl;
	//cerr<<"set view equal xyz"<<endl;
}

void plotcmd_frame(int sw[],double xc,double yc,double hc,double xf,double yf,double xh,double yh,double xfh[],double yfh[],double xfh0[],double yfh0[]){
	cerr<<"unset object"<<endl;
	cerr<<"unset arrow"<<endl;
	cerr<<"set xrange ["<<xc-20<<':'<<xc+20<<"]"<<endl;
	cerr<<"set yrange ["<<yc-20<<':'<<yc+20<<"]"<<endl;
	//cerr<<"set title \"V="<<int(vv)<<'\"'<<endl;
	for(int i=1;i<=4;i++) {cerr<<"set object "<<i<<" circle front at "<<xfh[i]<<','<<yfh[i]<<",0 size "<<.2
							<<" fc \""<<colors[i]<<"\" fs "<<(sw[i] ? "empty" : "solid")<<endl;}

	cerr<<"set object 11 circle front at "<<xc<<','<<yc<<','<<hc<<" size "<<.2<<" fc \"black\""<<endl;
	
	if(!sw[1] && !sw[2]) cerr<<"set arrow 1 from "<<xfh[1]<<','<<yfh[1]<<",0 to "<<xfh[2]<<','<<yfh[2]<<",0 nohead dt 2"<<endl;
	if(!sw[1] && !sw[3]) cerr<<"set arrow 2 from "<<xfh[1]<<','<<yfh[1]<<",0 to "<<xfh[3]<<','<<yfh[3]<<",0 nohead dt 2"<<endl;
	if(!sw[1] && !sw[4]) cerr<<"set arrow 3 from "<<xfh[1]<<','<<yfh[1]<<",0 to "<<xfh[4]<<','<<yfh[4]<<",0 nohead dt 2"<<endl;
	if(!sw[2] && !sw[3]) cerr<<"set arrow 4 from "<<xfh[2]<<','<<yfh[2]<<",0 to "<<xfh[3]<<','<<yfh[3]<<",0 nohead dt 2"<<endl;
	if(!sw[2] && !sw[4]) cerr<<"set arrow 5 from "<<xfh[2]<<','<<yfh[2]<<",0 to "<<xfh[4]<<','<<yfh[4]<<",0 nohead dt 2"<<endl;
	if(!sw[3] && !sw[4]) cerr<<"set arrow 6 from "<<xfh[3]<<','<<yfh[3]<<",0 to "<<xfh[4]<<','<<yfh[4]<<",0 nohead dt 2"<<endl;
	cerr<<"set arrow 7 from "<<xfh0[1]<<','<<yfh0[1]<<','<<hc<<" to "<<xfh[1]<<','<<yfh[1]<<",0 nohead"<<endl;
	cerr<<"set arrow 8 from "<<xfh0[2]<<','<<yfh0[2]<<','<<hc<<" to "<<xfh[2]<<','<<yfh[2]<<",0 nohead"<<endl;
	cerr<<"set arrow 9 from "<<xfh0[3]<<','<<yfh0[3]<<','<<hc<<" to "<<xfh[3]<<','<<yfh[3]<<",0 nohead"<<endl;
	cerr<<"set arrow 10 from "<<xfh0[4]<<','<<yfh0[4]<<','<<hc<<" to "<<xfh[4]<<','<<yfh[4]<<",0 nohead"<<endl;
	cerr<<"set arrow 11 from "<<xh<<','<<yh<<','<<hc<<" to "<<xf<<','<<yf<<','<<hc<<" lt -1 lw 5 nohead"<<endl;
	cerr<<"set arrow 12 from "<<xfh0[1]<<','<<yfh0[1]<<','<<hc<<" to "<<xfh0[2]<<','<<yfh0[2]<<','<<hc<<" nohead dt 3"<<endl;
	cerr<<"set arrow 13 from "<<xfh0[3]<<','<<yfh0[3]<<','<<hc<<" to "<<xfh0[4]<<','<<yfh0[4]<<','<<hc<<" nohead dt 3"<<endl;
	//cerr<<"set arrow 14 from "<<xfh0[3]<<','<<yfh0[3]<<','<<hc<<" to "<<xfh0[1]<<','<<yfh0[1]<<','<<hc<<" dt 3"<<endl;
	//cerr<<"set arrow 15 from "<<xfh0[4]<<','<<yfh0[4]<<','<<hc<<" to "<<xfh0[2]<<','<<yfh0[2]<<','<<hc<<" dt 3"<<endl;
	cerr<<"splot \"dat\" u 2:3:(0) w d"<<endl;
}


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
	int iniWALK=0;
	int iniTROT=0;
	int iniPACE=0;
	int iniRANDOM=0;
	char fdur[50]="dur";
	char ini_name[50]="ini";
	int nicepic=0;
	double start_out_time=1;
	srand (time(NULL));

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
	double load[5]={};		//loads of each foot
    double F[5][5]={};
	double total_load;		//total load
	static double Gpre=0; 	//total supporting force of the previous time step
	static double GP[5]={};
	static double tswpre[5]={0,0,0,0,0};
	double tsw[5]={};

	double kr=0.5;	//damp of rotation

	int swing_inhib_ipsi=0;	//indicator of ipsilateral inhibitions
	int swing_inhib_contra=0; //indicator of contralateral inhibitions
	int diff_kv=0;	//an interpretation of diagonal swing excitation using different balance parameters

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
		else if(strcmp(argv[i],"-ipsi_sw_inhib")==0) {swing_inhib_ipsi=1;}
		else if(strcmp(argv[i],"-contra_sw_inhib")==0) {swing_inhib_contra=1;}
		else if(strcmp(argv[i],"-diff_kv")==0) {diff_kv=1;}
		else if(strcmp(argv[i],"-cor")==0) {cor=1;}
		else if(strcmp(argv[i],"-i")==0) {ini=1;}
		else if(strcmp(argv[i],"-slowi")==0) {ini=2;}
		else if(strcmp(argv[i],"-fasti")==0) {ini=3;}
		else if(strcmp(argv[i],"-newi")==0) {newini=1;}
		else if(strcmp(argv[i],"-newslowi")==0) {newini=2;}
		else if(strcmp(argv[i],"-newfasti")==0) {newini=3;}
		else if(strcmp(argv[i],"-cookini")==0) {strcpy(inifiles[1],argv[++i]); newini=1;}
		else if(strcmp(argv[i],"-swpini")==0) {strcpy(inifiles[1],argv[++i]); ini=1;}
		else if(strcmp(argv[i],"-presetWALK")==0) {iniWALK=1;}
		else if(strcmp(argv[i],"-presetPACE")==0) {iniPACE=1;}
		else if(strcmp(argv[i],"-presetTROT")==0) {iniTROT=1;}
		else if(strcmp(argv[i],"-presetRANDOM")==0) {iniRANDOM=1;}
		else if(strcmp(argv[i],"-out_swp")==0) {strcpy(fdur,argv[++i]); stridedata_out=1;}
		else if(strcmp(argv[i],"-nicepic")==0) {nicepic=1;}
		else return 1;
	}

	//output settings
	ofstream out_timers(fdur);
	ofstream out_ini(inifiles[newini]);
	// ofstream out_timers("strideinfo");

	if(ani_out) plotcmd_head();
	if(nicepic) {start_out_time=5;}

	//initialize locomotion parameters
	F0=Fini;
	Gu=Guini;
	Gv=Gvini;
	kv=kvini;
	Tswc=Tswcini;

	if(ini==0)
	{
		xc=0; yc=0;
		theta=M_PI/4;
		vx=2*F0*cos(theta); vy=2*F0*sin(theta);
		omega=0;

		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);
		xfh0[1]=xf-lf*sin(theta), xfh0[2]=xf+lf*sin(theta);
		xfh0[3]=xh-lh*sin(theta), xfh0[4]=xh+lh*sin(theta);
		yfh0[1]=yf+lf*cos(theta), yfh0[2]=yf-lf*cos(theta);
		yfh0[3]=yh+lh*cos(theta), yfh0[4]=yh-lh*cos(theta);

		double qq=L/sqrt(2)/2;
		double r=0.3*qq;

		if(iniWALK==1)
		{
			sw[1]=0; sw[2]=1; sw[3]=0; sw[4]=0; 
			xfh[1]=xfh0[1], yfh[1]=yfh0[1];
			xfh[2]=xfh0[2]-qq, yfh[2]=yfh0[2]-qq;
			xfh[3]=xfh0[3], yfh[3]=yfh0[3];
			xfh[4]=xfh0[4]-qq, yfh[4]=yfh0[4]-qq;
		}
		else if(iniTROT==1)
		{
			sw[1]=1; sw[2]=0; sw[3]=0; sw[4]=1; 
			xfh[1]=xfh0[1], yfh[1]=yfh0[1];
			xfh[2]=xfh0[2]-qq, yfh[2]=yfh0[2]-qq;
			xfh[3]=xfh0[3]-qq, yfh[3]=yfh0[3]-qq+.001;
			xfh[4]=xfh0[4], yfh[4]=yfh0[4];
		}
		else if(iniPACE==1)
		{
			sw[1]=1; sw[2]=0; sw[3]=1; sw[4]=0;
			xfh[1]=xfh0[1], yfh[1]=yfh0[1];
			xfh[2]=xfh0[2]-qq, yfh[2]=yfh0[2]-qq;
			xfh[3]=xfh0[3], yfh[3]=yfh0[3];
			xfh[4]=xfh0[4]-qq, yfh[4]=yfh0[4]-qq;
		}
		else if(iniRANDOM=1)
		{
			sw[1]=1; sw[2]=0; sw[3]=0; sw[4]=1; 
			xfh[1]=xfh0[1]; yfh[1]=yfh0[1];
			xfh[2]=xfh0[2]-qq+r*rnd(); yfh[2]=yfh0[2]-qq+r*rnd();
			xfh[3]=xfh0[3]-qq+r*rnd(); yfh[3]=yfh0[3]-qq+r*rnd();
			xfh[4]=xfh0[4]; yfh[4]=yfh0[4];
		}
		else
		{
			xfh[1]=xfh0[1]-3*cos(theta); yfh[1]=yfh0[1]-3*sin(theta);
			xfh[2]=xfh0[2]; yfh[2]=yfh0[2];
			xfh[3]=xfh0[3]; yfh[3]=yfh0[3];
			xfh[4]=xfh0[4]-3*cos(theta); yfh[4]=yfh0[4]-3*sin(theta);
		}
	}
	else
	{
		ifstream(inifiles[ini])>>xc>>yc>>vx>>vy>>hc>>theta>>omega>>Gpre
		>>xfh[1]>>xfh[2]>>xfh[3]>>xfh[4]>>yfh[1]>>yfh[2]>>yfh[3]>>yfh[4]
		>>sw[1]>>sw[2]>>sw[3]>>sw[4]
		>>tswpre[1]>>tswpre[2]>>tswpre[3]>>tswpre[4];
	
		for(int i=1;i<=4;i++) {xfh[i]-=xc; yfh[i]-=yc;}
		xc=0; yc=0;	//reset to (0,0) position
		xf=xc+L1*cos(theta), yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta), yh=yc-L2*sin(theta);
		xfh0[1]=xf-lf*sin(theta), xfh0[2]=xf+lf*sin(theta);
		xfh0[3]=xh-lh*sin(theta), xfh0[4]=xh+lh*sin(theta);
		yfh0[1]=yf+lf*cos(theta), yfh0[2]=yf-lf*cos(theta);
		yfh0[3]=yh+lh*cos(theta), yfh0[4]=yh-lh*cos(theta);
	}
	vv=hypot(vx,vy);
	vexp=vv;


	//step cycles
	for(double t=0;t<=T;t+=dt)
	{
		swing_count=0;
		for(int i=1;i<=4;i++) if(sw[i]) swing_count++;	//swing legs count

		if(nicepic && t>=start_out_time)
		{
			F0=Fini+(Ffin-Fini)*(t-start_out_time)/(T-start_out_time);
			Gu=Guini+(Gufin-Guini)*(t-start_out_time)/(T-start_out_time);
			Gv=Gvini+(Gvfin-Gvini)*(t-start_out_time)/(T-start_out_time);
			kv=kvini+(kvfin-kvini)*(t-start_out_time)/(T-start_out_time);
			Tswc=Tswcini+(Tswcfin-Tswcini)*(t-start_out_time)/(T-start_out_time);
		}
		else if(!nicepic)
		{
			F0=Fini+(Ffin-Fini)*t/T;
			Gu=Guini+(Gufin-Guini)*t/T;
			Gv=Gvini+(Gvfin-Gvini)*t/T;
			kv=kvini+(kvfin-kvini)*t/T;
			Tswc=Tswcini+(Tswcfin-Tswcini)*t/T;
		}

		//load calculations
		for(int i=1;i<=4;i++) load[i]=0;
		total_load=0;
		update_F(xfh, yfh, xc, yc, F);

		Zx=0; Zy=0;
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
			//cerr<<"fall"<<endl;
			// break;
		}
		else if(total_load<0)
		{
			//cerr<<"fly"<<endl;
			break;
		}

		//swing limb movement
		double taus=.0001;
		double Vswc=(L+vv*Tswc)/Tswc;
		if(cor) update_corrected_swing_target(theta,vx,vy,vexp,xfh0,yfh0,xfh1,yfh1);

		for(int i=1;i<=4;i++) if(sw[i])
		{
			if(cor) {xfh[i]+=(xfh1[i]-xfh[i])*dt/taus; yfh[i]+=(yfh1[i]-yfh[i])*dt/taus;}
			else {xfh[i]+=(xfh0[i]-xfh[i])*dt/taus; yfh[i]+=(yfh0[i]-yfh[i])*dt/taus;}
		}

		//horizontal force calculations
		double Fx[5]={},Fy[5]={};
		double delta=0*omega;
		double FL=(1+delta)*F0, FR=(1-delta)*F0;
		
		for(int i=1;i<=4;i++) {d[i]=hypot(xfh0[i]-xfh[i],yfh0[i]-yfh[i]);}

		if(!sw[1] && d[1]>0.) { Fx[1]=FL*(xfh0[1]-xfh[1])/d[1]; Fy[1]=FL*(yfh0[1]-yfh[1])/d[1]; }
		if(!sw[2] && d[2]>0.) { Fx[2]=FR*(xfh0[2]-xfh[2])/d[2]; Fy[2]=FR*(yfh0[2]-yfh[2])/d[2]; }
		if(!sw[3] && d[3]>0.) { Fx[3]=FL*(xfh0[3]-xfh[3])/d[3]; Fy[3]=FL*(yfh0[3]-yfh[3])/d[3]; }
		if(!sw[4] && d[4]>0.) { Fx[4]=FR*(xfh0[4]-xfh[4])/d[4]; Fy[4]=FR*(yfh0[4]-yfh[4])/d[4]; }

		xc+=vx*dt;
		yc+=vy*dt;
		theta+=omega*dt;
		vx+=((Fx[1]+Fx[2]+Fx[3]+Fx[4])/tau-vx/tau+g*Zx/H)*dt;
		vy+=((Fy[1]+Fy[2]+Fy[3]+Fy[4])/tau-vy/tau+g*Zy/H)*dt;
		vv=hypot(vx,vy);

		double tauvv=1;
		vexp+=(vv-vexp)*dt/tauvv;

		xf=xc+L1*cos(theta); yf=yc+L1*sin(theta);
		xh=xc-L2*cos(theta); yh=yc-L2*sin(theta);
		xfh0[1]=xf-lf*sin(theta), xfh0[2]=xf+lf*sin(theta);
		xfh0[3]=xh-lh*sin(theta), xfh0[4]=xh+lh*sin(theta);
		yfh0[1]=yf+lf*cos(theta), yfh0[2]=yf-lf*cos(theta);
		yfh0[3]=yh+lh*cos(theta), yfh0[4]=yh-lh*cos(theta);

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

		
		for(int k=1;k<=4;k++) tsw[k]=0;

		for(int k=1;k<=4;k++) if(sw[k]) if(t-tswpre[k]>=Tswc) 
		{
			if(!positive_load_check(k,swing_count,sw,F)) continue;
			sw[k]=0;
			swing_count--;
		}	//stop swing a leg after specific time, and skip if it gets negative load after landing

		for(int k=1;k<=4;k++) if(!sw[k]) if(GP[k]>load[k] && load[k]<Gu)
		{
			sw[k]=1;
			swing_count++;
			tswpre[k]=t;
			if(swing_inhib_ipsi && k<=2 && sw[ipsilateral[k]] && positive_load_check(ipsilateral[k],swing_count,sw,F)) 
			{
				sw[ipsilateral[k]]=0; swing_count--;
			}
			if(swing_inhib_contra && sw[contralateral[k]] && positive_load_check(contralateral[k],swing_count,sw,F)) 
			{
				sw[contralateral[k]]=0; swing_count--;
			}
		}	//weak lifting conditions

		for(int k=1;k<=4;k++) if(!sw[k]) if(load[k]<0 || d[k]>DM[k])
		{
			sw[k]=1;
			swing_count++;
            tswpre[k]=t;
			if(swing_inhib_ipsi && k<=2 && sw[ipsilateral[k]] && positive_load_check(ipsilateral[k],swing_count,sw,F)) 
			{
				sw[ipsilateral[k]]=0; swing_count--;
			}
			if(swing_inhib_contra && sw[contralateral[k]] && positive_load_check(contralateral[k],swing_count,sw,F)) 
			{
				sw[contralateral[k]]=0; swing_count--;
			}
		}	//strong lifting conditions

		double balance=total_load>0? (Gpre-total_load)/dt:0;
		
		// for(int k=1;k<=4;k++) if(balance>(sw[diagnal[k]]? kv:0)) sw[k]=0;
		// for(int k=1;k<=4;k++) if(balance>kv) sw[k]=0;
		

		if(balance>0)
		{
			int kmax=0;
			for(int k=1;k<=4;k++) if(sw[k]) {tsw[k]=t-tswpre[k];}
			for(int k=1;k<=4;k++) if(tsw[k]>tsw[kmax]) {kmax=k;}
			
			int kmax_positiveload=positive_load_check(kmax,swing_count,sw,F);
			
			double losing_balance_threshold[5]={0,kv,kv,kv,kv};
			if(diff_kv) for(int k=1;k<=4;k++)
			{
				if(k>2) {losing_balance_threshold[k]=0;}
				else {losing_balance_threshold[k]=(sw[diagnal[k]]? kv:0);}
			}
			if(kmax_positiveload && (balance>losing_balance_threshold[kmax]))
			{
				sw[kmax]=0;
				swing_count--;
			}
		}	//stop swing a leg when lose balance
		// if(balance>kv)
		// {
		// 	int kmax=0;
		// 	for(int k=1;k<=4;k++) if(sw[k]) {tsw[k]=t-tswpre[k];}
		// 	for(int k=1;k<=4;k++) if(tsw[k]>tsw[kmax]) {kmax=k;}
		// 	int kmax_positiveload=1;
        //     if(swing_count==2)
		// 	{
		// 		int s=1; while(!sw[s] || s==kmax) s++;
		// 		if(F[kmax][s]<0.) kmax_positiveload=0;
		// 	}
		// 	else if(swing_count==1 && F[kmax][0]<0.) kmax_positiveload=0;
		// 	if(kmax_positiveload)
		// 	{
		// 		sw[kmax]=0;
		// 		swing_count--;
		// 	}
		// }	//stop swing a leg when lose balance
		

		if(swing_count==4)
		{
			int kmax=0;
			for(int k=1;k<=4;k++) if(sw[k]) {tsw[k]=t-tswpre[k];}
			for(int k=1;k<=4;k++) if(tsw[k]>tsw[kmax]) {kmax=k;}
			sw[kmax]=0;
			swing_count--;
		}	//stop swing a leg when flying (manage data for better figure)


		for(int k=1;k<=4;k++) GP[k]=load[k];
		Gpre=total_load;

		//falling indicator: if there is at least one leg having negative load for a long time (say 0.2s), it is falling
		// static double fallpre=0;
		// int iter_nonneg;
		// iter_nonneg=1;
		// while(load[iter_nonneg]>=0 && iter_nonneg<4) iter_nonneg++;
		// if(iter_nonneg<4) {if(t-fallpre>0.2) break;}
		// else fallpre=t;

		//output data
		if(ani_out && int(t/dt)%int(DT/dt)==0)
		{
			cout<<t<<'\t'<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<total_load;
			for(int i=1;i<=4;i++) cout<<'\t'<<load[i];
			for(int i=1;i<=4;i++) cout<<'\t'<<(sw[i]?i:0);
			// cout<<'\t'<<deltal<<'\t'<<deltar;
			cout<<endl;
			
			// animation frames commands using gnuplot
			plotcmd_frame(sw,xc,yc,hc,xf,yf,xh,yh,xfh,yfh,xfh0,yfh0);
		}

		//output timers
		static int swpre[5];
		static double ttswpre[5]={};
		double ttsw[5]={};
		static double ttstpre[5]={};
		double ttst[5]={};
		static double stridepre=0;
		static double stridetime=0;
		static double dutyfacor_pre=0;
		static double vvpre=0;
		for(int k=1;k<=4;k++) if(swpre[k] && !sw[k])	//stop swing tracker
		{
			ttsw[k]=t-ttswpre[k];
			ttstpre[k]=t; 
			if(k==1) { stridetime=t-stridepre; stridepre=t; }			
			//if(stridedata_out==1 && stridepre!=0) out_timers<<t<<'\t'<<vv<<'\t'<<k<<'\t'<<ttsw[k]<<'\t'<<-1<<'\t'<<1/stridetime<<endl;
		}
		for(int k=1;k<=4;k++) if(!swpre[k] && sw[k])	//start swing tracker
		{
			ttst[k]=t-ttstpre[k];
			ttswpre[k]=t;

			double dutyfactor=(stridetime==0? 0:(100*ttst[k]/stridetime));

			if(t>10 && abs(dutyfactor-dutyfacor_pre)>80) return 0;

			if(t>=1 && stridepre!=0) out_timers<<k<<'\t'<<t<<'\t'<<F0<<'\t'<<kv<<'\t'<<Gu<<'\t'<<Tswc<<'\t'<<vv<<'\t'<<ttst[k]<<'\t'<<(stridetime-ttst[k])<<'\t'<<dutyfactor<<endl;
			//if(stridepre!=0) out_timers<<k<<'\t'<<t<<'\t'<<F0<<'\t'<<kv<<'\t'<<Gu<<'\t'<<Tswc<<'\t'<<vv<<'\t'<<ttst[k]<<'\t'<<ttsw[k]<<'\t'<<dutyfactor<<endl;
			// Limb (k) is starts swing at time (t) with parameters (F0), (kv), (Gu) and (Tsw);
			// it has velocity (vv), stance time (ttst[k]), previous swing time (ttsw[k]) and duty factor (dutyfactor).
			
			//if(stridepre!=0) out_timers<<kv<<'\t'<<F0<<'\t'<<vv<<'\t'<<k<<'\t'<<-1<<'\t'<<ttst[k]<<'\t'<<dutyfactor<<endl;
			//if(stridepre!=0 && abs(dutyfactor-dutyfacor_pre)<0.15 && abs(vv-vvpre)<0.5) out_timers<<Tswc<<'\t'<<F0<<'\t'<<vv<<'\t'<<k<<'\t'<<-1<<'\t'<<ttst[k]<<'\t'<<dutyfactor<<endl;
			dutyfacor_pre=dutyfactor;
			vvpre=vv;
			//if(stridedata_out==1 && stridepre!=0) out_timers<<t<<'\t'<<vv<<'\t'<<k<<'\t'<<-1<<'\t'<<ttst[k]<<'\t'<<1/stridetime<<endl; && abs(Tswc-ttsw[k])<0.04
		}
		for(int k=1;k<=4;k++) swpre[k]=sw[k];
	}
	
	//output initial conditions for the future
	if(newini)
	{
		for(int k=1;k<=4;k++) tswpre[k]-=T;
		out_ini<<xc<<'\t'<<yc<<'\t'<<vx<<'\t'<<vy<<'\t'<<hc<<'\t'<<theta<<'\t'<<omega<<'\t'<<Gpre<<'\t'
		<<xfh[1]<<'\t'<<xfh[2]<<'\t'<<xfh[3]<<'\t'<<xfh[4]<<'\t'
		<<yfh[1]<<'\t'<<yfh[2]<<'\t'<<yfh[3]<<'\t'<<yfh[4]<<'\t'
		<<sw[1]<<'\t'<<sw[2]<<'\t'<<sw[3]<<'\t'<<sw[4]<<'\t'
		<<tswpre[1]<<'\t'<<tswpre[2]<<'\t'<<tswpre[3]<<'\t'<<tswpre[4]<<endl;
	}
	
	return 0;
}