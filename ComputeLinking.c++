 #include<iostream>
#include<stdlib.h>
#include<cmath>
#include<sstream>
#include<fstream>
#include<string>
#include<iomanip>
using namespace std;

const double PI  = 3.141592653589793238462;

long int k=0;
long int Kmax;
int Kstart = 0;
const int Nmax=160;
const int Nbeadsmax=2000;

int N,Nbeads;
int DeltaT=1000;
double dt=0.01;
double Lx,Ly,Lz;

//Observables
//t=0
double COMRing0[Nmax][3];
double position0[Nmax][Nbeadsmax][3];
//t>0
double COMRing[Nmax][3];
double COMRingWrapped[Nmax][3];
double position[Nmax][Nbeadsmax][3];
int ringl[Nmax];
double positionWrapped[Nmax][Nbeadsmax][3];
double COM[3];
double COM0[3];
//
string filename;
//
double gaussLnkNumb(double c1[][3], int n1, double c2[][3], int n2);
double min(double , double);
double Compdist(double, double, double, double,double,double);
double reduceInt(int v[], int s);

void computeLk(int n, int m);
double Dist(double p[6], double v[9]);
double distance(double p[6], double v[6]);
double* wedge(double v[3],double u[3]);
double norm(double u[3]);
double dot(double u[3], double v[3]);

void computelinking();
double TotalLinking=0;
double TotalLinkingIR=0;
int LinkingMap[Nmax][Nmax];

int main(int argc, char* argv[]){

cout << "Write argv[1]: datain; argv[2] = dataout; argv[3]=Tmax; argv[4]=N; argv[5]=Nbeads; argv[6]=Kstart argv[7]=DeltaT"<< endl;
cout << "Num of files to convert:";
Kmax=int(atoi(argv[3]));
cout << Kmax <<endl;
cout << "Number of Polymers:";
N=atoi(argv[4]);
cout << N <<endl;
cout<<"Number of Beads in each polymer:";
Nbeads=atoi(argv[5]);
cout << Nbeads<<endl; 
cout<<"Start frame:";
Kstart=int(atoi(argv[6]));
cout << Kstart <<endl;
cout<<"DeltaT:";
DeltaT=atoi(argv[7]);
cout << DeltaT <<endl;


for(int n=0;n<Nmax;n++)for(int i=0;i<3;i++) COMRing0[n][i]= COMRingWrapped[n][i]=0;

stringstream writeFileTS;
writeFileTS <<"LinkingMatrix_"<<argv[2];
ofstream writeTS(writeFileTS.str().c_str());
writeTS << "#time n m lk" <<endl;

filename=argv[2];
/////////////////////////////////////////////////////
//MAIN LOOP OVER TIME!!!
////////////////////////////////////////////////////
//Kmax-= Kstart;
for(k=0;k<Kmax;k++){
//
//
ifstream indata;
stringstream readFile;
readFile.str("");
readFile.clear();
long long int sum = int((k+Kstart));
if(sum==0) readFile << argv[1] <<sum;
if(DeltaT==100)if(sum>0)readFile << argv[1] <<sum << "00"; //n of zeros must match DeltaT
if(DeltaT==1000)if(sum>0)readFile << argv[1] <<sum << "000"; //n of zeros must match DeltaT
if(DeltaT==10000)if(sum>0)readFile << argv[1] <<sum << "0000"; //n of zeros must match DeltaT
indata.open(readFile.str().c_str());
cout << readFile.str().c_str()<<endl;
if(!indata){cout <<"file "<< readFile.str().c_str() << " is not there"<<endl; return 0;}
	
for(int n=0;n<N;n++)for(int i=0;i<3;i++)COMRing[n][i]=COM[i]=0;
for(int n=0;n<N;n++)for(int m=0;m<N;m++)LinkingMap[n][m]=0;

long int ring,ringlength;
double id,mol,type;
double num1,num2,num3;
double x,y,z;
string dummy;
long long int time;
long int Ntot;
double l1,l2;

//read 10 lines
for(int i=0;i<10;i++){
if(i==1) {
indata >> time;
time = k*DeltaT*dt;
cout << "time " << time <<endl;
}
if(i==3) {
indata >> Ntot;
cout << "Ntot " << Ntot<<endl;
}
if(i==5) {
indata >> l1 >> l2;
cout << "L " << l1<< " " << l2 <<endl;
Lx = l2-l1;
cout << "Lx " << Lx <<endl;
}
if(i==6) {
indata >> l1 >> l2;
cout << "L " << l1<< " " << l2 <<endl;
Ly = l2-l1;
cout << "Ly " << Ly <<endl;
}
if(i==7) {
indata >> l1 >> l2;
//cout << "L " << l1<< " " << l2 <<endl;
Lz = l2-l1;
//cout << "Lz " << Lz<<endl;
}

else {
      getline(indata,dummy);
//cout << dummy<<endl;
}
}
 
//////////////////////////
//READ FILES
////////////////////////////
int mem=1;
int old_mem=1;
int ncount=0;

for(int n=0; n<Ntot; n++){
	//indata >> id>>type>>  x>>y>>z>>num1>>num2>>num3;
	indata >> id >> mol >> ring >> ringlength >> type >> x >> y >> z >> num1 >> num2 >> num3;
	//cout << id << " " << x <<endl; cin.get();
        cout << "id " << id << " mol " << mol << " ring " << ring << " ringlength " << ringlength << " type " << type << " "  << x << " " << y << " " <<  z<< " "  << num1 << " " << num2 << " " << num3 <<endl;
	cout << Ntot;
        indata >>id >> mol ;
        //cout << id <<" "<< mol << " " <<endl;
        ringl[ring]=ringlength*Nbeads;
	mem = ring;
        //cout << "yo "  <<  ring << " mamma mia" <<endl;
	if(mem!=old_mem){
	ncount=0;
	mem=old_mem;
        
        }
	int Nring=ring; 
	//int Nring = floor((double)(id-1)/(1.0*ringlength*Nbeads));
	//int Nring = floor((double)(type-2));
	//int Nring = floor((double)(type-1));
	
	//cout << "hey " << id <<  " " << ring <<  " " <<  ncount <<endl; 
	//cin.get();
        //cout << "mamma mia" << Nring <<endl;
	if(Nring<N){
	position[Nring][ncount][0] = (x + Lx*num1);
	position[Nring][ncount][1] = (y + Ly*num2);
	position[Nring][ncount][2] = (z + Lz*num3);

	positionWrapped[Nring][ncount][0] = x;
	positionWrapped[Nring][ncount][1] = y;
	positionWrapped[Nring][ncount][2] = z;

	if(positionWrapped[Nring][ncount][0]-positionWrapped[Nring][ncount-1][0]>Lx/2.){
		positionWrapped[Nring][ncount][0]-Lx;
		}
	if(positionWrapped[Nring][ncount][0]-positionWrapped[Nring][ncount-1][0]<-Lx/2.)
		positionWrapped[Nring][ncount][0]+Lx;
	if(positionWrapped[Nring][ncount][1]-positionWrapped[Nring][ncount-1][1]>Ly/2.){
		positionWrapped[Nring][ncount][1]-Ly;
		//cout << "y pbc "<<endl;
		}
	if(positionWrapped[Nring][ncount][1]-positionWrapped[Nring][ncount-1][1]<-Ly/2.)
		positionWrapped[Nring][ncount][1]+Ly;
	if(positionWrapped[Nring][ncount][2]-positionWrapped[Nring][ncount-1][2]>Lz/2.)
		positionWrapped[Nring][ncount][2]-Lz;
	if(positionWrapped[Nring][ncount][2]-positionWrapped[Nring][ncount-1][2]<-Lz/2.)
		positionWrapped[Nring][ncount][2]+Lz;


	//if(Nring==1) cout << Nring << " " << ncount <<" " <<  positionWrapped[Nring][ncount][0] << " " << positionWrapped[Nring][ncount][1] << " " << positionWrapped[Nring][ncount][2] <<endl;
	
	if(k==0){
	position0[Nring][ncount][0] = position[Nring][ncount][0];
	position0[Nring][ncount][1] = position[Nring][ncount][1];
	position0[Nring][ncount][2] = position[Nring][ncount][2];

	COMRing0[Nring][0] += (position0[Nring][ncount][0])/(ringlength*Nbeads);
	COMRing0[Nring][1] += (position0[Nring][ncount][1])/(ringlength*Nbeads);
	COMRing0[Nring][2] += (position0[Nring][ncount][2])/(ringlength*Nbeads);

	COM0[0] += (x + Lx*num1)/(double(N*(ringlength*Nbeads)));
	COM0[1] += (y + Ly*num2)/(double(N*(ringlength*Nbeads)));	
	COM0[2] += (z + Lz*num3)/(double(N*(ringlength*Nbeads)));
	}
	
	COMRing[Nring][0] += (x + Lx*num1)/(ringlength*Nbeads);
	COMRing[Nring][1] += (y + Ly*num2)/(ringlength*Nbeads);
	COMRing[Nring][2] += (z + Lz*num3)/(ringlength*Nbeads);

	COM[0] += (x + Lx*num1)/(double(N*(ringlength*Nbeads)));
	COM[1] += (y + Ly*num2)/(double(N*(ringlength*Nbeads)));	
	COM[2] += (z + Lz*num3)/(double(N*(ringlength*Nbeads)));
	
	}	

ncount++;
}
 

//MAP Threadings/////////////////////////
//cout << "Checking Linking between rings ... " << endl;
cout << "ur maw" << endl;
computelinking();
//////////////////////////////
cout << "time t= " << k << " -> Nlinking = " <<  TotalLinking <<endl; // cin.get();
//cout << "check" << N <<endl;

//for(int n=0;n<N;n++)for(int m=0;m<N;m++)if(LinkingMap[n][m]!=0)writeTS << k+Kstart << " " << n << " " << m << " " << LinkingMap[n][m] << endl;
for(int n=0;n<N;n++)for(int m=0;m<N;m++)writeTS << k+Kstart << " " << n << " " << m << " " << LinkingMap[n][m] << endl;

}//closes time (loop over k)


return 0 ;
}
//
//cout << "ur maw" <<endl;
//////////////////////////////////////
//FUNCTIONS
////////////////////////////////////////
void computelinking(){

TotalLinking=0;
TotalLinkingIR=0;
double lk_pz=0;
double lk_mz=0;
double lk=0;
double Curve1[Nbeadsmax][3];
double Curve2[Nbeadsmax][3];
int n1,n2;


for(int n=0;n<N;n++){
	//indata >> id>>mol>>ring>>ringlength>>type>>x>>y>>z>>num1>>num2>>num3;

cout << "Ring " << n << " .. " <<endl;
//CREATE A CURVE1
n1=ringl[n];
for(int nb=0;nb<n1;nb++)for(int d=0;d<3;d++)Curve1[nb][d]=position[n][nb][d];

//for(int nb=0;nb<Nhead;nb++)write << nb<< " "<< n << " "  << Curve1[nb][0] << " " << Curve1[nb][1] << " " << Curve1[nb][2] <<endl;
	//now look at all other rings
	for(int m=n+1;m<N;m++){

	//FIRST- SHIFT RINGS TO REMOVE PBC EFFECTS (d<L/2)
	double d_ring_ring[3];
	int shift[3];
	double dd=0;
	for(int d=0;d<3;d++)d_ring_ring[d]=COMRing[m][d]-COMRing[n][d];
	shift[0]=round(d_ring_ring[0]*1.0/Lx);
	shift[1]=round(d_ring_ring[1]*1.0/Ly);
	shift[2]=round(d_ring_ring[2]*1.0/Lz);

    //CREATE CURVE 2
	n2=ringl[m];
	for(int mb=0;mb<n2;mb++){
	Curve2[mb][0]=position[m][mb][0]-shift[0]*Lx;
	Curve2[mb][1]=position[m][mb][1]-shift[1]*Ly;
	Curve2[mb][2]=position[m][mb][2]-shift[2]*Lz;
	}
	//for(int d=0;d<3;d++)cout <<n << " " << m << " " << (COMHead[n][d]-COMTailWrapped1[m][d])/Lx<<endl;
	//cin.get();
	//////////////////////////////////////

	lk_pz=0;lk_mz=0;lk=0;

//NOW COMPUTE LINKING BETWEEN THE 2 CURVES///////////
if(n!=m)lk=gaussLnkNumb(Curve1,n1,Curve2,n2);
/////////////////////////////////////////////////////////////

if(lk!=0)TotalLinking++;
if(lk!=0)LinkingMap[n][m]=1;
if(lk!=0 && n!=m)TotalLinkingIR++;
/////////////////////////////////////////////////

/////////////////////////////////////////////////////////
//IF LINKING !=0 WRITE CONFIGS
if(lk!=0 && k==0){

//DAT
//stringstream writeFile;
//writeFile <<"ThreadingConfigs_"<<n<<"-"<<m<<".dat";
//ofstream write(writeFile.str().c_str());
//for(int nb=0;nb<Nhead;nb++)write << nb<< " "<< n << " "  << Curve1[nb][0] << " " << Curve1[nb][1] << " " << Curve1[nb][2] <<endl;
//write<<endl;
//write<<endl;
//for(int mb=0;mb<Ntail+3*dclose;mb++)write << mb << " "<< m << " "  << Curve2[mb][0] << " " << Curve2[mb][1] << " " << Curve2[mb][2] <<endl;

//VMD (open with "vmd ThreadingConfig...")

stringstream writeFile;
writeFile <<"TLink/LinkedConfigs_"<<filename<<"-k"<<k<<"_"<<n<<"-"<<m<<".xyz";
ofstream write(writeFile.str().c_str());
write << 2*Nbeads<< endl;
write << "Atoms. Timestep : 0"<< endl;
for(int nb=0;nb<Nbeads;nb++)write << "C "  << Curve1[nb][0] << " " << Curve1[nb][1] << " " << Curve1[nb][2] <<endl;
for(int mb=0;mb<Nbeads;mb++)write  << "O "  << Curve2[mb][0] << " " << Curve2[mb][1] << " " << Curve2[mb][2] <<endl;
//cout << "STOOOP!!! " <<endl;cin.get();
}
///////////////////////////////////////////////////////////

	} //close loop over heads m
}//close loop over heads n

} //close function


double reduceInt(int v[], int s){
double r=0.;
	for(int n=0;n<s;n++) r+=v[n];
return r;	
}

double min(double a , double b){
	if(a<b) return a;
	else return b;
}
double Compdist(double cx,double cy, double cz, double c1x, double c1y, double c1z){
	double d=0.0;
	d = sqrt((cx-c1x)*(cx-c1x) + (cy-c1y)*(cy-c1y) + (cz-c1z)*(cz-c1z));
	return d;
}

double distance(double v1[6], double v2[6]){
double d=0;

d=(v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2]);

return sqrt(d);
}
double Dist(double v1[6], double v2[9]){
double d=0;

d=(v1[0]-v2[3])*(v1[0]-v2[3])+(v1[1]-v2[4])*(v1[1]-v2[4])+(v1[2]-v2[5])*(v1[2]-v2[5]);

return sqrt(d);
}

double* wedge(double v[3],double u[3]){
static double bnorm[3]={0.,0.,0.}; 
bnorm[0]=v[1]*u[2]-v[2]*u[1];
bnorm[1]=-v[0]*u[2]+v[2]*u[0];
bnorm[2]=v[0]*u[1]-v[1]*u[0];
return bnorm;
}

double norm(double u[3]){
	double n = u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
	return sqrt(n);
}
double dot(double u[3], double v[3]){
	double n = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
	return n;
}

double gaussLnkNumb(double c1[][3], int n1, double c2[][3], int n2){
double lk=0;
long double num=0;
double dr1[3],dr2[3],r[3],cross[3];
long double dist, dist3; 
//for(int n =0 ;n<n1; n++)cout << "1 "<<  c1[n][0] << " " << c1[n][1] << " " <<c1[n][2] << endl;
//for(int m=0;m<n2;m++)cout << "2 "<<  c2[m][0] << " " << c2[m][1] << " " <<c2[m][2] << endl;

for(int n =0 ;n<n1; n++){
	for(int m=0;m<n2;m++){
		r[0] = c2[m][0] - c1[n][0];
		r[1] = c2[m][1] - c1[n][1];
		r[2] = c2[m][2] - c1[n][2];
		dist = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
		dist3 = dist*dist*dist;
		
		dr1[0] = c1[(n+1)%n1][0] - c1[n][0];
		dr1[1] = c1[(n+1)%n1][1] - c1[n][1];
		dr1[2] = c1[(n+1)%n1][2] - c1[n][2];
		
		dr2[0] = c2[(m+1)%n2][0] - c2[m][0];
		dr2[1] = c2[(m+1)%n2][1] - c2[m][1];
		dr2[2] = c2[(m+1)%n2][2] - c2[m][2];
		
		cross[0] = dr1[1]*dr2[2] - dr1[2]*dr2[1];
		cross[1] = dr1[2]*dr2[0] - dr1[0]*dr2[2];
		cross[2] = dr1[0]*dr2[1] - dr1[1]*dr2[0];
		
		num = r[0]*cross[0] + r[1]*cross[1] + r[2]*cross[2];
		
		lk += num*1.0/dist3;
		
	//if(ch==1)cout << n << " " << m << " " << dist << " " << num << " " << lk*1.0/(4.0*M_PI)  <<endl;
	}
}	
	
return round(lk*1.0/(4.0*M_PI));
}


