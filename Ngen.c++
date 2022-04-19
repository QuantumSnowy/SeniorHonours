
using namespace std;
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<sstream>
#include<vector> 
#include <unistd.h>
#include<ctime>

int InitialConfg=0;

int bond1=2;
int ang1=3;
int c=1;
int nb=1;
int na=1;

const int npolymax=1000;
const int nbeadsmax=2048;

int bond[npolymax*nbeadsmax][2];
int angle[npolymax*nbeadsmax][3];
double position[npolymax][nbeadsmax][3];

int reps;
int nbond=0;
int nangle=0;


int main(int argc, char* argv[]){

cout << "Input argv: argv[1]:Npoly argv[2]=Nbeads; argv[3]=L; argv[4]=ntype"<<endl;
srand(time(NULL));

int npoly,nbeads, ntype;
npoly=atoi(argv[1]);
nbeads=atoi(argv[2]);
ntype=atoi(argv[4]);
double L=atof(argv[3]); 
     
ofstream writeW;
stringstream writeFileW;
writeFileW << "melt.ring.l"<<L<<"n"<<npoly<<"m"<<nbeads<<"t"<<ntype<<".0";
writeW.open(writeFileW.str().c_str());
    
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
for(reps=0;reps<npoly;reps++){
//sleep(0.5);
cout << reps <<endl;

//////////////////////////////////////////////////////////////
////////		MAKE POLYMER     /////////////////////
//////////////////////////////////////////////////////////////
cout << "BONDS" <<endl;
int cumulative=0;
for(int rr=0;rr<reps;rr++)cumulative+=nbeads;
int nbeadT=0;
for(int i=0;i<nbeads-1;i++){
    nbeadT++;
    //cout << " if " <<endl;
    bond[nbond][0] = i+cumulative;
    bond[nbond][1] = i+1+cumulative;
    nbond++;
    //cout << nbond << " " << i << endl;
}
    
for(int i=0;i<nbeads-2;i++){
    angle[nangle][0] = i+cumulative;
    angle[nangle][1] = i+1+cumulative;
    angle[nangle][2] = i+2+cumulative;
    nangle++;
}

cout << "MAKING THE POLYMER ..." <<endl;
///////////////////////////////////////////////////////////
//make polymer
double theta[2];
double phi[2];
double non1=((double)(rand())/((double)(RAND_MAX)));
double non=((double)(rand())/((double)(RAND_MAX)));
    
///////////////////////////////////////////////////////////////////////////
// RING
///////////////////////////////////////////////////////////////////////////
double zshift;
double yshift;
double xshift;
for(int m=0;m<nbeads;m++){
if(m==0){
position[reps][m][0]=((double)(rand())/((double)(RAND_MAX)))*L-L/2.;
position[reps][m][1]=((double)(rand())/((double)(RAND_MAX)))*L-L/2.;
position[reps][m][2]=((double)(rand())/((double)(RAND_MAX)))*L-L/2.;
}
else{
position[reps][m][0]=position[reps][m-1][0]+0.9*(((double)(rand())/((double)(RAND_MAX)))*2-1);
position[reps][m][1]=position[reps][m-1][1]+0.9*(((double)(rand())/((double)(RAND_MAX)))*2-1);
position[reps][m][2]=position[reps][m-1][2]+0.9*(((double)(rand())/((double)(RAND_MAX)))*2-1);
}
    
}
cout << "DONE Ring"<<reps << endl;
} //close loop over reps
////////////////////////////////////////////////////
////////		WRITE FILE     /////////////////////
////////////////////////////////////////////////////
int totb=npoly*nbeads;
writeW<< "LAMMPS data file from restart file: timestep = 0,\t procs = 1"<<endl;
writeW << totb << " atoms "<<endl;
writeW << nbond << " bonds "<<endl;
writeW << nangle << " angles "<<endl;
writeW << "\n";
writeW << ntype << " atom types "<<endl;
writeW << 1 << " bond types "<<endl;
writeW << 1 << " angle types "<<endl;
writeW << "\n";
writeW << -L/2.0 << " " << (L-L/2.0) << " xlo xhi"<<endl;
writeW << -L/2.0 << " " << (L-L/2.0) << " ylo yhi"<<endl;
writeW << -L/2.0 << " " << (L-L/2.0) << " zlo zhi"<<endl; ///TO BE CHANGED!!!
//
writeW << "\nMasses \n"<<endl;
for(int j=1; j<=ntype;j++) writeW << j << " " << 1 << endl;
 //
int cc=1;
writeW << "\nAtoms \n"<<endl;
int quant = npoly/(ntype-1);
for(int val=0; val<(ntype-1); val++){
    
    
    for(int nn=val*quant; nn< (val+1)*quant;nn++){
       
       for(int m=0;m<nbeads; m++){
         
         if(m==0 || m==nbeads-1) writeW << cc << " " << nn+1 << " " << val+2 << " " << position[nn][m][0]<<" " << position[nn][m][1] << " " << position[nn][m][2] << " " << 0 << " " << 0 << " " << 0 << endl;
        else writeW << cc << " " << nn+1 << " " << 1 << " " << position[nn][m][0]<<" " << position[nn][m][1] << " " << position[nn][m][2] << " " << 0 << " " << 0 << " " << 0 << endl;
        double distance=0;
        distance=sqrt(pow(position[nn][m][0]-position[nn][m-1][0],2.0)+pow(position[nn][m][1]-position[nn][m-1][1],2.0)+pow(position[nn][m][2]-position[nn][m-1][2],2.0));
        cout << nn << " " << m << " " << distance <<endl;
        cc++;

    }


}


 
}

 
/*
for(int nn=0; nn< floor(npoly/3); nn++){
    for(int m=0;m<nbeads; m++){
         if(m==0 || m==nbeads-1) writeW << cc << " " << nn+1 << " " << 2 << " " << position[nn][m][0]<<" " << position[nn][m][1] << " " << position[nn][m][2] << " " << 0 << " " << 0 << " " << 0 << endl;
        else writeW << cc << " " << nn+1 << " " << 1 << " " << position[nn][m][0]<<" " << position[nn][m][1] << " " << position[nn][m][2] << " " << 0 << " " << 0 << " " << 0 << endl;
        double distance=0;
        distance=sqrt(pow(position[nn][m][0]-position[nn][m-1][0],2.0)+pow(position[nn][m][1]-position[nn][m-1][1],2.0)+pow(position[nn][m][2]-position[nn][m-1][2],2.0));
        cout << nn << " " << m << " " << distance <<endl;
        cc++;
    }
}
for(int nn=int(npoly/3); nn< int(2*npoly/3); nn++){
    for(int m=0;m<nbeads; m++){
         if(m==0 || m==nbeads-1) writeW << cc << " " << nn+1 << " " << 3 << " " << position[nn][m][0]<<" " << position[nn][m][1] << " " << position[nn][m][2] << " " << 0 << " " << 0 << " " << 0 << endl;
        else writeW << cc << " " << nn+1 << " " << 1 << " " << position[nn][m][0]<<" " << position[nn][m][1] << " " << position[nn][m][2] << " " << 0 << " " << 0 << " " << 0 << endl;
        double distance=0;
        distance=sqrt(pow(position[nn][m][0]-position[nn][m-1][0],2.0)+pow(position[nn][m][1]-position[nn][m-1][1],2.0)+pow(position[nn][m][2]-position[nn][m-1][2],2.0));
        cout << nn << " " << m << " " << distance <<endl;
        cc++;
    }
}

for(int nn= int(2*(npoly/3)); nn< int(npoly); nn++){
    for(int m=0;m<nbeads; m++){
         if(m==0 || m==nbeads-1) writeW << cc << " " << nn+1 << " " << 4 << " " << position[nn][m][0]<<" " << position[nn][m][1] << " " << position[nn][m][2] << " " << 0 << " " << 0 << " " << 0 << endl;
        else writeW << cc << " " << nn+1 << " " << 1 << " " << position[nn][m][0]<<" " << position[nn][m][1] << " " << position[nn][m][2] << " " << 0 << " " << 0 << " " << 0 << endl;
        double distance=0;
        distance=sqrt(pow(position[nn][m][0]-position[nn][m-1][0],2.0)+pow(position[nn][m][1]-position[nn][m-1][1],2.0)+pow(position[nn][m][2]-position[nn][m-1][2],2.0));
        cout << nn << " " << m << " " << distance <<endl;
        cc++;
    }
}
*/
///////////////////////////////////////////////////////////
//FINISHED POLYMER
///////////////////////////////////////////////////////////

writeW << endl;
writeW << endl;
writeW << "\n Velocities \n" <<endl;
for(int j=0;j<totb; j++) writeW<<j+1 << " "<<0 << " "<<0<<" "<<0 <<endl;

writeW << endl;
writeW << "\n Bonds \n"<<endl; 
for(int i=0;i<nbond;i++) {
writeW << i+1 <<" "<< 1 <<" "<< bond[i][0]+1<<" "<<bond[i][1]+1 << endl;
}

writeW << "\n Angles \n"<<endl; 
for(int i=0;i<nangle;i++) writeW << i+1 <<" "<< 1 <<" "<< angle[i][0]+1<<" "<<angle[i][1]+1 <<" "<< angle[i][2]+1<< endl;

writeW.close();

return 0;
}


