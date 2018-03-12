/* DFT code. Copyright Nigel Wilding Nigel.Wilding@bristol.ac.uk

This code uses Classical Density Functional Theory to find the density profile of a fluid at a planar wall for a prescribed bulk density.
The fluid potential is written as the sum of a hard part and a Lennard-Jones-like attractive part with truncated interactions. The hard part is 
treated by Rosenfeld's functional measure theory with a choice of either the original form, or the White Bear version 
(See Roth J. Phys. Condensed Matter, 22, 063102 (2010) for more details)  The attractive part is treated in mean field via the random phase approximation. 
The code also allows calculation of the local compessibility profile, the adsorption  and the surface tension. 
Checking capability includes comparison with the pressure sum rule and the Gibbs adsoprtion theorem. NBW September 2017 */
 
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define MOD(n, N) ((n<0)? N+n : n)

#define PI 3.14159265359
#define PI_4 12.5663706144
#define N 10000  //Number of grid points
#define BARRIER 500 //Height of hard wall potential
#define R 0.5  // Particle diameter sigma=1.  sets relationship between rho and eta 
#define LJCUT 2.5  //Truncation radius for the Lennard-Jones potentials
#define dz 0.005   //Grid spacing
#define epsilon 1.0 // Sets the unit of energy
#define drho 0.00001 // Used for numerical compressibility derivative 
#define TOL 1e-15  //Tolerance on convergence of density distribution
#define MAXITER 40000  //Maximum number of iterations

//#define WHITEBEAR  //Switch to use WHite Bear Functional
#define ROSENFELD //Switch to use Rosenfeld Functional
//#define SHIFTEDWALL  //This will allow use of a 9-3 wall potential which is shifted so that its minimum is at the hard wall
//#define MUDIFF // This calculates two derivatives with respect to mu: the compressibility d\rho(z)/d\mu and the gibbs adsorption: -\d\gamma/\mu

//#define DIAG //Uncomment this line to get diagnostic information
#define LJ  //Turns on truncated Lennard-Jones-like fluid-fluid interactions
#define LR //Turns on the 9-3 wall-fluid potential

/* Function definitions */

void setVext(),setphiatt(),maken0(),initrho(),setwhts();
void make_dn0(),make_dn1(),make_dn2(),make_dn3(),make_dn1v(),make_dn2v();
void write_rho(),rs_convl();
void pres_profile();
double calcplanepot(double);
double pressure(), chempot(), sumrule(), omega(), adsorption();

//Declare global 1D storage.

double rho[N], rhonew[N], rhokeep[N], Vext[N];
double w0[N], w1[N], w2[N], w3[N], w1v[N], w2v[N], w1vn[N], w2vn[N];
double n0[N], n1[N], n2[N], n3[N], n1v[N], n2v[N];
double dn0[N], dn1[N], dn2[N], dn3[N], dn1v[N], dn2v[N];
double c0[N], c1[N], c2[N], c3[N], c1v[N], c2v[N], dcf[N];
double phiatt[N],cphiatt[N];
double phi[N],phiid[N];
double presprof[N],planepot[N];
double mu,new_mu,dmu,alpha,max,z,rhob,etab,ew,Rsqr;
double p,old_gamma,new_gamma;
double T,invT,rmin;

int isign=1,iter=1,isweep, iend;
int NiR=R/dz; //Number of grid points within particle radius
int NiW=R/dz; //Number of grid points within wall
int NiRCUT=LJCUT/dz; //Number of grid points within LJ cutoff

char rholive[120],runcode[120];
FILE *fpout, *fpVext, *fprholive, *fpwdens, *fpdirect, *fpdiag, *fpwhts, *fpfdivs, *fprho, *fprhonew, *fpdcf, *fptest;

int main(int argc,char *argv[])
{
  int i,j, converged = 0;
  double diff,dpz;
	
//{{{ Read in parameters
#ifdef LR
if(argc!=6) 
   {
     printf("\n Usage: %s [ T alpha ew rhob runcode\n", argv[0]);
     exit(1);
   }
T  = atof(argv[1]);
alpha = atof(argv[2]);
ew = atof(argv[3]);
rhob = atof(argv[4]); // Setting this also sets the chemical potential and pressure- see below
strcpy(runcode,argv[5]);  
strcat(runcode,".dat");
fpout = fopen(runcode,"w");
if(!strcmp(argv[5], "stderr")) *fpout=*stderr;
#else
if(argc!=5) 
   {
     printf("\n Usage: %s [T alpha rhob runcode\n", argv[0]);
     exit(1);
   }
T  = atof(argv[1]);
alpha = atof(argv[2]);
rhob = atof(argv[3]);
strcpy(runcode,argv[4]);  
strcat(runcode,".dat");
fpout = fopen(runcode,"w");
if(!strcmp(argv[4], "stderr")) *fpout=*stderr;
#endif
//}}}

mu=chempot(); //The user's choice of the bulk density sets the chemical potential and pressure.
p=pressure();
invT = 1.0/T;
Rsqr = R*R;
etab=rhob*PI/6.;
 
//{{{ Messages about run

printf("\nDFT for a fluid in planar geometry: NBW 2018\n");
printf("\nState parameters:\n  rho_b=%f\n  eta_b=%12.10f\n  mu=%12.10f\n  Pressure=%f\n  Temperature=%f\n  Inverse Temperature=%f\n\n",rhob,etab,mu,p,T,invT);
fprintf(fpout,"rho_b=%f\neta_b=%12.10f\nmu=%f\nPressure=%f\nTemperature=%f\ninverse Temperature=%f\n",rhob,etab,mu,p,T,invT);
printf("Model parameters:\n"); 
#ifdef ROSENFELD
  printf("  ROSENFELD Functional\n");fprintf(fpout,"ROSENFELD Functional\n");
#endif

#ifdef WHITEBEAR
  printf("  White Bear Functional\n");fprintf(fpout,"White Bear Functional\n");
#endif
#ifdef LJ
 printf("  LJ system\n"); fprintf(fpout,"LJ system\n");
#else
  printf("  HS system\n"); fprintf(fpout,"HS system\n");
#endif
#ifdef LR
 printf("  Wall-fluid potential switched ON: e_w=%f \n",ew); fprintf(fpout,"Wall-fluid potential switched ON: e_w=%f \n",ew);
#else 
 printf("  Hard wall potential\n");fprintf(fpout,"Hard wall potential\n");
#endif
#ifdef DIAG
 printf("  Diagnostics switched on\n");
#endif
printf("  dz=%f\n  N=%i\n  System Size(N*dz) =%f\n  NiR=%i\n  NiW=%i\n  mixing (alpha)=%f\n",dz,N,dz*N,NiR,NiW,alpha);
fprintf(fpout,"dz=%f\nN=%i\nSystem size=%f\nNiR=%i\nNiW=%i\nmixing (alpha)=%f\n",dz,N,dz*N,NiR,NiW,alpha);

//}}}

//{{{ Open diagnostic files

#ifdef DIAG
printf("Opening diagnostic files\n");
fpwdens = fopen("Diag/wdens","w"); if(fpwdens == NULL) {printf("Could not open diagnostic directory\n");exit(0);}
fpVext = fopen("Diag/Vext","w");
fpwhts = fopen("Diag/weights","w");
fpdirect = fopen("Diag/direct","w");
fprho = fopen("Diag/rho","w");
fprhonew = fopen("Diag/rhonew","w");
fpfdivs = fopen("Diag/fdivs","w");
fpdcf = fopen("Diag/dcf","w");
//printf("Fiished opening diagnostic files\n");
#endif
//}}}

#ifdef LJ
rmin=pow(2.0,1./6);
printf("  rmin=%f\n  NiRCUT=%i\n\n",rmin,NiRCUT);
fprintf(fpout,"rmin=%f\nNiRCUT=%i\n\n",rmin,NiRCUT);
#endif

setVext();  //Set the external (wall) potential

initrho(); //Set the initial density distribution and print it out

setwhts(); //Set the weights for planar geometry (see Roth 2010)

#ifdef MUDIFF  // If doing a derivative with respect to mu, change the density by a small amount and recaulate everything.
for(isweep=0;isweep<2;isweep++)
{
if(isweep>0)
{
	for(i=0;i<N;i++) rhokeep[i]=rho[i];	
    rhob+=drho;
    new_mu=chempot();
    dmu=new_mu-mu;
    mu=new_mu;
    p=pressure();
    printf("\nCalculating compressibility profile and performing Gibb's theorem check: drho=%f, dmu=%12.10f\n",drho,dmu);
    printf("\ew rhob=%f, New mu=%f, new p=%f\n",rhob,mu,p);
    converged=0;
    iter=0;
    initrho();
}
#endif

#ifdef LJ
//Form the attractive contribution. 
for(i=0;i<N;i++) planepot[i]=0;
for(i=0;i<N;i++)
  {
    if(i<=NiRCUT)
      {
      planepot[i]=calcplanepot(i*dz);
      if(i>0) planepot[N-i]=planepot[i];
      }
  }
#endif

#ifdef LJ // We stop calculating before we get to the end of the grid to avoid it acting like a second wall
  iend=N-3*NiRCUT;
#else
  iend=N-3*NiR;
#endif
 
printf("Starting iteration.....\n"); // Start picard iteration 

while(converged==0 && iter++ <MAXITER)  {

rs_convl(rho,w2,n2,NiR); //Get the weighted densities
rs_convl(rho,w3,n3,NiR);
rs_convl(rho,w2v,n2v,NiR);

for(i=0;i<N;i++)  { n0[i]=n2[i]/(PI_4*Rsqr); n1[i]=n2[i]/(PI_4*R); n1v[i]=n2v[i]/(PI_4*R); } //Others are simple related to n2,n3,n2v
 
//Now get the terms in the direct correlation function

make_dn0(); make_dn1(); make_dn2();
make_dn3(); make_dn1v(); make_dn2v();

for (i=N-2*NiR;i<N;i++) {  dn0[i]=0.0;  dn1[i]=0.0; dn2[i]=0.0; dn3[i]=0.0; dn1v[i]=0.0; dn2v[i]=0.0; } //Make sure convolutions don't see the far wall.

rs_convl(dn0,w0,c0,NiR);
rs_convl(dn1,w1,c1,NiR);
rs_convl(dn2,w2,c2,NiR);
rs_convl(dn3,w3,c3,NiR);
rs_convl(dn1v,w1vn,c1v,NiR); //NB the vector weights are odd, so use the negative- see Roth 2010
rs_convl(dn2v,w2vn,c2v,NiR);

// Now form the direct correlation function (dcf)

for(i=0;i<N;i++) dcf[i]=-(c0[i]+c1[i]+c2[i]+c3[i]+c1v[i]+c2v[i]);

#ifdef LJ
//Form the attractive contribution. 
rs_convl(rho,planepot,cphiatt,NiRCUT);

//use dcf to form the new density distribution with mixing to damp down the changes
 
 for(i=NiR;i<iend;i++) rhonew[i]=(1-alpha)*rho[i] + alpha*exp(invT*(mu-Vext[i])+dcf[i]-invT*cphiatt[i]);
#else
 for(i=NiR;i<iend;i++) rhonew[i]=(1-alpha)*rho[i] + alpha*exp(invT*(mu-Vext[i])+dcf[i]);
#endif

max=0.0; diff=0;                                                                                    
for(i=0;i<iend;i++) 
  {
    diff=fabs(rhonew[i]-rho[i]);
    if(diff>max) max=diff;
  }

if(max<TOL) converged=1;

if(iter%50==0) {printf("Iteration %i Max=%18.16f\n",iter++,max);write_rho();}

//{{{ Write out diagnostics
#ifdef DIAG
for(i=0;i<N;i++)  fprintf(fpwdens,"%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n",i*dz,n0[i],n1[i],n2[i],n3[i],n1v[i],n2v[i]);
for(i=0;i<N;i++)  fprintf(fpfdivs,"%f %f %f %f %f %f %f\n",i*dz,dn0[i],dn1[i],dn2[i],dn3[i],dn1v[i],dn2v[i]);
for(i=0;i<N;i++)  fprintf(fpdirect," %f %f %f %f %f %f %f\n",i*dz,c0[i],c1[i],c2[i],c3[i],c1v[i],c2v[i]);
for(i=0;i<N;i++)  fprintf(fpdcf," %f %f\n",i*dz,dcf[i]);
for(i=0;i<N;i++)  fprintf(fprho," %f %f\n",i*dz,rho[i]);
for(i=0;i<N;i++)  fprintf(fprhonew," %f %f %f %f\n",i*dz,rhonew[i],alpha*exp(invT*(mu-Vext[i])+dcf[i]),exp(invT*(mu-Vext[i])));
#endif //}}}

// Copy the new density profile to the old one.

for(i=0;i<iend;i++) rho[i]=rhonew[i];  // Note we fix the 'bulk' density far from the wall

} //end of picard loop

if(converged==1) printf("Converged to within tolerance in %i iterations\n",iter);
else printf("Failed to converge after %i iterations\n",MAXITER);

#ifdef MUDIFF
if(isweep==0) {
	old_gamma = omega();
	printf("gamma1= %12.10f\n",old_gamma);
}
else 
{
	new_gamma = omega();
    printf("gamma2= %12.10f\n",new_gamma);
}
}
printf("-d(gamma)/dmu= %f\nadsorption= %f\n",-(new_gamma-old_gamma)/dmu,adsorption());
for(i=0;i<iend;i++) {z=(i-NiW)*dz; if(rhokeep[i]>1e-8) fprintf(fpout,"A %f %12.10f %12.10f %12.10f %12.10f\n",z,rhokeep[i],rhokeep[i]/rhob,rhokeep[i]*PI/6,(rho[i]-rhokeep[i])/dmu);}
#else
printf("gamma= %12.10f\nadsorption= %12.10f\n",omega(),adsorption());
for(i=0;i<iend;i++) {z=(i-NiW)*dz; if(rho[i]>1e-8) fprintf(fpout,"A %f %12.10f %12.10f %12.10f \n",z,rho[i],rho[i]/rhob,rho[i]*PI/6);}
#endif

#ifdef LR
printf("Sum rule pressure: %f\n",sumrule());fprintf(fpout,"Sum rule pressure: %f\n",sumrule());
#else
printf("Hard wall k_BT*Contact density = %f \n",T*rho[NiW]);fprintf(fpout,"Hard wall k_BT*Contact density = %f \n",T*rho[NiW]);
#endif

  	  
} // End of main program

//{{{ setVext
void setVext() //Wall-fluid potential
{
int i;
double z,zwall,zt,zt3i,zt9i;

  for(i=0;i<N;i++)
    {
      Vext[i]=0.0;
      if(i<NiW) Vext[i]=BARRIER; 

#ifdef LR
if(i==NiW) Vext[i]=BARRIER; // We make V(0) "hard" for the LR case because the 9-3 potential diverges at this point
zwall=(i-NiW)*dz;
zt=zwall;
#ifdef SHIFTEDWALL // use this to have the minimum of the potential at the hard wall
zt=zwall+pow(0.4,1./6);
#endif

if(zwall>0) 
 {
 	 zt3i = 1./(zt*zt*zt);
 	 zt9i = zt3i*zt3i*zt3i;
     Vext[i] = epsilon*ew*((2.0/15.0)*zt9i-zt3i);
     if(Vext[i]>BARRIER) Vext[i]=BARRIER;
}
#endif
    }
    
#ifdef DIAG 
  for(i=0;i<N;i++)  fprintf(fpVext," %f %f\n",(i-NiW)*dz,Vext[i]);
  fclose(fpVext);
//  exit(0);
#endif
} //}}}

//{{{ initrho            
void initrho()  //set the density initially to be the bulk density and write it out
{
  int i;
  for(i=0;i<N;i++)  rho[i]=exp(-Vext[i]);
  for(i=0;i<N;i++) if(Vext[i]<10) rho[i]=rhob;
  write_rho(); 
} //}}}

//{{{ setwhts
void setwhts()
/* These are the weights for planar geometry (see Roth 2010)  They represent the 'response function' in the convolution.
Note that the negative parts appear in wrap-around order ie at the end of the array */

{
  int i,ind;
double z;
 for(i=0;i<N;i++) {w3[0]=0; w2[i]=0; w1[i]=0; w0[i]=0; w2v[i]=0;  w1v[i]=0; w1vn[i]=0; w2vn[i]=0;}

for(i=0;i<=NiR;i++) // The heaviside function is unity when it's argument is zero or greater. 
  {
       z = i*dz;  	  
	   w3[i]   =  PI*(R*R-z*z)*dz;      	 	 
	   w2[i]   =  2*PI*R*dz;              
	   w2v[i]  =  2*PI*z*dz;            
	   w2vn[i] =  -2*PI*z*dz;                      // Negative for the convolution to get direct correlation function, since the vectorial weight functions are odd
	   w1[i]   =  w2[i]/(4*PI*R);     
	   w1v[i]  =  w2v[i]/(4*PI*R);   
	   w1vn[i] =  -w2v[i]/(4*PI*R); 
	   w0[i]   =  w2[i]/(4*PI*R*R);   

	   if(i>0) 
	   {
	   	 z = -z;  
	   	 ind = N-i;
	   	 w3[ind]   =  PI*(R*R-z*z)*dz;
	   	 w2[ind]   =  2*PI*R*dz; 
	   	 w2v[ind]  =  2*PI*z*dz;
	   	 w2vn[ind] =  -2*PI*z*dz;
	   	 w1[ind]   =  w2[ind]/(4*PI*R);
	   	 w1v[ind]  =  w2v[ind]/(4*PI*R);
	   	 w1vn[ind] =  -w2v[ind]/(4*PI*R);
	   	 w0[ind]   =  w2[ind]/(4*PI*R*R);
	   }
	   
  }

#ifdef WHITEBEAR // This stops n3 becomng zero which blows up weighted densities
w3[NiR] = 1e-6; 
w3[N-NiR] = 1e-6;
#endif

#ifdef DIAG
for(i=0;i<N;i++)  fprintf(fpwhts,"%i %f %f %f %f %f %f %f %f %f\n",i,i*dz,w0[i],w1[i],w2[i],w3[i],w1v[i],w1vn[i],w2v[i],w2vn[i]);
fclose(fpwhts);
#endif

}
//}}}
	
//{{{ real space convolution
void rs_convl(const double *input, const double *response, double *output, int HALFWIDTH)   // real_space discrete convolution. 

{
	
  int i,j,k;
  int nterms;
  double this,next,term;
  double store[N];

  if(HALFWIDTH==NiR) //The weight functions have discontinities so use the extended closed formula from Numerical Recipes (eq 4.1.14).
  {	  
  	  for(i=0;i<N;i++)
  	  {
  	  	  k=0; 
  	  	  output[i]=0;
  	  	  nterms=2*HALFWIDTH+1;
  	  	  if(i<HALFWIDTH) nterms-=HALFWIDTH-i;
  	  	  for(j=i-HALFWIDTH;j<=i+HALFWIDTH;j++)
  	  	  {
  	  	  	  if(j>=0 && j<N-1)
  	  	  	  {
  	  	  	  	  term= input[j] * response[MOD(i-j, N)];
  	  	  	  	  store[k++]=term;              
  	  	  	  }
  	  	  }
  	  	  store[0]*=3./8; store[1]*=7./6; store[2]*=23./24;
  	  	  store[nterms-1]*=3./8; store[nterms-2]*=7./6; store[nterms-3]*=23./24; 
  	  	  for(k=0;k<nterms;k++) output[i]+=store[k];
  	  }
  	  
  }
  else
   {
     //This seems to work better for the smoother LJ convolution
  	  for(i=0;i<N;i++)
  	  {
  	  	  output[i]=0;   
  	  	  if(i-HALFWIDTH>=0) this=input[i-HALFWIDTH]*response[MOD(HALFWIDTH, N)];
  	  	  for(j=i-HALFWIDTH;j<i+HALFWIDTH;j++)
  	  	  {
  	  	  	  if(j>=0 && j<N-1)
  	  	  	  {
  	  	  	  	  next=   input[j+1] * response[MOD(i-j-1, N)];              
  	  	  	  	  output[i]+=(this+next)/2.;
  	  	  	  	  this=next;
  	  	  	  }
  	  	  }
  	  }
  }
 
 /* 
   //Rectangular rule.  This is a simpler but less accurate alternative - not used
  for(i=0;i<N;i++)
  {
   output[i]=0;
   for(j=i-HALFWIDTH;j<i+HALFWIDTH;j++)  
   if(j>=0 && j<N)  output[i] += input[j] * response[MOD(i-j, N)];

  }
*/
 
 
} //}}}

//{{{ make dn
void make_dn0() //Functional derivatives of the excess free energy
{
	int i;
    for(i=0;i<N;i++) dn0[i] = -log( 1 - n3[i] );
}
	
void make_dn1()
{
	int i;
    for(i=0;i<N;i++) dn1[i] = n2[i]/(1-n3[i]);
}	
	
void make_dn1v()
{
	int i;
        for(i=0;i<N;i++) dn1v[i] = -n2v[i]/(1-n3[i]);

}

void make_dn2()
{
	int i;
#ifdef WHITEBEAR
    for(i=0;i<N;i++) dn2[i] = n1[i]/(1-n3[i]) + \
    	(1/36.0) * (3*n2[i]*n2[i] -  3*n2v[i]*n2v[i] ) * ( n3[i] + (1-n3[i])*(1-n3[i])*log(1-n3[i]) ) / (PI*n3[i]*n3[i]*(1-n3[i])*(1-n3[i])); 
#endif

#ifdef ROSENFELD
    for(i=0;i<N;i++) dn2[i] = n1[i]/(1-n3[i]) + \
    	(1/24.0) * (3*n2[i]*n2[i] -  3*n2v[i]*n2v[i] )  / (PI*(1-n3[i])*(1-n3[i])); 

#endif
}

void make_dn2v()
{
	int i;

#ifdef WHITEBEAR
	for(i=0;i<N;i++) dn2v[i] = -n1v[i]/(1-n3[i]) -\
	      (1/6.0)*n2[i]*n2v[i]*(n3[i] + (1-n3[i])*(1-n3[i])*log(1-n3[i]) )/(PI*n3[i]*n3[i]*(1-n3[i])*(1-n3[i]));
#endif

#ifdef ROSENFELD
	for(i=0;i<N;i++) dn2v[i] = -n1v[i]/(1-n3[i]) -\
	      (1/4.0)*n2[i]*n2v[i]/(PI*(1-n3[i])*(1-n3[i]));

#endif
}

void make_dn3()
{
	int i;
	
#ifdef WHITEBEAR
	for(i=0;i<N;i++) dn3[i] = n0[i]/(1-n3[i])+\
	(n1[i]*n2[i]-n1v[i]*n2v[i])/((1-n3[i])*(1-n3[i]))+\
    (1/36.0)*(n2[i]*n2[i]*n2[i]-3*n2[i]*n2v[i]*n2v[i]) * (-(2*(1-n3[i]))*log(1-n3[i])+n3[i]) / (PI*n3[i]*n3[i]*(1-n3[i])*(1-n3[i]))-\
	(1/18.0)*(n2[i]*n2[i]*n2[i]-3*n2[i]*n2v[i]*n2v[i]) * (n3[i]+(1-n3[i])*(1-n3[i])*log(1-n3[i])) / (PI*n3[i]*n3[i]*n3[i]*(1-n3[i])*(1-n3[i]))+\
	(1/18.0)*(n2[i]*n2[i]*n2[i]-3*n2[i]*n2v[i]*n2v[i]) * (n3[i]+(1-n3[i])*(1-n3[i])*log(1-n3[i])) / (PI*n3[i]*n3[i]*(1-n3[i])*(1-n3[i])*(1-n3[i]));
#endif

#ifdef ROSENFELD
	for(i=0;i<N;i++) dn3[i] = n0[i]/(1-n3[i])+\
	(n1[i]*n2[i]-n1v[i]*n2v[i])/((1-n3[i])*(1-n3[i]))+\
    (1/12.0)*(n2[i]*n2[i]*n2[i]-3*n2[i]*n2v[i]*n2v[i])/ (PI*(1-n3[i])*(1-n3[i])*(1-n3[i]));

#endif

} //}}}

//{{{ write_rho
		void write_rho()
	{
		int i;
	    fprholive=fopen("rholive","w");
	    for(i=0;i<N;i++) fprintf(fprholive,"%f %f \n",i*dz,rho[i]);
	    fclose(fprholive);
    }


//}}}

//{{{ calcplanepot
double calcplanepot(double dp) //Calculates the potential at a point due to infinite plane of particles
{
  double pot;
  if(dp>rmin)                                            
	 pot=0.4/pow(dp,10.0)-1./pow(dp,4.0) -0.4/pow(LJCUT,10.0)+1./pow(LJCUT,4.0);  
  else
	 pot=0.5*(dp*dp-rmin*rmin)+0.4/pow(rmin,10.0)-1./pow(rmin,4.0) -0.4/pow(LJCUT,10.0)+1./pow(LJCUT,4.0);	   
  return(2*PI*dz*pot);
} //}}}
       
//{{{ pressure
double pressure()
/* The pressure for a given bulk rho; White Bear uses BMCSL expression for the pressure (see Roth 2010) 
   Rosenfeld uses the PY pressure */

{
    int i;
    double p,r,eta;

    eta = PI*rhob/6.;
    
#ifdef WHITEBEAR
        p=T*rhob*(1+eta+eta*eta-eta*eta*eta)/pow(1-eta,3.0);
#endif

#ifdef ROSENFELD
        p=T*rhob*(1+eta+eta*eta)/pow(1-eta,3.0);
#endif

#ifdef LJ
        p+=-2*PI*rhob*rhob*1.171861897*epsilon;
#endif

return(p);

} //}}}

//{{{ chempot
double chempot()
/* The pressure for a given bulk rho; White Bear uses BMCSL expression (see Roth 2010) 
   Rosenfeld uses the PY expression */

{
    int i;
    double c,r,eta;

    eta = PI*rhob/6.;
    
#ifdef WHITEBEAR
        c = T * ( log(rhob)+(8*eta - 9*eta*eta + 3*eta*eta*eta)/pow(1-eta,3.0) );
#endif

#ifdef ROSENFELD
        c=  T * (log(rhob) + (14*eta - 13*eta*eta + 5*eta*eta*eta)/(2*pow(1-eta,3.0)) -log(1-eta) );
        
#endif

#ifdef LJ
        c+=-4*PI*rhob*1.171861897*epsilon;
#endif

return(c);

} //}}}

//{{{ sumrule
double sumrule() //Estimates the pressure via the sumrule for a continuous potential, see Roth (2011) eq 26.
{
double sum=0,Vextdiff,zwall;
int i;
fptest = fopen("sumrule","w");
for (i=0;i<N-NiRCUT;i++) //Do a trapezium rule for the integration.
{
	zwall=(i-NiW)*dz;
    if(zwall>0) 
      {
    	Vextdiff=ew*epsilon*((-6.0/5.)*pow(zwall,-10.0)+3.0*pow(zwall,-4.0));
        sum+= (-1.0*rho[i]*Vextdiff-1.0*rho[i+1]*ew*epsilon*((-6.0/5.)*pow(zwall+dz,-10.0)+3.0*pow(zwall+dz,-4.0)))/2;
        if(zwall<20)  fprintf(fptest,"S %f %9.7f %9.7f %9.7f\n",zwall,rho[i],Vextdiff,rho[i]*Vextdiff);
      }
}
sum*=dz;
return(sum);
} //}}}

//{{{ omega

double omega()  

//Calculate the grand potential density omega and integrate it to get the surface tension, which is returned
 {
   	int i,end;
    double sumid,sumphi;

    for(i=0;i<N;i++) {phi[i]=0; phiid[i]=0.0;}
    for(i=0;i<N-2*NiR;i++)
       {
        phi[i] = -n0[i]*log( 1-n3[i] ) + ( n1[i]*n2[i]-n1v[i]*n2v[i] ) / ( 1-n3[i] );
#ifdef WHITEBEAR
        phi[i] += ( n2[i]*n2[i]*n2[i] -3*n2[i]*n2v[i]*n2v[i] )  * ( n3[i]+(1-n3[i])*(1-n3[i])*log(1-n3[i]) ) / ( 36*PI*n3[i]*n3[i]*(1-n3[i])*(1-n3[i]) );
#endif
        	
#ifdef ROSENFELD
       	phi[i] += ( n2[i]*n2[i]*n2[i] -3*n2[i]*n2v[i]*n2v[i] )   / ( 24*PI*(1-n3[i])*(1-n3[i]) );
#endif
        phi[i]*=T;
#ifdef LJ
phi[i] += 0.5*rho[i]*cphiatt[i];  
#endif
       phi[i]+=p;       
       if(i>=NiR) phiid[i-NiR] = T*rho[i]*(log(rho[i])-1.0) + rho[i]*(Vext[i]-mu);     // NB we have to shift the ideal gas contribution "inside the wall" before we add it in 
	   }

  sumid=0.0;sumphi=0.0;
  for (i=0;i<iend-1;i++) sumphi+=phi[i]+phi[i+1];        //Due to shift above, need to integrate over the same length. phi is not zero inside the wall.//
  for (i=0;i<iend-NiR-1;i++) sumid+=phiid[i]+phiid[i+1];
    
  return dz*(sumid+sumphi)/2.0;

 }
//}}}

//{{{ adsorption
double adsorption()
// Note this definition assumes that the zero of z starts at -R. This leads to a negative adsorption.
{
int i,end;
double a=0.0;

for (i=0;i<iend;i++) a+=rho[i]+rho[i+1]-2*rhob;
return (dz*a/2.);
} //}}}






