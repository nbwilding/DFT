/* DFT code. Copyright Nigel Wilding Nigel.Wilding@bristol.ac.uk

This code uses Classical Density Functional Theory to find the density
profile of a fluid at a planar wall for a prescribed bulk density.
The fluid potential is written as the sum of a hard part and a
Lennard-Jones-like attractive part with truncated interactions. The
hard part is treated by Rosenfeld's functional measure theory with a
choice of either the original form, or the White Bear version (See
Roth J. Phys. Condensed Matter, 22, 063102 (2010) for more details)
The attractive part is treated in mean field.  The code also allows calculation of the local
compessibility profile, the adsorption and the surface tension.
Checking capability includes comparison with the pressure sum rule and
the Gibbs adsoprtion theorem. NBW March 2018

*/
 
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MOD(n, N) ((n<0)? N+n : n)

#define PI 3.14159265359
#define PI_4 12.5663706144
#define N 10000            // Total number of grid points
#define BARRIER 500        // Height of hard wall potential
#define R 0.5              // Particle diameter sigma=1.  sets relationship between rho and eta 
#define LJCUT 2.5          // Truncation radius for the Lennard-Jones potentials
#define dz 0.005           // Grid spacing. Should be a factor of R above.
#define epsilon 1.0        // Sets the unit of energy. Don't change this.
#define drho 0.0000001     // Density step size for calculating derivatives w.r.t. mu 
#define TOL 1e-12          // Tolerance on convergence of density distribution
#define MAXITER 800000     // Maximum number of iterations
#define NGFREQ 3           // Ratio of Picard to Ng algorithm updates. Can be set as low as 3, but if in doubt set to 1000
//#define WHITEBEAR        // Switch to use White Bear Functional
#define ROSENFELD          // Switch to use Rosenfeld Functional
#define LJ                 // Turns on truncated Lennard-Jones-like fluid-fluid interactions
#define LR                 // Turns on the 9-3 wall-fluid potential
//#define MUDIFF           // Uncomment to calculate derivatives with respect to mu: the compressibility d\rho(z)/d\mu and the gibbs adsorption: -\d\gamma/\mu
//#define SHIFTEDWALL      // Use a 9-3 wall potential which is shifted so that its minimum is at the hard wall
//#define DIAG             // Uncomment this line to get diagnostic information written to files in a separate directory "Diag"
//#define READRHO          // Uncomment this line to read in an existing density profile as a starting guess

/* Function definitions */

void setVext(), setphiatt(), initrho(), setwhts();
void make_dn0(), make_dn1(), make_dn2(), make_dn3(), make_dn1v(), make_dn2v();
void write_rho(), rs_convl(), update();

double calcplanepot(double), omega(int);
double pressure(), chempot(), sumrule(), adsorption();

//Declare global 1D storage.

double rho[N], rhonew[N], rhokeep[N], Vext[N];
double g[N], g1[N], g2[N], d[N], d1[N], d2[N], d01[N], d02[N]; //Used by update function
double w0[N], w1[N], w2[N], w3[N], w1v[N], w2v[N], w1vn[N], w2vn[N];
double n0[N], n1[N], n2[N], n3[N], n1v[N], n2v[N];
double dn0[N], dn1[N], dn2[N], dn3[N], dn1v[N], dn2v[N];
double c0[N], c1[N], c2[N], c3[N], c1v[N], c2v[N], dcf[N];
double phiatt[N], cphiatt[N];
double phi[N], phiid[N], planepot[N];

//Other global variables

double mu,new_mu,dmu,alpha,z,rhob,etab,ew,Rsqr;
double p,old_gamma,new_gamma;
double T,invT,rmin,dev;

int iter=1, isweep, iend;
int NiR=R/dz;           // Number of grid points within particle radius
int NiW=R/dz;           // Number of grid points within wall
int NiRCUT=LJCUT/dz;    // Number of grid points within LJ cutoff

// Files. Mainly used for diagnostic output
char rholive[120],runcode[120];
FILE *fpout, *fpVext, *fprholive, *fpwdens, *fpdirect, *fpdiag, *fpwhts, *fpfdivs, *fprho, *fprhonew, *fpdcf, *fptest, *fprhostart;

int main(int argc,char *argv[])
{
  int i,j, converged = 0;

//{{{ Read in parameters
#ifdef WHITEBEAR
#ifdef ROSENFELD
printf("\nChoose EITHER WHITEBEAR or ROSENFELD DFT\n");
exit(0);
#endif
#endif
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

mu = chempot(); //The user's choice of the bulk density sets the chemical potential and pressure.
p = pressure();
invT = 1.0/T;
Rsqr = R*R;
etab = rhob*PI/6.;
 
//{{{ Messages about run.  Note "comments" of this form arise from the folding capabilities of some editors. I use jedit.

printf("\nDFT for a fluid in planar geometry: NBW 2018\n");
printf("\nState parameters:\n  rho_b= %12.10f\n  eta_b= %12.10f\n  mu= %12.10f\n  Pressure= %12.10f\n  Temperature= %10.8f\n  Inverse Temperature= %f\n\n",rhob,etab,mu,p,T,invT);
fprintf(fpout,"rho_b= %12.10f\neta_b= %12.10f\nmu= %12.10f\nPressure= %12.10f\nTemperature= %10.8f\nInverse Temperature= %10.8f\n",rhob,etab,mu,p,T,invT);
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
  printf("  Wall-fluid potential switched ON\n  e_w= %lg \n",ew); fprintf(fpout,"Wall-fluid potential switched ON\ne_w= %lg \n",ew);
#else 
  printf("  Hard wall potential\n");fprintf(fpout,"Hard wall potential\n");
#endif

#ifdef MUDIFF
printf("  Measuring mu derivatives for compressibility profile and Gibbs adsorption\n");fprintf(fpout,"Measuring mu derivatives for compressibility profile and Gibbs adsorption\n");
#endif
#ifdef DIAG
 printf("  Diagnostic file output switched on\n");
#endif
printf("  dz= %f\n  N= %i\n  System Size(N*dz)= %4.2f\n  NiR=%i\n  NiW=%i\n  mixing (alpha)=%4.2f\n  NGFREQ=%i\n  Tolerance=%lg\n",dz,N,dz*N,NiR,NiW,alpha,NGFREQ,TOL);
fprintf(fpout,"dz= %f\nN= %i\nSystem size= %4.2f\nNiR= %i\nNiW= %i\nmixing (alpha)= %4.2f\nNGFREQ=%i\nTolerance=%lg\n",dz,N,dz*N,NiR,NiW,alpha,NGFREQ,TOL);
//}}}

//{{{ Open diagnostic files. 

/* These write values of various functions to files in a sub directory "Diag" which the user should create. 
   Note these files get large so only run for a small number of iterations
*/

#ifdef DIAG
fpwdens = fopen("Diag/wdens","w"); if(fpwdens == NULL) {printf("Could not open diagnostic directory\n Pease create a subdirectory called 'Diag' ");exit(0);}
fpVext = fopen("Diag/Vext","w");
fpwhts = fopen("Diag/weights","w");
fpdirect = fopen("Diag/direct","w");
fprho = fopen("Diag/rho","w");
fprhonew = fopen("Diag/rhonew","w");
fpfdivs = fopen("Diag/fdivs","w");
fpdcf = fopen("Diag/dcf","w");
//printf("Finished opening diagnostic files\n");
#endif
//}}}

#ifdef LJ
 rmin = pow(2.0,1./6); //Location of LJ minimum
printf("  rmin=%f\n  NiRCUT=%i\n\n",rmin,NiRCUT);
fprintf(fpout,"rmin=%f\nNiRCUT=%i\n\n",rmin,NiRCUT);
#endif

setVext();  // Set the external (wall) potential

initrho();  // Set the initial density distribution and print it out

setwhts();  // Set the weights for planar geometry (see Roth 2010)

#ifdef MUDIFF  // If doing a derivative with respect to mu, change the density by a small amount and recalculate everything.

for(isweep=0; isweep<2; isweep++)
{
if(isweep>0)
{
    for(i=0;i<N;i++) rhokeep[i]=rho[i];	
    rhob += drho;
    new_mu = chempot();
    dmu = new_mu-mu;
    mu = new_mu;
    p = pressure();
    printf("\nCalculating compressibility profile and performing Gibb's theorem check: drho=%10.8f, dmu=%12.10f\n",drho,dmu);
    printf("New rhob=%12.10f, New mu=%12.10f, new p=%10.8f\n",rhob,mu,p);
    converged = 0;
    iter = 0;
    initrho();
}
#endif

#ifdef LJ
//Form the attractive contribution
for(i=0;i<N;i++) planepot[i]=0;
for(i=0;i<N;i++)
  {
    if(i<=NiRCUT)
      {
      planepot[i] = calcplanepot(i*dz);
      if(i>0) planepot[N-i] = planepot[i];
      }
  }
planepot[NiRCUT]*=3./8;   planepot[NiRCUT-1]*=7./6;    planepot[NiRCUT-2]*=23./24; // Apply the extended quadrature
planepot[N-NiRCUT]*=3./8; planepot[N-NiRCUT+1]*=7./6;  planepot[N-NiRCUT+2]*=23./24;
#endif

#ifdef LJ // We stop calculating before we get to the end of the grid to avoid it acting like a second wall
  iend=N-3*NiRCUT;
#else
  iend=N-3*NiR;
#endif
 
printf("Starting iteration.....\n"); // Start minimisation iteration 

while(converged==0 && iter++ <MAXITER)  {

rs_convl(rho,w2,n2,NiR); //Get the weighted densities via convolution
rs_convl(rho,w3,n3,NiR);
rs_convl(rho,w2v,n2v,NiR);

for(i=0;i<N;i++)  { n0[i]=n2[i]/(PI_4*Rsqr); n1[i]=n2[i]/(PI_4*R); n1v[i]=n2v[i]/(PI_4*R); } //Others are simple related to n2,n3,n2v
  
#ifdef WHITEBEAR //When we have a long ranged wall, we need to make sure that n3 is non-zero or the functional derivatives blow up
#ifdef LR
   for(i=0;i<N;i++) if(n3[i]<TOL) n3[i]=1e-12;
#endif 
#endif	

//Now get the terms in the direct correlation function

make_dn0(); make_dn1(); make_dn2();
make_dn3(); make_dn1v(); make_dn2v();

rs_convl(dn0,w0,c0,NiR);
rs_convl(dn1,w1,c1,NiR);
rs_convl(dn2,w2,c2,NiR);
rs_convl(dn3,w3,c3,NiR);
rs_convl(dn1v,w1vn,c1v,NiR); //NB the vector weights are odd, so use the negative- see Roth 2010
rs_convl(dn2v,w2vn,c2v,NiR);

// Now form the direct correlation function (dcf)

for(i=0;i<N;i++) dcf[i] = -( c0[i] + c1[i] + c2[i] + c3[i] + c1v[i] + c2v[i]);

#ifdef LJ
//Form the attractive contribution via a convolution (see notes file)
rs_convl(rho,planepot,cphiatt,NiRCUT);
#endif

update();   // Call the Picard-Ng update                                                                           

if(dev<TOL) converged = 1;

if(iter%50==0) {printf("Iteration %5i Deviation= %10.8lg\n",iter++,dev);write_rho();}

//{{{ Write out diagnostics
#ifdef DIAG
for(i=0;i<N;i++)  fprintf(fpwdens,"%f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n",i*dz,n0[i],n1[i],n2[i],n3[i],n1v[i],n2v[i]);
for(i=0;i<N;i++)  fprintf(fpfdivs,"%f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n",i*dz,dn0[i],dn1[i],dn2[i],dn3[i],dn1v[i],dn2v[i]);
for(i=0;i<N;i++)  fprintf(fpdirect,"%f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n",i*dz,c0[i],c1[i],c2[i],c3[i],c1v[i],c2v[i]);
for(i=0;i<N;i++)  fprintf(fpdcf,"%f %12.10f\n",i*dz,dcf[i]);
for(i=0;i<N;i++)  fprintf(fprho,"%f %12.10f\n",i*dz,rho[i]);
for(i=0;i<N;i++)  fprintf(fprhonew,"%f %12.10f %12.10f %12.10f\n",i*dz,rhonew[i],alpha*exp(invT*(mu-Vext[i])+dcf[i]),exp(invT*(mu-Vext[i])));
exit(0);  //Stops after first iteration to keep file sizes small. Remove this line to get data on all iterations
#endif //}}}

// Copy the new density profile to the old one.

for(i=0;i<iend;i++) rho[i]=rhonew[i];  // Note we fix the 'bulk' density far from the wall

} //end of iteration loop

if(converged==1) printf("Converged to within tolerance in %i iterations\n",iter);
else printf("Failed to converge after %i iterations\n",MAXITER);

#ifdef MUDIFF
if(isweep==0) {
	old_gamma = omega(1);
	printf("gamma_1= %12.10f\n",old_gamma);
}
else 
{
	new_gamma = omega(1);
    printf("gamma_2= %12.10f\n",new_gamma);
}
}
printf("-d(gamma)/dmu= %f\nadsorption= %f\n",-(new_gamma-old_gamma)/dmu,adsorption());
for(i=0;i<iend;i++) {z=(i-NiW)*dz; if(rhokeep[i]>1e-8) fprintf(fpout,"A %f %12.10f %12.10f %12.10f %12.10f\n",z,rhokeep[i],rhokeep[i]/rhob,rhokeep[i]*PI/6,(rho[i]-rhokeep[i])/dmu);}
#else
for(i=0;i<iend;i++) {z=(i-NiW)*dz; if(rho[i]>1e-8) fprintf(fpout,"B %f %12.10f %12.10f %12.10f %lg\n",z,rho[i],rho[i]/rhob,rho[i]*PI/6,d[i]);}
printf("gamma= %12.10f\nadsorption= %12.10f\n",omega(1),adsorption());fprintf(fpout,"gamma= %12.10f\nadsorption= %12.10f\n",omega(1),adsorption());
#endif


#ifdef LR
printf("Sum rule pressure: %10.8lg (%10.8lg)\n",sumrule(),p);fprintf(fpout,"Sum rule pressure: %10.8lg (%10.8lg)\n",sumrule(),p);
#else
printf("Hard wall k_BT*Contact density = %f (deviation = %10.8f)\n",T*rho[NiW],fabs(T*rho[NiW]-p));fprintf(fpout,"Hard wall k_BT*Contact density = %f (deviation = %10.8f)\n",T*rho[NiW],fabs(T*rho[NiW]-p));
#endif
  	
write_rho(); //This is written out to rholive incase we want to restart from the same profile

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

//{{{ write_rho
		void write_rho()
	{
		int i;
	    fprholive=fopen("rholive","w");
	    for(i=0;i<N;i++) fprintf(fprholive,"%f %f %lg\n",i*dz,rho[i],d[i]);
	    fclose(fprholive);
    }


//}}}


//{{{ initrho            
void initrho()  //set the density initially to be the bulk density and write it out
{
  int i;
  float dummy;
#ifdef READRHO
  printf("Reading in starting rho(z) from rholive\n");
  fprhostart=fopen("rholive","r");
  for(i=0;i<N;i++) fscanf(fprhostart,"%f %lg \n",&dummy,&rho[i]); 
  fclose(fprhostart);
#else
  for(i=0;i<N;i++)  rho[i]=exp(-Vext[i]);
  for(i=0;i<N;i++)  if(Vext[i]<10) rho[i]=rhob;
#endif
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

// Use the extended closed formula from Numerical Recipes (eq 4.1.14).

  w0[NiR]*=3./8;   w0[NiR-1]*=7./6;    w0[NiR-2]*=23./24;
  w1[NiR]*=3./8;   w1[NiR-1]*=7./6;    w1[NiR-2]*=23./24;
  w1v[NiR]*=3./8;  w1v[NiR-1]*=7./6;   w1v[NiR-2]*=23./24;
  w1vn[NiR]*=3./8; w1vn[NiR-1]*=7./6;  w1vn[NiR-2]*=23./24;
  w2[NiR]*=3./8;   w2[NiR-1]*=7./6;    w2[NiR-2]*=23./24;
  w2v[NiR]*=3./8;  w2v[NiR-1]*=7./6;   w2v[NiR-2]*=23./24;
  w2vn[NiR]*=3./8; w2vn[NiR-1]*=7./6;  w2vn[NiR-2]*=23./24;
  w3[NiR]*=3./8;   w3[NiR-1]*=7./6;    w3[NiR-2]*=23./24;

  w0[N-NiR]*=3./8;   w0[N-NiR+1]*=7./6;    w0[N-NiR+2]*=23./24;
  w1[N-NiR]*=3./8;   w1[N-NiR+1]*=7./6;    w1[N-NiR+2]*=23./24;
  w1v[N-NiR]*=3./8;  w1v[N-NiR+1]*=7./6;   w1v[N-NiR+2]*=23./24;
  w1vn[N-NiR]*=3./8; w1vn[N-NiR+1]*=7./6;  w1vn[N-NiR+2]*=23./24;
  w2[N-NiR]*=3./8;   w2[N-NiR+1]*=7./6;    w2[N-NiR+2]*=23./24;
  w2v[N-NiR]*=3./8;  w2v[N-NiR+1]*=7./6;   w2v[N-NiR+2]*=23./24;
  w2vn[N-NiR]*=3./8; w2vn[N-NiR+1]*=7./6;  w2vn[N-NiR+2]*=23./24;
  w3[N-NiR]*=3./8;   w3[N-NiR+1]*=7./6;    w3[N-NiR+2]*=23./24;
  
#ifdef WHITEBEAR // This stops n3 becoming zero which blows up weighted densities
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

   if(HALFWIDTH==NiR)
   {
  	  for(i=0;i<N;i++)
  	  {
  	  	  output[i]=0.0;
  	  	  for(j=i-HALFWIDTH;j<=i+HALFWIDTH;j++)
  	  	  {
  	  	  	  if(j>=0 && j<N) output[i]+= input[j] * response[MOD(i-j, N)];
  	  	  }
  	  }
  	  
   }  
  	    else
   {
     //Use a trapezoidal for the LJ convolution which seems to work better
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
        p+=-2*PI*rhob*rhob*1.171861897*epsilon;  //Van der Waals contribution
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
double sum=0, Vextdiff, zwall;
int i;

for (i=0;i<N-NiRCUT;i++) //Do a trapezium rule for the integration.
{
    zwall=(i-NiW)*dz;
    if(zwall>0) 
      {
    	Vextdiff=ew*epsilon*((-6.0/5.)*pow(zwall,-10.0)+3.0*pow(zwall,-4.0));
        sum+= (-1.0*rho[i]*Vextdiff-1.0*rho[i+1]*ew*epsilon*((-6.0/5.)*pow(zwall+dz,-10.0)+3.0*pow(zwall+dz,-4.0)))/2;
      }
}
sum*=dz;
return(sum);
} //}}}

//{{{ omega
//Set mode=0 to calculate the grand potential and set mode=1 to calculate the excess grand potential (ie surface tension)   

double omega(int mode)
 {
    int i,end;
    double sumid, sumphi;

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
         phi[i] += mode*p;       
         if(i>=NiR) phiid[i] = T*rho[i]*(log(rho[i])-1.0) + rho[i]*(Vext[i]-mu);  //Due to convolution phi has contributions over a arger range than the ideal part
       }

    sumid=0.0;sumphi=0.0;
    for (i=0;i<iend-1;i++) sumphi+=phi[i]+phi[i+1];        
    for (i=0;i<iend-1;i++) sumid+=phiid[i]+phiid[i+1];
    
  return dz*(sumid+sumphi)/2.0;

 }
//}}}

//{{{ adsorption
double adsorption()
// Note this definition assumes that the zero of z starts at -R. This can lead to a negative adsorption.
{
int i,end;
double a=0.0;

for (i=0;i<iend;i++) a+=rho[i]+rho[i+1]-2*rhob;

return (dz*a/2.);
} //}}}

//{{{ update
void update() 
/* Here we update the density distribution. 
Every NGFREQ iteration we substitute a Picard step with a step of the Ng algorithm, see J. Chem. Phys. 61, 2680 (1974) */
{	
double ip0101, ip0202, ip0102, ipn01, ipn02;	
double a1,a2,norm;
int i;

dev=0.0;
ip0101=0.0; ip0202=0.0; ip0102=0.0; ipn01=0.0; ipn02=0.0;

for(i=0;i<iend;i++) 
  {
#ifdef LJ
   g[i] = exp(invT*(mu-Vext[i]-cphiatt[i])+dcf[i]);
#else
   g[i] = exp(invT*(mu-Vext[i])+dcf[i]);
#endif
 
  d[i]   = g[i] - rho[i];
  dev += fabs(d[i]);
  if(iter==0) { g2[i]=g[i]; d2[i]=d[i];}
  if(iter==1) { g1[i]=g[i]; d1[i]=d[i];}
  if(iter>1)
  	  {
       d01[i] = d[i] - d1[i];
       d02[i] = d[i] - d2[i];
       ip0101 += d01[i] * d01[i]; //calculate the inner products
       ip0202 += d02[i] * d02[i];
       ip0102 += d01[i] * d02[i];
       ipn01  +=  d[i]  * d01[i];
       ipn02  +=  d[i]  * d02[i];
     }
  }
dev*=dz;

if(iter>3)
  {
   ip0101*=dz; ip0202*=dz; ip0102*=dz; ipn01*=dz; ipn02*=dz;
   norm = ip0102 * ip0102 - ip0101 * ip0202;
   a1 = (ip0102 * ipn02 - ipn01  * ip0202) / norm;
   a2 = (ip0102 * ipn01 - ip0101 * ipn02) / norm; 
  }

for(i=0;i<iend;i++) 
  {
    if(dev<0.01 && iter % NGFREQ==0) rhonew[i] = (1-a1-a2)*g[i] + a1*g1[i] + a2*g2[i]; //Only use Ng if first sufficiently converged by Picard
      else rhonew[i]=(1-alpha)*rho[i] + alpha*g[i];
 	  g2[i] = g1[i]; // Shuffle down
      g1[i] = g[i];
      d2[i] = d1[i];
      d1[i] = d[i];
   }
} //}}}
