/*
Compile as
$gcc -Wall -I/home/mario/gsl/include -c MXGamma.c
$gcc -L/home/mario/gsl/lib MXGamma.o -lgsl -lgslcblas -lm
$LD_LIBRARY_PATH=/home/mario/gsl/lib
$export LD_LIBRARY_PATH
Run
$./a.out
*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include<stdlib.h> 
#include<math.h> 
#include <time.h>
#define GNUPLOT "gnuplot -persist"
//2 state diffusivity model
double randn(double,double);
float ran2(long *idum);
int ProbSampleReplace(int n, double *p);
int main()
{
  
  double h,t,x,T,tauD,tau0,D,Dm,OT,c,x0,ta,delt;
  double mu=0,sigma=1;
  double meanp, meanm;
  int i,N,id,k,s;
  double ap,am,tp,tm;
  //long iseed;
  const gsl_rng_type * TT;
  gsl_rng * r;
  FILE *data;
  FILE *dat;
  char nom[50];
  int ini=1;
  int fin=300000;
  double pib[2]={0,0};//prob of plus=pib[0]; prob of minus=pib[1]
  //srand(time(NULL));
  gsl_rng_env_setup();

  TT = gsl_rng_default;
  r = gsl_rng_alloc (TT);
  x=0;
  x0=0;
  D=10;//D_{+}
  Dm=0;//D_{-}
  h=0.001;
  t=0;// t is t_N is Godreche as t_N=tau_1+tau_2+... +tau_N 
  OT=0;
  tauD=1;
  tau0=1;
  T=30;//final time, t in Godreche
  N=100;

  ap=0.5;//parameter a rv Gamma(ap,1) for tau_{+}
  am=0.5;//parameter a rv Gamma(am,1) for tau_{-}
  tp=8.0;//parameter theta_{+} rv Gamma(ap,thp) for tau_{+}
  tm=12.0;//parameter theta_{-} rv Gamma(am,thm) for tau_{-}
  meanp=ap*tp;//mean tau_{+}
  meanm=am*tp;//mean tau_{-}
  pib[0]=meanp/(meanp+meanm);//probability of IC plus
  pib[1]=meanm/(meanp+meanm);//probability of IC minus
  sprintf(nom,"Ppap4am6T30.dat");
  data=fopen(nom,"w");
  sprintf(nom,"Xpap4am6T30.dat");
  dat=fopen(nom,"w");
  for(id=ini;id<=fin;id++) {//loop over different systems
    //select randomly the initia; condition
    s=ProbSampleReplace(2,pib);//select initial condition
    if(s == 1){//initial condition plus state
      //printf("IC +\n");
      t=0;
      x=sqrt(2*D)*sqrt(1.0/N)*randn(mu,sigma);//starts walking;
      x0=x;
      OT=0;
      k=0;
      while(t<=T){
	//enter in plus state
	k++;//create one transition between states
	//iseed=-(long) rand();
	//tauD=gamdev(ap,&iseed);//tau_{+}
	tauD=gsl_ran_gamma(r,ap,tp);
	//sample a time shorter than T
	//while(tauD>T) tauD= ad + (bd-ad)*ran2(&iseed);
	//printf("step: k=%d; %lf\n",k,tauD);
	//start diffusion process at D_{+}
	h=tauD/N;//time step used in diffusion process +
	OT+=tauD;
	t+=tauD;
	ta=t-tauD;
	for(i=1;i<N;i++){//walk within tau_{+}
	  delt= ta + i*h;
	  if(delt  <= T)//walk until T
	    x+=sqrt(2*D)*sqrt(h)*randn(mu,sigma);//Lagevin eq.
	  //fprintf(data,"%f\t%f\n",t+i*h,x);
	}//end for
	//end diffusion process at D_{+}
	
		
	if(t>T) break;
	
	//enter in minus state
	k++;//create one transition between states
	//iseed=-(long) rand();
	//tau0=gamdev(am,&iseed);//tau_{-}
	tau0=gsl_ran_gamma(r,am,tm);
	//while(tau0>T) tau0=am + (bm-am)*ran2(&iseed);//tau_{-} less tah T
	//printf("step: k=%d; %lf\n",k,tau0);
	//start diffusion process at D_{-}
	h=tau0/N;//time step used in  process -
	t+=tau0;
	ta=t-tau0;
	for(i=1;i<N;i++) {//not move in - state
	  delt= ta + i*h;
	  if(delt  <= T)//walk until T
	    x+=sqrt(2*Dm)*sqrt(h)*randn(mu,sigma);//Lagevin eq.
	  //fprintf(data,"%f\t%f\n",t+i*h,x);
	}//end for
	//end diffusion process at D_{-}
	
	
      }//end while(t<T)

      if((k % 2) ==0) {//even number of jumps at T
	//OT for starting from + and \tilde{k} even
	//printf("tp=%g;  OT=%g; k=%d \n",t,OT,k);
	//printf("tp=%g;  p_{+}=%g; k=%d \n",t,OT/T,k);
	fprintf(data,"%.100g\t%d\n",OT/T,k);
      }
      else {//odd number of jumps at T
	//OT for starting from + and \tilde{k} odd
	//OT is the sum of tau_{+} until k minus the forward time (t-T)
	//printf("tp=%g;  OT=%g; k=%d \n",t,OT-(t-T),k);
	//printf("tp=%g;  p_{+}=%g; k=%d \n",t,(OT-(t-T))/T,k);
	fprintf(data,"%.100g\t%d\n",(OT-(t-T))/T,k);
      }
      fprintf(dat,"%.100g\n",x-x0);
    }//end if(s == 1)
    
    else if(s == 2){//initial condition minus state
      //printf("IC -\n");
      t=0;
      x=sqrt(2*Dm)*sqrt(1.0/N)*randn(mu,sigma);//starts walking at minus;
      x0=x;
      OT=0;
      k=0;
      while(t<=T){
	k++;
	//iseed=-(long) rand();
	//tau0=gamdev(am,&iseed);//tau_{-}
	tau0=gsl_ran_gamma(r,am,tm);
	//sample a time shorter than T
	//while(tau0>T) tau0=am + (bm-am)*ran2(&iseed);//tau_{-}
	//printf("step: k=%d; %lf\n",k,tau0);
	//start diffusion at D_{-}
	h=tau0/N;
	t+=tau0;
	ta=t-tau0;
	for(i=1;i<N;i++){//not move in - state
	  delt= ta + i*h;
	  if(delt  <= T)//walk until T
	    x+=sqrt(2*Dm)*sqrt(h)*randn(mu,sigma);//Lagevin eq.
	  //fprintf(data,"%f\t%f\n",t+i*h,x);
	}//end for
	//end diffusion at D_{-}
	
	if(t>T) break;
	k++;
	//iseed=-(long) rand();
	//tauD=gamdev(ap,&iseed);//tau_{+}
	tauD=gsl_ran_gamma(r,ap,tp);
	//while(tauD>T) tauD=gamdev(ap,&iseed);//tau_{+}
	//printf("step: k=%d; %lf\n",k,tauD);
	//start diffusion at D_{+}
	h=tauD/N;
	OT+=tauD;
	t+=tauD;
	ta=t-tauD;
	for(i=1;i<N;i++) {//walk within tau_{+}
	  delt= ta + i*h;
	  if(delt  <= T)//walk until T
	    x+=sqrt(2*D)*sqrt(h)*randn(mu,sigma);//Lagevin eq.
	  //fprintf(data,"%f\t%f\n",t+i*h,x);
	}//end for
	//end diffusion at D_{+}
	
      }//end while(t<T)

      /*saving data*/
      if((k % 2) ==0) {//even number of jumps at T
	//OT for starting from - and \tilde{k} even
	//OT is the sum of tau_{+} until k minus the forward time (t-T)
	//printf("tp=%g;  OT=%g; k=%d \n",t,OT-(t-T),k);
	//printf("tp=%g;  p_{+}=%g; k=%d \n",t,(OT-(t-T))/T,k);
	fprintf(data,"%.100g\t%d\n",(OT-(t-T))/T,k);
      }
      else {//odd number of jumps at T
	//OT for starting from - and \tilde{k} odd
	//printf("tp=%g;  OT=%g; k=%d \n",t,OT,k);
	//printf("tp=%g;  p_{+}=%g; k=%d \n",t,OT/T,k);
	fprintf(data,"%.100g\t%d\n",OT/T,k);
      }
      fprintf(dat,"%.100g\n",x-x0);
      
    }//end else IC minus
    
  }//end for id
  gsl_rng_free (r);
  fclose(data);
  fclose(dat);
  return 0;

}



#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX




int ProbSampleReplace(int n, double *p)
{
  
  float rU;
  int i, j,ans;
  int nm1 = n - 1;
  int perm[2]={0};//the number of elements in these 
  double q[2]={0};//two arrays must coincede with p[#N]
  long iseed;
  /* record element identities */
  for (i = 0; i < n; i++) {
    perm[i] = i + 1;
    q[i]=p[i];
  }
  /* compute cumulative probabilities */
  for (i = 1 ; i < n; i++)
    q[i] += q[i - 1];
  /* compute the sample */
  iseed=-(long) rand();
  rU = ran2(&iseed);
  for (j = 0; j < nm1; j++) {
    if (rU <= q[j])
      break;
  }
  //en C despues del for j=nm1
  ans = perm[j];
  return ans;
}


double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}
