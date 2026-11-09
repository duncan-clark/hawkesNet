#include <vector>
#include <cmath>
#include <algorithm>
// Positive batches: a Poisson new-family count and independent collapsed
// Poisson selections for new majors and existing major-line updates.
extern "C" void evaluate(int n,int nf,int ni,const double* initial,
 const double* changes,const long long* cp,const long long* selected,
 const long long* sp,const double* times,const long long* newborn,
 const double* par,double* ll,double* grad) {
 std::vector<double> state(nf*5,0),x(nf*6),weight(nf*2),eta(nf*2);
 std::vector<int> active;std::vector<bool> seen(nf,false);
 auto update=[&](const double* row){int u=(int)row[0];
  if(!seen[u]){seen[u]=true;active.push_back(u);}
  for(int k=0;k<5;k++)state[5*u+k]=row[k+1];};
 for(int j=0;j<ni;j++)update(initial+6*j);
 double m[3]={exp(par[0]),exp(par[1]),exp(par[2])};
 double total=m[0]+m[1]+m[2], logz=log(-expm1(-total));
 for(int b=0;b<n;b++) {
  double* g=grad+15*b;std::fill(g,g+15,0);
  ll[b]=newborn[b]*par[0]-lgamma(newborn[b]+1.)-total-logz;
  for(int c=0;c<3;c++)g[c]=-m[c]-m[c]/expm1(total);
  g[0]+=newborn[b];double maxeta[2]={-INFINITY,-INFINITY};
  for(int u:active){double age=times[b]-state[5*u];
   x[6*u]=age;x[6*u+1]=age*age;
   for(int k=2;k<6;k++)x[6*u+k]=state[5*u+k-1];
   for(int c=0;c<2;c++){double e=0;
    for(int k=0;k<6;k++)e+=par[3+6*c+k]*x[6*u+k];
    eta[2*u+c]=e;maxeta[c]=std::max(maxeta[c],e);}}
  double sum[2]={0,0},mean[12]={0};
  for(int u:active)for(int c=0;c<2;c++){
   double w=exp(eta[2*u+c]-maxeta[c]);weight[2*u+c]=w;sum[c]+=w;
   for(int k=0;k<6;k++)mean[6*c+k]+=w*x[6*u+k];}
  for(int c=0;c<2;c++)for(int k=0;k<6;k++)mean[6*c+k]/=sum[c];
  for(long long j=sp[b];j<sp[b+1];j++){
   int u=(int)selected[3*j],c=(int)selected[3*j+1];
   double lograte=par[1+c]+eta[2*u+c]-maxeta[c]-log(sum[c])-log((double)selected[3*j+2]);
   double r=exp(lograte),score=r<1e-8?1+r/2:r/(-expm1(-r));
   ll[b]+=lograte < -30?lograte:r+log(-expm1(-r));g[c+1]+=score;
   for(int k=0;k<6;k++)g[3+6*c+k]+=score*(x[6*u+k]-mean[6*c+k]);
  }
  for(long long j=cp[b];j<cp[b+1];j++)update(changes+6*j);
 }
}
