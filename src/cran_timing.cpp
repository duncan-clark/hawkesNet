// CRAN timing engine: exact calendar exposure and exponential recursions.
// The R interface restricts this engine to the two-exponential application model.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
List cran_timing_native(List data, NumericVector par, NumericVector beta,
                       IntegerVector shape, int model, bool separate,
                       bool last_only, double clip, double start, double end,
                       NumericVector baseline, NumericVector exposure,
                       bool details=false, double end_baseline=1.) {
  NumericVector times=data["times"], initial=data["initial"], birth=data["tz"],
    weight=data["tw"], logn=data["tx"], changed=data["cz"];
  NumericVector frozen=data.containsElementNamed("frozen")?as<NumericVector>(data["frozen"]):birth;
  if(frozen.size()!=birth.size()) stop("Invalid frozen topology dimension");
  IntegerVector tp=data["tp"], cp=data["cp"], tu=data["tu"], cu=data["cu"];
  int S=beta.size(), N=times.size(), P=initial.size(), G=separate?2:1;
  if(!std::isfinite(start)||!std::isfinite(end)||start<0||end<=start||
     !std::isfinite(clip)||clip<0)stop("Invalid timing window or clipping bound");
  auto finite=[](const NumericVector& x){for(double value:x)if(!std::isfinite(value))return false;return true;};
  if(!finite(times)||!finite(initial)||!finite(birth)||!finite(weight)||!finite(logn)||
     !finite(changed)||!finite(frozen)||!finite(par)||!finite(beta)||!finite(baseline)||!finite(exposure))
    stop("Timing inputs must be finite");
  for(int i=0;i<N;++i)if(times[i]<0||(i>0&&times[i]<=times[i-1]))stop("Timing times must be strictly increasing and nonnegative");
  if(tu.size()!=birth.size()||weight.size()!=birth.size()||logn.size()!=birth.size()||cu.size()!=changed.size())
    stop("Invalid timing trigger or change dimensions");
  auto pointers=[&](const IntegerVector& pointer,int size){
    if(pointer.size()!=N+1||pointer[0]!=0)stop("Invalid timing pointers");
    for(int i=0;i<pointer.size();++i)
      if(pointer[i]<0||pointer[i]>size||(i>0&&pointer[i]<pointer[i-1]))stop("Invalid timing pointers");
  };
  pointers(tp,birth.size());pointers(cp,changed.size());
  for(int id:tu)if(id<0||id>=P)stop("Timing trigger references an absent family");
  for(int id:cu)if(id<0||id>=P)stop("Timing change references an absent family");
  for(double value:weight)if(value<0)stop("Timing weights must be nonnegative");
  for(double value:baseline)if(value<0)stop("Timing baseline must be nonnegative");
  for(double value:exposure)if(value<0)stop("Timing exposure must be nonnegative");
  int base_H=G+(model==2)+1, ib=1+S, ic=G, id=base_H-1;
  int type_H=par.size()-ib-base_H;
  bool typed=type_H>0, type_by_kernel=type_H>2;
  if(type_H!=0 && type_H!=2 && type_H!=2*S) stop("Invalid release amplitude dimension");
  int H=base_H+type_H, D=ib+H;
  IntegerVector release_class;
  if(typed){
    if(!data.containsElementNamed("release_class")) stop("Missing timing release classes");
    release_class=as<IntegerVector>(data["release_class"]);
    if(release_class.size()!=birth.size()) stop("Invalid timing release class dimension");
    for(int c:release_class) if(c<0 || c>2) stop("Invalid timing release class");
  }
  if(par.size()!=D || shape.size()!=S || (model!=1 && model!=2) ||
     (last_only && model==2)) stop("Invalid model dimensions");
  if(baseline.size()!=N || exposure.size()!=N+1) stop("Invalid calendar dimensions");
  double mu=std::exp(par[0]);
  std::vector<double> alpha(S), Q(S,0.), R(S,0.), dQ(S*H,0.), dR(S*H,0.),
    grad(D,0.), z(initial.begin(),initial.end()), trigger_time(tu.size());
  std::vector<int> group(S,0);
  std::vector<std::vector<int>> history(P);
  for(int s=0;s<S;s++) {
    alpha[s]=std::exp(par[1+s]);
    group[s]=separate && shape[s]/beta[s]>=1.-1e-10;
    if(beta[s]<=0 || (shape[s]!=1 && shape[s]!=2)) stop("Invalid response");
  }
  auto amplitude=[&](int j,double current,int s,std::vector<double>& derivative) {
    int g=group[s];
    double change=model==2?std::max(-clip,std::min(clip,current-birth[j])):0.;
    double value=par[ib+g]*frozen[j]+par[ib+id]*logn[j];
    if(model==2) value+=par[ib+ic]*change;
    int c=typed?release_class[j]:0;
    int type_index=base_H+2*(type_by_kernel?s:0)+c-1;
    if(c>0) value+=par[ib+type_index];
    double a=weight[j]*std::exp(value);
    std::fill(derivative.begin(),derivative.end(),0.);
    derivative[g]=a*frozen[j]; derivative[id]=a*logn[j];
    if(model==2) derivative[ic]=a*change;
    if(c>0) derivative[type_index]=a;
    return a;
  };
  NumericVector contributions(details?N+1:0), comp(details?N+1:0), rates(details?N:0);
  double previous=0., ll=0., total_comp=0.;
  std::vector<double> da(H), db(H);
  for(int i=0;i<=N;i++) {
    double right=std::min(i<N?times[i]:end,end), left=std::max(previous,start);
    double mass=mu*exposure[i]; grad[0]-=mass;
    if(right>left) for(int s=0;s<S;s++) {
      double gap=left-previous, width=right-left, lag=beta[s]*gap;
      double decay=std::exp(-lag), q=Q[s]*decay,
        r=(R[s]+lag*Q[s])*decay, x=beta[s]*width;
      double J0=-std::expm1(-x), J1=J0-x*std::exp(-x);
      if(x<1e-4) J1=x*x/2-x*x*x/3+x*x*x*x/8-x*x*x*x*x/30;
      double a=alpha[s]*(shape[s]==1?q*J0:r*J0+q*J1);
      mass+=a; grad[1+s]-=a;
      for(int h=0;h<H;h++) {
        double dq=dQ[s*H+h]*decay,
          dr=(dR[s*H+h]+lag*dQ[s*H+h])*decay;
        grad[ib+h]-=alpha[s]*(shape[s]==1?dq*J0:dr*J0+dq*J1);
      }
    }
    total_comp+=mass; ll-=mass;
    if(details) { comp[i]=mass; contributions[i]=-mass; }
    if(i==N || times[i]>=end) break;
    double dt=times[i]-previous;
    for(int s=0;s<S;s++) {
      double x=beta[s]*dt, e=std::exp(-x);
      R[s]=(R[s]+x*Q[s])*e; Q[s]*=e;
      for(int h=0;h<H;h++) {
        dR[s*H+h]=(dR[s*H+h]+x*dQ[s*H+h])*e; dQ[s*H+h]*=e;
      }
    }
    if(times[i]>=start) {
      double rate=mu*baseline[i];
      for(int s=0;s<S;s++) rate+=alpha[s]*beta[s]*(shape[s]==1?Q[s]:R[s]);
      if(!(rate>0) || !R_finite(rate)) stop("Invalid intensity");
      ll+=std::log(rate); grad[0]+=mu*baseline[i]/rate;
      for(int s=0;s<S;s++) {
        double a=alpha[s]*beta[s];
        grad[1+s]+=a*(shape[s]==1?Q[s]:R[s])/rate;
        for(int h=0;h<H;h++) grad[ib+h]+=a*(shape[s]==1?dQ[s*H+h]:dR[s*H+h])/rate;
      }
      if(details) { rates[i]=rate; contributions[i]+=std::log(rate); }
    }
    if(last_only) {
      std::fill(Q.begin(),Q.end(),0.); std::fill(R.begin(),R.end(),0.);
      std::fill(dQ.begin(),dQ.end(),0.); std::fill(dR.begin(),dR.end(),0.);
    }
    if(model==2) for(int k=cp[i];k<cp[i+1];k++) {
      int u=cu[k];
      for(int j:history[u]) for(int s=0;s<S;s++) {
        double old=amplitude(j,z[u],s,da), now=amplitude(j,changed[k],s,db);
        double lag=beta[s]*(times[i]-trigger_time[j]), e=std::exp(-lag);
        Q[s]+=(now-old)*e; R[s]+=(now-old)*e*lag;
        for(int h=0;h<H;h++) {
          dQ[s*H+h]+=(db[h]-da[h])*e;
          dR[s*H+h]+=(db[h]-da[h])*e*lag;
        }
      }
      z[u]=changed[k];
    }
    for(int j=tp[i];j<tp[i+1];j++) {
      int u=tu[j]; trigger_time[j]=times[i];
      if(model==2 && std::abs(z[u]-birth[j])>1e-8) stop("Birth summary disagrees with graph state");
      for(int s=0;s<S;s++) {
        double a=amplitude(j,birth[j],s,da); Q[s]+=a;
        for(int h=0;h<H;h++) dQ[s*H+h]+=da[h];
      }
      if(model==2) history[u].push_back(j);
    }
    previous=times[i];
  }
  double end_rate=mu*end_baseline;for(int s=0;s<S;++s){double lag=beta[s]*(end-previous),decay=std::exp(-lag);end_rate+=alpha[s]*beta[s]*decay*(shape[s]==1?Q[s]:R[s]+lag*Q[s]);}
  List out=List::create(_["loglik"]=ll,_["gradient"]=wrap(grad),_["compensator"]=total_comp,_["end_intensity"]=end_rate);
  if(details) {out["contributions"]=contributions;out["interval_compensator"]=comp;out["intensity"]=rates;}
  return out;
}
