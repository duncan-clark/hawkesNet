#include <Rcpp.h>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include "cran_release_validation.h"

using namespace Rcpp;

// Exact release-selection probabilities and centered feature derivatives depend
// on the twelve shared preferences, but not on either component's count rates.
// [[Rcpp::export]]
List cran_release_mixture_cache(List data, NumericVector parameters, int families, int upto=0) {
  cran_input::release_data(data,parameters,families);
  NumericVector times=data["time"],newborn=data["new_count"];
  NumericVector change_pointer=data["change_ptr"],selected_pointer=data["selected_ptr"];
  NumericMatrix initial=data["initial"],changes=data["changes"],selected=data["selected"];
  if(upto<0)stop("Invalid release cache prefix");
  int full=times.size(),n=upto>0?upto:full;
  if(n<0||n>full||families<1)stop("Invalid release cache dimensions");
  if(initial.ncol()!=6||changes.ncol()!=6||selected.ncol()!=3||
     change_pointer.size()!=full+1||selected_pointer.size()!=full+1||newborn.size()!=full)
    stop("Invalid release cache arrays");
  int count=selected_pointer[n];
  NumericVector logprob(count);NumericMatrix centered(count,6);
  IntegerVector channel(count),batch(count);
  std::vector<std::array<double,5>> state(families);
  std::vector<std::array<double,6>> features(families);
  std::vector<std::array<double,2>> predictors(families);
  std::vector<int> active;std::vector<bool> seen(families,false);
  auto update=[&](const NumericMatrix& rows,int row){
    int family=rows(row,0);
    if(family<0||family>=families)stop("Release family index outside cache");
    if(!seen[family]){seen[family]=true;active.push_back(family);}
    for(int j=0;j<5;++j)state[family][j]=rows(row,j+1);
  };
  for(int row=0;row<initial.nrow();++row)update(initial,row);
  for(int event=0;event<n;++event){
    double maximum[2]={-INFINITY,-INFINITY};
    for(int family:active){
      double age=times[event]-state[family][0];
      auto& x=features[family];x[0]=age;x[1]=age*age;
      for(int j=2;j<6;++j)x[j]=state[family][j-1];
      for(int c=0;c<2;++c){
        double value=0.;for(int j=0;j<6;++j)value+=parameters[3+6*c+j]*x[j];
        predictors[family][c]=value;maximum[c]=std::max(maximum[c],value);
      }
    }
    double total[2]={0.,0.},mean[2][6]={{0.},{0.}};
    for(int family:active)for(int c=0;c<2;++c){
      double weight=std::exp(predictors[family][c]-maximum[c]);
      total[c]+=weight;
      for(int j=0;j<6;++j)mean[c][j]+=weight*features[family][j];
    }
    for(int c=0;c<2;++c){
      if(!(total[c]>0)||!std::isfinite(total[c]))stop("Invalid release selection normalizer");
      for(int j=0;j<6;++j)mean[c][j]/=total[c];
    }
    for(int row=selected_pointer[event];row<selected_pointer[event+1];++row){
      int family=selected(row,0),c=selected(row,1);double lines=selected(row,2);
      if(family<0||family>=families||!seen[family]||c<0||c>1||!(lines>0))
        stop("Invalid observed release selection");
      logprob[row]=predictors[family][c]-maximum[c]-std::log(total[c])-std::log(lines);
      for(int j=0;j<6;++j)centered(row,j)=features[family][j]-mean[c][j];
      channel[row]=c;batch[row]=event;
    }
    for(int row=change_pointer[event];row<change_pointer[event+1];++row)update(changes,row);
    if(event%100==0)checkUserInterrupt();
  }
  NumericVector selection(12);std::copy(parameters.begin()+3,parameters.end(),selection.begin());
  return List::create(_["logprob"]=logprob,_["centered"]=centered,_["channel"]=channel,
    _["batch"]=batch,_["newborn"]=newborn,_["fulln"]=full,_["upto"]=n,_["selection"]=selection);
}

// [[Rcpp::export]]
List cran_release_mixture_cached(List cache, NumericVector first, NumericVector second) {
  cran_input::release_parameters(first);cran_input::release_parameters(second);
  NumericVector selection=cache["selection"],logprob=cache["logprob"],newborn=cache["newborn"];
  NumericMatrix centered=cache["centered"];
  NumericVector channel=as<NumericVector>(cache["channel"]),batch=as<NumericVector>(cache["batch"]);
  double full_value=as<double>(cache["fulln"]),n_value=as<double>(cache["upto"]);
  if(!cran_input::integer(full_value)||full_value<0||full_value>INT_MAX/15||
     !cran_input::integer(n_value)||n_value<0||n_value>full_value)
    stop("Invalid cached release dimensions");
  int full=static_cast<int>(full_value),n=static_cast<int>(n_value);
  if(selection.size()!=12||newborn.size()!=full||centered.ncol()!=6||
     centered.nrow()!=logprob.size()||channel.size()!=logprob.size()||batch.size()!=logprob.size())
    stop("Invalid cached release arrays");
  for(double value:selection)if(!R_finite(value))stop("Non-finite cached release selection");
  for(double value:newborn)if(!cran_input::integer(value)||value<0||value>=std::ldexp(1.0,63))
    stop("Invalid cached release new count");
  for(int row=0;row<logprob.size();++row){
    if(!R_finite(logprob[row])||logprob[row]>0)stop("Invalid cached release log probability");
    cran_input::index(batch[row],n,"cached release batch");
    cran_input::index(channel[row],2,"cached release channel");
    for(int j=0;j<6;++j)if(!R_finite(centered(row,j)))stop("Non-finite cached release feature");
  }
  for(int j=0;j<12;++j)
    if(first[j+3]!=selection[j]||second[j+3]!=selection[j])stop("Stale release selection cache");
  List result(2);
  for(int component=0;component<2;++component){
    NumericVector parameters=component==0?first:second;
    double rates[3]={std::exp(parameters[0]),std::exp(parameters[1]),std::exp(parameters[2])};
    double total=rates[0]+rates[1]+rates[2],logz=std::log(-std::expm1(-total));
    NumericVector loglik(full);NumericMatrix gradient(full,15);
    for(int event=0;event<n;++event){
      loglik[event]=newborn[event]*parameters[0]-std::lgamma(newborn[event]+1.)-total-logz;
      for(int c=0;c<3;++c)gradient(event,c)=-rates[c]-rates[c]/std::expm1(total);
      gradient(event,0)+=newborn[event];
    }
    for(int row=0;row<logprob.size();++row){
      int event=batch[row],c=channel[row];
      double lograte=parameters[c+1]+logprob[row],rate=std::exp(lograte);
      double score=rate<1e-8?1+rate/2:rate/(-std::expm1(-rate));
      loglik[event]+=lograte< -30?lograte:rate+std::log(-std::expm1(-rate));
      gradient(event,c+1)+=score;
      for(int j=0;j<6;++j)gradient(event,3+6*c+j)+=score*centered(row,j);
    }
    result[component]=List::create(_["loglik"]=loglik,_["gradient"]=gradient);
  }
  return result;
}
