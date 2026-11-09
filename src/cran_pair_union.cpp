#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using namespace Rcpp;

// Positive frontier dynamic programming for the union of independent Poisson
// pair proposals. Vertices here are observed edits, not network vertices.
// Each selected pair has weight exp(rho)-1; a covered edit contributes f once.
// All incompatible proposals are accounted for by exp(-total rate) in R.
namespace {
struct PairEdge { int a,b; std::vector<int> rows; };
struct PairStep {
 int edge,apos,bpos,width;
 std::vector<int> keep;
};
struct PairComponent {
 int group,vertices,frontier=0;
 double work=0;
 std::vector<PairStep> steps;
};
struct PairPlan {
 int groups,edits,rows;
 std::vector<PairEdge> edges;
 std::vector<std::vector<int>> grouped_edges;
 std::vector<PairComponent> components;
};

PairComponent order_component(const PairPlan& p,int group,
                             const std::vector<int>& ids,bool degree_first) {
 PairComponent c;c.group=group;
 std::map<int,std::set<int>> incident;
 for(int e:ids){incident[p.edges[e].a].insert(e);incident[p.edges[e].b].insert(e);}
 c.vertices=incident.size();std::vector<int> frontier;
 while(!incident.empty()) {
  int chosen=-1;std::pair<int,int> best{INT_MAX,INT_MAX};
  for(auto& entry:incident) {
   std::set<int> expanded(frontier.begin(),frontier.end());
   for(int e:entry.second){expanded.insert(p.edges[e].a);expanded.insert(p.edges[e].b);}
   int degree=entry.second.size(),size=expanded.size();
   std::pair<int,int> score=degree_first?std::make_pair(degree,size):std::make_pair(size,degree);
   if(score<best){chosen=entry.first;best=score;}
  }
  std::vector<int> selected(incident[chosen].begin(),incident[chosen].end());
  std::stable_sort(selected.begin(),selected.end(),[&](int x,int y){
   auto rank=[&](int e){int v=p.edges[e].a==chosen?p.edges[e].b:p.edges[e].a;
    return std::make_pair(std::find(frontier.begin(),frontier.end(),v)==frontier.end(),incident[v].size());};
   return rank(x)<rank(y);
  });
  for(int e:selected) {
   auto edge=p.edges[e];
   for(int v:{edge.a,edge.b})if(std::find(frontier.begin(),frontier.end(),v)==frontier.end())frontier.push_back(v);
   PairStep s;s.edge=e;s.width=frontier.size();
   s.apos=std::find(frontier.begin(),frontier.end(),edge.a)-frontier.begin();
   s.bpos=std::find(frontier.begin(),frontier.end(),edge.b)-frontier.begin();
   incident[edge.a].erase(e);incident[edge.b].erase(e);
   std::vector<int> next;
   for(int j=0;j<(int)frontier.size();++j)if(!incident[frontier[j]].empty()){s.keep.push_back(j);next.push_back(frontier[j]);}
   if(incident[edge.a].empty())incident.erase(edge.a);
   if(incident[edge.b].empty())incident.erase(edge.b);
   c.frontier=std::max(c.frontier,s.width);c.work+=std::pow(2.,s.width);
   c.steps.push_back(s);frontier.swap(next);
  }
 }
 return c;
}

void check_budget(const PairPlan& p,double max_states) {
 if(!R_finite(max_states)||max_states<4)stop("pair_max_states must be finite and at least four");
 const PairComponent* worst=nullptr;
 for(auto& c:p.components)if(!worst||c.frontier>worst->frontier)worst=&c;
 if(worst&&(worst->frontier>=8*(int)sizeof(size_t)-1||std::pow(2.,worst->frontier)>max_states))
  stop("Exact Poisson-pair likelihood exceeds pair_max_states: releasing-line group %d has %d connected observed additions and a frontier of %d bits (%.0f states; budget %.0f). No likelihood approximation or candidate truncation was used.",
       worst->group+1,worst->vertices,worst->frontier,std::pow(2.,worst->frontier),max_states);
}
size_t compress(size_t mask,const std::vector<int>& keep) {
 size_t result=0;for(size_t j=0;j<keep.size();++j)if(mask&(size_t(1)<<keep[j]))result|=size_t(1)<<j;return result;
}
}

// [[Rcpp::export]]
SEXP cran_pair_union_plan(IntegerVector group,IntegerVector a,IntegerVector b,int groups,int edits) {
 if(group.size()!=a.size()||a.size()!=b.size()||groups<0||edits<0)stop("Invalid pair union dimensions");
 XPtr<PairPlan> p(new PairPlan(),true);p->groups=groups;p->edits=edits;p->rows=a.size();
 std::map<std::tuple<int,int,int>,int> index;
 std::vector<std::vector<int>> grouped(groups);
 for(int i=0;i<a.size();++i) {
  if(group[i]<0||group[i]>=groups||a[i]<0||a[i]>=edits||b[i]<0||b[i]>=edits||a[i]==b[i])stop("Invalid pair union endpoint");
  int u=std::min(a[i],b[i]),v=std::max(a[i],b[i]);auto key=std::make_tuple(group[i],u,v);
  auto found=index.find(key);int e;
  if(found==index.end()){e=p->edges.size();index[key]=e;p->edges.push_back({u,v,{}});grouped[group[i]].push_back(e);}else e=found->second;
  p->edges[e].rows.push_back(i);
 }
 p->grouped_edges=grouped;
 for(int g=0;g<groups;++g) {
  std::map<int,int> parent;
  for(int e:grouped[g]){parent.emplace(p->edges[e].a,p->edges[e].a);parent.emplace(p->edges[e].b,p->edges[e].b);}
  auto root=[&](int v){while(parent[v]!=v)v=parent[v];return v;};
  for(int e:grouped[g])parent[root(p->edges[e].a)]=root(p->edges[e].b);
  std::map<int,std::vector<int>> pieces;
  for(int e:grouped[g])pieces[root(p->edges[e].a)].push_back(e);
  for(auto& piece:pieces){
   auto first=order_component(*p,g,piece.second,true),second=order_component(*p,g,piece.second,false);
   p->components.push_back(std::make_pair(second.frontier,second.work)<std::make_pair(first.frontier,first.work)?second:first);
  }
  if(g%100==0)checkUserInterrupt();
 }
 return p;
}

// [[Rcpp::export]]
DataFrame cran_pair_union_complexity(SEXP plan) {
 XPtr<PairPlan> p(plan);int n=p->components.size();IntegerVector group(n),vertices(n),edges(n),frontier(n);NumericVector states(n),work(n);
 for(int i=0;i<n;++i){auto& c=p->components[i];group[i]=c.group+1;vertices[i]=c.vertices;edges[i]=c.steps.size();frontier[i]=c.frontier;states[i]=std::pow(2.,c.frontier);work[i]=c.work;}
 return DataFrame::create(_["group"]=group,_["vertices"]=vertices,_["pairs"]=edges,_["frontier"]=frontier,_["states"]=states,_["state_visits"]=work);
}

// [[Rcpp::export]]
void cran_pair_union_check(SEXP plan,double max_states) { XPtr<PairPlan> p(plan);check_budget(*p,max_states); }

// [[Rcpp::export]]
List cran_pair_union_evaluate(SEXP plan,NumericVector log_f,NumericMatrix dlog_f,
 NumericVector log_rate,NumericMatrix dlog_rate,double max_states=262144) {
 XPtr<PairPlan> p(plan);check_budget(*p,max_states);int dim=dlog_f.ncol(),stride=dim+1;
 if(log_f.size()!=p->edits||dlog_f.nrow()!=p->edits||log_rate.size()!=p->rows||dlog_rate.nrow()!=p->rows||dlog_rate.ncol()!=dim)stop("Invalid pair union values");
 for(double x:log_f)if(!R_finite(x))stop("Nonfinite pair coverage weight");
 for(double x:log_rate)if(!R_finite(x)&&x!=R_NegInf)stop("Nonfinite pair rate");
 for(double x:dlog_f)if(!R_finite(x))stop("Nonfinite pair coverage derivative");
 for(double x:dlog_rate)if(!R_finite(x))stop("Nonfinite pair rate derivative");
 NumericVector value(p->groups);NumericMatrix gradient(p->groups,dim);
 for(auto& c:p->components) {
  std::vector<double> current(stride,0.);current[0]=1.;double logz=0.;
  for(auto& s:c.steps) {
   auto& edge=p->edges[s.edge];double rho=0.;std::vector<double> score(dim,0.);
   for(int i:edge.rows){double r=std::exp(log_rate[i]);rho+=r;for(int j=0;j<dim;++j)score[j]+=r*dlog_rate(i,j);}
   if(!R_finite(rho))stop("Nonfinite summed pair rate");
   // d log(expm1(rho)) = d rho / (1-exp(-rho)).
   double included=-std::expm1(-rho),logw=rho>0?rho+std::log(included):R_NegInf;
   if(rho>0)for(double& v:score)v/=included;
   double shift=std::max(0.,logw+std::max(0.,log_f[edge.a])+std::max(0.,log_f[edge.b]));
   double absent=std::exp(-shift);
   size_t abit=size_t(1)<<s.apos,bbit=size_t(1)<<s.bpos,states=current.size()/stride;
   std::vector<double> next((size_t(1)<<s.keep.size())*stride,0.);
   for(size_t mask=0;mask<states;++mask) {
    const double* old=&current[mask*stride];if(old[0]==0.)continue;
    size_t no=compress(mask,s.keep),yes=compress(mask|abit|bbit,s.keep);
    bool add_a=!(mask&abit),add_b=!(mask&bbit);
    double selected=std::exp(logw+(add_a?log_f[edge.a]:0.)+(add_b?log_f[edge.b]:0.)-shift);
    next[no*stride]+=old[0]*absent;next[yes*stride]+=old[0]*selected;
    for(int j=0;j<dim;++j){
     double delta=score[j]+(add_a?dlog_f(edge.a,j):0.)+(add_b?dlog_f(edge.b,j):0.);
     next[no*stride+j+1]+=old[j+1]*absent;
     next[yes*stride+j+1]+=(old[j+1]+old[0]*delta)*selected;
    }
   }
   double norm=0.;for(size_t i=0;i<next.size();i+=stride)norm+=next[i];
   if(!(norm>0.)||!R_finite(norm))stop("Pair union numerical overflow or underflow");
   for(double& v:next)v/=norm;
   logz+=shift+std::log(norm);current.swap(next);
  }
  if(current.size()!=(size_t)stride)stop("Pair union frontier did not close");
  value[c.group]+=logz;for(int j=0;j<dim;++j)gradient(c.group,j)+=current[j+1];
  checkUserInterrupt();
 }
 return List::create(_["log_partition"]=value,_["gradient"]=gradient);
}

// Exact K=min(Poisson(M),2) mixture. Two pair draws are classified as equal,
// sharing one endpoint, or disjoint. Vertex prefix sums avoid enumerating all
// pairs of pairs; cancellation-prone disjoint sums use a positive direct sum.
// [[Rcpp::export]]
List cran_pair_capped2_evaluate(SEXP plan,NumericVector log_f,NumericMatrix dlog_f,
 NumericVector log_rate,NumericMatrix dlog_rate,NumericVector mean,NumericMatrix dmean) {
 XPtr<PairPlan> p(plan);int dim=dlog_f.ncol(),stride=dim+1;
 if(log_f.size()!=p->edits||dlog_f.nrow()!=p->edits||log_rate.size()!=p->rows||dlog_rate.nrow()!=p->rows||dlog_rate.ncol()!=dim||mean.size()!=p->groups||dmean.nrow()!=p->groups||dmean.ncol()!=dim)stop("Invalid capped-pair dimensions");
 NumericVector value(p->groups);NumericMatrix gradient(p->groups,dim);
 for(int g=0;g<p->groups;++g){
  double M=mean[g];if(!R_finite(M)||M<0)stop("Invalid addition-pair mean");
  value[g]=-M;for(int j=0;j<dim;++j)gradient(g,j)=-dmean(g,j);
  auto& ids=p->grouped_edges[g];if(ids.empty()||M==0)continue;
  int ne=ids.size();std::map<int,int> vertices;
  for(int e:ids){vertices.emplace(p->edges[e].a,vertices.size());vertices.emplace(p->edges[e].b,vertices.size());}
  // map insertion may reuse a tentative index; assign consecutive indices.
  int nv=0;for(auto& v:vertices)v.second=nv++;
  std::vector<int> a(ne),b(ne),count(nv,0);std::vector<long double> t(ne),score(ne*dim),lf(nv),df(nv*dim),s(nv*stride,0.),total(stride,0.),B(stride,0.);
  long double scale=-INFINITY;
  for(auto v:vertices){lf[v.second]=log_f[v.first];for(int j=0;j<dim;++j)df[v.second*dim+j]=dlog_f(v.first,j);}
  for(int k=0;k<ne;++k){
   auto& e=p->edges[ids[k]];a[k]=vertices[e.a];b[k]=vertices[e.b];long double rho=0;
   for(int i:e.rows){long double r=std::exp((long double)log_rate[i]);rho+=r;for(int j=0;j<dim;++j)score[k*dim+j]+=r*dlog_rate(i,j);}
   if(!(rho>0)) {t[k]=-INFINITY;continue;}
   t[k]=std::log(rho)+lf[a[k]]+lf[b[k]];scale=std::max(scale,t[k]);
   for(int j=0;j<dim;++j)score[k*dim+j]=score[k*dim+j]/rho+df[a[k]*dim+j]+df[b[k]*dim+j];
  }
  if(!std::isfinite(scale))continue;
  for(auto& x:t)x=std::exp(x-scale);
  for(int k=0;k<ne;++k){
   int u=a[k],v=b[k];long double x=t[k];if(x==0)continue;
   long double *su=&s[u*stride],*sv=&s[v*stride],*dx=&score[k*dim];
   // Same candidate drawn twice: the two edits are still retained just once.
   long double same=x*x*std::exp(-lf[u]-lf[v]);B[0]+=same;
   for(int j=0;j<dim;++j)B[j+1]+=same*(2*dx[j]-df[u*dim+j]-df[v*dim+j]);
   // Distinct candidates sharing one endpoint, with both draw orders.
   for(int z:{u,v}){
    long double factor=2*x*std::exp(-lf[z]);long double* sz=&s[z*stride];B[0]+=factor*sz[0];
    for(int j=0;j<dim;++j)B[j+1]+=factor*(sz[j+1]+sz[0]*(dx[j]-df[z*dim+j]));
   }
   int nd=k-count[u]-count[v];
   if(nd>0){
    std::vector<long double> disjoint(stride,0.);
    long double diff=total[0]-su[0]-sv[0];
    if(diff>1e-8L*total[0]){
     disjoint[0]=diff;for(int j=0;j<dim;++j)disjoint[j+1]=total[j+1]-su[j+1]-sv[j+1];
    }else{
     for(int h=0;h<k;++h)if(a[h]!=u&&a[h]!=v&&b[h]!=u&&b[h]!=v){disjoint[0]+=t[h];for(int j=0;j<dim;++j)disjoint[j+1]+=t[h]*score[h*dim+j];}
    }
    B[0]+=2*x*disjoint[0];for(int j=0;j<dim;++j)B[j+1]+=2*x*(disjoint[j+1]+disjoint[0]*dx[j]);
   }
   total[0]+=x;for(int j=0;j<dim;++j)total[j+1]+=x*dx[j];
   for(int z:{u,v}){s[z*stride]+=x;for(int j=0;j<dim;++j)s[z*stride+j+1]+=x*dx[j];count[z]++;}
  }
  long double c2,dc2;
  if(M<.001){
   c2=0;dc2=0;long double term=1;
   for(int k=0;k<12;++k){c2+=term/(k+2);if(k>0)dc2+=k*term/(M*(k+2));term*=-M/(k+1);}
  }else{long double mass=-std::expm1(-(long double)M)-M*std::exp(-(long double)M);c2=mass/(M*M);dc2=std::exp(-(long double)M)/M-2*c2/M;}
  long double l0=-M,l1=-M+scale+std::log(total[0]),l2=std::log(c2)+2*scale+std::log(B[0]);
  long double center=std::max(l0,std::max(l1,l2)),v0=std::exp(l0-center),v1=std::exp(l1-center),v2=std::exp(l2-center),sum=v0+v1+v2;
  value[g]=center+std::log(sum);
  for(int j=0;j<dim;++j)gradient(g,j)=(-(v0+v1)*dmean(g,j)+v1*total[j+1]/total[0]+v2*(B[j+1]/B[0]+dc2/c2*dmean(g,j)))/sum;
  if(!R_finite(value[g]))stop("Nonfinite capped-pair likelihood");
  if(g%100==0)checkUserInterrupt();
 }
 return List::create(_["loglik"]=value,_["gradient"]=gradient);
}
