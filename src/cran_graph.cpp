#include <Rcpp.h>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>
#include <random>
#include <queue>
#include <numeric>
#include <fstream>
#include <cstring>
using namespace Rcpp;
using Edges=std::set<std::pair<int,int>>;
static double choose2(double d){return d*(d-1)/2;}
static double choose3(double d){return d*(d-1)*(d-2)/6;}
static double sp(double z){return std::max(z,0.)+std::log1p(std::exp(-std::abs(z)));}
static double logistic(double z){return z>=0?1/(1+std::exp(-z)):std::exp(z)/(1+std::exp(z));}
static int intersect_size(const std::set<int>& a,const std::set<int>& b){const auto& small=a.size()<b.size()?a:b;const auto& large=a.size()<b.size()?b:a;int n=0;for(int v:small)n+=large.count(v);return n;}
struct Update {Edges add,del;std::vector<int> sources,born,release_class;int pairs=0;int component=1;};
struct CranGraph {
 std::vector<std::string> names,fnames;std::map<std::string,int> lookup,flookup;
 std::vector<int> family,major;std::vector<double> born,last;
 std::vector<std::set<int>> out,in,und,fo,fi;std::vector<std::vector<int>> fn;
 std::vector<std::vector<int>> personnel;
 std::unordered_map<int,std::vector<int>> personnel_families;
 bool has_personnel=false;
 std::map<std::pair<int,int>,int> arcs;std::vector<std::array<double,4>> reach;int synthetic=0;
 bool diamond_ready=false;
 std::vector<double> diamond_sum,diamond_sq,diamond_linked;
 double diamond_wedges=0.,diamond_closed=0.;
 bool consumer_ready=false;
 std::vector<double> consumer_sum,consumer_sq;
 double consumer_wedges=0.;
 int node(std::string name,double birth,double clock){auto it=lookup.find(name);if(it!=lookup.end())return it->second;size_t sep=name.rfind('@');if(sep==std::string::npos)stop("Major-line name requires @major");std::string base=name.substr(0,sep),version=name.substr(sep+1);
  if(base.empty()||version.empty()||version.find_first_not_of("0123456789")!=std::string::npos)
   stop("Major-line names must have the form package@nonnegative-integer");
  int maj;
  try{maj=std::stoi(version);}catch(const std::exception&){stop("Major version is outside the integer range");}
  int f;
  auto ff=flookup.find(base);if(ff==flookup.end()){f=fnames.size();flookup[base]=f;fnames.push_back(base);fn.emplace_back();fo.emplace_back();fi.emplace_back();last.push_back(-INFINITY);reach.push_back({0,0,0,0});if(has_personnel)personnel.emplace_back();}else f=ff->second;
  int u=names.size();lookup[name]=u;names.push_back(name);family.push_back(f);major.push_back(maj);born.push_back(birth);out.emplace_back();in.emplace_back();und.emplace_back();fn[f].push_back(u);last[f]=std::max(last[f],clock);
  if(diamond_ready){diamond_sum.push_back(0.);diamond_sq.push_back(0.);diamond_linked.push_back(0.);}
  if(consumer_ready){consumer_sum.push_back(0.);consumer_sq.push_back(0.);}
  return u;
 }
 double target_weight(int v) const {return 1.+in[v].size();}
 static double wedge_mass(double s,double sq){return (s*s-sq)/2.;}
 double linked_wedge_weight(int v,int w) const {
  return target_weight(v)*target_weight(w)*intersect_size(out[v],out[w]);
 }
 void linked_delta(const std::set<int>& affected,double sign){
  for(int v:affected)for(int w:und[v])if(v<w||!affected.count(w)){
   double d=sign*linked_wedge_weight(v,w);diamond_closed+=d;diamond_linked[v]+=d;diamond_linked[w]+=d;
  }
 }
 void enable_diamonds(){
  if(diamond_ready)return;
  diamond_ready=true;int n=names.size();diamond_sum.assign(n,0.);diamond_sq.assign(n,0.);diamond_linked.assign(n,0.);
  for(int v=0;v<n;++v){double a=target_weight(v);for(int z:out[v]){diamond_sum[z]+=a;diamond_sq[z]+=a*a;}}
  for(int z=0;z<n;++z)diamond_wedges+=wedge_mass(diamond_sum[z],diamond_sq[z]);
  for(int v=0;v<n;++v)for(int w:und[v])if(v<w){double x=linked_wedge_weight(v,w);diamond_closed+=x;diamond_linked[v]+=x;diamond_linked[w]+=x;}
 }
 std::set<int> diamond_before(const Edges& add,const Edges& del){
  std::set<int> affected;std::map<int,int> degree;
  for(auto e:add){affected.insert(e.first);affected.insert(e.second);++degree[e.second];}
  for(auto e:del){affected.insert(e.first);affected.insert(e.second);--degree[e.second];}
  linked_delta(affected,-1.);
  std::unordered_map<int,std::array<double,2>> delta;
  for(auto x:degree)if(x.second){double a=target_weight(x.first),b=a+x.second;
   for(int z:out[x.first]){auto& d=delta[z];d[0]+=b-a;d[1]+=b*b-a*a;}
  }
  auto edit=[&](std::pair<int,int> e,double sign){double b=target_weight(e.first)+degree[e.first];auto& d=delta[e.second];d[0]+=sign*b;d[1]+=sign*b*b;};
  for(auto e:del)edit(e,-1.);for(auto e:add)edit(e,1.);
  for(auto x:delta){int z=x.first;diamond_wedges-=wedge_mass(diamond_sum[z],diamond_sq[z]);
   diamond_sum[z]+=x.second[0];diamond_sq[z]+=x.second[1];diamond_wedges+=wedge_mass(diamond_sum[z],diamond_sq[z]);}
  return affected;
 }
 int open_diamond_partners(int u,int v,int w) const {
  if(und[v].count(w))return 0;
  return intersect_size(out[v],out[w])-(out[v].count(u)&&out[w].count(u));
 }
 double open_diamond_mass(int u,const std::set<int>& excluded) const {
  // Count each pair once per shared outgoing dependency. Remove linked
  // target pairs and pairs incident to excluded vertices.
  double wedges=diamond_wedges,closed=diamond_closed;
  std::unordered_map<int,double> removed;
  for(int v:excluded){double a=target_weight(v);
   for(int z:out[v]){double& r=removed[z];wedges-=a*(diamond_sum[z]-r-a);r+=a;}
   closed-=diamond_linked[v];
   for(int w:und[v])if(v<w&&excluded.count(w))closed+=linked_wedge_weight(v,w);
  }
  // The common dependency is a fourth vertex, distinct from the releaser.
  // Remove the open pairs in the incoming bucket of u.
  double source_mass=0.,sum=0.;
  for(int v:in[u])if(!excluded.count(v)){
   double a=target_weight(v);source_mass+=a*sum;sum+=a;
   for(int w:out[v])if(in[u].count(w)&&!excluded.count(w)&&(v<w||!out[w].count(v)))
    source_mass-=a*target_weight(w);
  }
  wedges-=source_mass;
  if(wedges<closed)stop("Negative open-diamond candidate mass");
  return wedges-closed;
 }
 std::pair<int,int> sample_open_diamond(int u,const std::set<int>& excluded) const {
  std::vector<double> s=diamond_sum,sq=diamond_sq,cum(s.size());
  for(int v:excluded){double a=target_weight(v);for(int z:out[v]){s[z]-=a;sq[z]-=a*a;}}
  s[u]=sq[u]=0.;
  double total=0.;for(size_t z=0;z<s.size();++z){total+=wedge_mass(s[z],sq[z]);cum[z]=total;}
  auto pick=[&](int z,int omit,bool first){
   double bound=first?2*wedge_mass(s[z],sq[z]):s[z]-target_weight(omit);
   double r=R::runif(0,bound);
   for(int v:in[z])if(v!=omit&&!excluded.count(v)){
    double a=target_weight(v);r-=first?a*(s[z]-a):a;if(r<0.)return v;
   }
   stop("Shared-dependency sampling mass mismatch");return -1;
  };
  for(int attempt=0;attempt<128;++attempt){
   int z=std::upper_bound(cum.begin(),cum.end(),R::runif(0,total))-cum.begin();
   int v=pick(z,-1,true),w=pick(z,v,false);if(!und[v].count(w))return {v,w};
  }
  // Exact fallback on small or unusually dense risk sets.
  double r=R::runif(0,open_diamond_mass(u,excluded));
  for(size_t z=0;z<s.size();++z)if(s[z]*s[z]>sq[z])
   for(int v:in[z])if(!excluded.count(v))for(auto it=in[z].upper_bound(v);it!=in[z].end();++it)
    if(!excluded.count(*it)&&!und[v].count(*it)){r-=target_weight(v)*target_weight(*it);if(r<0.)return {v,*it};}
  stop("Open-diamond sampling mass mismatch");return {-1,-1};
 }
 void enable_consumers(){
  if(consumer_ready)return;
  consumer_ready=true;
  consumer_sum.assign(names.size(),0.);
  consumer_sq.assign(names.size(),0.);
  for(size_t z=0;z<out.size();++z){
   for(int v:out[z]){
    double a=target_weight(v);
    consumer_sum[z]+=a;consumer_sq[z]+=a*a;
   }
   consumer_wedges+=wedge_mass(consumer_sum[z],consumer_sq[z]);
  }
 }
 void consumer_before(const Edges& add,const Edges& del){
  std::unordered_map<int,int> degree;
  for(auto e:add)++degree[e.second];
  for(auto e:del)--degree[e.second];
  std::unordered_map<int,std::array<double,2>> delta;
  for(auto change:degree)if(change.second){
   double a=target_weight(change.first),b=a+change.second;
   for(int z:in[change.first]){
    auto& d=delta[z];d[0]+=b-a;d[1]+=b*b-a*a;
   }
  }
  auto edit=[&](const std::pair<int,int>& e,double sign){
   double a=target_weight(e.second)+degree.at(e.second);
   auto& d=delta[e.first];d[0]+=sign*a;d[1]+=sign*a*a;
  };
  for(auto e:del)edit(e,-1.);
  for(auto e:add)edit(e,1.);
  for(auto change:delta){
   int z=change.first;
   consumer_wedges-=wedge_mass(consumer_sum[z],consumer_sq[z]);
   consumer_sum[z]+=change.second[0];consumer_sq[z]+=change.second[1];
   consumer_wedges+=wedge_mass(consumer_sum[z],consumer_sq[z]);
  }
 }
 double consumer_pair_mass(const std::set<int>& excluded) const {
  double mass=consumer_wedges;
  std::unordered_map<int,double> removed;
  for(int v:excluded){
   double a=target_weight(v);
   for(int z:in[v]){
    double& prior=removed[z];
    mass-=a*(consumer_sum[z]-prior-a);
    prior+=a;
   }
  }
  if(mass<0.)stop("Negative shared-consumer pair mass");
  return mass;
 }
 std::pair<int,int> sample_consumer_pair(const std::set<int>& excluded) const {
  // A pair receives one contribution for each pre-edit consumer that uses it.
  std::vector<double> s=consumer_sum,sq=consumer_sq,cum(s.size());
  for(int v:excluded){
   double a=target_weight(v);
   for(int z:in[v]){s[z]-=a;sq[z]-=a*a;}
  }
  double total=0.;
  for(size_t z=0;z<s.size();++z){total+=wedge_mass(s[z],sq[z]);cum[z]=total;}
  if(total<=0.)stop("Empty shared-consumer pair risk set");
  int z=std::upper_bound(cum.begin(),cum.end(),R::runif(0,total))-cum.begin();
  auto pick=[&](int omit){
   double r=R::runif(0,omit<0?2*wedge_mass(s[z],sq[z]):s[z]-target_weight(omit));
   for(int v:out[z])if(v!=omit&&!excluded.count(v)){
    double a=target_weight(v);r-=omit<0?a*(s[z]-a):a;
    if(r<0.)return v;
   }
   stop("Shared-consumer pair sampling mass mismatch");return -1;
  };
  int v=pick(-1);
  return {v,pick(v)};
 }
 void refresh(int f){auto counts=[&](const std::vector<std::set<int>>& adj){std::unordered_set<int> second;for(int v:adj[f])second.insert(adj[v].begin(),adj[v].end());for(int v:adj[f])second.erase(v);second.erase(f);return std::make_pair(adj[f].size(),second.size());};auto d=counts(fi),u=counts(fo);reach[f]={std::log1p(d.first),std::log1p(u.first),std::log1p(d.second),std::log1p(u.second)};}
 std::set<int> apply(const Edges& add,const Edges& del,bool refresh_reach=true){
  auto validate=[&](const Edges& edits){for(auto e:edits)
   if(e.first<0||e.second<0||e.first>=static_cast<int>(names.size())||e.second>=static_cast<int>(names.size()))
    stop("Graph edit references an absent vertex");};
  validate(add);validate(del);
  std::map<std::pair<int,int>,int> delta;std::set<int> touched;
  auto diamond_affected=diamond_ready?diamond_before(add,del):std::set<int>();
  if(consumer_ready)consumer_before(add,del);
  auto change=[&](int u,int v,int sign){if(u==v || bool(out[u].count(v))!=(sign<0))stop("Inadmissible graph edit");if(sign>0){out[u].insert(v);in[v].insert(u);und[u].insert(v);und[v].insert(u);}else{out[u].erase(v);in[v].erase(u);if(!out[v].count(u)){und[u].erase(v);und[v].erase(u);}}int f=family[u],g=family[v];if(f!=g)delta[{f,g}]+=sign;};
  for(auto e:del)change(e.first,e.second,-1);for(auto e:add)change(e.first,e.second,1);
  for(auto e:delta){int f=e.first.first,g=e.first.second,old=arcs[e.first];arcs[e.first]+=e.second;if(arcs[e.first]<0)stop("Negative family edge count");if(bool(old)==bool(arcs[e.first]))continue;touched.insert(f);touched.insert(g);touched.insert(fi[f].begin(),fi[f].end());touched.insert(fo[g].begin(),fo[g].end());if(arcs[e.first]){fo[f].insert(g);fi[g].insert(f);}else{fo[f].erase(g);fi[g].erase(f);}touched.insert(fi[f].begin(),fi[f].end());touched.insert(fo[g].begin(),fo[g].end());}
  if(diamond_ready)linked_delta(diamond_affected,1.);
  if(refresh_reach)for(int f:touched)refresh(f);return touched;
 }
 CranGraph(List state){
  CharacterVector n=state["names"];NumericVector b=state["birth"],clock=state["clock"];
  if(b.size()!=n.size()||clock.size()!=n.size())stop("Graph birth and clock vectors must match names");
  for(int i=0;i<n.size();++i){
   if(CharacterVector::is_na(n[i])||!std::isfinite(b[i])||!std::isfinite(clock[i]))stop("Graph names and clocks must be finite and nonmissing");
   std::string name=as<std::string>(n[i]);
   if(lookup.count(name))stop("Graph major-line names must be unique");
   node(name,b[i],clock[i]);
  }
  IntegerMatrix edges=state["edges"];Edges add;
  if(edges.ncol()!=2)stop("Graph edges must have two columns");
  for(int i=0;i<edges.nrow();++i)
   if(!add.insert({edges(i,0),edges(i,1)}).second)stop("Graph edges must be unique");
  apply(add,{});
  if(state.containsElementNamed("people")){
   List supplied=state["people"];
   CharacterVector packages=supplied.names();
   if(packages.size()!=supplied.size())stop("Personnel lists must be named by package family");
   has_personnel=true;
   personnel.resize(fn.size());
   std::unordered_set<std::string> seen_packages;
   for(int i=0;i<supplied.size();++i){
    if(CharacterVector::is_na(packages[i]))stop("Personnel package names must not be missing");
    if(!seen_packages.insert(as<std::string>(packages[i])).second)stop("Personnel package names must be unique");
    auto f=flookup.find(as<std::string>(packages[i]));
    if(f==flookup.end())stop("Unknown package family in personnel lists");
    auto& ids=personnel[f->second];
    ids=as<std::vector<int>>(supplied[i]);
    if(std::any_of(ids.begin(),ids.end(),[](int id){return id<=0;}))stop("Invalid person ID");
    std::sort(ids.begin(),ids.end());
    ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
    for(int id:ids)personnel_families[id].push_back(f->second);
   }
  }
 }
 List state(){
  std::vector<std::pair<int,int>> es;
  for(int u=0;u<(int)names.size();++u)for(int v:out[u])es.push_back({u,v});
  IntegerMatrix edges(es.size(),2);
  for(int i=0;i<(int)es.size();++i){edges(i,0)=es[i].first;edges(i,1)=es[i].second;}
  std::vector<double> clocks(names.size());
  for(int u=0;u<(int)names.size();++u)clocks[u]=last[family[u]];
  List result=List::create(_["names"]=names,_["birth"]=born,_["clock"]=clocks,_["edges"]=edges);
  if(has_personnel){
   List values=wrap(personnel);
   values.attr("names")=wrap(fnames);
   result["people"]=values;
  }
  return result;
 }
 std::vector<double> motif(int u,int v,bool deletion,bool directed,const std::unordered_map<int,int>& tr,const std::unordered_map<int,int>& cy){double sign=deletion?-1:1,di=(directed?in[v].size():und[v].size())-deletion,doo=(directed?out[u].size():und[u].size())-deletion;auto get=[&](const std::unordered_map<int,int>& m){auto it=m.find(v);return it==m.end()?0:it->second;};std::vector<double> x;
  if(directed)x={sign*di/100,sign*doo/100,sign*choose2(di)/10000,sign*choose2(doo)/10000,sign*get(tr)/10,sign*get(cy)/10};else{x={sign*(di+doo)/100,sign*(choose2(di)+choose2(doo))/10000,sign*get(tr)/10};if(out[v].count(u))std::fill(x.begin(),x.end(),0);}return x;}
 void counts(int u,bool directed,std::unordered_map<int,int>& tr,std::unordered_map<int,int>& cy){if(directed){for(int w:out[u]){for(int v:in[w])tr[v]++;for(int v:out[w])tr[v]++;}for(int w:in[u]){for(int v:out[w])tr[v]++;for(int v:in[w])cy[v]++;}}else for(int w:und[u])for(int v:und[w])tr[v]++;}
 std::unordered_map<int,int> diamond_completions(int u) const {
  // The focal line u is the top of the open diamond. Count each shared
  // dependency separately; all four vertices are distinct.
  std::unordered_map<int,int> ans;
  for(int w:out[u])for(int z:out[w])if(z!=u)
   for(int v:in[z])if(v!=u&&v!=w&&!und[w].count(v))++ans[v];
  return ans;
 }
 int common_consumers(int u,int v,int w) const {
  return intersect_size(in[v],in[w])-(in[v].count(u)&&in[w].count(u));
 }
 std::vector<int> consumer_completions(int u) const {
  std::vector<int> counts(names.size()),seen(names.size(),-1);
  for(int w:out[u])for(int z:in[w])if(z!=u)
   for(int v:out[z])if(v!=u&&v!=w&&seen[v]!=w){++counts[v];seen[v]=w;}
  return counts;
 }
 std::vector<double> pair_delta(int u,int old,int v){double du=und[u].size(),doo=und[old].size(),dn=und[v].size(),d=!in[u].count(old),a=!in[u].count(v),du1=du-d;return {(-d*(du+doo-2)+a*(du1+dn))/100,(-d*(choose2(du-1)+choose2(doo-1))+a*(choose2(du1)+choose2(dn)))/10000,(-d*intersect_size(und[u],und[old])+a*(intersect_size(und[u],und[v])-d*und[v].count(old)))/10};}
 std::array<double,4> directed_pair_delta(int u,int old,int v,const std::unordered_map<int,int>& tr,const std::unordered_map<int,int>& cy){
  auto get=[](const std::unordered_map<int,int>& x,int j){auto it=x.find(j);return it==x.end()?0:it->second;};
  // A migration preserves its sender's out-degree, so both out-star changes
  // vanish. The two transitive configurations containing both candidate edges
  // must be removed from the addition score after deleting the old edge.
  return {((double)in[v].size()-in[old].size()+1)/100,
          (choose2(in[v].size())-choose2(in[old].size()-1))/10000,
          (get(tr,v)-get(tr,old)-(int)out[old].count(v)-(int)out[v].count(old))/10.,
          (get(cy,v)-get(cy,old))/10.};
 }
 static int categorical(const std::vector<double>& cumulative,double sum){double q=R::runif(0,sum);return std::lower_bound(cumulative.begin(),cumulative.end(),q)-cumulative.begin();}
 static std::vector<double> cumweights(const std::vector<double>& logw){if(logw.empty())return {};double mx=*std::max_element(logw.begin(),logw.end()),s=0;std::vector<double> c;for(double w:logw){s+=std::exp(w-mx);c.push_back(s);}return c;}
 double triangle_weight(int v,int w) const {return (1.+in[v].size())*(1.+in[w].size());}
 double triangle_incident_mass(const std::set<int>& vertices) const {
  double mass=0.;for(int v:vertices)for(int w:und[v])if(v<w||!vertices.count(w))mass+=triangle_weight(v,w);
  return mass;
 }
 template<typename F> void triangle_pairs(int u,F visit) const {
  for(int v=0;v<(int)und.size();++v)if(v!=u&&!out[u].count(v))
   for(int w:und[v])if(v<w&&w!=u&&!out[u].count(w))if(visit(v,w))return;
 }
 int previous_major(int u,int oldn) const {
  int p=-1;for(int v:fn[family[u]])if(v<oldn&&major[v]<major[u]&&(p<0||major[v]>major[p]))p=v;
  return p;
 }
 Update draw(double t,const std::string& requested,const std::vector<double>& parameters,const std::vector<double>& rp){
  std::string requested_law=requested;std::vector<double> supplied=parameters;
  std::vector<double> release=rp;
  std::array<double,4> activity_shift={0.,0.,0.,0.};
  int component=1;
  const std::string mixture_suffix="-mixture2";
  if(requested_law.size()>=mixture_suffix.size()&&requested_law.compare(requested_law.size()-mixture_suffix.size(),mixture_suffix.size(),mixture_suffix)==0){
   requested_law.resize(requested_law.size()-mixture_suffix.size());
   if(release.size()!=19||supplied.size()<4)stop("Invalid two-component mixture parameters");
   for(double value:release)if(!std::isfinite(value))stop("Nonfinite mixture release parameter");
   for(size_t j=supplied.size()-4;j<supplied.size();++j)if(!std::isfinite(supplied[j]))stop("Nonfinite mixture edit-rate shift");
   // The component is drawn once for the complete update. Each component's
   // release and edit laws retain their own nonempty conditioning below.
   component=R::runif(0,1)<logistic(release[18])?2:1;
   if(component==2){
    for(int j=0;j<3;++j)release[j]+=release[15+j];
    std::copy(supplied.end()-4,supplied.end(),activity_shift.begin());
   }
   supplied.resize(supplied.size()-4);release.resize(15);
  }
  if(release.size()!=15)stop("Supply 15 release parameters or a complete two-component mixture");
  if(!std::isfinite(t))stop("Update time must be finite");
  for(double value:release)if(!std::isfinite(value))stop("Release parameters must be finite");
  for(double value:supplied)if(!std::isfinite(value))stop("Edit parameters must be finite");
  double release_total=std::exp(release[0])+std::exp(release[1])+std::exp(release[2]);
  if(!std::isfinite(release_total)||release_total<=0)stop("Invalid total release proposal rate");
  if(fn.empty())stop("CRAN release selection requires an initial package family");
  const std::string capped_suffix="-capped2";
  bool capped_pairs=requested_law.size()>=capped_suffix.size()&&requested_law.compare(requested_law.size()-capped_suffix.size(),capped_suffix.size(),capped_suffix)==0;
  if(capped_pairs)requested_law.resize(requested_law.size()-capped_suffix.size());
  const std::string poisson_suffix="-poisson";
  bool poisson_pairs=requested_law.size()>=poisson_suffix.size()&&requested_law.compare(requested_law.size()-poisson_suffix.size(),poisson_suffix.size(),poisson_suffix)==0;
  if(poisson_pairs)requested_law.resize(requested_law.size()-poisson_suffix.size());
  const std::string degree_suffix="-olddegree";
  bool olddegree=requested_law.size()>=degree_suffix.size()&&requested_law.compare(requested_law.size()-degree_suffix.size(),degree_suffix.size(),degree_suffix)==0;
  double degree_effect=0;
  if(olddegree){requested_law.resize(requested_law.size()-degree_suffix.size());degree_effect=supplied.back();supplied.pop_back();}
  const std::string people_suffix="-people";
  bool people_model=requested_law.size()>=people_suffix.size()&&
   requested_law.compare(requested_law.size()-people_suffix.size(),people_suffix.size(),people_suffix)==0;
  std::array<double,3> people_effects{};
  if(people_model){
   if(!has_personnel)stop("Shared-person forecasts require personnel at the forecast origin");
   if(supplied.size()<3)stop("Missing shared-person parameters");
   std::copy(supplied.end()-3,supplied.end(),people_effects.begin());
   supplied.resize(supplied.size()-3);
   requested_law.resize(requested_law.size()-people_suffix.size());
  }
  const std::string consumer_suffix="-consumers";
  bool consumers=requested_law.size()>=consumer_suffix.size()&&requested_law.compare(requested_law.size()-consumer_suffix.size(),consumer_suffix.size(),consumer_suffix)==0;
  const std::string suffix=consumers?consumer_suffix:"-patterns";
  bool patterns=requested_law.size()>=suffix.size()&&requested_law.compare(requested_law.size()-suffix.size(),suffix.size(),suffix)==0;
  std::string law=patterns?requested_law.substr(0,requested_law.size()-suffix.size()):requested_law;
  std::vector<double> input=supplied,pair_preferences,consumer_rates;
  if(consumers&&law=="CS-pair-migration-add"){
   if(input.size()!=37)stop("Invalid shared-consumer paired-addition parameter vector");
   consumer_rates.assign(input.end()-3,input.end());input.resize(input.size()-3);
   enable_consumers();
  }
  if(patterns&&(law=="CS-pair"||law=="CS-pair-migration-add")){
   if(input.size()<3)stop("Missing independent pair preferences");
   pair_preferences.assign(input.end()-3,input.end());input.resize(input.size()-3);
  }
  bool additions_only=law=="CS-pair-add";
  bool migration_add=law=="CS-pair-migration-add";
  bool diamond_model=law=="CS-pair-diamond"||additions_only||migration_add,triangle_model=law=="CS-pair-triangle"||diamond_model;
  bool paired_model=law=="CS-pair"||triangle_model,cs=law=="CS-1"||paired_model;
  if((poisson_pairs||capped_pairs)&&!triangle_model)stop("Poisson addition pairs require a paired-addition law");
  if(triangle_model&&input.size()!=(diamond_model?31:28))stop("Invalid paired-addition parameter vector");
  if(diamond_model)enable_diamonds();
  int expbase=paired_model?20:15,legacybase=law=="PA"?8:paired_model?27:law=="CS-1"?21:16;
  bool exponential=triangle_model||(cs&&(input.size()==(size_t)expbase||input.size()==(size_t)expbase+4));
  bool tergm_logdegree=law=="TERGM-logdegree";
  bool tergm=law=="TERGM"||tergm_logdegree;
  if(tergm && input.size()!=(tergm_logdegree?23:19))stop("Invalid TERGM parameter vector");
  int base=exponential?expbase:legacybase;bool varying=triangle_model||(!tergm&&input.size()==(size_t)base+4);
  if(!tergm&&!triangle_model&&input.size()!=static_cast<size_t>(base)&&input.size()!=static_cast<size_t>(base+4))
   stop("Invalid edge parameter vector");
  if(law!="PA"&&law!="CS-1"&&law!="CS-1U"&&!paired_model&&!tergm)stop("Unknown CRAN update law");
  input[0]+=activity_shift[0];input[tergm?10:1]+=activity_shift[1];
  if(paired_model&&!additions_only)input[exponential?19:26]+=activity_shift[2];
  if(triangle_model)for(int j=24;j<27;++j)input[j]+=activity_shift[3];
  if(diamond_model)for(int j=28;j<31;++j)input[j]+=activity_shift[3];
  for(double& rate:consumer_rates)rate+=activity_shift[3];
  std::vector<double> p=input;
  if(exponential){p.assign(paired_model?27:21,0);p[0]=input[0];p[1]=input[1];int a[]={3,5,7,8,9,10,11},d[]={13,15,17,18};for(int j=0;j<7;++j)p[a[j]]=input[2+j];for(int j=0;j<4;++j)p[d[j]]=input[9+j];p[19]=input[13];p[20]=input[14];if(paired_model){for(int j=0;j<4;++j)p[22+j]=input[15+j];p[26]=input[19];}}
  Update ans;ans.component=component;int oldn=names.size(),nf=fn.size();std::vector<double> lm(nf),lu(nf);for(int f=0;f<nf;++f){double age=(t-last[f])/365.25;double x[6]={age,age*age,reach[f][0],reach[f][1],reach[f][2],reach[f][3]};for(int j=0;j<6;++j){lm[f]+=release[3+j]*x[j];lu[f]+=release[9+j]*x[j];}}auto cm=cumweights(lm),cu=cumweights(lu);
  int nnew,nm,nu;do{nnew=R::rpois(std::exp(release[0]));nm=R::rpois(std::exp(release[1]));nu=R::rpois(std::exp(release[2]));}while(nnew+nm+nu==0);
  std::set<int> majors,selected;for(int i=0;i<nm;++i)majors.insert(categorical(cm,cm.back()));for(int i=0;i<nu;++i){int f=categorical(cu,cu.back());selected.insert(fn[f][(int)R::runif(0,fn[f].size())]);}
  std::map<int,int> pred;for(int f:majors){int v=*std::max_element(fn[f].begin(),fn[f].end(),[&](int a,int b){return major[a]<major[b];});int u=node(fnames[f]+"@"+std::to_string(major[v]+1),t,t);pred[u]=v;ans.born.push_back(u);selected.insert(u);}
  for(int i=0;i<nnew;++i){std::string name;do{name="__hawkesNet_sim_"+std::to_string(++synthetic)+"@1";}while(lookup.count(name));int u=node(name,t,t);ans.born.push_back(u);selected.insert(u);}
  ans.sources.assign(selected.begin(),selected.end());
  for(int u:ans.sources)ans.release_class.push_back(u<oldn?0:family[u]>=nf?1:2);
  bool pa=law=="PA",directed=cs||tergm;int na=pa?3:directed?10:7,nd=pa?1:directed?7:5,ageA=2+na+nd,ageD=ageA+1;
  for(int u:ans.sources){int predecessor=-1;if(pred.count(u))predecessor=pred[u];else for(int v:fn[family[u]])if(v<oldn&&major[v]<major[u]&&(predecessor<0||major[v]>major[predecessor]))predecessor=v;
   std::vector<bool> shares(people_model?fn.size():0,false);
   if(people_model){
    const auto& source=personnel[family[u]];
    for(int id:source)for(int f:personnel_families.at(id))shares[f]=true;
   }
   std::set<int> mig,copy;for(int v:out[u])for(int w:fn[family[v]])if(major[w]>major[v])mig.insert(w);if(out[u].empty()&&predecessor>=0)copy=out[predecessor];
   double add_shift=0.,delete_shift=0.;
   if(tergm_logdegree){
    double log_degree=std::log1p(out[u].size());
    add_shift=p[19]*(u>=oldn&&family[u]>=nf)+p[20]*(u>=oldn&&family[u]<nf)+p[21]*log_degree;
    delete_shift=p[22]*log_degree;
   }
   std::unordered_map<int,int> tr,cy;if(!pa&&!people_model)counts(u,directed,tr,cy);auto dm=patterns&&!people_model?diamond_completions(u):std::unordered_map<int,int>();std::vector<int> targetsA,targetsD;std::vector<double> logA,logD;
   auto cm=consumers&&!people_model?consumer_completions(u):std::vector<int>();
   for(int v=0;v<(int)names.size();++v){if(v==u)continue;bool del=out[u].count(v);double d=in[v].size(),b=(born[v]-365.)/365.25; // days since 2019 -> birth years since 2020
    double z=0;int offset=tergm?(del?10:0):(del?2+na:2);
    if(people_model){
     if(del)z=(tergm?p[10]:0.)-p[tergm?16:18]*std::log(std::max(1.,d));
     else z=(tergm?p[0]:0.)+p[tergm?7:9]*std::log1p(d)+
      p[tergm?8:10]*mig.count(v)+p[tergm?9:11]*copy.count(v);
    }
    else if(pa){z=p[offset]*std::log1p(d);if(!del)z+=p[offset+1]*mig.count(v)+p[offset+2]*copy.count(v);}
    else{double x[10]={1,0,0,0,0,0,0,0,0,0},sign=del?-1:1;auto ti=tr.find(v),ci=cy.find(v);double tri=ti==tr.end()?0:ti->second,cycle=patterns?dm[v]:(ci==cy.end()?0:ci->second);
      if(directed){double di=d-del,doo=out[u].size()-del;x[1]=sign*di/100;x[2]=sign*doo/100;x[3]=sign*choose2(di)/10000;x[4]=sign*choose2(doo)/10000;x[5]=sign*tri/10;if(del)x[6]=-std::log(std::max(1.,d));else{x[6]=cycle/10;x[7]=std::log1p(d);x[8]=mig.count(v);x[9]=copy.count(v);}}
      else{double di=und[v].size()-del,doo=und[u].size()-del;if(!out[v].count(u)){x[1]=sign*(di+doo)/100;x[2]=sign*(choose2(di)+choose2(doo))/10000;x[3]=sign*tri/10;}x[4]=del?-std::log(std::max(1.,d)):std::log1p(d);if(!del){x[5]=mig.count(v);x[6]=copy.count(v);}}
      if(consumers&&!del)x[3]=cm[v]/10.;
      for(int j=0;j<(del?nd:na);++j)z+=p[offset+j]*x[j];}
    double lw;if(tergm)lw=z+(del?delete_shift:add_shift)-p[del?18:17]*(t-born[v])/365.25;else lw=(pa||exponential?z:-sp(-z))+p[del?ageD:ageA]*b;
    if(olddegree&&!paired_model&&!del&&mig.count(v)){
     double degree=0;for(int old:out[u])if(family[old]==family[v]&&major[old]<major[v])degree+=in[old].size();
     lw+=degree_effect*std::log1p(degree);
    }
    if(people_model&&shares[family[v]]){
     lw+=people_effects[del?1:0];
     if(!del&&!paired_model&&mig.count(v))lw+=people_effects[2];
    }
    (del?targetsD:targetsA).push_back(v);(del?logD:logA).push_back(lw);
   }
   if(tergm){std::vector<double> logits=logA;logits.insert(logits.end(),logD.begin(),logD.end());double H=0;for(double z:logits)H+=sp(z);bool condition=u<oldn;double threshold=condition?-std::log1p(-R::runif(0,1)*(-std::expm1(-H))):INFINITY;double running=0;bool first=false;
    for(int j=0;j<(int)logits.size();++j){bool yes=false;if(condition&&!first){running+=sp(logits[j]);if(running>=threshold){yes=true;first=true;}}else yes=R::runif(0,1)<logistic(logits[j]);if(yes){if(j<(int)targetsA.size())ans.add.insert({u,targetsA[j]});else ans.del.insert({u,targetsD[j-targetsA.size()]});}}
   }else{
    auto ca=cumweights(logA),cd=cumweights(logD);std::vector<std::pair<int,int>> pairs;std::vector<double> lp;bool paired=paired_model&&!additions_only;
    if(paired)for(int old:out[u])for(int v:fn[family[old]])if(v!=u&&!out[u].count(v)&&major[v]>major[old]){
     double z=p[21];
     if(!people_model){
      auto x=directed_pair_delta(u,old,v,tr,cy);
      if(patterns)x[3]=(dm[v]-dm[old]-(und[old].count(v)?0:open_diamond_partners(u,old,v)))/10.;
      if(consumers)x[1]=(cm[v]-cm[old]-(common_consumers(u,old,v)>0))/10.;
      for(int j=0;j<4;++j)z+=p[22+j]*x[j];
     }
     double popA=patterns?pair_preferences[0]:p[9],popD=patterns?pair_preferences[1]:p[18],age=patterns?pair_preferences[2]:p[19];
     z+=popA*std::log1p(in[v].size())-popD*std::log((double)in[old].size());
     lp.push_back((exponential?z:-sp(-z))+age*(born[v]-365.)/365.25);pairs.push_back({old,v});
    }
    double logma=p[0],logmd=p[1];if(varying){double degree=std::log1p(out[u].size());logma+=input[base]*(u>=oldn&&family[u]>=nf)+input[base+1]*(u>=oldn&&family[u]<nf)+input[base+2]*degree;logmd+=input[base+3]*degree;}
    double degree_multiplier=1.;
    if((olddegree||people_model)&&!lp.empty()){
     double center=*std::max_element(lp.begin(),lp.end()),before=0,after=0;
     for(size_t j=0;j<lp.size();++j){
      before+=std::exp(lp[j]-center);
      if(olddegree)lp[j]+=degree_effect*std::log1p(in[pairs[j].first].size());
      if(people_model&&shares[family[pairs[j].second]])lp[j]+=people_effects[2];
      after+=std::exp(lp[j]-center);
     }
     degree_multiplier=after/before;
    }
    auto cp=cumweights(lp);double ma=ca.empty()?0:std::exp(logma),md=cd.empty()?0:std::exp(logmd),mp=cp.empty()?0:std::exp(p[26])*degree_multiplier,total=ma+md+mp;if(total<=0||!std::isfinite(total))stop("Invalid edit proposal rate");
    bool triangle_draw=false;
    if(triangle_model){
     bool active=false;triangle_pairs(u,[&](int,int){active=true;return true;});
     int release_class=u<oldn?0:family[u]>=nf?1:2;
     double mt=active?std::exp(input[24+release_class]):0.,mdia=0.,mc=0.;
     std::set<int> excluded=out[u];excluded.insert(u);
     if(diamond_model&&std::exp(input[28+release_class])>0){
      if(open_diamond_mass(u,excluded)>0.)mdia=std::exp(input[28+release_class]);
     }
     if(!consumer_rates.empty()&&std::exp(consumer_rates[release_class])>0.&&consumer_pair_mass(excluded)>0.)
      mc=std::exp(consumer_rates[release_class]);
     double pair_rate=mt+mdia+mc;
     if(!std::isfinite(pair_rate))stop("Invalid addition-pair proposal rate");
     if(pair_rate>0){
      double probability=-std::expm1(-pair_rate)/(u<oldn?-std::expm1(-total-pair_rate):1.);
      triangle_draw=probability>0&&R::runif(0,1)<probability;
     }
     // Conditional on activation, an unrestricted channel has a zero-truncated
     // Poisson count. Inverse survival sampling remains efficient for tiny rates.
     double proposed=triangle_draw?(poisson_pairs?R::qpois(R::runif(0,1)*(-std::expm1(-pair_rate)),pair_rate,false,false):1):0;
     if(triangle_draw&&capped_pairs)proposed=R::runif(0,1)<pair_rate/std::expm1(pair_rate)?1:2;
     if(!std::isfinite(proposed)||proposed>INT_MAX)stop("Addition-pair proposal count exceeds numerical capacity");
     int pair_draws=(int)proposed;
     std::vector<std::pair<int,int>> candidate;std::vector<double> mass;
     if(pair_draws>0&&mt>0.){
      double accum=0.,copy_boost=(additions_only||migration_add)?1.:std::exp(input[27]);
      triangle_pairs(u,[&](int v,int w){
       accum+=triangle_weight(v,w)*(copy.count(v)&&copy.count(w)?copy_boost:1.);
       candidate.emplace_back(v,w);mass.push_back(accum);return false;
      });
     }
     for(int draw=0;draw<pair_draws;++draw){
      double type=mc>0.?R::runif(0,pair_rate):0.;
      if(mc>0.&&type<mc){
       auto edge=sample_consumer_pair(excluded);ans.add.insert({u,edge.first});ans.add.insert({u,edge.second});
      }else if(mdia>0.&&(mt==0.||(mc>0.?type<mc+mdia:R::runif(0,mt+mdia)<mdia))){
       auto edge=sample_open_diamond(u,excluded);ans.add.insert({u,edge.first});ans.add.insert({u,edge.second});
      }else{
      auto edge=candidate[categorical(mass,mass.back())];ans.add.insert({u,edge.first});ans.add.insert({u,edge.second});
      }
      if(poisson_pairs||capped_pairs)ans.pairs++;
     }
    }
    if(total>1000){
     // Exact Poisson splitting avoids iterating an arbitrarily large number
     // of repeated proposals when the covariate-dependent means are large.
     bool have;do{have=false;
      auto singles=[&](const std::vector<double>& logw,const std::vector<double>& cumulative,const std::vector<int>& targets,double mean,Edges& edits){if(logw.empty())return;double mx=*std::max_element(logw.begin(),logw.end()),lz=mx+std::log(cumulative.back());for(int j=0;j<(int)logw.size();++j)if(R::runif(0,1)<-std::expm1(-mean*std::exp(logw[j]-lz))){edits.insert({u,targets[j]});have=true;}};
      singles(logA,ca,targetsA,ma,ans.add);singles(logD,cd,targetsD,md,ans.del);
      if(!lp.empty()){double mx=*std::max_element(lp.begin(),lp.end()),lz=mx+std::log(cp.back());for(int j=0;j<(int)lp.size();++j){int k=R::rpois(mp*std::exp(lp[j]-lz));if(k){ans.del.insert({u,pairs[j].first});ans.add.insert({u,pairs[j].second});ans.pairs+=k;have=true;}}}
     }while(u<oldn&&!triangle_draw&&!have);continue;
    }
    int count;do{count=R::rpois(total);}while(u<oldn&&!triangle_draw&&count==0);
    for(int i=0;i<count;++i){double which=R::runif(0,total);if(which<ma)ans.add.insert({u,targetsA[categorical(ca,ca.back())]});else if(which<ma+md)ans.del.insert({u,targetsD[categorical(cd,cd.back())]});else{auto pair=pairs[categorical(cp,cp.back())];ans.del.insert({u,pair.first});ans.add.insert({u,pair.second});ans.pairs++;}}
   }
  }return ans;
 }
 List summary(){int n=names.size();std::map<int,int> ih,oh,uh,esp;double ne=0,triples=0,espsum=0;int maxin=0,maxout=0;std::vector<double> pr(n,1./n);
  for(int u=0;u<n;++u){ih[in[u].size()]++;oh[out[u].size()]++;uh[und[u].size()]++;ne+=out[u].size();triples+=choose2(und[u].size());maxin=std::max(maxin,(int)in[u].size());maxout=std::max(maxout,(int)out[u].size());for(int v:und[u])if(u<v){int s=intersect_size(und[u],und[v]);esp[s]++;espsum+=s;}}
  auto histogram=[](const std::map<int,int>& h){NumericVector v(h.size());CharacterVector names(h.size());int k=0;for(auto x:h){v[k]=x.second;names[k++]=std::to_string(x.first);}v.attr("names")=names;return v;};
  for(int it=0;it<1000;++it){double dangling=0;for(int u=0;u<n;++u)if(out[u].empty())dangling+=pr[u];std::vector<double> next(n,.15/n+.85*dangling/n);for(int u=0;u<n;++u)if(!out[u].empty())for(int v:out[u])next[v]+=.85*pr[u]/out[u].size();double error=0;for(int u=0;u<n;++u)error+=std::abs(next[u]-pr[u]);pr.swap(next);if(error<1e-11)break;if(it==999)stop("PageRank did not converge");}
  auto gini=[](std::vector<double> x){std::sort(x.begin(),x.end());double total=0,weighted=0;for(int i=0;i<(int)x.size();++i){total+=x[i];weighted+=(i+1)*x[i];}return total>0?2*weighted/(x.size()*total)-(x.size()+1.)/x.size():0.;};
  std::vector<double> indeg,outdeg;for(int u=0;u<n;++u){indeg.push_back(in[u].size());outdeg.push_back(out[u].size());}double pr2=0;for(double q:pr)pr2+=q*q;
  std::vector<bool> seen(n,false);int giant=0;for(int s=0;s<n;++s)if(!seen[s]){std::vector<int> stack(1,s);seen[s]=true;int count=0;while(!stack.empty()){int u=stack.back();stack.pop_back();count++;for(int v:und[u])if(!seen[v]){seen[v]=true;stack.push_back(v);}}giant=std::max(giant,count);}
  std::vector<int> sources(n);std::iota(sources.begin(),sources.end(),0);std::mt19937 rng(427);std::shuffle(sources.begin(),sources.end(),rng);int landmarks=std::min(32,n);std::vector<double> bc(n,0);
  for(int j=0;j<landmarks;++j){int s=sources[j];std::vector<int> distance(n,-1),order;std::vector<double> sigma(n,0),dependency(n,0);std::vector<std::vector<int>> pred(n);std::queue<int> queue;queue.push(s);distance[s]=0;sigma[s]=1;
    while(!queue.empty()){int u=queue.front();queue.pop();order.push_back(u);for(int v:out[u]){if(distance[v]<0){distance[v]=distance[u]+1;queue.push(v);}if(distance[v]==distance[u]+1){sigma[v]+=sigma[u];pred[v].push_back(u);}}}
    for(auto it=order.rbegin();it!=order.rend();++it){int w=*it;for(int v:pred[w])dependency[v]+=sigma[v]/sigma[w]*(1+dependency[w]);if(w!=s)bc[w]+=dependency[w];}}
  double bcscale=(double)n/landmarks/std::max(1.,(n-1.)*(n-2.)),bcmean=0,bcmax=0;for(double b:bc){bcmean+=b*bcscale/n;bcmax=std::max(bcmax,b*bcscale);}
  NumericVector metrics=NumericVector::create(_["nodes"]=n,_["edges"]=ne,_["max_indegree"]=maxin,_["max_outdegree"]=maxout,_["transitivity"]=triples>0?espsum/triples:0,_["triangles"]=espsum/3,_["indegree_centralization"]=(n*(double)maxin-ne)/std::max(1.,(n-1.)*(n-1.)),_["indegree_gini"]=gini(indeg),_["outdegree_gini"]=gini(outdeg),_["pagerank_gini"]=gini(pr),_["pagerank_effective_fraction"]=1/(n*pr2));
  metrics.push_back(ne/n,"mean_indegree");metrics.push_back(ne/n,"mean_outdegree");metrics.push_back((double)giant/n,"weak_giant_fraction");metrics.push_back((double)uh[0]/n,"isolated_fraction");metrics.push_back(bcmean,"mean_betweenness");metrics.push_back(bcmax,"max_betweenness");metrics.push_back(landmarks,"betweenness_landmarks");
  return List::create(_["metrics"]=metrics,_["distributions"]=List::create(_["indegree"]=histogram(ih),_["outdegree"]=histogram(oh),_["undirected_degree"]=histogram(uh),_["esp"]=histogram(esp)));
 }
};
// [[Rcpp::export]]
SEXP cran_graph_create(List state){return XPtr<CranGraph>(new CranGraph(state),true);}
// [[Rcpp::export]]
NumericVector cran_directed_pair_change(SEXP state,int sender,int old,int target){
 XPtr<CranGraph> g(state);if(!g->out[sender].count(old)||g->out[sender].count(target)||sender==target||old==target)stop("Inadmissible migration pair");
 std::unordered_map<int,int> tr,cy;g->counts(sender,true,tr,cy);auto d=g->directed_pair_delta(sender,old,target,tr,cy);
 return NumericVector::create(d[0],0.,d[1],0.,d[2],d[3]);
}
// [[Rcpp::export]]
List cran_directed_pair_cache(List state,List updates,CharacterVector sources,IntegerVector batches){
 CranGraph g(state);std::vector<std::array<double,13>> values;std::vector<std::array<double,12>> observed;int group=0;
 for(int b=0;b<updates.size();++b){
  List rec=updates[b];double t=rec["time"];CharacterVector births=rec["births"];NumericVector clocks=rec["clocks"];
  for(int j=0;j<births.size();++j)g.node(as<std::string>(births[j]),t,clocks[j]);
  CharacterMatrix aa=rec["add"],dd=rec["delete"];Edges add,del;
  for(int j=0;j<aa.nrow();++j)add.insert({g.lookup.at(as<std::string>(aa(j,0))),g.lookup.at(as<std::string>(aa(j,1)))});
  for(int j=0;j<dd.nrow();++j)del.insert({g.lookup.at(as<std::string>(dd(j,0))),g.lookup.at(as<std::string>(dd(j,1)))});
  while(group<sources.size()&&batches[group]==b){
   int u=g.lookup.at(as<std::string>(sources[group]));std::unordered_map<int,int> tr,cy;g.counts(u,true,tr,cy);
   auto lexical=[&](int a,int c){return g.names[a]<g.names[c];};std::vector<int> old(g.out[u].begin(),g.out[u].end());std::sort(old.begin(),old.end(),lexical);
   for(int v:old){auto candidates=g.fn[g.family[v]];std::sort(candidates.begin(),candidates.end(),lexical);
    for(int w:candidates)if(w!=u&&!g.out[u].count(w)&&g.major[w]>g.major[v]){
     auto d=g.directed_pair_delta(u,v,w,tr,cy);auto und=g.pair_delta(u,v,w);int row=values.size();double birth=(g.born[w]-365.)/365.25;
     values.push_back({(double)group,d[0],0.,d[1],0.,d[2],d[3],birth,std::log1p(g.in[w].size()),-std::log((double)g.in[v].size()),und[0],und[1],und[2]});
     if(add.count({u,w})&&del.count({u,v}))observed.push_back({(double)row,(double)group,(double)g.in[w].size(),(double)g.out[u].size(),(double)tr[w],(double)cy[w],birth,(double)g.in[v].size(),(double)g.out[u].size(),(double)tr[v],(double)cy[v],(g.born[v]-365.)/365.25});
    }
   }++group;
  }
  g.apply(add,del);if(b%100==0)checkUserInterrupt();
 }
 if(group!=sources.size())stop("Incomplete migration candidate replay");
 NumericMatrix v(values.size(),13),o(observed.size(),12);for(int i=0;i<v.nrow();++i)for(int j=0;j<13;++j)v(i,j)=values[i][j];for(int i=0;i<o.nrow();++i)for(int j=0;j<12;++j)o(i,j)=observed[i][j];
 return List::create(_["values"]=v,_["observed"]=o);
}
// Exact focal open-diamond scores for single additions. The new cache uses
// the vacated cyclic-addition coordinate; archived caches are never changed.
// [[Rcpp::export]]
NumericVector cran_diamond_completion(SEXP state,int sender) {
 XPtr<CranGraph> g(state);auto scores=g->diamond_completions(sender);
 NumericVector ans(g->names.size());for(auto x:scores)ans[x.first]=x.second;return ans;
}
// [[Rcpp::export]]
List cran_pattern_preference_cache(List state,List updates,CharacterVector sources,IntegerVector batches,std::string rawpath,bool consumers=false) {
 CranGraph g(state);std::ofstream raw(rawpath,std::ios::binary|std::ios::trunc);
 if(!raw)stop("Cannot create pattern preference cache");
 if(consumers)g.enable_consumers();
 std::map<uint64_t,uint32_t> ids;std::vector<uint64_t> codes;std::map<double,uint32_t> birthids;std::vector<double> births;
 std::vector<double> offset(sources.size()+1);uint64_t n=0;int group=0;double maxdiamond=0;
 std::vector<std::array<double,11>> observed;
 std::vector<std::array<double,8>> history;
 std::vector<std::array<double,2>> pair_observed;
 std::vector<std::array<double,5>> triangles,diamonds,consumer_pairs;
 std::vector<double> pairdiamond,pairconsumer;
 NumericVector consumer_mass(consumers?sources.size():0);
 int consumer_checks=0;
 for(int b=0;b<updates.size();++b){
  List rec=updates[b];double t=rec["time"];CharacterVector newborn=rec["births"];NumericVector clocks=rec["clocks"];
  int oldn=g.names.size();for(int j=0;j<newborn.size();++j)g.node(as<std::string>(newborn[j]),t,clocks[j]);
  auto read=[&](CharacterMatrix x){Edges e;for(int j=0;j<x.nrow();++j)e.emplace(g.lookup.at(as<std::string>(x(j,0))),g.lookup.at(as<std::string>(x(j,1))));return e;};
  Edges add=read(as<CharacterMatrix>(rec["add"])),del=read(as<CharacterMatrix>(rec["delete"]));
  while(group<sources.size()&&batches[group]==b){
   int u=g.lookup.at(as<std::string>(sources[group]));std::unordered_map<int,int> tr,cy;g.counts(u,true,tr,cy);auto dm=g.diamond_completions(u);
   auto cm=consumers?g.consumer_completions(u):std::vector<int>();
   if(consumers){
    std::set<int> excluded=g.out[u];excluded.insert(u);
    consumer_mass[group]=g.consumer_pair_mass(excluded);
    if(group%5000==0){
     double direct=0.;
     for(const auto& targets:g.out){
      double sum=0.;
      for(int v:targets)if(!excluded.count(v)){
       double a=g.target_weight(v);direct+=a*sum;sum+=a;
      }
     }
     if(direct!=consumer_mass[group])stop("Shared-consumer normaliser disagrees with direct sum");
     ++consumer_checks;
    }
   }
   std::set<int> mig,copy,candidates;for(int v:g.out[u])for(int w:g.fn[g.family[v]])if(g.major[w]>g.major[v]&&w!=u&&!g.out[u].count(w))mig.insert(w);
   if(g.out[u].empty()){int p=g.previous_major(u,oldn);if(p>=0)copy=g.out[p];copy.erase(u);}
   if(!consumers){for(auto x:tr)candidates.insert(x.first);for(auto x:dm)candidates.insert(x.first);}
   auto birth=[&](int v){return (g.born[v]-365.)/365.25;};
   auto score=[](const std::unordered_map<int,int>& counts,int v){
    auto it=counts.find(v);return it==counts.end()?0:it->second;
   };
   auto write_candidate=[&](int v){
    int triangle=score(tr,v),diamond=score(dm,v),co=consumers?cm[v]:0;
    if(!(triangle||diamond||co)||v==u||g.out[u].count(v))return;
    int degree=g.in[v].size(),field=consumers?co:g.out[u].size();
    maxdiamond=std::max(maxdiamond,(double)diamond);
    if(degree>=65535||field>=65535||triangle>=65535||diamond>=65535)stop("Pattern score exceeds exact cache encoding; no truncation performed");
    uint64_t code=((uint64_t)degree<<48)|((uint64_t)field<<32)|((uint64_t)triangle<<16)|(uint64_t)diamond;
    auto si=ids.emplace(code,ids.size());if(si.second)codes.push_back(code);
    auto bi=birthids.emplace(birth(v),birthids.size());if(bi.second)births.push_back(birth(v));
    uint32_t row[2]={si.first->second,bi.first->second};raw.write((char*)row,sizeof(row));++n;
   };
   if(consumers){for(int v=0;v<(int)cm.size();++v)write_candidate(v);}
   else for(int v:candidates)write_candidate(v);
   std::set<int> special=mig;special.insert(copy.begin(),copy.end());
   for(int v:special)if(v!=u&&!g.out[u].count(v))history.push_back({(double)group,(double)g.in[v].size(),birth(v),consumers?(double)cm[v]:(double)g.out[u].size(),(double)tr[v],(double)dm[v],(double)mig.count(v),(double)copy.count(v)});
   std::map<int,int> observed_index;
   for(auto it=add.lower_bound({u,-1});it!=add.end()&&it->first==u;++it){int v=it->second;
    observed_index[v]=observed.size();
    observed.push_back({(double)group,(double)g.in[v].size(),(double)g.out[u].size(),(double)tr[v],(double)cy[v],birth(v),(double)mig.count(v),(double)copy.count(v),(double)dm[v],(double)v,consumers?(double)cm[v]:0.});
   }
   auto lexical=[&](int v,int w){return g.names[v]<g.names[w];};std::vector<int> old(g.out[u].begin(),g.out[u].end());std::sort(old.begin(),old.end(),lexical);
   for(int v:old){auto newer=g.fn[g.family[v]];std::sort(newer.begin(),newer.end(),lexical);
    for(int w:newer)if(w!=u&&!g.out[u].count(w)&&g.major[w]>g.major[v]){
     int interaction=g.und[v].count(w)?0:g.open_diamond_partners(u,v,w);
     if(add.count({u,w})&&del.count({u,v}))pair_observed.push_back({(double)pairdiamond.size(),(double)observed_index.at(w)});
     pairdiamond.push_back((dm[w]-dm[v]-interaction)/10.);
     if(consumers)pairconsumer.push_back((cm[w]-cm[v]-(g.common_consumers(u,v,w)>0))/10.);
    }
   }
   for(auto a=observed_index.begin();a!=observed_index.end();++a)for(auto c=std::next(a);c!=observed_index.end();++c){
    int v=a->first,w=c->first;double weight=g.triangle_weight(v,w);
    if(g.und[v].count(w))triangles.push_back({(double)group,(double)a->second,(double)c->second,weight,(double)(copy.count(v)&&copy.count(w))});
    else {int count=g.open_diamond_partners(u,v,w);if(count)diamonds.push_back({(double)group,(double)a->second,(double)c->second,weight*count,0.});}
    if(consumers){
     int count=g.common_consumers(u,v,w);
     if(count)consumer_pairs.push_back({(double)group,(double)a->second,(double)c->second,weight*count,0.});
    }
   }
   offset[++group]=n;
  }
  g.apply(add,del,false);if(b%100==0)checkUserInterrupt();
 }
 if(group!=sources.size()||!raw)stop("Incomplete pattern preference cache");raw.close();
 RawVector packed(codes.size()*8);std::memcpy(packed.begin(),codes.data(),packed.size());
 NumericMatrix o(observed.size(),consumers?11:10),h(history.size(),5),flags(history.size(),2);NumericVector hg(history.size());
 for(int i=0;i<o.nrow();++i)for(int j=0;j<o.ncol();++j)o(i,j)=observed[i][j];
 for(int i=0;i<h.nrow();++i){hg[i]=history[i][0];for(int j=0;j<5;++j)h(i,j)=history[i][j+1];flags(i,0)=history[i][6];flags(i,1)=history[i][7];}
 NumericMatrix po(pair_observed.size(),2),tt(triangles.size(),5),dd(diamonds.size(),5);
 for(int i=0;i<po.nrow();++i)for(int j=0;j<2;++j)po(i,j)=pair_observed[i][j];
 for(int i=0;i<tt.nrow();++i)for(int j=0;j<5;++j)tt(i,j)=triangles[i][j];
 for(int i=0;i<dd.nrow();++i)for(int j=0;j<5;++j)dd(i,j)=diamonds[i][j];
 List result=List::create(_["index"]=List::create(_["offsets"]=offset,_["birth"]=births),_["codes"]=packed,
   _["history"]=List::create(_["group"]=hg,_["features"]=h,_["flags"]=flags),_["observed"]=o,
   _["pair_diamond"]=pairdiamond,_["pair_observed"]=po,_["triangles"]=tt,_["diamonds"]=dd,
   _["max_diamond"]=maxdiamond,_["candidate_records"]=(double)n,
   _["schema"]=consumers?"shared-consumer-v1":"single-pattern-preferences-v2");
 if(consumers){
  List index=result["index"];index["shared_consumer"]=true;result["index"]=index;
  NumericMatrix cc(consumer_pairs.size(),5);
  for(int i=0;i<cc.nrow();++i)for(int j=0;j<5;++j)cc(i,j)=consumer_pairs[i][j];
  result["consumer_pairs"]=cc;
  result["consumer_mass"]=consumer_mass;
  result["consumer_checks"]=consumer_checks;
  result["pair_consumer"]=pairconsumer;
 }
 return result;
}

// Exact candidate masses for the paired-addition channel. Only pairs fully
// contained in the observed additions need individual records.
// [[Rcpp::export]]
List cran_triangle_pair_cache(List state,List updates,CharacterVector sources,IntegerVector batches,bool diamond=false){
 CranGraph g(state);double graph_mass=0.;
 if(diamond)g.enable_diamonds();
 for(int v=0;v<(int)g.und.size();++v)for(int w:g.und[v])if(v<w)graph_mass+=g.triangle_weight(v,w);
 int group=0,migration=0,checks=0;NumericMatrix mass(sources.size(),2);
 std::vector<std::array<double,9>> edits;
 std::vector<std::array<double,5>> observed;
 NumericVector diamond_mass(diamond?sources.size():0);
 std::vector<std::array<double,4>> observed_diamonds;
 for(int b=0;b<updates.size();++b){
  List rec=updates[b];double t=rec["time"];CharacterVector births=rec["births"];NumericVector clocks=rec["clocks"];
  int oldn=g.names.size();for(int j=0;j<births.size();++j)g.node(as<std::string>(births[j]),t,clocks[j]);
  auto read=[&](CharacterMatrix m){Edges a;for(int j=0;j<m.nrow();++j)a.emplace(g.lookup.at(as<std::string>(m(j,0))),g.lookup.at(as<std::string>(m(j,1))));return a;};
  Edges add=read(as<CharacterMatrix>(rec["add"])),del=read(as<CharacterMatrix>(rec["delete"]));
  while(group<sources.size()&&batches[group]==b){
   int u=g.lookup.at(as<std::string>(sources[group]));
   std::set<int> excluded=g.out[u],copy;excluded.insert(u);
   mass(group,0)=graph_mass-g.triangle_incident_mass(excluded);
   if(mass(group,0)<0.)stop("Negative eligible pair mass");
   if(g.out[u].empty()){
    int p=g.previous_major(u,oldn);if(p>=0)copy=g.out[p];copy.erase(u);
    for(int v:copy)for(int w:g.und[v])if(v<w&&copy.count(w))mass(group,1)+=g.triangle_weight(v,w);
   }
   if(mass(group,1)>mass(group,0))stop("Invalid copied pair mass");
   if(diamond)diamond_mass[group]=g.open_diamond_mass(u,excluded);
   if(group%5000==0){double direct=0.,both=0.;g.triangle_pairs(u,[&](int v,int w){double x=g.triangle_weight(v,w);direct+=x;if(copy.count(v)&&copy.count(w))both+=x;return false;});
    if(direct!=mass(group,0)||both!=mass(group,1))stop("Triangle normaliser disagrees with direct sum");
    if(diamond){double raw=0.,closed=0.;for(int z=0;z<(int)g.in.size();++z)if(z!=u){double sum=0.;for(int v:g.in[z])if(!excluded.count(v)){double a=g.target_weight(v);raw+=a*sum;sum+=a;}}
     g.triangle_pairs(u,[&](int v,int w){closed+=g.linked_wedge_weight(v,w)-g.triangle_weight(v,w)*(g.out[v].count(u)&&g.out[w].count(u));return false;});
     if(raw-closed!=diamond_mass[group])stop("Open-diamond normaliser disagrees with direct sum");
    }
    ++checks;
   }
   std::unordered_map<int,int> eligible_migration;
   auto lex=[&](int v,int w){return g.names[v]<g.names[w];};
   std::vector<int> removed,added;
   for(auto it=del.lower_bound({u,-1});it!=del.end()&&it->first==u;++it)removed.push_back(it->second);
   for(auto it=add.lower_bound({u,-1});it!=add.end()&&it->first==u;++it)added.push_back(it->second);
   std::sort(removed.begin(),removed.end(),lex);std::sort(added.begin(),added.end(),lex);
   for(int v:removed)for(int w:added)if(g.family[v]==g.family[w]&&g.major[w]>g.major[v]){
    if(!eligible_migration.emplace(w,migration++).second)stop("Overlapping observed migration pairs");
   }
   if(!added.empty()){
    std::unordered_map<int,int> tr,cy;g.counts(u,true,tr,cy);
    int first=edits.size();
    for(int w:added){
     bool mig=false;for(int v:g.out[u])if(g.family[v]==g.family[w]&&g.major[v]<g.major[w]){mig=true;break;}
     edits.push_back({(double)group,(double)g.in[w].size(),(double)g.out[u].size(),(double)tr[w],(double)cy[w],
        (g.born[w]-365.)/365.25,(double)mig,(double)copy.count(w),eligible_migration.count(w)?(double)eligible_migration.at(w):-1.});
    }
    for(int j=0;j<(int)added.size();++j)for(int k=j+1;k<(int)added.size();++k){
     int v=added[j],w=added[k];
     if(g.und[v].count(w))observed.push_back({(double)group,(double)(first+j),(double)(first+k),g.triangle_weight(v,w),(double)(copy.count(v)&&copy.count(w))});
     else if(diamond){int common=g.open_diamond_partners(u,v,w);if(common)observed_diamonds.push_back({(double)group,(double)(first+j),(double)(first+k),g.triangle_weight(v,w)*common});}
    }
   }
   ++group;
  }
  std::set<int> targets;for(auto e:add)targets.insert(e.second);for(auto e:del)targets.insert(e.second);
  graph_mass-=g.triangle_incident_mass(targets);g.apply(add,del,false);graph_mass+=g.triangle_incident_mass(targets);
  if(b%100==0)checkUserInterrupt();
 }
 if(group!=sources.size())stop("Incomplete paired-addition replay");
 NumericMatrix e(edits.size(),9),p(observed.size(),5);
 for(int i=0;i<e.nrow();++i)for(int j=0;j<9;++j)e(i,j)=edits[i][j];
 for(int i=0;i<p.nrow();++i)for(int j=0;j<5;++j)p(i,j)=observed[i][j];
 List result=List::create(_["mass"]=mass,_["edits"]=e,_["observed"]=p,_["normaliser_checks"]=checks);
 if(diamond){NumericMatrix d(observed_diamonds.size(),4);
  for(int i=0;i<d.nrow();++i)for(int j=0;j<4;++j)d(i,j)=observed_diamonds[i][j];
  result["diamond_mass"]=diamond_mass;result["observed_diamonds"]=d;
 }
 return result;
}
// [[Rcpp::export]]
List cran_graph_snapshot(SEXP state){XPtr<CranGraph> g(state);return g->state();}
// [[Rcpp::export]]
List cran_graph_summary(SEXP state){XPtr<CranGraph> g(state);return g->summary();}
// [[Rcpp::export]]
List cran_graph_step(SEXP state,double time,std::string law,NumericVector parameters,NumericVector release_parameters,bool apply=true){XPtr<CranGraph> g(state);auto u=g->draw(time,law,as<std::vector<double>>(parameters),as<std::vector<double>>(release_parameters));if(apply)g->apply(u.add,u.del);auto matrix=[](const Edges& e){IntegerMatrix x(e.size(),2);int i=0;for(auto p:e){x(i,0)=p.first;x(i++,1)=p.second;}return x;};return List::create(_["sources"]=u.sources,_["born"]=u.born,_["add"]=matrix(u.add),_["delete"]=matrix(u.del),_["pair_proposals"]=u.pairs,_["mixture_component"]=u.component,_["release_class"]=u.release_class);}
// Topology of an earlier release's own major line or family.
enum class CranTopology { FamilyIn, MajorIn, FamilyOut };
static int timing_size(const CranGraph& g,CranTopology type){
 return type==CranTopology::MajorIn?g.names.size():g.fn.size();
}
static int timing_anchor(const CranGraph& g,int u,CranTopology type){
 return type==CranTopology::MajorIn?u:g.family[u];
}
static double timing_summary(const CranGraph& g,int id,CranTopology type){
 if(type==CranTopology::MajorIn)return std::log1p(g.in[id].size());
 return std::log1p(type==CranTopology::FamilyOut?g.fo[id].size():g.fi[id].size());
}
static std::map<int,double> timing_before(const CranGraph& g,const Edges& add,const Edges& del,CranTopology type){
 std::map<int,double> before;
 for(const Edges* edges:{&add,&del})for(auto e:*edges){
  int id=timing_anchor(g,type==CranTopology::FamilyOut?e.first:e.second,type);
  before.emplace(id,timing_summary(g,id,type));
 }
 return before;
}
// [[Rcpp::export]]
List cran_graph_replay(SEXP state,List updates,bool timing_details=false){
 XPtr<CranGraph> g(state);List results(updates.size());
 std::array<CranTopology,2> types={CranTopology::MajorIn,CranTopology::FamilyOut};
 std::set<int> dirty;
 if(timing_details){
  List initial(2);
  for(int k=0;k<2;++k){
   NumericVector x(timing_size(*g,types[k]));
   for(int j=0;j<x.size();++j)x[j]=timing_summary(*g,j,types[k]);
   initial[k]=x;
  }
  results.attr("timing_initial")=initial;
 }
 for(int b=0;b<updates.size();++b){
  List rec=updates[b];double t=rec["time"];CharacterVector born=rec["births"];
  int oldn=g->names.size(),oldf=g->fn.size();
  NumericVector clocks=rec.containsElementNamed("clocks")?as<NumericVector>(rec["clocks"]):NumericVector(born.size(),t);
  for(int j=0;j<born.size();++j)g->node(as<std::string>(born[j]),t,clocks[j]);
  CharacterMatrix aa=rec["add"],dd=rec["delete"];Edges add,del;
  for(int i=0;i<aa.nrow();++i)add.insert({g->lookup.at(as<std::string>(aa(i,0))),g->lookup.at(as<std::string>(aa(i,1)))});
  for(int i=0;i<dd.nrow();++i)del.insert({g->lookup.at(as<std::string>(dd(i,0))),g->lookup.at(as<std::string>(dd(i,1)))});
  std::array<std::map<int,double>,2> before;
  if(timing_details)for(int k=0;k<2;++k)before[k]=timing_before(*g,add,del,types[k]);
  // The observed replay reads degrees, not release-selection reach covariates.
  // Refresh those once before returning the final graph to its caller.
  auto touched=g->apply(add,del,!timing_details);
  if(timing_details)dirty.insert(touched.begin(),touched.end());
  CharacterVector sources=rec["sources"];
  IntegerVector family(sources.size()),lines(sources.size()),release_class(sources.size());NumericVector s(sources.size());
  for(int j=0;j<sources.size();++j){
   int u=g->lookup.at(as<std::string>(sources[j]));lines[j]=u;
   family[j]=g->family[u];release_class[j]=u<oldn?0:family[j]>=oldf?1:2;s[j]=timing_summary(*g,family[j],CranTopology::FamilyIn);
  }
  List result=List::create(_["time"]=t,_["families"]=family,_["summary"]=s,_["release_class"]=release_class);
  if(timing_details){
   List topo(2);
   for(int k=0;k<2;++k){
    IntegerVector ids=k==0?lines:family;NumericVector v(ids.size());
    std::vector<int> cu;std::vector<double> cz;
    for(int j=0;j<ids.size();++j)v[j]=timing_summary(*g,ids[j],types[k]);
    for(auto prev:before[k]){
     double z=timing_summary(*g,prev.first,types[k]);
     if(z!=prev.second){cu.push_back(prev.first);cz.push_back(z);}
    }
    topo[k]=List::create(_["id"]=ids,_["summary"]=v,_["cu"]=cu,_["cz"]=cz);
   }
   result["topology"]=topo;
  }
  results[b]=result;if(b%100==0)checkUserInterrupt();
 }
 for(int f:dirty)g->refresh(f);
 return results;
}
// Evaluate directed first creations on the post-update graph.
static std::array<bool,4> cran_directed_hits(const CranGraph& g,const Edges& add,const Edges& del){
 std::array<bool,4> hit={false,false,false,false};
 for(auto a:add){
  int u=a.first,v=a.second;
  for(auto d:del)if(d.first==u && g.family[d.second]==g.family[v] && g.major[d.second]<g.major[v])hit[0]=true;
  if(!hit[1] && (intersect_size(g.out[u],g.in[v]) || intersect_size(g.out[u],g.out[v]) || intersect_size(g.in[u],g.in[v])))hit[1]=true;
  // The new edge is an upper edge of u -> {v,w} -> x.
  if(!hit[2])for(int w:g.out[u])if(w!=v && !g.und[v].count(w)){
   for(int x:g.out[v])if(x!=u && g.out[w].count(x)){hit[2]=true;break;}
   if(hit[2])break;
  }
  // The new edge is a lower edge: {u,w} -> v, with common root z.
  if(!hit[2])for(int w:g.in[v])if(w!=u && !g.und[u].count(w)){
   for(int z:g.in[u])if(z!=v && g.in[w].count(z)){hit[2]=true;break;}
   if(hit[2])break;
  }
  // u and z both use v and w: all four vertices are distinct.
  if(!hit[3])for(int z:g.in[v])if(z!=u){
   for(int w:g.out[u])if(w!=v && g.out[z].count(w)){hit[3]=true;break;}
   if(hit[3])break;
  }
 }
 // Removing the last link between middle vertices can open a diamond.
 if(!hit[2])for(auto d:del){
  int v=d.first,w=d.second;if(g.und[v].count(w))continue;
  for(int u:g.in[v])if(g.in[w].count(u)){
   for(int x:g.out[v])if(x!=u && g.out[w].count(x)){hit[2]=true;break;}
   if(hit[2])break;
  }
  if(hit[2])break;
 }
 return hit;
}

// [[Rcpp::export]]
List cran_joint_path(List state,List history,double origin,double horizon,std::string law,NumericVector parameters,NumericVector release_parameters,NumericVector timing,NumericVector calendar,std::string model="M0",int event_limit=10000,bool directed_waits=false,Nullable<List> migration_cohort=R_NilValue,double clip=1.3862943611198906){
 if(timing.size()!=8 && timing.size()!=10 && timing.size()!=12)stop("Invalid joint timing dimensions");
 if(!std::isfinite(clip)||clip<0)stop("Clipping bound must be finite and nonnegative");
 bool typed=timing.size()>8, type_by_kernel=timing.size()==12;
 if(typed && model=="M0")stop("M0 uses unit amplitudes");
 auto class_offset=[&](int c,int k){
  if(c<0 || c>2)stop("Invalid forecast release class");
  return typed && c>0?timing[7+2*(type_by_kernel?k:0)+c]:0.;
 };
 CranGraph g(state);auto p=as<std::vector<double>>(parameters),rp=as<std::vector<double>>(release_parameters);double mu=timing[0],alpha[2]={timing[1],timing[2]},tau[2]={timing[3],timing[4]},delta=timing[5],gb=timing[6],gd=timing[7];
 bool dynamic=model=="M2"||model=="M2_major_in"||model=="M2_family_out";
 if(model!="M0"&&model!="M1"&&!dynamic)stop("Invalid CRAN timing model");
 CranTopology type=model=="M2_major_in"?CranTopology::MajorIn:model=="M2_family_out"?CranTopology::FamilyOut:CranTopology::FamilyIn;
 int feature=type==CranTopology::MajorIn?0:type==CranTopology::FamilyOut?1:-1;
 struct Term{double time,summary;std::array<double,2> base;};
 auto term_for=[&](double t,double summary,double initial,double scale,int c){
  return Term{t,summary,{scale*std::exp(gb*initial+class_offset(c,0)),
                         scale*std::exp(gb*initial+class_offset(c,1))}};
 };
 std::vector<std::vector<Term>> terms(timing_size(g,type));double q[2]={0,0};
 auto multiplier=[&](const Term& term,double current){
  double change=std::max(-clip,std::min(clip,current-term.summary));
  return std::exp(gd*change);
 };
 for(int j=0;j<history.size();++j){
  List rec=history[j];double t=rec["time"];if(t>origin || (t==origin && migration_cohort.isNull()))stop("Future event supplied to forecast history");
  if(model=="M0"){
   for(int k=0;k<2;++k)q[k]+=std::exp(-(origin-t)/tau[k]);
  }else{
   IntegerVector ids=rec["families"];NumericVector s=rec["summary"],ref=s;
   if(feature>=0){
    List all=rec["topology"],x=all[feature];
    ids=as<IntegerVector>(x["id"]);ref=as<NumericVector>(x["summary"]);
    if(ids.size()!=s.size()||ref.size()!=s.size())stop("Invalid forecast topology history");
   }
   if(ids.size()!=s.size())stop("Invalid forecast history dimension");
   IntegerVector release_class;
   if(typed){
    if(!rec.containsElementNamed("release_class"))stop("Missing forecast release classes");
    release_class=as<IntegerVector>(rec["release_class"]);
    if(release_class.size()!=s.size())stop("Invalid forecast release class dimension");
   }
   double scale=std::pow((double)s.size(),delta-1);
   for(int l=0;l<ids.size();++l){
    int f=ids[l];Term term=term_for(t,ref[l],s[l],scale,typed?release_class[l]:0);double a=1.;
    if(dynamic){
     if(f<0||f>=static_cast<int>(terms.size()))stop("Forecast topology history references an absent anchor");
     terms[f].push_back(term);a=multiplier(term,timing_summary(g,f,type));
    }
    for(int k=0;k<2;++k)q[k]+=std::exp(-(origin-t)/tau[k])*(term.base[k]*a);
   }
  }
 }
 int target=-1;
 std::set<int> exposed,migrated;
 std::vector<double> migration_times;
 std::vector<std::string> migration_families;
 if(migration_cohort.isNotNull()){
  List cohort(migration_cohort);
  target=g.lookup.at(as<std::string>(cohort["target"]));
  CharacterVector names=cohort["exposed"];
  for(auto name:names)exposed.insert(g.flookup.at(as<std::string>(name)));
 }
 int origin_day=std::floor(origin);if(origin_day<0||origin_day>=calendar.size())stop("Calendar does not cover origin");double initial_rate=mu*calendar[origin_day];for(int k=0;k<2;++k)initial_rate+=alpha[k]*q[k]/tau[k];
 double t=origin,end=origin+horizon;NumericVector waits(4,horizon);LogicalVector found(4,false),directed_found(4,false);NumericVector directed_wait(4,horizon);
 std::vector<double> events;std::vector<std::array<int,4>> sizes;
 auto decay=[&](double dt){for(int k=0;k<2;++k)q[k]*=std::exp(-dt/tau[k]);};
 while(t<end){int day=std::floor(t);if(day<0||day>=calendar.size())stop("Calendar does not cover forecast");double boundary=std::min(end,day+1.),base=mu*calendar[day],bound=base;for(int k=0;k<2;++k)bound+=alpha[k]*q[k]/tau[k];if(!std::isfinite(bound)||bound<0)stop("Invalid thinning bound");if(bound==0){decay(boundary-t);t=boundary;continue;}double wait=R::rexp(1/bound);if(t+wait>=boundary){decay(boundary-t);t=boundary;continue;}decay(wait);t+=wait;double rate=base;for(int k=0;k<2;++k)rate+=alpha[k]*q[k]/tau[k];if(R::runif(0,1)>rate/bound)continue;
  Update u=g.draw(t,law,p,rp);std::set<int> touched;for(auto e:u.add){touched.insert(e.first);touched.insert(e.second);}for(auto e:u.del){touched.insert(e.first);touched.insert(e.second);}
  std::map<int,std::set<int>> before;for(int v:touched)before[v]=g.und[v];auto old=dynamic?timing_before(g,u.add,u.del,type):std::map<int,double>();
  bool migration=false;
  for(auto d:u.del)for(auto a:u.add)if(d.first==a.first&&g.family[d.second]==g.family[a.second]&&g.major[a.second]>g.major[d.second]){
   migration=true;
   int f=g.family[a.first];
   if(a.second==target && exposed.count(f) && migrated.insert(f).second){
    migration_times.push_back(t-origin);
    migration_families.push_back(g.fnames[f]);
   }
  }
  g.apply(u.add,u.del);bool triangle=false,star2=false,star3=false;for(auto a:u.add)if(!before[a.first].count(a.second)&&intersect_size(g.und[a.first],g.und[a.second])>0)triangle=true;
  for(int v:touched){double survivors=intersect_size(before[v],g.und[v]);if(choose2(g.und[v].size())>choose2(survivors))star2=true;if(choose3(g.und[v].size())>choose3(survivors))star3=true;}
  bool hits[4]={migration,triangle,star2,star3};for(int j=0;j<4;++j)if(!found[j]&&hits[j]){found[j]=true;waits[j]=t-origin;}
  if(directed_waits){
   auto directed=cran_directed_hits(g,u.add,u.del);
   for(int j=0;j<4;++j)if(!directed_found[j]&&directed[j]){directed_found[j]=true;directed_wait[j]=t-origin;}
  }
  if(dynamic){
   terms.resize(timing_size(g,type));
   for(auto prev:old){
    double now=timing_summary(g,prev.first,type);if(now==prev.second)continue;
    for(const auto& term:terms[prev.first]){
     double change=multiplier(term,now)-multiplier(term,prev.second);
     for(int k=0;k<2;++k)q[k]+=std::exp(-(t-term.time)/tau[k])*term.base[k]*change;
    }
   }
  }
  if(model=="M0"){for(int k=0;k<2;++k)q[k]+=1;}
  else{
   double scale=std::pow((double)u.sources.size(),delta-1);
   for(size_t l=0;l<u.sources.size();++l){
    int v=u.sources[l];
    int f=timing_anchor(g,v,type);double s=timing_summary(g,g.family[v],CranTopology::FamilyIn);
    Term term=term_for(t,timing_summary(g,f,type),s,scale,u.release_class[l]);
    if(dynamic)terms[f].push_back(term);
    for(int k=0;k<2;++k)q[k]+=term.base[k];
   }
  }
  events.push_back(t);sizes.push_back({(int)u.born.size(),(int)u.add.size(),(int)u.del.size(),(int)u.sources.size()});if((int)events.size()>=event_limit)stop("Forecast event guard reached; no censoring substituted");if(events.size()%100==0)checkUserInterrupt();
  if(target<0 && is_true(all(found)) && (!directed_waits || is_true(all(directed_found))))break;
 }
 NumericMatrix counts(sizes.size(),4);for(int i=0;i<(int)sizes.size();++i)for(int j=0;j<4;++j)counts(i,j)=sizes[i][j];waits.attr("names")=CharacterVector::create("Migration","Triangle","Two-star","Three-star");List result=List::create(_["initial_intensity"]=initial_rate,_["wait"]=waits,_["found"]=found,_["times"]=events,_["sizes"]=counts);
 if(directed_waits){
  directed_wait.attr("names")=CharacterVector::create("Migration","Transitive","Diamond","Consumer");
  result["directed_wait"]=directed_wait;
 }
 if(target>=0)result["migration_cohort"]=List::create(_["time"]=migration_times,_["family"]=migration_families);
 return result;
}

// Evaluate edit laws on a common draw from an autonomous M0 clock. The input
// times are simulated independently of the marks and never use future data.
// [[Rcpp::export]]
List cran_shared_clock_waits(List state,NumericVector times,double origin,double horizon,
                            std::string law,NumericVector parameters,NumericVector release_parameters){
 CranGraph g(state);auto p=as<std::vector<double>>(parameters),rp=as<std::vector<double>>(release_parameters);
 NumericVector waits(3,horizon);IntegerVector ranks(3,NA_INTEGER);
 waits.attr("names")=ranks.attr("names")=CharacterVector::create("Migration","Triangle","Transitive");
 double previous=origin;int first_diamond=0;
 for(int i=0;i<times.size();++i){
  double t=times[i];if(t<=previous||t>=origin+horizon)stop("Invalid common forecast times");previous=t;
  Update u=g.draw(t,law,p,rp);bool hits[3]={false,false,false};
  if(i==0)for(auto a=u.add.begin();a!=u.add.end();++a)
   for(auto b=std::next(a);b!=u.add.end()&&b->first==a->first;++b)
    if(g.open_diamond_partners(a->first,a->second,b->second))++first_diamond;
  if(ranks[0]==NA_INTEGER)for(auto d:u.del)for(auto a:u.add)
   if(d.first==a.first&&g.family[d.second]==g.family[a.second]&&g.major[a.second]>g.major[d.second])hits[0]=true;
  Edges new_und;
  if(ranks[1]==NA_INTEGER)for(auto a:u.add)if(!g.und[a.first].count(a.second))new_und.insert(a);
  g.apply(u.add,u.del);
  for(auto a:new_und)if(intersect_size(g.und[a.first],g.und[a.second])){hits[1]=true;break;}
  if(ranks[2]==NA_INTEGER)for(auto a:u.add){int v=a.first,w=a.second;
   if(intersect_size(g.out[v],g.in[w])||intersect_size(g.out[v],g.out[w])||intersect_size(g.in[v],g.in[w])){hits[2]=true;break;}
  }
  for(int k=0;k<3;++k)if(ranks[k]==NA_INTEGER&&hits[k]){ranks[k]=i+1;waits[k]=t-origin;}
  if(std::all_of(ranks.begin(),ranks.end(),[](int r){return r!=NA_INTEGER;}))break;
  if(i%100==0)checkUserInterrupt();
 }
 return List::create(_["wait"]=waits,_["rank"]=ranks,_["first_diamond"]=first_diamond);
}

// Old-target degrees in precisely the history/observed order of the pattern cache.
// [[Rcpp::export]]
List cran_migration_degree_cache(List state,List updates,CharacterVector sources,IntegerVector batches){
 CranGraph g(state);std::vector<double> history,observed;int group=0;
 for(int b=0;b<updates.size();++b){
  List rec=updates[b];double t=rec["time"];CharacterVector newborn=rec["births"];NumericVector clocks=rec["clocks"];
  int oldn=g.names.size();for(int j=0;j<newborn.size();++j)g.node(as<std::string>(newborn[j]),t,clocks[j]);
  auto read=[&](CharacterMatrix x){Edges e;for(int j=0;j<x.nrow();++j)e.emplace(g.lookup.at(as<std::string>(x(j,0))),g.lookup.at(as<std::string>(x(j,1))));return e;};
  Edges add=read(as<CharacterMatrix>(rec["add"])),del=read(as<CharacterMatrix>(rec["delete"]));
  while(group<sources.size()&&batches[group]==b){
   int u=g.lookup.at(as<std::string>(sources[group]));std::map<int,double> degree;
   for(int v:g.out[u])for(int w:g.fn[g.family[v]])if(g.major[w]>g.major[v]&&w!=u&&!g.out[u].count(w))degree[w]+=g.in[v].size();
   std::set<int> special;for(auto x:degree)special.insert(x.first);
   if(g.out[u].empty()){int p=g.previous_major(u,oldn);if(p>=0)special.insert(g.out[p].begin(),g.out[p].end());special.erase(u);}
   for(int v:special)if(v!=u&&!g.out[u].count(v))history.push_back(std::log1p(degree[v]));
   for(auto it=add.lower_bound({u,-1});it!=add.end()&&it->first==u;++it)observed.push_back(std::log1p(degree[it->second]));
   ++group;
  }
  g.apply(add,del,false);if(b%100==0)checkUserInterrupt();
 }
 if(group!=sources.size())stop("Incomplete migration degree replay");
 return List::create(_["history"]=history,_["observed"]=observed,_["schema"]="migration-old-degree-v1");
}

// Shared-person candidates are sparse. Replay all dated changes to personnel,
// including publications that do not change the dependency graph.
// [[Rcpp::export]]
List cran_shared_people_cache(List state,List updates,CharacterVector sources,
                              IntegerVector batches,List history){
 CranGraph g(state);
 IntegerVector bins=history["bin"];
 CharacterVector packages=history["package"];
 List people=history["people"];
 std::map<std::string,std::vector<int>> active;
 std::unordered_map<int,std::set<std::string>> by_person;
 std::vector<std::array<double,5>> shared_add,deletion;
 std::vector<std::array<double,6>> observed_add;
 std::vector<std::array<double,7>> pairs;
 int h=0,group=0;
 for(int b=0;b<updates.size();++b){
  List record=updates[b];double time=record["time"];
  const int bin=2577168+std::llround(time*144.); // 1 January 2019 in 600-second bins
  while(h<bins.size() && bins[h]<=bin){
   auto package=as<std::string>(packages[h]);
   auto& old=active[package];
   for(int p:old)by_person[p].erase(package);
   old=as<std::vector<int>>(people[h]);
   for(int p:old)by_person[p].insert(package);
   ++h;
  }
  CharacterVector newborn=record["births"];NumericVector clocks=record["clocks"];
  int old_size=g.names.size();
  for(int j=0;j<newborn.size();++j)g.node(as<std::string>(newborn[j]),time,clocks[j]);
  auto edges=[&](CharacterMatrix x){
   Edges result;
   for(int j=0;j<x.nrow();++j)
    result.emplace(g.lookup.at(as<std::string>(x(j,0))),g.lookup.at(as<std::string>(x(j,1))));
   return result;
  };
  Edges add=edges(as<CharacterMatrix>(record["add"])),del=edges(as<CharacterMatrix>(record["delete"]));
  while(group<sources.size() && batches[group]==b){
   int u=g.lookup.at(as<std::string>(sources[group]));
   std::set<int> shared_families,migration,copy;
   for(int p:active[g.fnames[g.family[u]]])for(const auto& name:by_person[p]){
    auto family=g.flookup.find(name);
    if(family!=g.flookup.end())shared_families.insert(family->second);
   }
   auto shared=[&](int v){return (double)shared_families.count(g.family[v]);};
   auto birth=[&](int v){return (g.born[v]-365.)/365.25;};
   auto popularity=[&](int v){return std::log1p(g.in[v].size());};
   auto eligible=[&](int old,int v){return v!=u && !g.out[u].count(v) && g.major[v]>g.major[old];};
   for(int old:g.out[u])for(int v:g.fn[g.family[old]])if(eligible(old,v))migration.insert(v);
   if(g.out[u].empty()){
    int predecessor=g.previous_major(u,old_size);
    if(predecessor>=0)copy=g.out[predecessor];
   }
   for(int f:shared_families)for(int v:g.fn[f])if(v!=u && !g.out[u].count(v))
    shared_add.push_back({(double)group,popularity(v),birth(v),
                          (double)migration.count(v),(double)copy.count(v)});
   std::map<int,int> add_rows,delete_rows;
   for(auto it=add.lower_bound({u,-1});it!=add.end() && it->first==u;++it){
    int v=it->second;
    add_rows[v]=observed_add.size();
    observed_add.push_back({(double)group,popularity(v),birth(v),
                            (double)migration.count(v),(double)copy.count(v),shared(v)});
   }
   auto lexical=[&](int v,int w){return g.names[v]<g.names[w];};
   std::vector<int> old(g.out[u].begin(),g.out[u].end());
   std::sort(old.begin(),old.end(),lexical);
   for(int v:old){
    double removed=del.count({u,v});
    if(removed)delete_rows[v]=deletion.size();
    deletion.push_back({(double)group,(double)g.in[v].size(),birth(v),removed,shared(v)});
   }
   for(int v:old){
    auto newer=g.fn[g.family[v]];
    std::sort(newer.begin(),newer.end(),lexical);
    for(int w:newer)if(eligible(v,w)){
     bool observed=delete_rows.count(v) && add_rows.count(w);
     pairs.push_back({(double)group,birth(w),popularity(w),-std::log((double)g.in[v].size()),
                      shared(w),observed?(double)delete_rows.at(v):-1.,
                      observed?(double)add_rows.at(w):-1.});
    }
   }
   ++group;
  }
  g.apply(add,del,false);
  if(b%100==0)checkUserInterrupt();
 }
 if(group!=sources.size())stop("Incomplete shared-person replay");
 auto matrix=[](const auto& rows,int columns){
  NumericMatrix output(rows.size(),columns);
  for(size_t i=0;i<rows.size();++i)for(int j=0;j<columns;++j)output(i,j)=rows[i][j];
  return output;
 };
 return List::create(_["shared_add"]=matrix(shared_add,5),_["observed_add"]=matrix(observed_add,6),
                     _["deletion"]=matrix(deletion,5),_["pairs"]=matrix(pairs,7));
}
