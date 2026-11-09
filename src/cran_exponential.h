// Exact positive-sum normalisation for exponential directed CS weights.
// A segment tree holds target degree/birth cohorts. Removing excluded and
// specially scored candidates uses integer counts, avoiding cancellation of
// nearly equal masses. All weights and derivatives are combined in log space.
struct CranExpCache {
 struct Cohort {int degree;double birth;};
 struct Extra {int cohort,tr,cy,mig,copy,count,consumer;};
 std::vector<Cohort> cohorts;
 std::map<std::pair<int,double>,int> lookup;
 std::vector<int> initial,oldid,newid,update_batch,group_batch;
 std::vector<std::vector<std::pair<int,int>>> removed;
 std::vector<std::vector<Extra>> extras;
 CranCache* source;std::vector<bool> prepared;
 int id(int d,double b){auto key=std::make_pair(d,b);auto it=lookup.find(key);if(it!=lookup.end())return it->second;int i=cohorts.size();lookup[key]=i;cohorts.push_back({d,b});return i;}
 CranExpCache(CranCache& c):source(&c){
  for(int i=0;i<c.ni;++i)initial.push_back(id(c.id[i],c.ib[i]));
  for(int i=0;i<c.nu;++i){update_batch.push_back(c.ub[i]);oldid.push_back(c.uo[i]<0?-1:id(c.uo[i],c.bu[i]));newid.push_back(id(c.un[i],c.bu[i]));}
  group_batch.assign(c.gb.begin(),c.gb.end());removed.resize(c.ng);extras.resize(c.ng);prepared.resize(c.ng,false);
 }
 void prepare(int g){if(prepared[g])return;auto& c=*source;
   std::map<int,int> remove;
   remove[id(c.sd[g],c.sb[g])]++;
   for(int k=c.eoff[g];k<c.eoff[g+1];++k)remove[id(c.ed[k],c.eb[k])]++;
   using Key=std::tuple<int,double,int,int,int>;
   std::map<Key,int> ordinary;
   auto motif=[&](int sid,int bid){uint64_t code=c.codes[sid];int d=code>>48,tr=(code>>16)&65535,cy=code&65535;
    if(cy==65535)Rcpp::stop("Exponential directed cache received an undirected score");
    if(c.packed[2*sid+1]!=(uint64_t)d)Rcpp::stop("Exponential cache requires directed target in-degree");
    int consumer=c.consumers?(code>>32)&65535:0;
    double b=c.birth[bid];remove[id(d,b)]++;ordinary[Key(d,b,tr,cy,consumer)]++;
   };
   if(c.birth_ids){for(int k=c.grouped_off[g];k<c.grouped_off[g+1];++k)for(int j=c.birth_off[k];j<c.birth_off[k+1];++j)motif(c.grouped_sid[k],c.birth_ids[j]);}
   else for(int k=c.off[g];k<c.off[g+1];++k)motif(c.raw[2*k],c.raw[2*k+1]);
   for(int k=c.hoff[g];k<c.hoff[g+1];++k){int d=c.hf[5*k],tr=c.hf[5*k+3],cy=c.hf[5*k+4];double b=c.hf[5*k+1];
    int consumer=c.consumers?c.hf[5*k+2]:0;
    if(tr||cy||consumer){auto it=ordinary.find(Key(d,b,tr,cy,consumer));if(it==ordinary.end()||it->second<=0)Rcpp::stop("History/motif candidate mismatch");it->second--;}
    else remove[id(d,b)]++;
    extras[g].push_back({id(d,b),tr,cy,(int)c.hflags[2*k],(int)c.hflags[2*k+1],1,consumer});
   }
   for(auto& kv:ordinary)if(kv.second){int d,tr,cy,consumer;double b;std::tie(d,b,tr,cy,consumer)=kv.first;extras[g].push_back({id(d,b),tr,cy,0,0,kv.second,consumer});}
   for(auto kv:remove)removed[g].push_back(kv);
  prepared[g]=true;
 }
};
template<int N> struct CranLogMass {
 double logmass=-INFINITY;std::array<double,N> mean{};
 void add(double z,const std::array<double,N>& x){if(z==-INFINITY)return;if(logmass==-INFINITY){logmass=z;mean=x;return;}
  double high=std::max(z,logmass),a=std::exp(logmass-high),b=std::exp(z-high),s=a+b;
  for(int j=0;j<N;++j)mean[j]=(a*mean[j]+b*x[j])/s;logmass=high+std::log(s);
 }
};
struct CranCohortTree {
 CranExpCache& c;const double* p;int size=1;std::vector<int> count;std::vector<CranLogMass<4>> tree;
 CranCohortTree(CranExpCache& cache,const double* par):c(cache),p(par),count(c.cohorts.size(),0){while(size<(int)count.size())size*=2;tree.resize(2*size);for(int i:c.initial)count[i]++;for(int i=0;i<(int)count.size();++i)leaf(i);for(int i=size-1;i;--i)merge(i);}
 std::array<double,4> feature(int i){double d=c.cohorts[i].degree;return {d/100,c.source->consumers?0.:d*(d-1)/20000,std::log1p(d),c.cohorts[i].birth};}
 void leaf(int i){auto& v=tree[size+i];v=CranLogMass<4>();if(count[i]){v.mean=feature(i);v.logmass=std::log((double)count[i])+p[0]*v.mean[0]+p[1]*v.mean[1]+p[4]*v.mean[2]+p[7]*v.mean[3];}}
 void merge(int i){tree[i]=tree[2*i];tree[i].add(tree[2*i+1].logmass,tree[2*i+1].mean);}
 void change(int i,int delta){count[i]+=delta;if(count[i]<0)Rcpp::stop("Negative candidate cohort count");leaf(i);for(int k=(size+i)/2;k;k/=2)merge(k);}
};
