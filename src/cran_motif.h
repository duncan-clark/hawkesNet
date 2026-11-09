#include <cmath>
#include <cstdint>
#include <vector>
#include <thread>
#include <unordered_map>
#include <set>
#include <array>
static double sigmoid(double x){return 1/(1+std::exp(-x));}
static std::array<double,9> weight(double di,double dout,double tr,double cy,double pop,const double* p,bool consumers=false){
 double x[8]={1,di/100.,dout/100.,di*(di-1)/20000.,dout*(dout-1)/20000.,tr/10.,cy/10.,std::log1p(pop)};
 if(consumers){x[2]=x[4]=0.;x[3]=dout/10.;}
 if(cy==65535){for(int j=1;j<7;++j)x[j]=0;}
 double z=0;for(int j=0;j<8;++j)z+=p[j]*x[j];double f=p[9]?std::exp(z):sigmoid(z);std::array<double,9> v;v[0]=f;
 for(int j=0;j<8;++j)v[j+1]=f*(p[9]?1:1-f)*x[j];return v;
}
extern "C" void motif_sum(int64_t ng,int64_t ns,int64_t nb,const int64_t* offsets,const uint32_t* rows,const uint64_t* codes,const double* birth,const double* p,double* output,bool consumers=false){
 std::vector<std::array<double,9>> table(ns);std::vector<double> age(nb);
 for(int64_t b=0;b<nb;++b)age[b]=std::exp(p[8]*birth[b]);
 for(int64_t s=0;s<ns;++s){auto c=codes[2*s];double pop=codes[2*s+1];double di=c>>48,dout=(c>>32)&65535,tr=(c>>16)&65535,cy=c&65535;
 auto v=weight(di,dout,tr,cy,pop,p,consumers),v0=weight(di,consumers?0:dout,0,0,pop,p,consumers);for(int j=0;j<9;++j)table[s][j]=v[j]-v0[j];}
 auto work=[&](int64_t lo,int64_t hi){for(int64_t g=lo;g<hi;++g){double* out=output+g*10;
  for(int64_t i=offsets[g];i<offsets[g+1];++i){auto sid=rows[2*i],bid=rows[2*i+1];double w=age[bid];auto& v=table[sid];
   for(int j=0;j<9;++j)out[j]+=w*v[j];out[9]+=w*v[0]*birth[bid];}}};
 std::vector<std::thread> threads;for(int k=0;k<4;++k)threads.emplace_back(work,ng*k/4,ng*(k+1)/4);for(auto& t:threads)t.join();
}
extern "C" void motif_prefix(int64_t ng,int64_t ni,int64_t nu,const int64_t* id,const double* ib,
 const int64_t* ub,const int64_t* uk,const int64_t* uo,const int64_t* un,const double* birth,
 const int64_t* gb,const int64_t* gd,const double* p,double* output,bool consumers=false){
 std::unordered_map<int,std::array<double,3>> hist[2];std::set<int> active[2];
 auto change=[&](int k,int d,double b,int sign){if(d<0)return;auto& h=hist[k][d];double w=std::exp(p[8]*b);h[0]+=sign*w;h[1]+=sign*w*b;h[2]+=sign;
 if(h[2]>0)active[k].insert(d);else {active[k].erase(d);h[0]=h[1]=0;}};
 for(int64_t i=0;i<ni;++i)change(i%1,id[i],ib[i],1);
 std::unordered_map<uint64_t,std::array<double,9>> cache;int64_t u=0;
 for(int64_t g=0;g<ng;++g){while(u<nu&&ub[u]<=gb[g]){for(int k=0;k<1;++k)if(uk[u]<0||uk[u]==k){change(k,uo[u],birth[u],-1);change(k,un[u],birth[u],1);}++u;}
 double* out=output+g*10;
 for(int k=0;k<1;++k){int dout=consumers?0:gd[g+k];for(int d:active[k]){uint64_t key=((uint64_t)dout<<32)|d;
 auto it=cache.find(key);if(it==cache.end())it=cache.emplace(key,weight(d>>16,dout,0,0,d&65535,p,consumers)).first;
 auto& v=it->second;auto& h=hist[k][d];for(int j=0;j<9;++j)out[j]+=h[0]*v[j];out[9]+=h[1]*v[0];}}
 }
}
// At eta=0, aggregate repeated statistic IDs exactly; retain sum(birth) for its derivative.
extern "C" void motif_compact(int64_t ng,int64_t ns,const int64_t* offsets,const uint32_t* sid,const uint32_t* count,const double* birthsum,const uint64_t* codes,const double* p,double* output){
 std::vector<std::array<double,9>> table(ns);
 for(int64_t s=0;s<ns;++s){auto c=codes[2*s];double pop=codes[2*s+1];double di=c>>48,dout=(c>>32)&65535,tr=(c>>16)&65535,cy=c&65535;
 auto v=weight(di,dout,tr,cy,pop,p),v0=weight(di,dout,0,0,pop,p);for(int j=0;j<9;++j)table[s][j]=v[j]-v0[j];}
 for(int64_t g=0;g<ng;++g){double* out=output+g*10;for(int64_t i=offsets[g];i<offsets[g+1];++i){auto& v=table[sid[i]];
 for(int j=0;j<9;++j)out[j]+=count[i]*v[j];out[9]+=birthsum[i]*v[0];}}
}

// Exact grouping by candidate feature vector, retaining every individual birth time.
extern "C" void motif_grouped(int64_t ng,int64_t ns,int64_t nb,const int64_t* offsets,const uint32_t* sid,const int64_t* bo,const uint32_t* bi,const uint64_t* codes,const double* birth,const double* p,double* output){
 std::vector<std::array<double,9>> table(ns);std::vector<double> age(nb),agebirth(nb);
 for(int64_t b=0;b<nb;++b){age[b]=std::exp(p[8]*birth[b]);agebirth[b]=age[b]*birth[b];}
 for(int64_t s=0;s<ns;++s){auto c=codes[2*s];double pop=codes[2*s+1],di=c>>48,dout=(c>>32)&65535,tr=(c>>16)&65535,cy=c&65535;
 auto v=weight(di,dout,tr,cy,pop,p),v0=weight(di,dout,0,0,pop,p);for(int j=0;j<9;++j)table[s][j]=v[j]-v0[j];}
 auto work=[&](int64_t lo,int64_t hi){for(int64_t g=lo;g<hi;++g){double* out=output+g*10;
  for(int64_t i=offsets[g];i<offsets[g+1];++i){double w=0,wb=0;for(int64_t k=bo[i];k<bo[i+1];++k){w+=age[bi[k]];wb+=agebirth[bi[k]];}
   auto& v=table[sid[i]];for(int j=0;j<9;++j)out[j]+=w*v[j];out[9]+=wb*v[0];}}};
 std::vector<std::thread> threads;for(int k=0;k<4;++k)threads.emplace_back(work,ng*k/4,ng*(k+1)/4);for(auto& t:threads)t.join();
}
