#include <cmath>
#include <cstdint>
#include <vector>
#include <thread>
#include <unordered_map>
#include <array>
#include <algorithm>

// Exact sums over all target candidates. Arrays are read-only and the output
// rows are disjoint between workers. No risk-set or history subsampling.
struct Cell { int d,n; double b,x1,x3,pop,z,ez; };
static inline void link(double z,double& h,double& s) {
  if(z>=0){double x=std::exp(-z);h=z+std::log1p(x);s=1/(1+x);}
  else {double x=std::exp(z);h=std::log1p(x);s=x/(1+x);}
}
static inline void fastlink(double z,double ez,double q,double offset,double& h,double& s){
  if(z>-600 && z<600 && offset>-100 && offset<100){double x=ez*q;h=std::log1p(x);s=x/(1+x);}
  else link(z+offset,h,s);
}
extern "C" void normalizer(int64_t ng,int64_t ni,int64_t nu,int64_t ns,
    int64_t upto,int64_t nth,void** aa,const double* p,double* output,const double* sender_offsets,bool consumers=false,bool motif_free=false){
  auto id=(int64_t*)aa[0];auto ib=(double*)aa[1];auto ub=(int64_t*)aa[2];
  auto uo=(int64_t*)aa[3];auto un=(int64_t*)aa[4];auto bu=(double*)aa[5];
  auto gb=(int64_t*)aa[6];auto go=(int64_t*)aa[7];auto gt=(double*)aa[8];
  auto off=(int64_t*)aa[9];auto raw=(uint32_t*)aa[10];auto codes=(uint64_t*)aa[11];auto birth=(double*)aa[12];
  auto ed=(int64_t*)aa[13];auto eb=(double*)aa[14];auto et=(int64_t*)aa[15];auto eo=(int64_t*)aa[16];
  auto hf=(double*)aa[17];auto flags=(int64_t*)aa[18];auto sd=(int64_t*)aa[19];auto sb=(double*)aa[20];
  auto eoff=(int64_t*)aa[21];auto hoff=(int64_t*)aa[22];
  double eta=p[17];
  // The integer IDs identify exact birth times, not rounded age categories.
  std::unordered_map<double,int> bid;std::vector<double> bvalues;
  auto bidx=[&](double b){auto it=bid.find(b);if(it!=bid.end())return it->second;
    int k=bvalues.size();bid.emplace(b,k);bvalues.push_back(b);return k;};
  std::vector<int> ibi(ni),ubi(nu);
  for(int64_t i=0;i<ni;++i)ibi[i]=bidx(ib[i]);
  for(int64_t i=0;i<nu;++i)ubi[i]=bidx(bu[i]);
  auto work=[&](int64_t lo,int64_t hi){
    std::unordered_map<uint64_t,int> positions;std::vector<Cell> cells;cells.reserve(45000);
    auto change=[&](int d,int bi,int delta){if(d<0)return;
      uint64_t key=((uint64_t)bi<<32)|(uint32_t)d;auto it=positions.find(key);
      if(it==positions.end()){
        double x1=d/100.,x3=consumers?0.:d*(d-1.)/20000.,pop=std::log1p((double)d),b=bvalues[bi];
        double z=p[1]*x1+p[3]*x3+p[7]*pop+eta*b;
        positions.emplace(key,cells.size());cells.push_back({d,delta,b,x1,x3,pop,z,std::exp(std::min(700.,z))});
      }else{
        int j=it->second;cells[j].n+=delta;
        if(cells[j].n==0){positions.erase(it);if(j!=(int)cells.size()-1){cells[j]=cells.back();
          uint64_t other=((uint64_t)bid.at(cells[j].b)<<32)|(uint32_t)cells[j].d;positions[other]=j;}
          cells.pop_back();}
      }
    };
    for(int64_t i=0;i<ni;++i)change(id[i],ibi[i],1);
    int64_t update=0;
    auto add=[&](double* out,int d,int o,int tr,int cy,double b,int mig,int copy,double sign,double t,double sender_offset,int co=0){
      double x[10]={1.,d/100.,o/100.,consumers?co/10.:d*(d-1.)/20000.,o*(o-1.)/20000.,tr/10.,cy/10.,std::log1p((double)d),(double)mig,(double)copy};
      double z=sender_offset-eta*(t-b);for(int j=0;j<10;++j)z+=p[j]*x[j];double h,s;link(z,h,s);
      out[0]+=sign*h;for(int j=0;j<10;++j)out[1+j]+=sign*s*x[j];out[18]-=sign*s*(t-b);
    };
    for(int64_t g=lo;g<hi;++g){
      while(update<nu&&ub[update]<=gb[g]){change(uo[update],ubi[update],-1);change(un[update],ubi[update],1);++update;}
      double* out=output+20*g;double t=gt[g],o=go[g],ox1=o/100.,ox3=o*(o-1.)/20000.;
      double addition_offset=sender_offsets?sender_offsets[2*g]:0.;
      double deletion_offset=sender_offsets?sender_offsets[2*g+1]:0.;
      double offset=p[0]+p[2]*ox1+p[4]*ox3-eta*t+addition_offset,q=std::exp(std::min(700.,offset));
      double total=0,btotal=0;
      for(const auto& c:cells){double h,s;fastlink(c.z,c.ez,q,offset,h,s);h*=c.n;s*=c.n;
        out[0]+=h;total+=s;btotal+=s*c.b;out[2]+=s*c.x1;out[4]+=s*c.x3;out[8]+=s*c.pop;}
      out[1]=total;out[3]=total*ox1;out[5]=total*ox3;out[18]=btotal-t*total;
      // Prefix includes self-loops and present edges; remove them exactly.
      add(out,sd[g],o,0,0,sb[g],0,0,-1,t,addition_offset);
      for(int64_t i=eoff[g];i<eoff[g+1];++i){
        add(out,ed[i],o,0,0,eb[i],0,0,-1,t,addition_offset);
        double d=ed[i]-1.,od=eo[i]-1.;
        double x[7]={1.,-d/100.,-od/100.,-d*(d-1.)/20000.,-od*(od-1.)/20000.,-et[i]/10.,-std::log((double)std::max((int64_t)1,ed[i]))};
        double z=deletion_offset-p[18]*(t-eb[i]);for(int j=0;j<7;++j)z+=p[10+j]*x[j];double h,s;link(z,h,s);
        out[0]+=h;for(int j=0;j<7;++j)out[11+j]+=s*x[j];out[19]-=s*(t-eb[i]);
      }
      // Correct all candidates with nonzero motif changes, then history flags.
      if(!motif_free)for(int64_t i=off[g];i<off[g+1];++i){auto code=codes[raw[2*i]];
        int d=code>>48,oo=(code>>32)&65535,tr=(code>>16)&65535,cy=code&65535;double b=birth[raw[2*i+1]];
        int co=consumers?oo:0;if(consumers)oo=go[g];
        add(out,d,oo,tr,cy,b,0,0,1,t,addition_offset,co);add(out,d,oo,0,0,b,0,0,-1,t,addition_offset);
      }
      for(int64_t i=hoff[g];i<hoff[g+1];++i){auto f=hf+5*i;
        int oo=consumers?go[g]:f[2],co=consumers?f[2]:0;
        add(out,f[0],oo,f[3],f[4],f[1],flags[2*i],flags[2*i+1],1,t,addition_offset,co);
        add(out,f[0],oo,f[3],f[4],f[1],0,0,-1,t,addition_offset,co);
      }
    }
  };
  nth=std::min(nth,upto);std::vector<std::thread> threads;
  for(int k=0;k<nth;++k)threads.emplace_back(work,upto*k/nth,upto*(k+1)/nth);
  for(auto& thread:threads)thread.join();
}
