#include <Rcpp.h>
#include <fstream>
#include <regex>
#include <cstring>
#include <numeric>
#include <memory>
#include <thread>
#include <unordered_map>
#include <set>
#include <array>
#include <map>
#include <tuple>
#include <type_traits>
#include <limits>
#include "cran_release_validation.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
namespace proposal {
#include "cran_motif.h"
}
namespace transition {
#include "cran_tergm.h"
}
namespace release {
#include "cran_release.h"
}
using namespace Rcpp;
// The import reads the legacy numerical cache only. No Python runtime is used.
// [[Rcpp::export]]
SEXP cran_read_npy(std::string path) {
 std::ifstream f(path, std::ios::binary);
 unsigned char magic[8]; f.read(reinterpret_cast<char*>(magic), 8);
 if (!f || std::memcmp(magic, "\x93NUMPY", 6)) stop("Invalid NPY file: %s", path);
 if (magic[6] < 1 || magic[6] > 3 || magic[7] != 0) stop("Unsupported NPY version");
 const uint16_t endian = 1;
 if (*reinterpret_cast<const unsigned char*>(&endian) != 1) stop("NPY cache import requires a little-endian host");
 unsigned char length_bytes[4] = {0,0,0,0};
 f.read(reinterpret_cast<char*>(length_bytes), magic[6] == 1 ? 2 : 4);
 uint32_t hlen = 0;
 for (int j = 0; j < 4; ++j) hlen |= static_cast<uint32_t>(length_bytes[j]) << (8*j);
 auto header_start = f.tellg(); f.seekg(0, std::ios::end); auto file_end = f.tellg();
 if (!f || header_start < 0 || file_end < header_start ||
     static_cast<uint64_t>(file_end - header_start) < hlen) stop("Truncated NPY header");
 f.seekg(header_start); std::string header(hlen, ' ');
 if (hlen) f.read(&header[0], hlen);
 if (!f) stop("Truncated NPY header");
 std::smatch match;
 if (!std::regex_search(header, match, std::regex("fortran_order['\"]? *: *(True|False)")))
   stop("Missing or invalid NPY array order");
 if (match[1] == "True") stop("Fortran-order cache is unsupported");
 if (!std::regex_search(header, match, std::regex("descr['\"]? *: *['\"]([^'\"]+)"))) stop("Missing NPY dtype");
 std::string type = match[1];
 const std::set<std::string> supported = {"f8","f4","i8","i4","i2","i1","u8","u4","u2","u1","b1"};
 if (type.size() < 3 || (type[0] != '<' && type[0] != '|' && type[0] != '=') ||
     !supported.count(type.substr(1))) stop("Unsupported NPY dtype/endian: %s", type);
 int size = std::stoi(type.substr(2));
 if (!std::regex_search(header, match, std::regex("shape['\"]? *: *\\(([^)]*)\\)"))) stop("Missing NPY shape");
 std::string dims = match[1];
 if (!std::regex_match(dims, std::regex("[ \\t]*(?:[0-9]+[ \\t]*(?:,[ \\t]*[0-9]+[ \\t]*)*,?[ \\t]*)?"))) stop("Invalid NPY shape");
 std::vector<int> shape; std::regex digits("[0-9]+");
 size_t n = 1;
 for (auto it = std::sregex_iterator(dims.begin(), dims.end(), digits); it != std::sregex_iterator(); ++it) {
   double value = std::stod(it->str());
   if (!R_finite(value) || value > INT_MAX) stop("NPY dimension exceeds supported range");
   int d = static_cast<int>(value); shape.push_back(d);
   if (d && n > static_cast<size_t>(R_XLEN_T_MAX) / d) stop("NPY shape is too large");
   n *= d;
 }
 auto payload = f.tellg();
 if (payload < 0 || n > static_cast<uint64_t>(file_end - payload) / size) stop("Truncated NPY cache");
 std::vector<unsigned char> bytes(n * size);
 if (!bytes.empty()) f.read(reinterpret_cast<char*>(bytes.data()), bytes.size());
 if (!f) stop("Truncated NPY cache");
 if (type.substr(1) == "u8") { RawVector ans(bytes.begin(), bytes.end()); ans.attr("npy_shape")=shape; return ans; }
 NumericVector ans(n);
 for (size_t i=0; i<n; ++i) {
   const void* b=bytes.data()+size*i; double value=0;
   if (type[1]=='f') { if(size==8) std::memcpy(&value,b,8); else {float q;std::memcpy(&q,b,4);value=q;} }
   else if(type[1]=='i') { if(size==8){int64_t q;std::memcpy(&q,b,8);value=q;} else if(size==4){int32_t q;std::memcpy(&q,b,4);value=q;} else if(size==2){int16_t q;std::memcpy(&q,b,2);value=q;} else value=*static_cast<const int8_t*>(b); }
   else { if(size==4){uint32_t q;std::memcpy(&q,b,4);value=q;} else if(size==2){uint16_t q;std::memcpy(&q,b,2);value=q;} else value=type[1]=='b'?(*static_cast<const uint8_t*>(b)!=0):*static_cast<const uint8_t*>(b); }
   ans[i]=value;
 }
 if(shape.size()==2) {NumericMatrix out(shape[0],shape[1]);for(int i=0;i<shape[0];++i)for(int j=0;j<shape[1];++j)out(i,j)=ans[(size_t)i*shape[1]+j];return out;}
 return ans;
}

// A mapping owns its resources even when a later constructor check throws.
std::shared_ptr<uint32_t> map_uint32_cache(const std::string& path, size_t& bytes) {
 int fd = open(path.c_str(), O_RDONLY);
 if (fd < 0) stop("Cannot open candidate cache: %s", path);
 struct stat info;
 if (fstat(fd, &info) != 0 || info.st_size < 0 || info.st_size % sizeof(uint32_t)) {
   close(fd); stop("Invalid candidate cache size: %s", path);
 }
 bytes = static_cast<size_t>(info.st_size);
 if (!bytes) {close(fd); return {};}
 void* address = mmap(nullptr, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
 close(fd);
 if (address == MAP_FAILED) stop("Cannot map candidate cache: %s", path);
 return std::shared_ptr<uint32_t>(static_cast<uint32_t*>(address),
   [bytes](uint32_t* p) {munmap(p, bytes);});
}

template<class T> void check_offsets(const std::vector<T>& offsets, size_t groups, size_t available, const char* name) {
 if (offsets.size() != groups+1 || offsets.front() != 0 || offsets.back() < 0 ||
     static_cast<size_t>(offsets.back()) > available || !std::is_sorted(offsets.begin(), offsets.end()))
   stop("Invalid %s offsets", name);
}

template<class T> std::vector<T> vec(SEXP s) {
 if (TYPEOF(s) != REALSXP && TYPEOF(s) != INTSXP) stop("Cache arrays must be numeric");
 NumericVector x = as<NumericVector>(s);
 auto checked = [](double value) -> T {
   if (!R_finite(value)) stop("Cache arrays must contain finite values");
   if constexpr (std::is_integral<T>::value) {
     double upper = std::ldexp(1.0, std::numeric_limits<T>::digits);
     double lower = std::is_signed<T>::value ? -upper : 0;
     if (value != std::floor(value) || value < lower || value >= upper)
       stop("Cache integer array contains a fractional or out-of-range value");
   }
   return static_cast<T>(value);
 };
 std::vector<T> v(x.size());
 if (Rf_isMatrix(s)) {
   NumericMatrix m(x);
   for (int i=0;i<m.nrow();++i) for (int j=0;j<m.ncol();++j)
     v[static_cast<size_t>(i)*m.ncol()+j] = checked(m(i,j));
 } else for (R_xlen_t i=0;i<x.size();++i) v[i] = checked(x[i]);
 return v;
}
std::vector<double> rowmajor(NumericMatrix x){std::vector<double> v(x.size());for(int i=0;i<x.nrow();++i)for(int j=0;j<x.ncol();++j)v[(size_t)i*x.ncol()+j]=x(i,j);return v;}
struct CranCache {
 std::vector<int64_t> id,ub,uk,uo,un,gb,gd,off,ed,et,eo,sd,eoff,hoff,hflags;
 std::vector<double> ib,bu,gt,birth,eb,sb,hf;
 std::vector<uint64_t> codes,packed;
 std::vector<int64_t> grouped_off,birth_off;
 std::vector<uint32_t> grouped_sid;
 size_t birth_bytes=0, bytes=0; uint32_t *birth_ids=nullptr, *raw=nullptr;
 std::shared_ptr<uint32_t> raw_mapping, birth_mapping; int ng,ni,nu;
 bool consumers=false;
 CranCache(List a,List ix,RawVector code, std::string rawpath,List hist,NumericVector time,bool undir,Nullable<List> groups,std::string birthpath){
  for(double t:time)if(!R_finite(t))stop("Candidate times must be finite");
  id=vec<int64_t>(a["initial_degree"]);ib=vec<double>(a["initial_birth"]);ni=id.size();gb=vec<int64_t>(a["group_batch"]);ng=gb.size();gd=vec<int64_t>(a["group_outdegree"]);
  ub=vec<int64_t>(a["update_batch"]);uk=vec<int64_t>(a["update_kind"]);uo=vec<int64_t>(a["update_old_degree"]);un=vec<int64_t>(a["update_new_degree"]);bu=vec<double>(a["update_birth"]);nu=ub.size();
  if(ib.size()!=id.size() || gd.size()!=gb.size() || uk.size()!=ub.size() || uo.size()!=ub.size() || un.size()!=ub.size() || bu.size()!=ub.size())stop("Inconsistent candidate-cache array lengths");
  if(!std::is_sorted(gb.begin(),gb.end()))stop("Candidate groups must be ordered by batch");
  for(auto b:gb)if(b<0||b>=time.size())stop("Invalid candidate batch index");
  // Changes after the final batch use the terminal sentinel n; they are
  // retained in archived caches but never enter a pre-event normalizer.
  for(auto b:ub)if(b<0||b>time.size())stop("Invalid degree-update batch index");
  std::vector<size_t> order(nu);std::iota(order.begin(),order.end(),0);std::stable_sort(order.begin(),order.end(),[&](size_t i,size_t j){return ub[i]<ub[j];});
  auto sortby=[&](auto& x){auto y=x;for(int i=0;i<nu;++i)x[i]=y[order[i]];};sortby(ub);sortby(uk);sortby(uo);sortby(un);sortby(bu);
  off=vec<int64_t>(ix["offsets"]);birth=vec<double>(ix["birth"]);
  consumers=ix.containsElementNamed("shared_consumer")&&as<bool>(ix["shared_consumer"]);
  if(code.size()%8)stop("Invalid packed candidate codes");
  size_t nc=code.size()/8;codes.resize(nc);if(nc)std::memcpy(codes.data(),code.begin(),code.size());
  if(!undir){packed.resize(2*nc);for(size_t i=0;i<nc;++i){packed[2*i]=codes[i];packed[2*i+1]=codes[i]>>48;}}
  else packed=codes;
  raw_mapping=map_uint32_cache(rawpath,bytes);raw=raw_mapping.get();
  check_offsets(off,ng,bytes/(2*sizeof(uint32_t)),"candidate");
  for(int64_t i=0;i<off.back();++i)if(raw[2*i]>=packed.size()/2||raw[2*i+1]>=birth.size())stop("Invalid candidate statistic or birth index");
  ed=vec<int64_t>(a["existing_degree"]);et=vec<int64_t>(a["existing_triangle"]);eo=vec<int64_t>(a["existing_outdegree"]);eb=vec<double>(a["existing_birth"]);sd=vec<int64_t>(a["group_source_degree"]);sb=vec<double>(a["group_source_birth"]);
  NumericVector eg=a["existing_group"],hg=hist["group"];eoff.assign(ng+1,0);hoff.assign(ng+1,0);
  if(ed.size()!=eb.size()||ed.size()!=et.size()||ed.size()!=eo.size()||ed.size()!=static_cast<size_t>(eg.size())||sd.size()!=static_cast<size_t>(ng)||sb.size()!=static_cast<size_t>(ng))stop("Inconsistent existing-candidate arrays");
  if(!std::is_sorted(eg.begin(),eg.end()))stop("Existing candidates must be ordered by group");
  for(double g:eg){if(!R_finite(g)||g!=std::floor(g)||g<0||g>=ng)stop("Invalid existing-candidate group");eoff[g+1]++;}
  for(double g:hg){if(!R_finite(g)||g!=std::floor(g)||g<0||g>=ng)stop("Invalid history group");hoff[g+1]++;}for(int g=0;g<ng;++g){eoff[g+1]+=eoff[g];hoff[g+1]+=hoff[g];}
  NumericMatrix h=hist["features"],flags=hist["flags"];if((h.ncol()!=5&&!(undir&&h.ncol()==7))||flags.ncol()!=2||h.nrow()!=flags.nrow()||h.nrow()!=hg.size())stop("Invalid history feature dimensions");std::vector<int> horder(h.nrow());std::iota(horder.begin(),horder.end(),0);std::stable_sort(horder.begin(),horder.end(),[&](int i,int j){return hg[i]<hg[j];});
  auto history_values=vec<double>(h);auto flag_values=vec<int64_t>(flags);
  hf.resize(h.size());hflags.resize(flags.size());for(int i=0;i<h.nrow();++i){for(int j=0;j<h.ncol();++j)hf[(size_t)i*h.ncol()+j]=history_values[(size_t)horder[i]*h.ncol()+j];for(int j=0;j<2;++j)hflags[2*i+j]=flag_values[2*horder[i]+j];}
  gt.resize(ng);for(int g=0;g<ng;++g)gt[g]=time[gb[g]];
  if(groups.isNotNull()){
   List gr(groups);grouped_off=vec<int64_t>(gr["offsets"]);birth_off=vec<int64_t>(gr["birth_offsets"]);grouped_sid=vec<uint32_t>(gr["sid"]);
   birth_mapping=map_uint32_cache(birthpath,birth_bytes);birth_ids=birth_mapping.get();
   check_offsets(grouped_off,ng,grouped_sid.size(),"grouped candidate");
   check_offsets(birth_off,grouped_sid.size(),birth_bytes/sizeof(uint32_t),"grouped birth");
   for(auto sid:grouped_sid)if(sid>=packed.size()/2)stop("Invalid grouped statistic index");
   for(int64_t j=0;j<birth_off.back();++j)if(birth_ids[j]>=birth.size())stop("Invalid grouped birth index");
  }
 }
};
#include "cran_exponential.h"
// [[Rcpp::export]]
SEXP cran_exponential_cache(SEXP cache){XPtr<CranCache> c(cache);if(!c)stop("Invalid candidate cache pointer");XPtr<CranExpCache> ans(new CranExpCache(*c),true);ans.attr("parent")=cache;return ans;}
// [[Rcpp::export]]
NumericMatrix cran_exponential_sum(SEXP cache,NumericVector p,Nullable<IntegerVector> groups=R_NilValue){
 XPtr<CranExpCache> c(cache);if(!c)stop("Invalid exponential cache pointer");if(p.size()!=8)stop("Eight exponential addition coefficients required");for(double value:p)if(!R_finite(value))stop("Non-finite exponential coefficient");
 CranCohortTree tree(*c,p.begin());NumericMatrix ans(c->group_batch.size(),9);int u=0;std::vector<bool> take(c->group_batch.size(),groups.isNull());if(groups.isNotNull()){IntegerVector ix(groups);for(int g:ix){if(g<0||g>=(int)take.size())stop("Invalid candidate group");take[g]=true;}}
 for(int g=0;g<(int)c->group_batch.size();++g){
  if(!take[g])continue;c->prepare(g);if(c->cohorts.size()!=tree.count.size())stop("Incomplete target cohort history");
  while(u<(int)c->update_batch.size()&&c->update_batch[u]<=c->group_batch[g]){if(c->oldid[u]>=0)tree.change(c->oldid[u],-1);tree.change(c->newid[u],1);++u;}
  for(auto r:c->removed[g])tree.change(r.first,-r.second);
  CranLogMass<8> total;std::array<double,8> x{};auto base=tree.tree[1];x[0]=base.mean[0];x[1]=base.mean[1];x[4]=base.mean[2];x[7]=base.mean[3];total.add(base.logmass,x);
  for(auto e:c->extras[g]){auto v=tree.feature(e.cohort);std::array<double,8> z={v[0],c->source->consumers?e.consumer/10.:v[1],e.tr/10.,e.cy/10.,v[2],(double)e.mig,(double)e.copy,v[3]};double score=std::log((double)e.count);for(int j=0;j<8;++j)score+=p[j]*z[j];total.add(score,z);}
  ans(g,0)=total.logmass;for(int j=0;j<8;++j)ans(g,j+1)=total.mean[j];
  for(auto r:c->removed[g])tree.change(r.first,r.second);
  if(g%100==0)checkUserInterrupt();
 }
 return ans;
}
// [[Rcpp::export]]
SEXP cran_native_cache(List a,List index,RawVector codes,std::string rawpath,List history,NumericVector time,bool undirected,Nullable<List> groups=R_NilValue,std::string birthpath=""){return XPtr<CranCache>(new CranCache(a,index,codes,rawpath,history,time,undirected,groups,birthpath),true);}
// [[Rcpp::export]]
NumericMatrix cran_candidate_sum(SEXP cache,NumericVector p,bool undirected=false,bool tergm=false,int upto=0,int threads=4,bool exponential=false,Nullable<NumericMatrix> sender_offsets=R_NilValue,bool motif_free=false){
 XPtr<CranCache> c(cache);if(!c)stop("Invalid candidate cache pointer");
 if(upto<0||upto>c->ng||threads<1)stop("Invalid candidate prefix or thread count");
 if(!tergm&&p.size()!=9)stop("Nine proposal coefficients required");
 for(double value:p)if(!R_finite(value))stop("Non-finite candidate coefficient");
 int n=upto>0?upto:c->ng;int stride=tergm?20:10;std::vector<double> out(c->ng*stride,0),baseline(c->ng);std::vector<double> pp(p.begin(),p.end());if(!tergm){pp.resize(10);pp[9]=exponential;}
 if(tergm){
  if(undirected||p.size()!=19)stop("TERGM requires 19 directed edge coefficients");
  if(motif_free)for(int j:{1,2,3,4,5,6,11,12,13,14,15})
   if(p[j]!=0)stop("Motif-free TERGM normalizer requires fixed motif coefficients");
  std::vector<double> offsets;
  if(sender_offsets.isNotNull()){
   NumericMatrix shifts(sender_offsets);
   if(shifts.nrow()!=c->ng||shifts.ncol()!=2)stop("Invalid TERGM sender offsets");
   offsets=rowmajor(shifts);
   for(double value:offsets)if(!R_finite(value))stop("Non-finite TERGM sender offset");
  }
  void* arrays[]={c->id.data(),c->ib.data(),c->ub.data(),c->uo.data(),c->un.data(),c->bu.data(),c->gb.data(),c->gd.data(),c->gt.data(),c->off.data(),c->raw,c->codes.data(),c->birth.data(),c->ed.data(),c->eb.data(),c->et.data(),c->eo.data(),c->hf.data(),c->hflags.data(),c->sd.data(),c->sb.data(),c->eoff.data(),c->hoff.data()};transition::normalizer(c->ng,c->ni,c->nu,c->codes.size(),n,threads,arrays,p.begin(),out.data(),offsets.empty()?nullptr:offsets.data(),c->consumers,motif_free);
 }
 else {auto id=c->id,uo=c->uo,un=c->un;if(!undirected){for(auto& d:id)if(d>=0)d*=65537;for(auto& d:uo)if(d>=0)d*=65537;for(auto& d:un)if(d>=0)d*=65537;}
  proposal::motif_prefix(c->ng,c->ni,c->nu,id.data(),c->ib.data(),c->ub.data(),c->uk.data(),uo.data(),un.data(),c->bu.data(),c->gb.data(),c->gd.data(),pp.data(),out.data(),c->consumers);
  for(int i=0;i<c->ng;++i)baseline[i]=out[i*stride];
  if(c->birth_ids)proposal::motif_grouped(c->ng,c->packed.size()/2,c->birth.size(),c->grouped_off.data(),c->grouped_sid.data(),c->birth_off.data(),c->birth_ids,c->packed.data(),c->birth.data(),pp.data(),out.data());
  else proposal::motif_sum(c->ng,c->packed.size()/2,c->birth.size(),c->off.data(),c->raw,c->packed.data(),c->birth.data(),pp.data(),out.data(),c->consumers);}
 NumericMatrix ans(n,stride+(exponential?1:0));for(int i=0;i<n;++i){for(int j=0;j<stride;++j)ans(i,j)=out[(size_t)i*stride+j];if(exponential)ans(i,stride)=baseline[i];}return ans;
}
// [[Rcpp::export]]
NumericMatrix cran_group_sum(NumericVector group,NumericMatrix x,int n){
 if(n<0||group.size()!=x.nrow())stop("Invalid group-sum dimensions");
 for(double g:group)if(!R_finite(g)||g!=std::floor(g)||g<0||g>=n)stop("Invalid group-sum index");
 NumericMatrix ans(n,x.ncol());for(int i=0;i<x.nrow();++i)for(int j=0;j<x.ncol();++j)ans(group[i],j)+=x(i,j);return ans;}
// [[Rcpp::export]]
List cran_release_evaluate(List a,NumericVector p,int families){cran_input::release_data(a,p,families);NumericVector time=a["time"];int n=time.size();auto init=rowmajor(as<NumericMatrix>(a["initial"])),changes=rowmajor(as<NumericMatrix>(a["changes"]));auto ptr=vec<long long>(a["change_ptr"]),sel=vec<long long>(a["selected"]),sptr=vec<long long>(a["selected_ptr"]),nc=vec<long long>(a["new_count"]);std::vector<double> ll(n),g(n*15);release::evaluate(n,families,init.size()/6,init.data(),changes.data(),ptr.data(),sel.data(),sptr.data(),time.begin(),nc.data(),p.begin(),ll.data(),g.data());NumericMatrix grad(n,15);for(int i=0;i<n;++i)for(int j=0;j<15;++j)grad(i,j)=g[i*15+j];return List::create(_["loglik"]=ll,_["gradient"]=grad);}
