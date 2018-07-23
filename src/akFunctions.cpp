#include "akFunctions.h"


//******************************************************************************
double CLkk(int L, int ka, int kb)
{
  int la = ATI_l_k(ka);
  int lb = ATI_l_k(kb);
  int two_ja = ATI_twoj_k(ka);
  int two_jb = ATI_twoj_k(kb);
  double ja = 0.5*two_ja;
  double jb = 0.5*two_jb;

  if((la+lb+L)%2!=0) return 0; //Parity rule
  if((la+lb<L)||(abs(la-lb)>L)) return 0; //triangle rule (l)

  double tjB = WIG_3j(jb,ja,L,-0.5,0.5,0);
  return (2*ja+1)*(2*jb+1)*(2*L+1)*pow(tjB,2);
}



//******************************************************************************
double CLkk_OLD(int L, int ka, int kb)
/*
Calculates the angular coeficient (averaged over all m)
B. M. Roberts, V. A. Dzuba, V. V. Flambaum, M. Pospelov, and Y. V. Stadnik,
Phys. Rev. D 93, 115037 (2016). [arXiv:1604.04559]
*/
{
  int two_ja = ATI_twoj_k(ka);
  int two_jb = ATI_twoj_k(kb);
  double ja = 0.5*two_ja;
  double jb = 0.5*two_jb;
  int la = ATI_l_k(ka);
  int lb = ATI_l_k(kb);

  double tjB = WIG_3j(jb,L,ja,-0.5,0,0.5);
  if(fabs(tjB)==0) return 0;
  double B = 1./pow(tjB,2);

  //(-1)^(ja etc) -> calc sign
  int s1 = -1;
  if((two_ja+two_jb-2*(la+lb))%4==0) s1=1;

  double tj1 = WIG_3j(lb,la,L,0,0,0);
  double A = (1./4)*s1*(2*L+1)*pow(tj1,2);
  double X = s1*(two_ja+1)*(two_jb+1)*pow(tj1,2);
  double tj2 = WIG_3j(lb,la,L,-1,1,0);
  double Y = 8*sqrt(la*(la+1)*lb*(lb+1))*tj1*tj2;
  double Z = -4*(ka+1)*(kb+1)*pow(tj2,2);

  return (A*B)*(X+Y+Z);
}


//******************************************************************************
/*
Helper functions for the binary read in/out
*/
template<typename T>
int binary_rw(std::fstream& stream, T& value, bool write){
  if(write) stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
  else stream.read(reinterpret_cast<char*>(&value), sizeof(T));
  return 0;
}
template<typename T>
int binary_str_rw(std::fstream& stream, T& value, bool write){
  if(write){
    size_t temp_len = value.length();
    stream.write(reinterpret_cast<const char*>(&temp_len), sizeof(size_t));
    stream.write(value.c_str(), value.length());
  }else{
    size_t temp_len;
    stream.read(reinterpret_cast<char*>(&temp_len), sizeof(size_t));
    char* tvalue = new char[temp_len+1];
    stream.read(tvalue, temp_len);
    tvalue[temp_len] = '\0'; //null 'end of string' character
    value = tvalue;
    delete [] tvalue;
    return 0;
  }
  return 0;
}


//******************************************************************************
int akReadWrite(std::string fname, bool write,
  std::vector< std::vector< std::vector<float> > > &AK,
  std::vector<float> &dElst,
  std::vector<std::string> &nklst,
  std::vector<float> &qlst)
{
  std::fstream iof;
  if(write) iof.open(fname,std::ios_base::out|std::ios_base::binary);
  else      iof.open(fname,std::ios_base::in |std::ios_base::binary);

  if(write){
    int tmp1 = AK.size();
    int tmp2 = AK[0].size();
    int tmp3 = AK[0][0].size();
    binary_rw(iof,tmp1,write);
    binary_rw(iof,tmp2,write);
    binary_rw(iof,tmp3,write);
  }else{
    int nq,ns,nde;
    binary_rw(iof,nde,write);
    binary_rw(iof,ns,write);
    binary_rw(iof,nq,write);
    AK.resize(nde,std::vector< std::vector<float> >(ns,std::vector<float>(nq)));
    dElst.resize(nde);
    nklst.resize(ns);
    qlst.resize(nq);
  }
  for(size_t ie=0; ie<AK.size(); ie++){
    binary_rw(iof,dElst[ie],write);
    for(size_t in=0; in<AK[0].size(); in++){
      if(ie==0) binary_str_rw(iof,nklst[in],write);
      for(size_t iq=0; iq<AK[0][0].size(); iq++){
        if(ie==0 && in==0) binary_rw(iof,qlst[iq],write);
        binary_rw(iof,AK[ie][in][iq],write);
      }
    }
  }

  return 0;
}




//******************************************************************************
int calculateKpw_nk(ElectronOrbitals &wf, int nk, double dE,
  std::vector< std::vector<float> > &jl_q_r,
  std::vector< std::vector<float> > &K_nk;
)
/*
For plane-wave final state.
Only has f-part....Can restore g-part, but need to be sure of plane-wave!
Chi(q) - Int[ f_nk*j_l(qr)*r , {r,0,inf}]
Should be called once per initial state
*/
{
  if (nk>=(int)wf.p.size()) return 1; //should never occur XXX

  int kappa = wf.klist[nk];
  int l = ATI_l_k(kappa);
  int twoj = ATI_twoj_k(kappa);

  int qsteps = (int)jl_q_r.size();
  std::vector<double> tmpK_q(qsteps);

  double eps = dE - wf.en[nk];
  int maxir = wf.pinflist[nk]; //don't bother going further

  for(int iq=0; iq<qsteps; iq++){
    if(eps<=0) break;
    double chi_q=0;
    for(int ir=0; ir<maxir; ir++){
      chi_q += wf.p[nk][ir]*jl_qr[iq][ir]*wf.r[ir];
    }
    chi_q *= wf.h;
    tmpK_q[iq] = (2./M_PI)*(twoj+1)*pow(chi_q,2)*sqrt(2.*eps);
  }

  K_nk.push_back(tmpK_q);
  return 0;
}