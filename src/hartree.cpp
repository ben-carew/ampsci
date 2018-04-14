#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "ContinuumOrbitals.h"

int main(void){

  clock_t ti,tf;
  ti = clock();

  double varalpha=1; //need same number as used for the fitting!

  //Input options
  std::string Z_str;
  int A;
  double r0,rmax;
  int ngp;

  double Tf,Tt,Tg;  //Teitz potential parameters
  double Gf,Gh,Gd;  //Green potential parameters

  int n_max,l_max;
  std::vector<std::string> str_core;

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("hartree.in");
    std::string jnk;
    // read in the input parameters:
    ifs >> Z_str >> A;            getline(ifs,jnk);
    while(true){
      std::string str;
      ifs >> str;
      if(str=="."||str=="|"||str=="!") break;
      str_core.push_back(str);
    }
    getline(ifs,jnk);
    ifs >> r0 >> rmax >> ngp;     getline(ifs,jnk);
    ifs >> Gf >> Gh >> Gd;        getline(ifs,jnk);
    ifs >> Tf >> Tt >> Tg;        getline(ifs,jnk);
    ifs >> n_max >> l_max;        getline(ifs,jnk);
    ifs.close();
  }

  int Z = ATI_get_z(Z_str);
  if(Z==0) return 2;
  if(A==0) A=ATI_a[Z]; //if none given, get default A

  //Normalise the Teitz/Green weights:
  if(Gf!=0 || Tf!=0){
    double TG_norm = Gf + Tf;
    Gf /= TG_norm;
    Tf /= TG_norm;
  }

  //If H,d etc are zero, use default values
  if(Gf!=0 && Gh==0) PRM_defaultGreen(Z,Gh,Gd);
  if(Tf!=0 && Tt==0) PRM_defaultTietz(Z,Tt,Tg);


  printf("\nRunning parametric potential for %s, Z=%i A=%i\n",
    Z_str.c_str(),Z,A);
  printf("*************************************************\n");
  if(Gf!=0) printf("%3.0f%% Green potential: H=%.4f  d=%.4f\n",Gf*100.,Gh,Gd);
  if(Tf!=0) printf("%3.0f%% Tietz potential: T=%.4f  g=%.4f\n",Tf*100.,Tt,Tg);

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  //if(A!=0) wf.sphericalNucleus();

  printf("Grid: pts=%i h=%7.5f Rmax=%5.1f\n",wf.ngp,wf.h,wf.r[wf.ngp-1]);

  //Determine which states are in the core:
  std::vector<int> core_list; //should be in the class!
  int core_ok = wf.determineCore(str_core,core_list);
  if(core_ok==2){
    std::cout<<"Problem with core: ";
    for(size_t i=0; i<str_core.size(); i++) std::cout<<str_core[i]<<" ";
    std::cout<<"\n";
    return 1;
  }





  //Fill the electron part of the potential
  wf.vdir.resize(wf.ngp);
  for(int i=0; i<wf.ngp; i++){
    double tmp = 0;
    if(Gf!=0) tmp += Gf*PRM_green(Z,wf.r[i],Gh,Gd);
    if(Tf!=0) tmp += Tf*PRM_tietz(Z,wf.r[i],Tt,Tg);
    wf.vdir[i] = tmp;
  }

  // Solve for each core state:
  int tot_el=0; // for working out Z_eff
  for(size_t i=0; i<core_list.size(); i++){
    int num = core_list[i];
    if(num==0) continue;

    int n = ATI_core_n[i];
    int l = ATI_core_l[i];

    //effective Z (for energy guess) -- not perfect!
    double Zeff =  double(Z - tot_el - num);
    if(l==1) Zeff = 1. + double(Z - tot_el - 0.5*num);
    if(l==2) Zeff = 1. + double(Z - tot_el - 0.5*num);
    if(Zeff<1.) Zeff=1.;
    tot_el+=num;

    double en_a = -0.5 * pow(Zeff/n,2);
    if(n>1) en_a *= 0.5;

    int k1 = l; //j = l-1/2
    if(k1!=0) {
      wf.solveLocalDirac(n,k1,en_a);
      en_a = 0.95*wf.en[wf.nlist.size()-1]; //update guess for next same l
    }
    int k2 = -(l+1); //j=l+1/2
    if(num>2*l) wf.solveLocalDirac(n,k2,en_a);

  }


  //make list of energy indices in sorted order:
  std::vector<int> sort_list;
  wf.sortedEnergyList(sort_list);

  //Output results:
  printf("\n n l_j    k Rinf its    eps      En (au)        En (/cm)\n");
  for(size_t m=0; m<sort_list.size(); m++){
    int i = sort_list[m];
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    printf("%2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %11.5f %15.3f\n",
        n,ATI_l(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni, eni*HARTREE_ICM);
  }



  tf = clock();
  double total_time = 1000.*double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\nt=%.3f ms.\n",total_time);

  return 0;
}