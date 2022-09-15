#include "/home/fr/fr_fr/fr_rs1121/PhD/code/itensor/itensor/all.h"
#include "/home/fr/fr_fr/fr_rs1121/PhD/code/itensor/itensor/util/print_macro.h"
#include "bits/stdc++.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include <sstream>
#include "/home/fr/fr_fr/fr_rs1121/PhD/code/headers/randutils.hpp"
using std::vector;
using std::move;

using namespace itensor;
using namespace std;

itensor::AutoMPO HamiltonianMPO(auto sites, int L, double t0, double Eu, double Ed)
{
	auto ampo = AutoMPO(sites);
	for(int i = 1; i < L; ++i)
    	{
		// NN hopping
	    	ampo += -t0,"Cdagup",i,"Cup",i+1;
		ampo += -t0,"Cdagup",i+1,"Cup",i;
		ampo += -t0,"Cdagdn",i,"Cdn",i+1;
		ampo += -t0,"Cdagdn",i+1,"Cdn",i;
	}
	for(int i = 1; i <= L; ++i)
	{
		// Onsite
		ampo += Eu,"Nup",i;
		ampo += Ed,"Ndn",i;
	}
	return ampo;
}


int main()
{
	// Parameters of run
	int N = 3;
	double t0 = 1.0;
	double Eu = 1.0;
	double Ed = -1.0;
	int indEE = 2;
	int mdim = 200;
	// 1D TB Lattice 
	auto sites = Electron(N);
	auto state = InitState(sites);
	for(int i = 1; i <= N; ++i)
	{	
		if (i==1) state.set(i,"UpDn");
		else state.set(i,"Emp");
	}
	auto psi0 = MPS(state);
	auto psi = MPS(state);
	// Initial MPO
	auto ampo = HamiltonianMPO(sites,N,t0,Eu,Ed);
	
	// Not necessary, but FYI
	auto H = toMPO(ampo);
	printfln("Maximum bond dimension of H is %d",maxLinkDim(H));
	// Generate unitary for time-evolution
	//	auto expH = toExpH(ampo,dt*Cplx_i);
//	auto args = Args("Method=","Fit","Cutoff=",1E-9,"MaxDim=",3000);
	auto args = Args("Method=","DensityMatrix","Cutoff=",1E-9,"MaxDim=",mdim);
	auto dt = 0.001;
	auto tend = 25;
	auto nt = int(tend/dt+(1e-9*(tend/dt)));
	complex<double> pu_l[nt][N];
	complex<double> pd_l[nt][N];
		
	double EE_l[nt];

	// Time-evolution by applying MPO exp(-i*H*dt) to MPS
	for(int n = 1; n <= nt; ++n)
	{
//		auto expH = toExpH(ampo,dt*Cplx_i);
		auto dexpH= toExpH(ampo,-dt*Cplx_i);
			
		// Time-evolution of normal psi
		psi = applyMPO(dexpH,psi,args);
		psi.noPrime().normalize();

		// Measure local electrons
		for (int j = 1; j <= N; j++)
		{
			auto ndu = op(sites,"Nup",j);
			auto ndd = op(sites,"Ndn",j);
			psi.position(j);
			auto ket = psi(j);
			auto bra = dag(prime(ket,"Site"));
			auto pu = eltC(bra*ndu*ket);
			auto pd = eltC(bra*ndd*ket);
			pu_l[n-1][j-1] = pu;
			pd_l[n-1][j-1] = pd;
		}

		//"Gauge" the MPS to site b
		psi.position(indEE);

		//SVD this wavefunction to get the spectrum
		//of density-matrix eigenvalues
		auto l = leftLinkIndex(psi,indEE);
		auto s = siteIndex(psi,indEE);
		auto [U,S,V] = svd(psi(indEE),{l,s});
		auto u = commonIndex(U,S);

		//Apply von Neumann formula
		//to the squares of the singular values
		Real SvN = 0.;
		for(auto j : range1(dim(u)))
		{
			auto Sn = elt(S,j,j);
		    	auto p = sqr(Sn);
		    	if(p > 1E-12) SvN += -p*log(p);
			{
				EE_l[n-1] = SvN;
			}
		}
	}
	// Save data to file
	stringstream streamT0,streamDT,streamMD;
	streamT0 << fixed << setprecision(2) << t0;
	streamDT << fixed << setprecision(3) << dt;
	streamMD << fixed << setprecision(1) << mdim;
	string filname = "N"+to_string(N)+"_t0"+streamT0.str()+"_mdim"+streamMD.str()+"_dt"+streamDT.str()+".dat";
	ofstream fil;
	fil.open(filname);
	fil << "# ---------- PARAMETERS -----------#\n" ;
	fil << "#N = " + to_string(N) + "\n";
	fil << "#t0 = " + to_string(t0) + "\n";
	fil << "#meth=MPO\n";
	fil << "# -------- THIS IS THE END --------#\n" ;	
	fil << "#(1)t " ;
	for (int j = 1; j<=N; j++)
	{
		fil << "(" << to_string(2*j) << ")pu"+to_string(j)+" "+"(" << to_string(2*j+1) << ")pd"+to_string(j)+" " ;
	}
	fil << " ("+to_string(2*N+2)+")EE"+to_string(indEE) << endl;
	for(int i = 1; i<=nt; i++)
	{
		fil << i*dt << " "; 
		for (int j = 0; j<N; j++)
		{
			fil << real(pu_l[i-1][j]) << " " << real(pd_l[i-1][j]) << " " ;
		}
		fil << EE_l[i-1] << endl;
	}
	fil.close();	
	return 0;
}
