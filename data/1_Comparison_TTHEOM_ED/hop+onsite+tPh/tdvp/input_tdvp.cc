#include "/home/fr/fr_fr/fr_rs1121/PhD/code/itensor/itensor/all.h"
#include "/home/fr/fr_fr/fr_rs1121/PhD/code/itensor/itensor/util/print_macro.h"
#include "bits/stdc++.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include <complex>
#include <sstream>
#include "/home/fr/fr_fr/fr_rs1121/PhD/code/headers/randutils.hpp"
#include "/home/fr/fr_fr/fr_rs1121/PhD/code/TDVP/tdvp.h"
#include "/home/fr/fr_fr/fr_rs1121/PhD/code/TDVP/basisextension.h"

using std::vector;
using std::move;

using namespace itensor;
using namespace std;

void crossProduct(double v_A[], double v_B[], double c_P[])
{
   c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

void vDot(double A[], double tso, complex<double> O[])
{
   std::complex<double> im(0.0,1.0);
   O[0] = im*tso*A[2];
   O[1] = im*tso*(A[0] - A[1]*im);
   O[2] = im*tso*(A[0] + A[1]*im);
   O[3] = im*tso*(-A[2]);
}
void Diff(double v1[], double v2[], double vo[], int n)
{
        for (int i = 0; i<n; i++)
        {
                vo[i] = v1[i]-v2[i];
        }
}
void UnitVectorHelix(double ahel, double chel, double i, double num_laps, double L, double rj[])
{
   const double PI  =3.141592653589793238463;
   auto phij = 2*PI*(i-1)*num_laps/(L-1);
   auto cj = chel*(i-1)/(L-1);
   rj[0] = ahel*cos(phij);
   rj[1] = ahel*sin(phij);
   rj[2] = cj;
}
double absv(double v[], int n)
{
        double a = 0.0;
        for (int i = 0; i < n; i++)
        {
                a += v[i]*v[i];
        }
        return sqrt(a);
}

itensor::AutoMPO HamiltonianMPO(auto sites, int L, int num_laps, double ahel, double chel, double t0, double Eu, double Ed, double tso, double wPh, double tPh)
{
	auto ampo = AutoMPO(sites);
	std::complex<double> im(0.0,1.0);
	for(int i = 1; i <= L-2; i+=2)
    	{
		// NN hopping

	    	ampo += -t0,"Cdagup",i,"Cup",i+2;
		ampo += -t0,"Cdagup",i+2,"Cup",i;
		ampo += -t0,"Cdagdn",i,"Cdn",i+2;
		ampo += -t0,"Cdagdn",i+2,"Cdn",i;
	}
	for(int i = 1; i <= L-1; i+=2)
	{
		// Onsite
		ampo += Eu,"Cdagup",i,"Cup",i;
		ampo += Ed,"Cdagdn",i,"Cdn",i;
	}
	for (int i = 2; i <= L; i += 2)
	{
		// Phonon harmonic energy
		ampo += wPh,"Adag",i,"A",i;
	}
	for (int i = 1; i <= L-3; i+=2)
	{
		// NN el-ph coupling

		ampo += -tPh,"Cdagup",i,"Cup",i+2,"A",i+1;
		ampo += -tPh,"Cdagup",i,"Cup",i+2,"Adag",i+1;
		ampo += -tPh,"Cdagup",i+2,"Cup",i,"A",i+1;
                ampo += -tPh,"Cdagup",i+2,"Cup",i,"Adag",i+1;
		ampo += -tPh,"Cdagdn",i,"Cdn",i+2,"A",i+1;
                ampo += -tPh,"Cdagdn",i,"Cdn",i+2,"Adag",i+1;
                ampo += -tPh,"Cdagdn",i+2,"Cdn",i,"A",i+1;
                ampo += -tPh,"Cdagdn",i+2,"Cdn",i,"Adag",i+1;
	}
	return ampo;
}


int main()
{
	// Parameters of run
	int fS = 3;
	int N = 2*fS;
	auto dt = 0.005;
	auto tend = 50;
	int num_laps = 1;
	double ahel = 1.0;
	double chel = 1.0;
	double t0 = 1.0;
	double Eu = 1.0;
	double Ed = -1.0;
	double tso = 0.0;//t0/40;
	double wPh = 1.4*t0; 
	double tPh = 0.1;//0.1*t0;
	int indEE = 2;
	// 1D TB Lattice 	
	using Holstein = MixedSiteSet<ElectronSite,BosonSite>;	
	auto sites = Holstein(N,{"ConserveQNs=",false,"MaxOcc=",5});
	auto state = InitState(sites);

	for(int i = 1; i <= N-1; i+=2)
	{	
		if (i==1) state.set(i,"UpDn");
		else state.set(i,"Emp");
	}
	
	for(int i = 2; i<=N-2; i+=2)
	{
		state.set(i,"0");
	}

	//auto psi0 = MPS(state);
	auto psi = MPS(state);

	// Initial MPO
	auto ampo = HamiltonianMPO(sites,N,num_laps,ahel,chel,t0,Eu,Ed,tso,wPh,tPh);
	cout << ampo << endl;
	// Not necessary, but FYI
	auto H = toMPO(ampo);
	printfln("Maximum bond dimension of H is %d",maxLinkDim(H));
	auto nt = int(tend/dt+(1e-9*(tend/dt)));
	double pu_l[nt][N];
	double pd_l[nt][N];
		
	double EE_l[nt];

	auto sweeps = Sweeps(1);
	sweeps.maxdim() = 3000;
	sweeps.cutoff() = 1E-11;
	sweeps.niter() = 50;
	// Time-evolution by using tdvp 
	for(int n = 1; n <= nt; ++n)
	{
		//if (n<3)
		//{
		// Global subspace expansion
		std::vector<Real> epsilonK = {1E-10, 1E-10};
		addBasis(psi,H,epsilonK,{"Cutoff",1E-8,
					 "Method","DensityMatrix",
					 "KrylovOrd",3,
					 "DoNormalize",true,
					 "Quiet",true});
		//}
		tdvp(psi,H,-dt*Complex_i,sweeps,{"Truncate",true,
                                        "DoNormalize",true,
                                        "Quiet",true,
                                        "NumCenter",1,
                                        "ErrGoal",1E-9});
		
		// psi.noPrime().normalize();
		// Measure local electrons
		int k = 0;
		for (int j = 1; j <= N-1; j+=2)
		{
			auto ndu = op(sites,"Nup",j);
			auto ndd = op(sites,"Ndn",j);
			psi.position(j);
			auto ket = psi(j);
			auto bra = dag(prime(ket,"Site"));
			auto pu = real(eltC(bra*ndu*ket));
			auto pd = real(eltC(bra*ndd*ket));
			pu_l[n-1][k] = pu;
			pd_l[n-1][k] = pd;
			k += 1;
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
	stringstream streamT0,streamTSO,streamDT,streamMD,streamTPH;
	streamT0 << fixed << setprecision(2) << t0;
	streamDT << fixed << setprecision(5) << dt;
	streamTSO << fixed << setprecision(3) << tso;
	streamTPH << fixed << setprecision(3) << tPh; 

	string filname = "N"+to_string(fS)+"_t0"+streamT0.str()+"_tso"+streamTSO.str()+"_tPh"+streamTPH.str()+"_dt"+streamDT.str()+".dat";
	ofstream fil;
	fil.open(filname);
	fil << "# ---------- PARAMETERS -----------#\n" ;
	fil << "#N = " + to_string(N) + "\n";
	fil << "#Eu = " + to_string(Eu) + "\n";
	fil << "#Ed = " + to_string(Ed) + "\n";
	fil << "#t0 = " + to_string(t0) + "\n";
	fil << "#tso = " + to_string(tso) + "\n";
	fil << "#wPh = " + to_string(wPh) + "\n";
        fil << "#tPh = " + to_string(tPh) + "\n";	
	fil << "#meth=MPO\n";
	fil << "# -------- THIS IS THE END --------#\n" ;	
	fil << "#(1)t " ;
	for (int j = 1; j<=fS; j++)
	{
		fil << "(" << to_string(2*j) << ")pu"+to_string(j)+" "+"(" << to_string(2*j+1) << ")pd"+to_string(j)+" " ;
	}
	fil << " ("+to_string(2*N+2)+")EE"+to_string(indEE) << endl;
	for(int i = 1; i<=nt; i++)
	{
		fil << i*dt << " "; 
		for (int j = 0; j<fS; j++)
		{
			fil << real(pu_l[i-1][j]) << " " << real(pd_l[i-1][j]) << " " ;
		}
		fil << EE_l[i-1] << endl;
	}
	fil.close();	
	return 0;
}
