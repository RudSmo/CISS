#include "/home/fr/fr_fr/fr_rs1121/PhD/code/itensor/itensor/all.h"
#include "/home/fr/fr_fr/fr_rs1121/PhD/code/itensor/itensor/util/print_macro.h"
#include "bits/stdc++.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include <sstream>
#include "/home/fr/fr_fr/fr_rs1121/PhD/code/headers/randutils.hpp"
#include <complex>
using std::vector;
using std::move;

using namespace itensor;
using namespace std;
using namespace std::complex_literals;

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

int main()
{
   int n = 3;
   double L = 9;
   double ahel = 1.0;
   double chel = 1.0;
   double tso = 1.0;
   int i = 4;
   double num_laps = 1;

   double rj[n];
   double rj1[n];
   double rj2[n];

   UnitVectorHelix(ahel,chel,i,num_laps,L,rj);   
   UnitVectorHelix(ahel,chel,i+1,num_laps,L,rj1);
   UnitVectorHelix(ahel,chel,i+2,num_laps,L,rj2);

   double rjrj1[n];
   Diff(rj,rj1,rjrj1,n);
   double rjrj2[n];
   Diff(rj,rj2,rjrj2,n);
   
   double dj1[n],dj2[n];
   double n1,n2;
   n1 = absv(rjrj1,n);
   n2 = absv(rjrj2,n);
   for (int i = 0; i<n; i++)
   {
	   dj1[i] = rjrj1[i]/n1;
	   dj2[i] = rjrj2[i]/n2;
   }   

   double vj[n];
   crossProduct(dj1,dj2,vj);
  
   complex<double> Vj[n];
   vDot(vj,tso,Vj);    

   for (int i = 0; i<=0; i++)
   {
	   cout << Vj[i] << endl;
   }
   return 0;
}
