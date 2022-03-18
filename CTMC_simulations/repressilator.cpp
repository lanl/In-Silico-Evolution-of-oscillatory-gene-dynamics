/*

Continuous-time Markov Chain (Gillespie) simulation for the repressilator gene circuit.
Creator: Yen Ting Lin, CCS-3, LANL
Note: For the manuscript "Gene expression noise accelerates the evolution of a biological oscillator", co-authored by Nicolas E. Buchler, NCSU
The code has been reviewed by Richard P. Feynman Center for Innovation at the Los Alamos National Laboratory, with a C number C21109

This script generates a sample path of the repressilator model.
The script takes three inputs:
argv[1]: r_0 of the gene circuit
argv[2]: r_1 of the gene circuit
argv[3]: a unique ID for specifying the sample path

*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>

using namespace std;

struct parameter
{

	double h;	// Hill coefficient
	double r0;  // Basal transcription rate
	double r1;	// Fully suppresed transcription rate
	double k;  // Transcription rate of gene B with activation

	double gamma;	// Degradation rate

	double Omega; // system size

	int n_species;	// number of species
	int n_reactions;	// number of reactions

	double ** stoi;		// stoichiometry

	int ens_N;			// ensemble number

	double endT;
	double dt;
};

struct state
{

	int A; 			// number of transcription factor A
	int B; 			// number of transcription factor B
	int C; 			// number of transcription factor B

	double t;

};

void evolve_until_T(state *sta, parameter par, double tend);
double rndexclusive();
double exponential(double rate);

int main(int argc, char *argv[])
{

	srand (time(NULL));

	parameter par;

	par.Omega = 500;

	par.r0=atof(argv[1]);
	par.r1=atof(argv[2]);
	par.k=0.5;
	par.h=3.0;
	par.gamma=1.0;

	par.r0 *= par.Omega;
	par.r1 *= par.Omega;
	par.k *= par.Omega;

	par.n_species = 3;
	par.n_reactions = 6;

	par.endT = 30;
	par.dt = 0.01;

	ofstream output;
	stringstream filename;
	filename << "./repressilator-"<< (int)atoi(argv[3])<< ".txt";
	output.open(filename.str().c_str());

	par.n_species = 3;
	par.n_reactions = 6;
	par.ens_N = 1;

	par.stoi = new double * [par.n_reactions+1];

	for (int i=0;i<=par.n_reactions;i++)
	{
		par.stoi[i] = new double [par.n_species];

		for (int j=0;j<par.n_species;j++)
		{
			par.stoi[i][j]=0;
		}

    }


	// 0: A
	// 1: B
	// 2: C

	par.stoi[1][0]  = +1.0;										// synthesis of A
	par.stoi[2][0]  = -1.0;										// degradation of A
	par.stoi[3][1]  = +1.0;										// synthesis of B
	par.stoi[4][1]  = -1.0;										// degradation of B
	par.stoi[5][2]  = +1.0;										// synthesis of C
	par.stoi[6][2]  = -1.0;										// degradation of C

	state *sta = new state [par.ens_N];

	for (int i=0;i<par.ens_N;i++)
	{
		sta[i].A =0.1*par.Omega;
		sta[i].B =0;
		sta[i].C =0;
		sta[i].t = 0.0;
	}

	state sta_temp;

	for (double t=0; t<=par.endT; t+=par.dt)
	{

		for (int i=0;i<par.ens_N;i++)
		{
			evolve_until_T(sta+i, par, t);
		}

		output << sta[0].t << "\t";

		for (int i=0;i<par.ens_N;i++)
		{
			output << sta[i].A/par.Omega << "\t";
		}

		for (int i=0;i<par.ens_N;i++)
		{
			output << sta[i].B/par.Omega << "\t";
		}

		for (int i=0;i<par.ens_N;i++)
		{
			output << sta[i].C/par.Omega << "\t";
		}

		output << endl;

	}

	output.close();

  return 0;

}

void evolve_until_T(state *sta, parameter par, double tend)
{
	// compute all the possible reactions
	double *temp_rates = new double [par.n_reactions+1];

	while ((*sta).t < tend)
	{
		temp_rates[0]  =  0.0;

		temp_rates[1] = par.r0 + (par.r1-par.r0) * pow((*sta).C, par.h)/(pow((*sta).C, par.h)+pow(par.k, par.h));
		temp_rates[2] = par.gamma*(*sta).A;

		temp_rates[3] = par.r0 + (par.r1-par.r0) * pow((*sta).A, par.h)/(pow((*sta).A, par.h)+pow(par.k, par.h));
		temp_rates[4] = par.gamma*(*sta).B;

		temp_rates[5] = par.r0 + (par.r1-par.r0) * pow((*sta).B, par.h)/(pow((*sta).B, par.h)+pow(par.k, par.h));
		temp_rates[6] = par.gamma*(*sta).C;

    for (int i=1;i<par.n_reactions+1;i++)
		{
			temp_rates[i] = temp_rates[i] + temp_rates[i-1];
		}

		double dt = exponential(temp_rates[par.n_reactions]);

		if ((*sta).t + dt < tend)
		{
			(*sta).t += dt;

			double r = rndexclusive()*temp_rates[par.n_reactions];

			int i=0;

			while (temp_rates[i]<r)
				i++;

			(*sta).A += par.stoi[i][0];
			(*sta).B += par.stoi[i][1];
			(*sta).C += par.stoi[i][2];

		}else{
			(*sta).t = tend;
		}

	}


}

double exponential(double rate)
{
	double out = -1.0/(rate) * log(rndexclusive());
	return out;
}

double rndexclusive()
{
	double i=0;

	while ((i==0)||(i==1))
	{
		i = double(rand())/RAND_MAX;
	}

	return i;
}
