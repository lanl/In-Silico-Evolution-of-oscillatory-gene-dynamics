/*

Continuous-time Markov Chain (Gillespie) simulation for the titration gene circuit.
Creator: Yen Ting Lin, CCS-3, LANL
Note: For the manuscript "Gene expression noise accelerates the evolution of a biological oscillator", co-authored by Nicolas E. Buchler, NCSU
The code has been reviewed by Richard P. Feynman Center for Innovation at the Los Alamos National Laboratory, with a C number C21109

This script generates a sample path of the repressilator model.
The script takes three inputs:
argv[1]: \beta^F_X of the gene circuit
argv[2]: \beta^B_X of the gene circuit
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

	int nA; 		// Number of binding sites on gene A
	int nB; 		// Number of binding sites on gene B

	double betaAf;	// Transcription rate of gene A without activation
	double betaAb;  // Transcription rate of gene A with activation
	double betaBf;	// Transcription rate of gene B without activation
	double betaBb;  // Transcription rate of gene B with activation

	double kappaA;	// Binding affinity to gene A
	double kappaB;	// Binding affinity to gene B

	double alpha;	// rate of forming heterodimers
	double gammaA;	// Degradation of the transcription factor A
	double gammaB;	// Degradation of the transcription factor B

	double delta; 	// dissociation rate

	int n_species;	// number of species
	int n_reactions;	// number of reactions

	double ** stoi;		// stoichiometry

	int ens_N;			// ensemble number

	double endT;
	double dt;
};

struct state
{

	int GA; 		// number of bound molecules on gene A
	int GB; 		// number of bound molecules on gene Binding

	int A; 			// number of transcription factor A
	int B; 			// number of transcription factor B

	double t;

};

void evolve_until_T(state *sta, parameter par, double tend);
double rndexclusive();
double exponential(double rate);

int main(int argc, char *argv[])
{

	srand (time(NULL));

	parameter par;

	double kappaX=1.2;
	double kappaY=0.9;
	double theta=1;
	double nA=3;
	double nB=3;
	double alpha=10;


	double betaFX=atof(argv[1]);
	double betaBX=atof(argv[2]);
	par.endT = 50;
	par.dt = 0.001;

	double betaBY=400;
	double betaFY=10;
	double gammaX=1;
	double gammaY=0.05;

	ofstream output;
	stringstream filename;
	filename << "./titration-"<< (int)atoi(argv[3])<< ".txt";
	output.open(filename.str().c_str());

	double Omega = 1000;

	par.n_species = 4;
	par.n_reactions = 9;

	par.ens_N = 1;

	double fastSwitching=50*Omega;

	par.nA=nA;
	par.nB=nB;

	par.betaAf = Omega*betaFX;	// Transcription rate of gene A without activation
	par.betaAb = Omega*betaBX;   // Transcription rate of gene A with activation
	par.betaBf = Omega*betaFY;	// Transcription rate of gene B without activation
	par.betaBb = Omega*betaBY;   // Transcription rate of gene B with activation

	par.kappaA=fastSwitching*kappaX/Omega;	// Binding affinity to gene A
	par.kappaB=fastSwitching*kappaY/Omega;	// Binding affinity to gene B
	par.delta = fastSwitching*theta;	// dissociation rate

	par.alpha=alpha/Omega;  // rate of forming heterodimers
	par.gammaA=gammaX;	// Degradation of the transcription factor
	par.gammaB=gammaY;

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
	// 2: GA
	// 3: GB

	par.stoi[1][0]  = +1.0;										// synthesis of A
	par.stoi[2][0]  = -1.0;										// degradation of A
	par.stoi[3][1]  = +1.0;										// synthesis of B
	par.stoi[4][1]  = -1.0;										// degradation of B


	par.stoi[5][0]=-1; par.stoi[5][1]=-1;						// forming heterodimer AB and removed from the system

	par.stoi[6][0]  = -1.0; par.stoi[6][2] = +1.0;				// an A binds to a site on gene A
	par.stoi[7][0]  = +1.0; par.stoi[7][2] = -1.0;				// dissociation of an A (bound to gene A)
	par.stoi[8][0]  = -1.0; par.stoi[8][3] = +1.0;				// an A binds to a site on gene A
	par.stoi[9][0]  = +1.0; par.stoi[9][3] = -1.0;				// dissociation of an B (bound to gene B)

	state *sta = new state [par.ens_N];

	for (int i=0;i<par.ens_N;i++)
	{
		sta[i].A =0;
		sta[i].B =0;
		sta[i].GA =0;
		sta[i].GB =0;
		sta[i].t = 0.0;
	}

	state sta_temp;


	for (double t=0; t<=par.endT; t+=par.dt)
	{

		//cout << t << endl;

		for (int i=0;i<par.ens_N;i++)
		{
			evolve_until_T(sta+i, par, t);
		}

		//output << t << "\t";

		output << sta[0].t << "\t";

		for (int i=0;i<par.ens_N;i++)
		{
			output << sta[i].GA << "\t";
		}

		for (int i=0;i<par.ens_N;i++)
		{
			output << sta[i].GB << "\t";
		}

		for (int i=0;i<par.ens_N;i++)
		{
			output << sta[i].A/Omega << "\t";
		}

		for (int i=0;i<par.ens_N;i++)
		{
			output << sta[i].B/Omega << "\t";
		}
		output << endl;

	}


	output.close();



}

void evolve_until_T(state *sta, parameter par, double tend)
{
	// compute all the possible reactions
	double *temp_rates = new double [par.n_reactions+1];

	while ((*sta).t < tend)
	{
		temp_rates[0]  =  0.0;

		if ((*sta).GA==par.nA)
			temp_rates[1]  =  par.betaAb;
		else
			temp_rates[1]  =  par.betaAf;

		temp_rates[2]  =  par.gammaA*(*sta).A;

		if ((*sta).GB==par.nB)
			temp_rates[3]  =  par.betaBb;
		else
			temp_rates[3]  =  par.betaBf;

		temp_rates[4]  =  par.gammaB*(*sta).B;




		temp_rates[5]  =  par.alpha   * (*sta).A   * (*sta).B;

		if ((*sta).GA==par.nA)
			temp_rates[6]  =  0;
		else
			temp_rates[6]  =  par.kappaA * (*sta).A;

		if ((*sta).GA==0)
			temp_rates[7]  =  0;
		else
			temp_rates[7]  =  par.delta;

		if ((*sta).GB==par.nB)
			temp_rates[8]  = 0;
		else
			temp_rates[8]  =  par.kappaB * (*sta).A;

		if ((*sta).GB==0)
			temp_rates[9] = 0;
		else
			temp_rates[9]  =  par.delta;


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
			(*sta).GA += par.stoi[i][2];
			(*sta).GB += par.stoi[i][3];


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
