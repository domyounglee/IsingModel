#include<stdio.h>
#include<math.h>

/* Ising model Hyperparameters */
#define X_max	16 // number of spin in x axis 
#define Y_max	16 // number of spin in y axis
#define N_spin  (X_max * Y_max) //totoa number of spin

#define T_init	4.0
#define T_final		1.0
#define dT		-0.2 // temperature of each step 

#define KB		1.0  // Boltzman Constant 1.0J/K
#define Beta(T)		(1.0/(KB*T))	//Define the Beta

#define MCS		5000	//maximum epoch 
#define MCS_init		2000	//get the ensembles after 2000 epoch
#define Max_step	(MCS*N_spin)	//real step 
#define MAX_init_step (MCS_init*N_spin)

/* Random number Hyperparameter */
#define R_A		16807
#define R_M		((int)pow(2,31)-1)
#define R_Q		127773
#define R_R		2836
#define RANDOM_MAX	(R_M-1) 	// max Random number 	

/* Global variables */
int spin[X_max][Y_max];		//spin as value of 1 or -1
float	T=T_init;					//Temperature
long 	N_accept;			//the number of changes by Metropolis algor
long RANDOM_NUMBER =1;

/* Generate random number 0~RANDOM_MAX */ 
long make_random_number(void){
	long U,V;

	U = RANDOM_NUMBER/R_Q;
	V = RANDOM_NUMBER%R_Q;
	RANDOM_NUMBER=R_A*V - R_R*U;
	if(RANDOM_NUMBER<0) RANDOM_NUMBER = RANDOM_NUMBER + R_M;
	return (RANDOM_NUMBER);
}

/* generate 0~1 float random number */
float make_unit_real_random_num(void){

	return ((float)make_random_number()/RANDOM_MAX);

}

/* generate integer random number from imin to imax*/
int make_int_random_number(int imin, int imax){

	return ((int)((float)make_random_number()/(RANDOM_MAX+1.)*(imax-imin+1))+imin);

} 

void set_random_number_seed(long seed){
	RANDOM_NUMBER=seed;
}

/*initialize spin */
void initialize_spin_configuration(void){
	int x,y;
	for(x=0;x<X_max;x++)
		for(y=0;y<Y_max;y++)
			spin[x][y]=make_int_random_number(0,1)*2 -1 ;
}

/*update using Metropolis algorithm*/
void update_spin_configuration(void){

	int x,y,x_left,x_right,y_up,y_down;
	float dE,A,P;
	
	//choose random spin
	x =  make_int_random_number(0,X_max-1);
	y =  make_int_random_number(0,Y_max-1);

	//nearest neighborhood
	x_left = x-1;
	x_right= x+1;
	y_up =   y+ 1 ;
	y_down = y-1 ;

	//consider boundary condition
	if(x==0){ x_left = X_max-1;}
	if(x==(X_max-1)) x_right = 0;
	if(y==0) y_up=Y_max-1;
	if(y==(Y_max-1)) y_down=0;

	//dE
	dE = 2*spin[x][y]*(spin[x_left][y]+spin[x_right][y]+spin[x][y_up]+spin[x][y_down]);
	
	/* Update with  Metropolis algorithm */
	A = exp(-1 * Beta(T) * dE);

	if(A>=1){

		spin[x][y] = -1*spin[x][y];
		N_accept++;
	} 
	else{
		P = make_unit_real_random_num();

		if(A>=P){

			spin[x][y]= -1*spin[x][y];
			N_accept++;
		}
	}
}

float get_energy(void){
	float energy =0.;
	int x,y,x_right,y_down;
	for(x=0;x<X_max;x++)
		for(y=0;y<Y_max;y++){
			x_right=x+1;
			y_down=y+1;
			if(x==(X_max-1)) x_right =0;
			if(y==(Y_max-1)) y_down =0;
			energy = energy - spin[x][y]*(spin[x][y_down] + spin[x_right][y]);
				
		}

	return energy;
}

float get_magnetization(void){
	float mag =0.;
	int x,y;
	
	for(x=0;x<X_max;x++)
		for(y=0;y<Y_max;y++){
			mag+=spin[x][y];		
		}

	return mag;
}


void main(){

	long i;
	float energy, energy2, energy_per_spin, E,Eavg2,E2avg;
	float magnetization, magnetization_per_spin,M;
	float specific_heat,specific_heat_per_spin;
	float accept_rate;

	set_random_number_seed(1);
	/*initialize spin*/
	initialize_spin_configuration();

	printf("# %d x %d Ising model \n",X_max,Y_max);
	printf("From %4.1fK to %4.1fK, Step : %4.1fk) \n",T_init, T_final,dT);
	printf("# Monte Carlo Step (MCS) : %d, Discarded MCS : %d\n", MCS, MCS_init);
	
	while(T>=T_final){
		
		energy=energy2=0;
		magnetization=0;
		specific_heat=0;
		Eavg2=E2avg=0;

		/*until MAX_init_step just update*/
		for(i=0;i<MAX_init_step;i++)
			update_spin_configuration();




		N_accept=0;

		/* after MCS_init step*/
		for(i=0;i<Max_step;i++){
			update_spin_configuration();
			E=get_energy();
			//printf("%.1f",E);
			energy=energy+E;
			energy2=energy2+E*E;
			M=get_magnetization();
			magnetization=magnetization+M;

		} 

		/*measure each physical state*/
		energy= energy/Max_step;
		energy_per_spin = energy/N_spin;
		Eavg2=energy*energy;
		E2avg= energy2/Max_step;
		specific_heat = KB*Beta(T)*Beta(T)*(E2avg-Eavg2);
		specific_heat_per_spin = specific_heat/N_spin;
		magnetization = magnetization/Max_step;
		magnetization_per_spin = magnetization/N_spin;
		accept_rate = (float)N_accept/Max_step;
		for(int x=0;x<X_max;x++){
			for(int y=0;y<Y_max;y++){
				if(spin[x][y]==-1){
				printf("0  ");	
				continue;}
				printf("%d  ",spin[x][y]);		
			}
			printf("\n");
		}


		printf("%.1f %f %f %f %f\n",
			T,energy_per_spin,specific_heat_per_spin,magnetization_per_spin,accept_rate);

		T+=dT;
	}
}

