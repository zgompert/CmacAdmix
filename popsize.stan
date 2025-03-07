data{
	int N; /* # of observations */
	int Np; /* # of source populations, including hybrid */
	int Nr; /* # of replicates */
	int<lower=0> Y[N]; /* vector of cumulative counts */
	int<lower=0> reps[N]; /* vector of rep number */
	vector[N] day; /*  vector of days */
	vector[N] day2; /* vector of days2 */
	matrix<lower=0, upper=1>[N,Np] pops; /* design matrix of pop contribution weights */
	
}

parameters{
	vector[Np] BpopD; /* population day effects */
	vector[Np] BpopD2; /* population day2 effects */
	vector[Nr-1] AlphaD_raw; /* rep random effect for slope */
	vector[Nr-1] AlphaD2_raw; /* rep random effect for slope2 */
	real sigma; /* SD for residual error */
	real sigmaA; /* SD for rep effects on slope */
	real sigmaA2; /* SD for rep effects on slope2 */
}

transformed parameters{
	vector[Nr] AlphaD = append_row(AlphaD_raw, -sum(AlphaD_raw)); /* for s2z constraint */
	vector[Nr] AlphaD2 = append_row(AlphaD2_raw, -sum(AlphaD2_raw)); /* for s2z constraint */

}
model{

	/* GLM with Poisson likelihood */
	for(i in 1:N){
		target += normal_lpdf(Y[i] | (pops[i,] * BpopD) * day[i] + (pops[i,] * BpopD2) * day2[i] + 
			AlphaD[reps[i]] * day[i] + AlphaD2[reps[i]] * day2[i], sigma);
	} 


	/* increment prior on pop-specific betas */
	for(j in 1:Np){
		target += normal_lpdf(BpopD[j] | 0, 100);
		target += normal_lpdf(BpopD2[j] | 0, 100);

	}
	
	/* increment prior on replicate effects */
	for(j in 1:(Nr-1)){
		target += normal_lpdf(AlphaD_raw[j] | 0, sigmaA);
		target += normal_lpdf(AlphaD2_raw[j] | 0, sigmaA2);
	}

	/* increment priors on SDs */
	target += gamma_lpdf(sigma | 0.1, 0.01);
	target += gamma_lpdf(sigmaA | 0.1, 0.01);
	target += gamma_lpdf(sigmaA2 | 0.1, 0.01);
}

generated quantities {
        /* expected values for each observation */
        real mu[N]; /* vector of expected counts */
        for(i in 1:N){
		mu[i] = (pops[i,] * BpopD) * day[i] + (pops[i,] * BpopD2) * day2[i] +
			AlphaD[reps[i]] * day[i] + AlphaD2[reps[i]] * day2[i];
	}
}

