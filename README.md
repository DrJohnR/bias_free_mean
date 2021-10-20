# bias_free_mean
Implements the Metropolis algorithm to generate a Markov chain data set which conforms to a Gaussian distribution. Statistical thermalisation is (optionally) verified, and then
bootstrapping with binning is used to compute the bias-free mean of the data set.

A brief descriptive overview of the computational and statistical methods used within the program is presented below in (1), and then a guide on how to actually use the 
program is provided in (2).

-------------------------------- (1) --------------------------------


------ Monte Carlo methods ------

Monte Carlo methods are a varied class of computational algorithms which may be used to obtain numerical information from certain mathematical and physical models, where the use 
of other methods is impractical; their implementation is often most useful for systems which cannot be solved analytically. A defining feature of Monte Carlo algorithms is their 
reliance on repeated random sampling; a given problem may be 'solved' (whereby a desired result is extracted) by generating random numbers within an appropriate domain, and
observing what fraction of these data points satisfy some relevant property or condition. The accuracy of results obtained using a Monte Carlo algorithm is dependent on 
the number of random values generated in the process, with larger sample sizes typically inducing convergence towards a more accurate result.


------ Metropolis algorithm ------ 

The Metropolis algorithm is one example of a Monte Carlo method whose generated output is a sequence of random samples from a desired "target" distribution T(x); this 
discrete list (a Markov chain) will asymptotically approximate the target distribution for a sufficiently high number of samples. The algorithm is implemented by starting 
with a "proposal" distribution P(x) which is assumed to be proportional to the desired target distribution, and an arbitrary initialisation sample x_0 taken from the 
domain of P.

For each iteration ("sweep") of the algorithm:
1) Given the inital sample x_0, generate a candidate sample y = g(x_0) according to some (arbitrarily chosen, but unchanging) probability density function g. The idea is to 
	determine whether this candidate sample should be appended to the output sequence, based on its occurence probability relative to the previous sample.
2) Determine acceptance of candidate sample by computing acc = P(y)/P(x_0). If acc > 1 (i.e. the candidate is more probable) then y is immediately accepted and appended 
to the output sequence, is assigned as the new "initial" sample x_0 = y, and the algorithm returns to step 1). Otherwise, the algorithm proceeds to step 3).
3) A random value r is generated from the domain [0,1]. If acc >= r then the candidate sample y is accepted and appended to the output sequence, is assigned as the new "initial" 
	sample x_0 = y, and the algorithm returns to step 1). Otherwise, the candidate is rejected and the algorithm returns to step 1) using the same initial sample x_0.

Acceptance is a measure of how successfully the algorithm generates accepted sample points; it is defined as the proportion of the total number of sweeps that produces 
an accepted candidate sample. In the case of the Metropolis algorithm, a suitably tuned acceptance parameter ensures that the system doesn’t get "stuck" at any particular 
value, but doesn’t reach an equilibrium (the approximation of the target distribution) too slowly either. An estimate of acceptance should only be obtained from 
a thermalised system.   


------ Thermalisation ------

Thermalisation is the process whereby a data set is brought to a form of "statistical equilibrium"; this means that the mean of the data set varies only within an 
acceptable range (suitably defined by the user) as the number of resamples is changed. If a condition proposed to detect thermalisation is not satisfied, then one might 
(for example) discard a certain fraction of the total data set and check the condition again; the process can be reiterated until the condition is satisfied. The precise
implementation of a thermalisation condition for this program is described below in (2).


------ Bootstrapping with binning -------

Bias elimination is the process by which bias is removed (or minimised) for some computed quantity or observable such as the mean or standared deviation of a data set; bias
elimination is necessary within Monte Carlo simulations (such as the Metropolis algorithm) as it ensures that the expected value of a given statistical quantity is close to 
its true value, and is hence a more accurate representation of the data set. One example of such a bias elimination tool is "bootstrapping" with binned data. 

To apply the bootstrap method to binned data: 
1) Randomly choose N bins with replacement (i.e. the same bin may be chosen multiple times) to obtain a bootstrap sample, and then calculate the mean (referred to as the 
  sample mean) of this bootstrap sample. 
2) Repeat step 1) a for a sufficiently large number of iterations to obtain a set of separate sample means.
3) Compute the mean and standard deviation of these sample means in order to obtain the bias-free mean and bias-free error of the data set.



-------------------------------- (2) --------------------------------

The complete program is the file 'BiasFreeMean.c', and is comprised of four stages. Each of these stages has also been included as its own standalone file to facilitate
modification and to make the program modular; these four files should be executed in sequential order to replicate the complete program. The four stages are as follows:

1) ('Gaussian.c') 

Employs the Metropolis algorithm to generate a data set corresponding to random samples of a target Gaussian distribution. The parameters alpha and beta control the generator 
of candidate samples y, and the proportional distribution P, respectively; suggested ranges for these parameters are provided in comments. The user is requested to provide the 
number of Metropolis algorithm iterations ("sweeps") and the number of data points to discard (see step 3 below for details). By default the output data is written to a file 
named 'xnData.dat'. An acceptance ratio is returned to the user, and should ideally be between 0.5 --> 0.8 for optimal performance (refer to 'Metropolis algorithm' description
above). 

2) (DivideData.c) 

Simply imports the output data contained within 'xnData.dat' and divides it into two subsets of equal size, exporting each subset to its own .dat file; by default these files
are named 'HalfSample1.dat' and 'HalfSample2.dat' (the program handles an odd number of data samples by appending the additional point to the second of these files). The sole
purpose of this stage is to prepare the data for a thermalisation check based on the two data subsets, but this could easily be modified if the user wanted to impose a
thermalisation condition based on three or more subsets instead.

3) (Thermalisation.c)
 
Separately imports the two data subsets from 'HalfSample1.dat' and 'HalfSample2.dat', and computes their thermalisation statistics (mean, variance, standard deviation). The
complete data set is considered to be thermalised if both mean values lie within the range given by the other mean +/- its standard deviation. Explicitly:

thermalised if( mean1 - stand_dev1 <= mean2 <= mean1 + stand_dev1  &&  mean2 - stand_dev2 <= mean1 <= mean2 + stand_dev2 )

Note that, strictly speaking, the acceptance value returned by 'Gaussian.c' should only be computed after the data has been confirmed as thermalised. The complete program
'BiasFreeMean.c' instead gives the user the option to verify and impose the thermalisation condition, after a number of samples have already been discarded by an earlier prompt. 
If the condition is requested and fails then the program terminates, and the user should consider discarding more samples when prompted (dicarding approximately 10% of the 
total sweeps is typically sufficient to ensure thermalisation). 

4) (Binning.c)
The separate program 'binning.c' can be executed immediately after running 'Gaussian.c' (skipping the thermalisation checks) if required. This is also true for the complete
program 'BiasFreeMean.c' if the user inputs 'n' when prompted, though in both cases no thermalisation statistics will be computed / displayed.
 
The data set is imported from 'xnData.dat' and the user is prompted to provide the desired number of data bins (ensure that this is an exact divisor of the number of retained
data values), and the number of resamples (bootstrap samples) to generate; increasing the number of resamples will increase the computation time more significantly than
generating more bins. The bias-free mean of the data set is computed and returned to the user, and the program terminates.
