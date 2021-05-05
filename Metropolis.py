import mc3
import numpy as np

"""
Based on: https://mc3.readthedocs.io/en/latest/mcmc_tutorial.html, last viewed in November 2020

Required arguments for mc3.sample():

data: input data as array
uncert: corresponding array of uncertainty
func: parameterized modeling function for fitting the data
params: initial-guess values for the model fitting parameters as array
sampler: Sampler algorithm, three choices: 'snooker', 'demc' or 'mrw'. Snooker is recommended.
    Snooker uses the demc-zs algorithm with snooker proposals. mrw uses the classical Metropolis-Hastings algorithm.
nsamples: total number of MCMC samples as integer

Optional arguments for mc3.sample():

indparams: contains any additional argument required by func.
pmin: Default: -np.inf. Lower exploration boundary of MCMC as float array. 
pmax: Default: np.inf. Upper exploration boundary of MCMC as float array. 
pstep: Default values: 1. Positive values in the array denote step size. The value 0 keeps a parameter fixed.
    Negative integer force the parameter to use the value of the parameter in position of the value of the negative int
prior: Default values: 0. Prior estimate as float array. 
priorlow: Default values: 0. Prior lower uncertainty as float array. 
priorup: Default values: 0. Prior upper uncertainty as float array. 
pnames: defines names of model parameters up to 11 characters for screen output. Array
texnames: enables names using LaTeX syntax. Array. If both pnames and texnames are None defaults to [Param 1, ...] 
burnin: Default: 0. Sets the number of burned-in / removed iterations at the beginning of each chain. 
nchains: Default: 7. Sets the number of parallel chains to use. 
ncpu: Default: nchains. Sets the number of CPUs to use for the chains. The central MCMC hub uses an extra CPU. 
thinning: Default: 1. Sets the thinning factor, all but every thinning-th sample gets discarded. Reduces memory usage. 
kickoff: Sets starting point of MCMC chains. 'normal' (Default) for normal-distribution draw at params with 
    standard deviation pstep or set 'uniform' for a uniform draw between pmin and pmax.  
hsize: Default: 10. Initial sample size that doesn't count towards posterior-distribution statistics. Needed by snooker.
leastsq: Default: None. If needed, a least-squares optimization can be run before the MCMC.
    Choose 'lm' for the Levenberg-Marquardt algorithm or 'trf' for the Trust Region Refelctive algorithm.
chisqscale: Default: False. If enabled, multiplies all uncertainties by a common scale factor.
grtest: Default: False. If enabled runs the Gelman-Rubin convergence test.
grbreak: Default: 0.0. Sets a convergence threshold to stop an MCMC when GR drops below grbreak.
grnmin: Default: 0.5.  
wlike: Default: False. Enables the Wavelet-based method to account for time-correlated noise. Three additional
    fitting parameters have to be appended to params, and added to pmin, pmax, stepsize, prior, priorlow and priorup.
fgamma: Default: 1.0. Scale factor for DEMC's gamma jump.
fepsilon: Default: 0.0. Jump scale factor for DEMC's "e" distribution.
log: Default: None. A filename or mc3.utils.Log object can be given so the screen output will be saved.
savefile: Default: None. A filename can be given to store the MCMC outputs in .npz format.
plots: Default: False. Enables the generation of a data plot, the MCMC-chain trace plot for each parameter 
    and the marginalized and pair-wise posterior  plots.
ioff: Default: False. True will turn the display interactive mode off.
rms: Default: False. Enables the computing and plotting of a time-averaging test for time-correlated noise. 
"""

# Beginning Tutorial for MC3:

def quad(p, x):
    """

    :param p: Polynomial constant, linear and quadratic coefficients.
    :param x: Array of dependent variables where to evaluate the polynomial.

    :return: y: Polinomial evaluated at x: y(x) = p0 + p1*x + p2*x²)

    """

    y = p[0] + p[1] * x + p[2] * x ** 2.0
    return y


# Creating synthetic dataset

# 1000 evenly spaced numbers between 0 and 10 on quadratic function
x = np.linspace(0, 10, 1000)
p = [3.14, -2.4, 1.8]
y = quad(p, x)

# Determining array of values on the basis of y values for shifting y values
uncert = np.sqrt(np.abs(y))

# Calculating normally distributed values on the basis of y
error = np.random.normal(0, uncert)

# Adding dependent normally distributed values to quadratic function y
data = y + error

# Defining the modeling function as a callable:
func = quad

# List of additional arguments of func:
indparams = [x]

# Creating an array with initial-guess values of fitting parameters:
params = np.array([10.0, -3.0, 0.5])

# Lower and upper boundaries for the MCMC exploration.
# Proposed steps outside this boundary are automatically rejected.
pmin = np.array([-10.0, -20.0, -10.0])
pmax = np.array([ 40.0,  20.0,  10.0])

# Parameter names
pnames = ['y0', 'alpha', 'beta']

# Sampler algorithm, here Metropolis-Hastings
sampler = 'mrw'

# MCMC setup
nsamples = 1e6
burnin = 200
nchains = 2
ncpu = 1

# Optimization
leastsq = 'lm'

# MCMC Convergence
grtest = True
grbreak = 1.01
grnmin = 0.5

# Logging
log = 'MCMC_tutorial.log'

# File outputs
savefile = 'MCMC_tutorial.npz'
plots = True
rms = True

mc3_output = mc3.sample(data=data, uncert=uncert, func=func, params=params, indparams=indparams, pmin=pmin, pmax=pmax,
                        pnames=pnames, sampler=sampler, nsamples=nsamples, nchains=nchains, ncpu=ncpu,
                        burnin=burnin, leastsq=leastsq, grtest=grtest, grbreak=grbreak, grnmin=grnmin, log=log,
                        plots=plots, savefile=savefile, rms=rms)


def helpmc3():
    help(mc3.sample)


