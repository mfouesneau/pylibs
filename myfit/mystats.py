"""
Statistical package

CONTENT:
	calc_min_interval	Internal method to determine the minimum interval of a given width
	hpd			Calculate HPD (minimum width BCI) of array for given alpha
	make_indices 		Generates complete set of indices for given dimensions
	MAP			Returns a MAP estimate, using error as default resolution
	optbins			Determine the optimal binning of the data based on common estimators
	quantiles		Computes quantiles from an array
	stats			Returns common statistics of a sample

REQUIREMENTS:
	numpy		mandatory

HISTORY:
	0.1   -- 15/09/2011, MF, Initial
			First clollection of statistical functions 
	0.2   -- 10/10/2011, MF, added functions, added full documentation
		 	MAP     - estimate the maximum a posteriori
			optBins - Estimate the optimal binning of a sample
	0.2.1 -- 10/10/2011, MF, corrected hpd for N-D data estimation
	
"""
__version__ = '0.2.1'

import numpy
from math import *

def quantiles(x, qlist=[2.5, 25, 50, 75, 97.5]):
    """computes quantiles from an array

	Quantiles :=  points taken at regular intervals from the cumulative
	distribution function (CDF) of a random variable. Dividing ordered data
	into q essentially equal-sized data subsets is the motivation for
	q-quantiles; the quantiles are the data values marking the boundaries
	between consecutive subsets.

	The quantile with a fraction 50 is called the median 
	(50% of the distribution)	
	
	Inputs: 
		x     - variable to evaluate from
		qlist - quantiles fraction to estimate (in %)

	Outputs: 
		Returns a dictionary of requested quantiles from array
    """
    # Make a copy of trace
    x = x.copy()

    # For multivariate node
    if x.ndim>1:
        # Transpose first, then sort, then transpose back
        sx = numpy.transpose(numpy.sort(numpy.transpose(x)))
    else:
        # Sort univariate node
        sx = numpy.sort(x)

    try:
        # Generate specified quantiles
        quants = [sx[int(len(sx)*q/100.0)] for q in qlist]

        return dict(zip(qlist, quants))

    except IndexError:
        print "Too few elements for quantile calculation"



def calc_min_interval(x, alpha):
    """Internal method to determine the minimum interval of
    a given width"""

    # Initialize interval
    min_int = [None,None]

    try:

        # Number of elements in trace
        n = len(x)

        # Start at far left
        start, end = 0, int(n*(1-alpha))

        # Initialize minimum width to large value
        min_width = numpy.inf

        while end < n:

            # Endpoints of interval
            hi, lo = x[end], x[start]

            # Width of interval
            width = hi - lo

            # Check to see if width is narrower than minimum
            if width < min_width:
                min_width = width
                min_int = [lo, hi]

            # Increment endpoints
            start +=1
            end += 1

        return min_int

    except IndexError:
        print 'Too few elements for interval calculation'
        return [None,None]

def make_indices(dimensions):
    """ Generates complete set of indices for given dimensions """

    level = len(dimensions)

    if level==1: return range(dimensions[0])

    indices = [[]]

    while level:

        _indices = []

        for j in range(dimensions[level-1]):

            _indices += [[j]+i for i in indices]

        indices = _indices

        level -= 1

    try:
        return [tuple(i) for i in indices]
    except TypeError:
        return indices


def hpd(x, alpha):
    """Calculate HPD (minimum width BCI) of array for given alpha"""

    # Make a copy of trace
    x = x.copy()

    # For multivariate node
    if x.ndim>1:

        # Transpose first, then sort
        tx = numpy.transpose(x, range(x.ndim)[1:]+[0])
        dims = numpy.shape(tx)

        # Container list for intervals
        intervals = numpy.resize(0.0, dims[:-1]+(2,))

        for index in make_indices(dims[:-1]):

            try:
                index = tuple(index)
            except TypeError:
                pass

            # Sort trace
            sx = numpy.sort(tx[index])

            # Append to list
            intervals[index] = calc_min_interval(sx, alpha)

        # Transpose back before returning
        return numpy.array(intervals)

    else:
        # Sort univariate node
        sx = numpy.sort(x)

        return numpy.array(calc_min_interval(sx, alpha))

def MAP(x, error=0.01, *args, **kwargs):
	""" Return a MAP estimate, using error as default resolution """
	if not ('bins' in kwargs):
		bins = numpy.linspace(numpy.min(x), numpy.max(x), numpy.round((numpy.max(x)-numpy.min(x))/error))
		n, b = numpy.histogram(x, bins=bins, *args, **kwargs)	
	else:
		n, b = numpy.histogram(x, *args, **kwargs)	
	
	c = 0.5*(b[:-1]+b[1:])
	return c[n.argmax()]

def optbins(data, method='freedman', ret='N'):
	""" Determine the optimal binning of the data based on common estimators
	and returns either the number of bins of the width to use.

	input:
		data	1d dataset to estimate from
	keywords:
		method	the method to use: str in {sturge, scott, freedman}
		ret	set to N will return the number of bins
			set to W will return the width
	refs:
		Sturges, H. A. (1926)."The choice of a class interval". J. American Statistical Association, 65-66
		Scott, David W. (1979), "On optimal and data-based histograms". Biometrika, 66, 605-610
		Freedman, D.; Diaconis, P. (1981). "On the histogram as a density estimator: L2 theory". 
				Zeitschrift fur Wahrscheinlichkeitstheorie und verwandte Gebiete, 57, 453-476
	"""
	x = numpy.asarray(data)
	n = x.size
	r = x.max()-x.min()

	def sturge():
		if (n<=30):
			print "Warning: Sturge estimator can perform poorly for small samples"
		k = int(numpy.log(n)+1)
		h = r/k
		return h,k
	
	def scott():
		h = 3.5*numpy.std(x)*float(n)**(-1./3.)
		k = int(r/h)
		return h,k
	
	def freedman():
		q = quantiles(x, [25, 75])
		h = 2*(q[75]-q[25])*float(n)**(-1./3.)
		k = int(r/h)
		return h,k
	
	m = {'sturge':sturge, 'scott':scott, 'freedman': freedman}

	if method.lower() in m:
		s = m[method.lower()]()
		if ret.lower() == 'n':
			return s[1]
		elif ret.lower() == 'w':
			return s[0]
	else:
		return None
	
	
	
	

def stats(var, alpha=0.05, qlist=[2.5, 25, 50, 75, 97.5]):
	assert(len(var) > 0), 'Variable is empty'
	
	return {
                'n': len(var),
                'standard deviation': numpy.std(var),
                'mean': numpy.mean(var),
		'variation range': [numpy.min(var), numpy.max(var)],
                '%s%s HPD interval' % (int(100*(1-alpha)),'%'): hpd(var, alpha),
                'quantiles': quantiles(var, qlist=qlist)
		}

def simpleFit(x0,y0, yerr0=None, p0=None, expression=None):
	from scipy.optimize import leastsq
	if expression == None:
		return None
	fitfunc  = lambda p,x: eval(expression)
	if yerr0 == None:
		errfunc  = lambda p,x,y,yerr: (fitfunc(p,x)-y) 
	else:
		errfunc  = lambda p,x,y,yerr: (fitfunc(p,x)-y) / yerr
	errfuncsumsq = lambda p,x,y,yerr: numpy.sum(errfunc(p,x,y,yerr)**2)

	optim = leastsq( errfunc, p0[:], args=(x0,y0,yerr0), full_output=1 )

	return optim

def mcmcfit(x0,y0, yerr0=None,p0=None,expression=None, bounds=None, steps=None):
	from PyAstronomy import funcFit as fuf
	
	pnames = [ 'p%d' % k for k in range(len(p0)) ]
	guess0 = {}
	for k in range(len(p0)):
		guess0['p%d' % k] = p0[k]
	if yerr0 == None:
		yerr0 = numpy.asarray(y0)*0.1
	if p0 == None:
		return
	if expression == None:
		return

	class mymodel(fuf.OneDFit):
		def __init__(self):
			self.expression = expression
			self.pnames = pnames
			fuf.OneDFit.__init__(self, pnames)
		
		def func(self,p,x):
			return eval(expression)

		def evaluate(self, x):
			return self.func([self[k] for k in pnames],x)

	cMod = mymodel()
	for k in range(len(p0)):
		cMod['p%d'%k] = p0[k]

	if bounds != None:
		for k in bounds:
			cMod.setRestriction({k:bounds[k]})

	cMod.thaw(pnames)
	print '--direct Fit'
	cMod.fit(x0,y0)
	cMod.parameterSummary()

	X0 = cMod.parameters()
	print "-- MCMC sampling"
	_steps = {}
	for k in pnames:
		if steps != None:
			if k in steps:
				_steps[k]=steps[k]
		else:
			_steps[k] = 0.1
	Lims = {}
	for k in pnames:
		if bounds != None:
			if k in bounds:
				m, M = bounds[k]	
				if m != None:
					m = numpy.min([X0[k],m])
				if M != None:
					M = numpy.min([X0[k],M])
				Lims[k] = (m,M)
			else:
				Lims[k] = (None,None)
		else:
			Lims[k] = (None,None)

	ppa = {}
	print _steps
	cMod.fitMCMC(x0, y0, guess0, Lims, _steps, yerr=yerr0, \
		   pymcPars=ppa, iter=10000, burn=0, thin=1, \
		   dbfile="mcmcExample.tmp")
	
	db = pymc.database.pickle.load('mcmcExample.tmp')
	return db
