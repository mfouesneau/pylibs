
def example():
	""" Usage example of this package 
		It consider data following y=6-4*x**2, with a centered Gaussian
		noise of 1. as dispersion.

		This example considers a Quadric model to fit the data
		And plots some results

	"""
	from .models import Quadric
	from .likelihoods import LeastSq
	from .quickfit import fit
	from .plot import quickPlot
	import numpy
	print example.__doc__

	_model = lambda x, p0,p1: p0+p1*x**2 
	def _data(p0=6., p1=-4, noise=1., npts=5.):
		x = numpy.linspace(0.,5., npts)
		y = _model(x,p0,p1)+numpy.random.normal(0,noise, len(x))
		yerr = noise*numpy.ones_like(x)
		return x, y, yerr

	emcee_cfg = dict( nwalkers=16, threads=4, burn=250, nsample=50, guess_sig=0.005)

	x,y,yerr   = _data(npts=10)	
	lnp        = LeastSq( x,y,yerr, mod = Quadric() )
	lnp.priors = [ (-5, 20), (-10, 10), (-10,10) ]
	p0         = [ 3., 2., 1. ]
	p0, p1, p2, stats, t, s = fit(lnp, p0, verbose=True, emcee_cfg=emcee_cfg, skip_scipy=False, full_output=True)
	quickPlot(lnp, s, p0, p1, p2)

