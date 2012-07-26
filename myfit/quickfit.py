""" Provide a quick and almost dirty fitting code
	Using scipy optimize and emcee
"""

import numpy
from scipy import optimize
import emcee
import mytables

default_emcee_cfg = dict( nwalkers=256, threads=4, burn=250, nsample=50, guess_sig=0.005)


def fit(lnp, p0, verbose=True, emcee_cfg=default_emcee_cfg, skip_scipy=False, full_output=False):
	""" Do the actual fit
		Run scipy leastsq first if the likelihood is based on leastsq in
		order to find a better place to start (if skip_scipy is set,
		this step is skipped) 
		Then run emcee to do the full analysis
	"""


	if hasattr(lnp, '__leastsq__') & (not skip_scipy):
		if verbose:
			print "Running scipy.optimize.leastsq"
		p1, success = optimize.leastsq(lnp.__leastsq__, p0[:])
		if verbose:
			for k in range( len(lnp.names) ):
				print "%10s, %g" %(lnp.names[k], p1[k])
	else:
		p1 = p0[:]
	
	if verbose:
		print "Running Emcee"
	sampler = emcee.EnsembleSampler(emcee_cfg['nwalkers'], len(lnp), lnp, threads=emcee_cfg['threads'])
	pos = [	lnp.getRandomStart(p1, emcee_cfg['guess_sig']) for k in range(emcee_cfg['nwalkers']) ]
	pos,prob,state = sampler.run_mcmc(pos, emcee_cfg['burn'])
	sampler.reset()
	sampler.run_mcmc(pos, emcee_cfg['nsample'], rstate0=state)

	p2 = [  numpy.median(sampler.flatchain[:,k]) for k in range(len(lnp)) ]
	
	d = { lnp.names[k]: sampler.flatchain[:,k] for k in range(len(lnp)) }
	d['lnp'] = sampler.lnprobability.flatten()
	t = mytables.Table(d)
	if verbose:
		stats = t.stats(val='q')
		for k in range(stats.nrows):
			nk = stats['Field'][k]
			m  = stats['q50'][k]
			mm = m - stats['q25'][k]
			mp = stats['q75'][k] - m
			print '%s: %g [+ %g / - %g ]' % (nk, m, mm, mp)

	if not full_output:
		return p0, p1, p2, t.stats(val='q')
	else:
		return p0, p1, p2, t.stats(val='q'), t, sampler
