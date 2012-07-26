import numpy 

class Base(object):
	def __init__(self, x, y, yerr=None, mod=None, priors = None):
		self.mod = mod
		if yerr != None:
			self.yerr = yerr 
		else:
			self.yerr = numpy.ones(len(y)) #or y[:]*0.1
		self.names = self.mod.names
		self.x = x
		self.y = y
		self.priors = priors
		self.__childInit__()

	def __childInit__(self):
		pass

	def __len__(self):
		return len(self.mod)

	def ym(self, theta):
		return self.mod(self.x, theta)
		
	def isInside(self, v, interval):
		""" Check if v within interval"""
		_interv = numpy.sort(interval)
		return ( _interv[0] <= v <= _interv[1] )
	
	def evalPriors(self, theta):
		""" evaluate if theta satisfies the priors """
		if self.priors != None:
			return [self.isInside(theta[k], self.priors[k]) for k in range(len(self))]
		else:
			return [True]

	def getRandomStart(self, guess, sig=0.005):
		return guess*(1.+ numpy.random.normal(0,sig,size=len(self)))


	def __lnp__(self, theta):
		""" Return lnp"""
		raise Exception("Not implemented")

	def __call__(self, theta):
		if False in self.evalPriors(theta):
			return -numpy.inf
		else:
			return self.__lnp__(theta)

class LeastSq(Base):
	def __childInit__(self):
		self._D0 = +len(self.x)*0.5*numpy.log(2*numpy.pi)
		self._D1 = 0.5*numpy.log(self.yerr**2).sum() 

	def __leastsq__(self, theta):
		""" Return a normal chi2 lnp"""
		return (self.y - self.ym(theta) )/self.yerr

	def __lnp__(self, theta):
		""" Return a normal chi2 lnp"""
		return -self._D0 - self._D1 -0.5 * (self.__leastsq__(theta)**2).sum()
