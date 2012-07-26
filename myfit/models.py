import numpy
class Model(object):
	"""defines the general scheme of a model"""
	def __init__(self, *args, **kwargs):
		self.name = 'General Model'
		self.len  = 0
		self.names = None
	def __len__(self):
		return self.len
	def __call__(self, x, *args, **kwargs):
		pass
	def __repr__(self):
		return self.name+'\t'+self.expr+'\n'+object.__repr__(self)
	def __str__(self):
		return self.name+'\t'+self.expr
	def __add__(self, a):
		assert(issubclass(a,Model))
		m      = Model()
		m.name = self.name + ' + ' + a.name
		m.expr = self.expr + ' + ' + a.expr
		m.len  = len(self) + len(a)
		m.__call__ = lambda x, p: self.__call__(x, p[:len(self)]) + a.__call__(x, p[len(self):])
		return m
	def __sub__(self, a):
		assert(issubclass(a,Model))
		m      = Model()
		m.name = self.name + ' - ' + a.name
		m.expr = self.expr + ' - ' + a.expr
		m.len  = len(self) + len(a)
		m.__call__ = lambda x, p: self.__call__(x, p[:len(self)]) - a.__call__(x, p[len(self):])
		return m
	def __mul__(self, a):
		assert(issubclass(a,Model))
		m      = Model()
		m.name = self.name + ' * ' + a.name
		m.expr = self.expr + ' * ' + a.expr
		m.len  = len(self) + len(a)
		m.__call__ = lambda x, p: self.__call__(x, p[:len(self)]) * a.__call__(x, p[len(self):])
		return m
	def __div__(self, a):
		assert(issubclass(a,Model))
		m      = Model()
		m.name = self.name + ' / ' + a.name
		m.expr = self.expr + ' / ' + a.expr
		m.len  = len(self) + len(a)
		m.__call__ = lambda x, p: self.__call__(x, p[:len(self)]) / a.__call__(x, p[len(self):])
		return m
		
class Constant(Model):
	""" Constant model, lambda x: p0 """
	def __init__(self, *args, **kwargs):
		self.name='Constant'
		self.expr='lambda x: p0'
		self.len =1 
	@property
	def names(self):
		return [ 'p%d' %k for k in range(len(self)) ]
	def __call__(self, x, p):
		return p

class Linear(Model):
	""" Linear relation, lambda x: p0 + p1*x"""
	def __init__(self, *args, **kwargs):
		self.name='Linear'
		self.expr='lambda x: p0 + p1*x'
		self.len = 2 
	@property
	def names(self):
		return [ 'p%d' %k for k in range(len(self)) ]
	def __call__(self, x, p):
		return p[0]+p[1]*x

class Quadric(Model):
	""" 2nd order polynom, lambda x: p0 + p1*x + p2*x**2"""
	def __init__(self, *args, **kwargs):
		self.name='Quadric'
		self.expr='lambda x: p0 + p1*x + p2*x**2'
		self.len = 3 
	@property
	def names(self):
		return [ 'p%d' %k for k in range(len(self)) ]
	def __call__(self, x, p):
		return p[0] + p[1]*x + p[2]*x**2

class Polynomial(Model):
	"""
	Creates a polynomial model of the form
	a0 + a1*x + a2*x**2 + ... + an*x**n
	where the order parameter is n
	"""
	def __init__(self, order, *args, **kwargs):
		self.order = order
		self.name  = 'Polynomial'
		self.expr  = 'lambda x:'
		self.len   = order + 1  
		for o in range(order):
			param = 'p%d' % o
			if o > 0:
				self.expr += ' + x**%d' %o
			else:
				self.expr += ' p0'
			def base(x, p, order=o):
				return p*(x**order)
			self.basis[param] = base

		self.params=self.basis.keys()
	@property
	def names(self):
		return [ 'p%d' %k for k in range(len(self)) ]
	def __call__(self, x, *params, **keyword_params):

		for key, val in zip(self.params, params):
			keyword_params[key] = val

		sum = 0.
		for key, val in keyword_params.iteritems():
			sum += self.basis[key](x, val)

		return sum

class Gaussian(Model):
	""" Gaussian-law, lambda x: p0*exp( -0.5 * [(x-p1)/p2]**2 ) """
	def __init__(self, *args, **kwargs):
		self.name = 'Gaussian'
		self.expr = 'lambda x: p0*exp( -0.5 * [(x-p1)/p2]**2 )'
		self.len  = 3 
		self.names=['mean', 'sigma']
	def __call__(self, x, p):
		z = (x - p[1]) / p[2]
		return p[0]*numpy.exp(-0.5*z**2)

class PowerLaw(Model):
        """ power-law relation, 
		lambda M: N * M**alpha
	"""
        def __init__(self, *args, **kwargs):
                self.name='Power law'
                self.expr='lambda x: p0 * (x)**p1'
		self.names=['N', 'alpha']

	def __len__(self):
		return len(self.names)

        def __call__(self, x, p):
                return p[0] * x**p[1]

class logPowerLaw(Model):
        """ power-law relation, in log-space 
		lambda M: log10( p[0] +  p[1]*log(x) )
	"""
        def __init__(self, *args, **kwargs):
                self.name='Power law'
                self.expr='lambda x: p0 * (x)**p1'
		self.names=['lnN', 'alpha']

	def __len__(self):
		return len(self.names)

        def __call__(self, x, p):
                return numpy.log10(numpy.exp( p[0] + p[1] * numpy.log(x) ))


class logSchechter(Model):
        """ Shechter relation, 
		lambda M: N * (M/Mc)**alpha * exp (-M/Mc)
	"""
        def __init__(self, *args, **kwargs):
                self.name='log space Schechter (Schechter 1974)'
                self.expr='lambda x: logN + alpha * log(x/Mc) - x/Mc'
		self.names=['lnMc', 'alpha', 'lnN']

	def __len__(self):
		return len(self.names)
        
	def __call__(self, x, p):
                return numpy.log10(numpy.exp(p[2] + p[1] * numpy.log(x/10**p[0]) - x/10**p[0]))
