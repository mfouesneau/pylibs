import figure
import numpy

def plotData(lnp, ax=None, **kwargs):
		""" Do the actual plot of the imput data """
		if ax == None:
			_ax = figure.gca()
		else:
			_ax = ax
		_ax.errorbar(lnp.x, lnp.y, yerr=lnp.yerr, **kwargs)

		if ax == None:
			_ax.set_xlabel('X')
			_ax.set_ylabel('Y')
			figure.theme(ax=_ax)

		figure.draw_if_interactive()

def plotFit(lnp, p0, p1, p2, ax=None, **kwargs):

	if ax == None:
		ax = figure.gca()

	mody = lnp.mod(lnp.x, p2)
	ax.plot(lnp.x, mody, label='Best fit', **kwargs)

	figure.draw_if_interactive()

def plotCI(lnp, sampler, ax = None, **kwargs):
	""" Plot confidence interval (fill_between) """
	if ax == None:
		_ax = figure.gca()
	else:
		_ax = ax

	for k in range(min(len(sampler.flatchain), 100)):
		mody = lnp.mod(lnp.x, sampler.flatchain[k,:])
		_ax.plot(lnp.x, mody, color='0.', ls='solid', alpha=0.2)
	
	figure.draw_if_interactive()

def plotResiduals(lnp, sampler, p0,p1,p2, ax=None, **kwargs):
	""" Plot the residuals of the fit """
	if ax == None:
		_ax = figure.gca()
	else:
		_ax = ax
	_ax.errorbar(lnp.x, lnp.y-lnp.mod(lnp.x, p2), yerr=lnp.yerr, **kwargs)
	if ax == None:
		_ax.set_xlabel('X')
		_ax.set_ylabel('Residuals [Data-Model]')
		figure.theme(ax=_ax)
	figure.draw_if_interactive()

def quickPlot(lnp, sampler, p0, p1, p2):
	""" Generates a figure showing the data, the best fit and 3 sigmas
	intervals (if possible), and the residuals """
	left = 0.125
	width = 1.-2*0.125
	rect1 = [left, 0.34, width, 1.-0.1-0.34]
	rect2 = [left, 0.1, width, 0.34-0.1]

	fig = figure.figure()
	ax0 = fig.add_axes(rect1)	
	ax1 = fig.add_axes(rect2, sharex=ax0)	

	plotData(lnp, ax=ax0, lw=2., label='Data', zorder=5, color='#ff0000')
	plotFit(lnp, p0,p1,p2, ax=ax0, zorder=0, lw=1, color='0.0')
	#plotCI(lnp, sampler, ax=ax0, zorder=-5, alpha=0.3)
	ax0.set_xlabel('X')
	ax0.set_ylabel('Y')
	l = ax0.legend(numpoints=1, scatterpoints=1)
	l.draw_frame(False)
	l.draggable(True)
	figure.setp(ax0.get_xticklabels(), visible=False)
	figure.theme(ax=ax0)
	ax0.yaxis.set_major_locator(figure.MaxNLocator(5, prune='lower'))

	plotResiduals(lnp, sampler, p0,p1,p2, ax=ax1)
	ax1.set_xlabel('X')
	ax1.set_ylabel('Residuals')
	figure.theme(ax=ax1)
	ax0.yaxis.set_major_locator(figure.MaxNLocator(4, prune='both'))

	return ax0, ax1

