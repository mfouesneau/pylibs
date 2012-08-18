import aplpy
import numpy
from matplotlib.pylab import arrow

class FITSFigure(aplpy.FITSFigure):

    "A class for plotting FITS files."

    def show_compass(self, xw, yw, size=10.,color='red', 
    	head_width=20, **kwargs):
	
	if type(xw) == str:
		xw = self.hms2deg(xw)
	if type(yw) == str:
		yw = self.dms2deg(yw)

	r = xw + size/3600.*numpy.array([0,1,0])
    	d = yw + size/3600.*numpy.array([0,0,1])
	rp, dp = self.world2pixel (r,d)
	arrow(rp[0], dp[0],rp[1]-rp[0], dp[1]-dp[0],
		color=color, head_width=head_width, **kwargs)
	arrow(rp[0], dp[0],rp[2]-rp[0], dp[2]-dp[0],
		color=color, head_width=head_width, **kwargs)

    def hms2deg(self, str):
	if str[0] == '-': 
		neg = -1
		str = str[1:]
	else:
		neg = 1
	str = numpy.str.split(str,':')
	return neg*((((float(str[-1])/60.+float(str[1]))/60.+float(str[0]))/24.*360.))

    def dms2deg(self, str):
	if str[0] == '-': 
		neg = -1
		str = str[1:]
	else:
		neg = 1
	str = numpy.str.split(str,':')
	return (neg*((float(str[-1])/60.+float(str[1]))/60.+float(str[0])))

