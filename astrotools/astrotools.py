import sys
import urllib2,StringIO
import matplotlib.pyplot  as plt
import aplpy,pyfits
from .suds.client import Client
from lxml import etree

def get_dss(ra, dec, survey='all', radius = 15, debug=False):
	"""
	Get SDSS image at a given position
	INPUTS:
		ra, dec position in degrees
		radius is in arcmins
		survey: 'all', 'poss2ukstu', 'dss1'.. 
	OUTPUTS:
		returns an aplpy FITSfigure
	"""
	url='http://archive.stsci.edu/cgi-bin/dss_search?v=%s&r=%fd&d=%fd&e=J2000&h=%f&w=%f&f=fits&c=none&fov=NONE&v3='%(survey,ra,dec,radius,radius)
	
	response = urllib2.urlopen(url)
	if debug:
		print url
	html = response.read()
	f=StringIO.StringIO(html)
	
	dat=pyfits.open(f)
	plt.clf()
	gc=aplpy.FITSFigure(dat,figure=plt.gcf())
	#gc.set_tick_labels_format('ddd.dddd','ddd.dddd')
	gc.show_grayscale()
	return gc


def resolve(objectName):
	"""
	Resolve the object by name using CDS
	Inputs: 
		objectName	Name to resolve
	Outputs:
		return ra,dec in degrees
	Example:
	>> resolve.resolve('M31')
	(10.684708329999999, 41.268749999999997)

	Requires the following modules:
		suds, lxml
	"""
	url = 'http://cdsws.u-strasbg.fr/axis/services/Sesame?wsdl'
	client = Client(url)    
	xml=client.service.SesameXML(objectName)           
	tree = etree.fromstring(xml.encode('utf-8'))
	# take the first resolver
	pathRa = tree.xpath('/Sesame/Target/Resolver[1]/jradeg')
	pathDec = tree.xpath('/Sesame/Target/Resolver[1]/jdedeg')
	if len(pathRa)==0:
		return []
	ra=float(pathRa[0].text)
	dec=float(pathDec[0].text)
	return ra,dec
from numpy import deg2rad,rad2deg,sin,cos,sqrt,arcsin

def sphdist (ra1, dec1, ra2, dec2):
	"""measures the spherical distance between 2 points 
	Inputs:
		(ra1, dec1)	in degrees
		(ra2, dec2)	in degrees
	Outputs:
		returns a distance in degrees
	"""
	dec1_r = deg2rad(dec1)
	dec2_r = deg2rad(dec2)
	return 2 * rad2deg ( arcsin ( sqrt ( ( sin((dec1_r - dec2_r) / 2))**2 +  cos(dec1_r) * cos(dec2_r) * ( sin((deg2rad(ra1 - ra2)) / 2))**2)))	
