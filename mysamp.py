"""
	SAMPy Client package
	
	This class aimed at providing a very small VO interactivity using the
	SAMP protocol. This allows anyone to easily send and receive data to VO
	applications such as Aladin or Topcat.

	It provides 2 Classes:
		
		* Hub    -- Samp hub that is required to manage the communications
			    between all the VO applications

		* Client -- Python object that is a proxy to send and receive
			    data from/to applications
	REQUIREMENTS:
		mytables	to manage any table types at least FITS
			 	(this depends on pyfits)
		sampy		samp access
	
	see the function demo() for a usage example. (demo??)
"""
__version__ = '1.0'
import atexit, os, urllib
import sampy
import mytables

#==========================================================================
# HUB -- Generate a SAMP hub to manage communications
#==========================================================================
class Hub(object):
    """ 
    This class is a very minimalistic class that provides a working SAMP hub
    Usage:
	> h = Hub()
	> h.stop()
    """
    def __init__(self,addr, *args, **kwargs):
        self.SAMPHubServer=sampy.SAMPHubServer(addr=addr, *args, **kwargs)
        atexit.register(self.SAMPHubServer.stop)

    def __del__(self):
        self.SAMPHubServer.stop()

#==========================================================================
# Client -- Generate a SAMP client able to send and receive data
#==========================================================================
class Client(object):
    """
	This class implenent an interface to SAMP applications like Topcat and
	Aladin using Sampy module.
	It allows you to exchange tables with topcat using the SAMP system.

	To instanciate a connection:
	> client = Client(addr='localhost', hub=True)
	This could create a local hub to connect to a local session
	Nb: Topcat need to be running when sending messages, however this
		could be done later.

	To send a table to Topcat:
	> client['tblName'] = { <table content> }

	To receive a table from Topcat:
	(broadcast the table from Topcat)
	> table = client['tblName']

    """
    #Destructor ===============================================================
    def __del__(self):
    	self.client.disconnect()
	if self.hub != None: 
		self.hub.stop()

    #Constructor ==============================================================
    def __init__(self, addr='localhost', hub=True):
        # Before we start, let's kill off any zombies
        
        if hub:
        # sampy seems to fall over sometimes if 'localhost' isn't specified,
	# even though it shouldn't
		self.hub = sampy.SAMPHubServer(addr=addr)
        	self.hub.start()
	else: 
		self.hub = None
        
        self.metadata = {
                            'samp.name' : 'Python Session',
                            'samp.icon.url' :'http://docs.scipy.org/doc/_static/scipyshiny_small.png',
                            'samp.description.text' : 'Python Samp Module',
                            'client.version' : '0.1a'
                        }
        
        self.client = sampy.SAMPIntegratedClient(metadata=self.metadata, addr=addr)
        self.client.connect()
        atexit.register(self.client.disconnect)
        
	# Bind interaction functions - we will register that we want to listen
	# for table.highlight.row (the SAMP highlight protocol), and all the
	# typical SAMP-stuff We could include this to listen for other sorts of
	# subscriptions like pointAt(), etc.
        self.client.bindReceiveNotification('table.highlight.row', self.highlightRow)
        self.client.bindReceiveNotification('table.select.rowList', self.highlightRow)
        self.client.bindReceiveNotification('table.load.fits', self.receiveNotification)
        self.client.bindReceiveNotification('samp.app.*', self.receiveNotification)

        self.client.bindReceiveCall('samp.app.*', self.receiveCall)
        self.client.bindReceiveCall('table.load.fits', self.receiveCall)

        self.client.bindReceiveNotification('image.load.fits', self.receiveCall)
        self.client.bindReceiveCall('image.load.fits', self.receiveCall)
        
        self.tables= {}
        self.images= {}
	self.lastMessage = None
        

    # Generic response protocols ===============================================
    def receiveNotification(self, privateKey, senderID, mType, params, extra):
        print '[SAMP] Notification ', privateKey, senderID, mType, params, extra
       	self.lastMessage = {'label':'Notification',
			 'privateKey':privateKey, 
			 'senderID': senderID, 
			 'mType': mType, 
			 'params': params, 
			 'extra': extra }
	print mType,  mType == 'image.load.fits'
	if mType == 'image.load.fits':
		self.receiveCall(privateKey, senderID, msgID, mType, params, extra)
		
    def receiveResponse(self, privateKey, senderID, msgID, response):
        print '[SAMP] Response ', privateKey, senderID, msgID, response
        
    def receiveCall(self, privateKey, senderID, msgID, mType, params, extra):
	print 'CAlling'
	data = { 'privateKey':privateKey, 'senderID':senderID, 'msgID':msgID, 'mType':mType, 'params':params}
        print '[SAMP] Call'
	if mType == 'table.load.fits':
		print "[SAMP] Table received."
		self.tables[params['name']]=data['params']
		self.tables[params['name']]['data'] = None
	elif mType == 'image.load.fits':
		print "[SAMP] Image received."
		self.images[params['name']]=data['params']
		self.images[params['name']]['data'] = None
	else:
        	print '[SAMP] Call'
		print mType
	self.client.ereply(msgID, sampy.SAMP_STATUS_OK, result = {'txt' : 'printed' })

    def __getitem__(self, k):
	return self.get(k)	

    def __setitem__(self, k, data):
	""" Broadcast data to all """
	return self.send(k, data, to='all')

    def __call__(self):
    	""" print detailed info about the current client """
	self.info()	
	neighbours = self.getNeighbours()
	if len(neighbours) > 0:
		print "%d detected client(s):" % len(neighbours)
		for n in neighbours:
			print '\t'+self.client.getMetadata(n)['samp.name']
	
	print "Registered tables: %d" % len(self.tables)
	for k in self.tables.keys(): 
		print '      %s' % k

    # Application functions ===================================================
    def get(self, k):
    	""" get table """
	if k in self.tables:
		cTab = self.tables[k]
		if cTab['data'] == None:
			u = urllib.urlretrieve(cTab['url'])
			cTab['data'] = data  = mytables.load(u[0])
			return data
		else:	
			return cTab['data']
	if k in self.images:
		cTab = self.images[k]
		if cTab['data'] == None:
			cTab['data'] = data  = mytables.pyfits.getdata(cTab['image-id'])
			return data
		else:	
			return cTab['data']
	
    def send(self, k, data, to='All'):
	""" Broadcast data to one or all applications """
	
	if k[-5:]!='.fits':
		k+='.fits'

	if isinstance(data, mytables.Table):
		data.write(k, append=False, clobber=True)
	else:
    		mytables.Table(data).write(k, append=False, clobber=True)

	self.tables[k] = { 'name':k, 
		    	   'url': 'file://' + os.getcwd() + '/' + k,
			   'data': mytables.load(k) 
			   }
	return self._broadcastTable(k)

    def info(self):
    	""" print information about the current client """
	for k in self.metadata.keys():
		print k, self.metadata[k]
    
    def getNeighbours(self):
    	""" returns print information about the current client """
	return self.client.getRegisteredClients()
   
       
    def getAppId(self, name):
	"""
	Returns the registered Id of a given application name
	"""
        neighbours = self.client.getRegisteredClients()
        for neighbour in neighbours:
            	metadata = self.client.getMetadata(neighbour)
		try:
			if (metadata['samp.name'] == name):
				return neighbour
		except KeyError:
			continue

    def isAppRunning(self, name):
    	""" 
	Check if a given application is running and registered to the SAMP
	server
	"""
    
        neighbours = self.client.getRegisteredClients()
        
        for neighbour in neighbours:
            metadata = self.client.getMetadata(neighbour)
            
            try:
                if (metadata['samp.name'] == name):
                    self.topcat = neighbour
                    return True
                
            except KeyError:
                continue

        return False
        
    # Broadcast a table file to TOPCAT
    def _broadcastTable(self, table, to='All'):
   	""" Broadcast message to given recipient """ 
        metadata = {
                    'samp.mtype' : 'table.load.fits',
                    'samp.params' : {
                                     'name' : table,
                                     'table-id' : table,
                                     'url' : 'file://' + os.getcwd() + '/' + table
                                    }
                   }
                   
	if len(self.client.getRegisteredClients()) > 0:
		if to.lower() == 'all':
			return self.client.notifyAll(metadata)
		else:
        		if self.isAppRunning(name): 
				return self.client.notify(self.getAppId(to), metadata)
        else: 
		return False
    
    def highlightRow(self, privateKey, senderID, mType, params, extra):
    	""" Not working yet... Only receiving the selected row """
        print '[SAMP] Highlighted row', privateKey, senderID, mType, params, extra
       	self.lastMessage = {'label':'Highlighted Row',
			 'privateKey':privateKey, 
			 'senderID': senderID, 
			 'mType': mType, 
			 'params': params, 
			 'extra': extra }

        if(senderID == self.topcat):
           
            try:
                filename, row = [params['url'], int(params['row'])]
            except KeyError:
                print '[SAMP] Highlighted row was missing vital information'
            else:
                print 'TOPCAT tells us that row %s of file %s was highlighted!' % (row, filename)
        
#==========================================================================
#==========================================================================

def demo(client=None):
	""" 
	This function is a short example of how to use this package

	It creates a minimal table. 
	However it is prefered to use the mytables module to do so, in order to
	provide more data descriptions such as units or comments.
	"""
	import numpy
	# first start a topcat session. (Topcat application can be launched after)
	if client == None:
		c  = Client()
	else:
		c = client

	# Some data is generated from my program. Let's use some mock data.
	x = numpy.arange(0, 100)
	y = x**2

	# broadcast a table to topcat
	c['t0'] = {'x':x, 'y':y }

	print """	
	 to obtain a table from topcat,
	 broadcast a message from topcat and then
	     mytable = topcat['mytable'] """
	if client == None:	
		return c
