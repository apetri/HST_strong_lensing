import os
import urllib.request

import numpy as np
from astropy.io import fits
import astropy.units as u

import lenstools as lt

models = {
	"cats":"v4.1",
	"williams":"v4",
	"bradac":"v2",
	"glafic":"v4"
	}

class Abell2744(object):

	_local_root = os.path.join('abell2744','models')
	_url_format = "https://archive.stsci.edu/pub/hlsp/frontier/abell2744/models/{method}/{version}/range/hlsp_frontier_model_abell2744_{method}-map{n:03d}_{version}_kappa.fits" 

	def __init__(self,method,version):
		self.method = method
		self.version = version

	def path(self,n):
		return os.path.join(self._local_root,self.method,self.version,"range","hlsp_frontier_model_abell2744_{method}-map{n:03d}_{version}_kappa.fits".format(method=self.method,version=self.version,n=n))

	def isDownoaded(self,n):
		return os.path.exists(self.path(n))

	def download(self,n,overwrite=False):
		
		# Construct map url
		map_url =  self._url_format.format(method=self.method,version=self.version,n=n)
		if not(overwrite) and self.isDownoaded(n):
			print("[+] Map file {0} already downloaded".format(map_url))
			return

		# Create directory structure
		pth = self.path(n)

		pth_parts = pth.split(os.path.sep)
		cur = pth_parts[0]
		for p in pth_parts[1:]:
			try:
				os.mkdir(cur)
			except OSError:
				pass
			finally:
				cur = os.path.join(cur,p)

		# Download
		print("[+] Downloading: "+map_url)
		urllib.request.urlretrieve(map_url,pth)

	###########################################################################

	# FITS loader
	@staticmethod
	def loadFits(fn):
		with fits.open(fn) as hdu:
			head = hdu[0].header
			data = hdu[0].data

		# Return
		angle = (data.shape[0]*np.abs(head['CDELT1'])*u.deg).to(u.arcsec)
		return angle,data


	# Read one map into a ConvergenceMap instance
	def loadKappaMap(self,n):
		fn = self.path(n)
		return lt.image.convergence.ConvergenceMap.load(fn,format=self.loadFits)
