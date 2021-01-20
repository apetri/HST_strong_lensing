import os
import urllib.request

from bs4 import BeautifulSoup

import numpy as np
import pandas as pd
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
	_url_format = "https://archive.stsci.edu/pub/hlsp/frontier/abell2744/models/{method}/{version}/range/"
	_map_filename_format = "hlsp_frontier_model_abell2744_{method}-map{n:03d}_{version}_kappa.fits" 

	def __init__(self,method,version,metadata=True):
		self.method = method
		self.version = version
		if metadata:
			self.getMapMetadata()

	def path(self,n):
		return os.path.join(self._local_root,self.method,self.version,"range",self._map_filename_format.format(method=self.method,version=self.version,n=n))

	@property
	def datapath(self):
		return os.path.join(self._local_root,self.method,self.version)

	def getMapMetadata(self):
		with urllib.request.urlopen(self._url_format.format(method=self.method,version=self.version)) as resp:
			html = BeautifulSoup(resp.read(),"html.parser")

		map_files = []
		for tag in html.find_all("a"):
			fn = tag.attrs.get("href","")
			if fn.startswith("hlsp"):
				map_files.append(fn)

		df = pd.DataFrame({"map_filename":map_files})
		df["kind"] = df.map_filename.apply(lambda f:f.split(".fits")[0].split("_")[-1])
		df["realization"] = df.map_filename.apply(lambda f:int(f.split("map")[1].split("_")[0]))

		# Done
		self.maps_metadata = df
		return df

	@property
	def realizations(self):
		return self.maps_metadata[self.maps_metadata.kind=="kappa"].realization.values

	def isDownoaded(self,n):
		return os.path.exists(self.path(n))

	def download(self,n,overwrite=False):
		
		# Construct map url
		map_url =  self._url_format.format(method=self.method,version=self.version) + self._map_filename_format.format(method=self.method,version=self.version,n=n)
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

		# Fill NaN with 0
		data[np.isnan(data)] = 0.0

		# Cut map into a square
		nside = min(data.shape)
		data = data[:nside,:nside]

		# Return
		angle = (data.shape[0]*np.abs(head['CDELT1'])*u.deg).to(u.arcsec)
		return angle,data


	# Read one map into a ConvergenceMap instance
	def loadKappaMap(self,n,callback=None):
		fn = self.path(n)
		img = lt.image.convergence.ConvergenceMap.load(fn,format=self.loadFits)
		if callback is not None:
			return callback(img)
		else:
			return img
