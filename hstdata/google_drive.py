from __future__ import print_function

import io,os,pickle
import itertools

from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request

import numpy as np
import pandas as pd

from astropy.cosmology import WMAP9 as cosmo
from astropy.io import fits
import astropy.units as u
import astropy.constants as cnst

import lenstools as lt

#If modifying these scopes, delete the file token.pickle.
SCOPES = ['https://www.googleapis.com/auth/drive']

# Credentials
def getCredentials():

	creds = None
	# The file token.pickle stores the user's access and refresh tokens, and is
	# created automatically when the authorization flow completes for the first
	# time.
	if os.path.exists('token.pickle'):
		with open('token.pickle', 'rb') as token:
			creds = pickle.load(token)
			return creds

	# If there are no (valid) credentials available, let the user log in.
	if not creds or not creds.valid:
		if creds and creds.expired and creds.refresh_token:
			creds.refresh(Request())
		else:
			flow = InstalledAppFlow.from_client_secrets_file(
				'credentials.json', SCOPES)
			creds = flow.run_local_server(port=0)
		# Save the credentials for the next run
		with open('token.pickle', 'wb') as token:
			pickle.dump(creds, token)

		# Return credentials to caller
		return creds

# Handle to service
def getServiceHandle(creds):
	return build('drive', 'v2', credentials=creds)

########################################################################
########################################################################

class SimulatedMapsManager(object):

	# Constructor
	def __init__(self,local_root,maps_root_id='1ZjN46EvDirPeBVv3Ucsfn289VUiPh2e6'):
		self.local_root = local_root
		self.maps_root_id = maps_root_id
		if os.path.exists(maps_root_id+".csv"):
			self.maps_metadata = pd.read_csv(maps_root_id+".csv")
		else:
			self.maps_metadata = None

	# Available realizations
	@property
	def realizations(self):
		return self.maps_metadata.realization.unique()

	# Path to a specific realization/redshift/projection
	def path(self,r,z,p):
		return os.path.join(self.local_root,"D"+str(r),"map_{0:03d}_{1}_{2}_sph.fits".format(int(100*(z+1e-16)),r,p))

	# Check if map is downloaded
	def isDownloaded(self,r,z,p):
		return os.path.exists(self.path(r,z,p))

	# Download map
	def download(self,srv,r,z,p,overwrite=False):
		if self.maps_metadata is None:
			raise OSError("No maps metadata present, need to generate that first")

		fn = os.path.basename(self.path(r,z,p))
		if not(overwrite) and self.isDownloaded(r,z,p):
			print("[+] Map file {0} already downloaded".format(fn))
			return

		# Find google drive ID of the map
		try:
			fid = self.maps_metadata[self.maps_metadata.map_filename==fn].map_id.values[0]
		except IndexError:
			raise ValueError("Map file {0} does not exist on drive".format(fn))

		# Download
		print("[+] Downloading map: {0}, id={1}".format(fn,fid))
		buf = self.downloadFile(srv,fid)

		# Write the map to disk
		try:
			os.mkdir(self.local_root)
		except OSError:
			pass

		try:
			os.mkdir(os.path.join(self.local_root,'D'+str(r)))
		except OSError:
			pass

		print("[+] Saving downloaded map to: "+self.path(r,z,p))
		with open(self.path(r,z,p),"wb") as fp:
			buf.seek(0)
			fp.write(buf.read())

	# Download all maps in range
	def downloadRange(self,srv,r=[],z=[],p=[]):
		for c in itertools.product(r,z,p):
			ri,zi,pi = c
			self.download(srv,ri,zi,pi)

	# Build maps directory/name/metadata dataset
	def getMapMetadata(self,srv):

		if self.maps_metadata is not None:
			return self.maps_metadata

		df_rows = []

		# Realizations
		realizations = srv.children().list(folderId=self.maps_root_id).execute().get('items',[])
		for r in realizations:
			sub = srv.files().get(fileId=r['id']).execute()
			if not sub['title'].startswith('D'):
				break

			print("[+] Found realization: "+sub['title'])
			maps = srv.children().list(folderId=r['id']).execute().get('items',[])
			for mi in maps:
				m = srv.files().get(fileId=mi['id']).execute()
				if not m['title'].startswith('map_'):
					break

				print("[+][{0}] Found map file: {1}".format(sub['title'],m['title']))
				df_rows.append([sub['id'],sub['title'],m['id'],m['title']])

		# Parse into DataFrame
		df = pd.DataFrame.from_records(df_rows,columns=['realization_id','realization_name','map_id','map_filename'])

		# Compute additional columns
		df['redshift'] = df.apply(lambda x:float(x['map_filename'].split("_")[1])/100,axis=1)
		df['realization'] = df.apply(lambda x:int(x['map_filename'].split("_")[2]),axis=1)
		df['projection'] = df.apply(lambda x:int(x['map_filename'].split("_")[3]),axis=1)

		# Done
		self.maps_metadata = df.sort_values(["realization","redshift","projection"]).reset_index().drop('index',axis=1)
		return self.maps_metadata

	# Download file with a specific ID to a BytesIO object
	@staticmethod
	def downloadFile(srv,file_id):
		request = srv.files().get_media(fileId=file_id)
		fh = io.BytesIO()
		downloader = MediaIoBaseDownload(fh, request)
		done = False
		while done is False:
			status, done = downloader.next_chunk()
			print("Download %d%%." % int(status.progress() * 100))
		return fh

	###########################################################################

	# FITS loader
	@staticmethod
	def loadFits(fn):
		with fits.open(fn) as hdu:
			head = hdu[0].header
			data = hdu[0].data

		# Compute critical density scaling
		crit = 0.35/cosmo.angular_diameter_distance(head['ZL']).to(u.Gpc).value
		data = (data*cnst.M_sun/(u.kpc**2)).to(u.g/(u.cm**2)).value
		data = data/crit

		# Return
		angle = (data.shape[0]*np.abs(head['CDELT1'])*u.deg).to(u.arcsec)
		return angle,data


	# Read one map into a ConvergenceMap instance
	def loadKappaMap(self,r,z,p):
		fn = self.path(r,z,p)
		return lt.image.convergence.ConvergenceMap.load(fn,format=self.loadFits)

