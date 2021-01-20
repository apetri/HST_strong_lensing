import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

def plotPowerSpectrum(ell,p,label="default",theta_ticks=np.array([1.5,2,3,4,5,6])*u.arcsec,ax=None,figsize=(12,8)):

	# Init axes
	if ax is None:
		init_ax = True
		fig,ax = plt.subplots(figsize=figsize)
	else:
		fig = None
		init_ax = False
	
	# Plot
	ax.plot(ell,ell*ell*p/(2*np.pi),label=label)

	if init_ax:
	
		# Overlay theta
		ax1 = ax.twiny()
		ax1.set_xticks(2*np.pi/theta_ticks.to(u.rad).value)
		ax1.set_xlim(ax.get_xlim())
		ax1.set_xticklabels([str(x.value) for x in theta_ticks])

		# Legends
		ax.set_xlabel(r'$\ell$',fontsize=18)
		ax1.set_xlabel(r'$\theta$ (arcsec)',fontsize=18)
		ax.set_ylabel(r'$\ell^2P_{\kappa\kappa}(\ell)/2\pi$',fontsize=18)
	
	# Return
	return fig,ax