#########################
#########################
# Need to account for limit in input period
#########################
#########################

# Colossus GC analyse script -- NO crowding
# White Dwarf (WD) version of analyse script
# New script copied from quest - want to take p and ecc from each population (all, obs, rec) and put them into separate file
# Doing this so we don't have to run analyse each time

import pandas as pd
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units, constants
from astropy.modeling import models, fitting
import scipy.stats
from scipy.integrate import quad

#for Quest
import matplotlib
matplotlib.use('Agg')

doIndividualPlots = True

from matplotlib import pyplot as plt

def file_len(fname): 
	i = 0
	with open(fname, encoding = 'us-ascii') as f:
		for i, l in enumerate(f):
			# print('index, line: ', i,l)
			pass
	return i + 1

def getPhs(sigma, m1=1*units.solMass, m2=1*units.solMass, m3=0.5*units.solMass):
	Phs = np.pi*constants.G/np.sqrt(2.)*(m1*m2/m3)**(3./2.)*(m1 + m2)**(-0.5)*sigma**(-3.)
	return Phs.decompose().to(units.day)

#similar to field, but limiting by the hard-soft boundary
def fitRagfb():
	x = [0.05, 0.1, 1, 8, 15]  #estimates of midpoints in bins, and using this: https://sites.uni.edu/morgans/astro/course/Notes/section2/spectralmasses.html
	y = [0.20, 0.35, 0.50, 0.70, 0.75]
	init = models.PowerLaw1D(amplitude=0.5, x_0=1, alpha=-1.)
	fitter = fitting.LevMarLSQFitter()
	fit = fitter(init, x, y)

	return fit

def RagNormal(x, cdf = False):
	mean = 5.03
	std = 2.28
	if (cdf):
		return scipy.stats.norm.cdf(x,mean,std)

	return scipy.stats.norm.pdf(x,mean,std)


def saveHist(histAll, histObs, histRec, bin_edges, xtitle, fname, filters = ['u_', 'g_', 'r_', 'i_', 'z_', 'y_','all']):
	c1 = '#5687A6' #Dali Blue (Andrew's AAS Poster)
	c2 = '#A62B1F' #Dai Red 
	c3 = '#BF8A26' #Dali Beige
	fig,ax1 = plt.subplots(figsize=(8,6), sharex=True)#can change to include cdf with ax1, ax2

	histAll = np.insert(histAll,0,0)
	histObs = np.insert(histObs,0,0)
	for f in filters:
		histRec[f] = np.insert(histRec[f],0,0)

	#PDF
	ax1.step(bin_edges, histAll/np.sum(histAll), color=c1)
	ax1.step(bin_edges, histObs/np.sum(histObs), color=c2)
	for f in filters:
		lw = 1
		if (f == 'all'):
			lw = 0.5
	ax1.step(bin_edges, histRec[f]/np.sum(histRec[f]), color=c3, linewidth=lw)
	ax1.set_ylabel('PDF')
	ax1.set_yscale('log')
	ax1.set_title('Globular Clusters - Colossus', fontsize = 16)
	ax1.set_xlabel(xtitle)
	#CDF
	#cdfAll = []
	#cdfObs = []
	#cdfRec = dict()
	#for f in filters:
#		cdfRec[f] = []

#	for i in range(len(histAll)):
#		cdfAll.append(np.sum(histAll[:i])/np.sum(histAll))
#	for i in range(len(histObs)):
#		cdfObs.append(np.sum(histObs[:i])/np.sum(histObs))
#	for f in filters:
	#	for i in range(len(histRec[f])):
	#		cdfRec[f].append(np.sum(histRec[f][:i])/np.sum(histRec[f]))
	#ax2.step(bin_edges, cdfAll, color=c1)
	#ax2.step(bin_edges, cdfObs, color=c2)
	#for f in filters:
	#	lw = 1
	#	if (f == 'all'):
	#		lw = 0.5
	#	ax2.step(bin_edges, cdfRec[f], color=c3, linewidth=lw)
	#ax2.set_ylabel('CDF')

	#ax2.set_xlabel(xtitle)
	fig.subplots_adjust(hspace=0)
	fig.savefig('./plots/'  + fname+'.pdf',format='pdf', bbox_inches = 'tight')

	#write to a text file
	with open('./eblsst_files/' + fname+'.csv','w') as fl:
		outline = 'binEdges,histAll,histObs'
		for f in filters:
			outline += ','+f+'histRec'
		outline += '\n'
		fl.write(outline)
		for i in range(len(bin_edges)):
			outline = str(bin_edges[i])+','+str(histAll[i])+','+str(histObs[i])
			for f in filters:
				outline += ','+str(histRec[f][i])
			outline += '\n'
			fl.write(outline)

# ######################################## Testing functions ########################################
def writeCornerFiles(df, binParams, population, clusterType, strategy, crowding):
	'''
	Function to write certain binary parameters to csv data files - will feed white dwarf dataframes (loc statement on all, obs, rec dfs) 
	df - input dataframe
	binParams - listlike, needs to be strings included in df.columns
	Below are naming parameters for output hist files
	population : either all, obs, rec
	clustertype : Globular (G), Open (O), M10, M67
	strategy : baseline (B) or colossus (C)
	crowding : either no-crowding scenario (N) or with crowding (C)
	'''

	newDF = pd.DataFrame(columns = binParams)
	param_list = []
	for param in binParams:
		param_list.append(df[param].values)
		params = np.concatenate(param_list)


		newDF[param] = params#_list
	print(newDF)

	# newDF.to_csv(f'./data/{population}-{clusterType}{strategy}{crowding}-histData.csv')
	# print('dataframes 2 write: ', newDF)

	# return newDF

def wdMRrelation(mass):
	'''
	Function for mass-radius relationship for MS stars, 
	WDs will have radii smaller than the predicted R from this relation
	Source: http://personal.psu.edu/rbc3/A534/lec18.pdf
	'''
	xi = 0
	for m in mass:
		if (m < 1.0):
			xi = 0.8
			radius = m ** xi
			# print('radius: ', radius)
			return radius

		elif m >= 1.0:
			xi = 0.57
			radius = m ** xi
			# will not have negative radii values
			return 0.0

#####################################################################################################


if __name__ == "__main__":

	filters = ['u_', 'g_', 'r_', 'i_', 'z_', 'y_', 'all']

	#get the Raghavan binary fraction fit
	fbFit= fitRagfb()
	print(fbFit)
		
	#to normalize
	intAll, err = quad(RagNormal, -20, 20)
	intCut, err = quad(RagNormal, -20, np.log10(365*10.))
	intNorm = intCut/intAll

	#cutoff in percent error for "recovered"
	Pcut = 0.1

	#assumed mean stellar mass
	mMean = 0.5

	#minimum number of lines to consider in file
	Nlim = 3

	if (doIndividualPlots):
		fmass, axmass = plt.subplots()
		fqrat, axqrat = plt.subplots()
		fecc, axecc = plt.subplots()
		flper, axlper = plt.subplots()
		fdist, axdist = plt.subplots()
		fmag, axmag = plt.subplots()
		frad, axrad = plt.subplots()

	#bins for all the histograms
	Nbins = 25
	mbins = np.arange(0,10, 0.1, dtype='float')
	qbins = np.arange(0,1, 0.1, dtype='float')
	ebins = np.arange(0, 1.05, 0.05, dtype='float')
	lpbins = np.arange(-2, 10, 0.5, dtype='float')
	dbins = np.arange(0, 40, 1, dtype='float')
	magbins = np.arange(11, 25, 1, dtype='float')
	rbins = np.arange(0, 100, 0.2, dtype='float')

	#blanks for the histograms
	#All
	m1hAll = np.zeros_like(mbins)[1:]
	qhAll = np.zeros_like(qbins)[1:]
	ehAll = np.zeros_like(ebins)[1:]
	lphAll = np.zeros_like(lpbins)[1:]
	dhAll = np.zeros_like(dbins)[1:]
	maghAll = np.zeros_like(magbins)[1:]
	rhAll = np.zeros_like(rbins)[1:]
	#Observable
	m1hObs = np.zeros_like(mbins)[1:]
	qhObs = np.zeros_like(qbins)[1:]
	ehObs = np.zeros_like(ebins)[1:]
	lphObs = np.zeros_like(lpbins)[1:]
	dhObs = np.zeros_like(dbins)[1:]
	maghObs = np.zeros_like(magbins)[1:]
	rhObs = np.zeros_like(rbins)[1:]
	#Recovered
	m1hRec = dict()
	qhRec = dict()
	ehRec = dict()
	lphRec = dict()
	dhRec = dict()
	maghRec = dict()
	rhRec = dict()
	for f in filters:
		m1hRec[f] = np.zeros_like(mbins)[1:]
		qhRec[f] = np.zeros_like(qbins)[1:]
		ehRec[f] = np.zeros_like(ebins)[1:]
		lphRec[f] = np.zeros_like(lpbins)[1:]
		dhRec[f] = np.zeros_like(dbins)[1:]
		maghRec[f] = np.zeros_like(magbins)[1:]
		rhRec[f] = np.zeros_like(rbins)[1:]

	RA = []
	Dec = []
	recFrac = []
	recN = []
	rawN = []
	obsN = []
	fileN = []
	fileObsN = []
	fileRecN = []

	allNPrsa = []
	obsNPrsa = []
	recNPrsa = []


	# Lists for different params
	eccAll = []
	eccObs = []
	eccRec = []

	pAll = []
	pObs = []
	pRec = []

	iAll = []
	iObs = []
	iRec = []

	m1All = []
	m1Obs = []
	m1Rec = []

	m2All = []
	m2Obs = []
	m2Rec = []

	r1All = []
	r1Obs = []
	r1Rec = []	

	r2All = []
	r2Obs = []
	r2Rec = []

	magAll = []
	magObs = []
	magRec = []

	# Lists to append WD dataframes to, will concatenate after looping thru files
	all_WD = []
	obs_WD = []
	rec_WD = []
	# Using prsa dataframes for these lists because of period cutoff at 1000 days

	# Dataframes to write to files later; 3 files for each sub-population - append everything to these
	All = pd.DataFrame(columns = ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r'])
	Obs = pd.DataFrame(columns = ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r'])
	Rec = pd.DataFrame(columns = ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r'])

	#Read in all the data and make the histograms
	d = "/projects/p30137/ageller/testing/EBLSST/clusters/GlobularClusters/colossus/output_files/"
	files = os.listdir(d)
	IDs = []
	for i, f in enumerate(files):
		print(round(i/len(files),4), f)
		fl = file_len(d+f)
		if (fl >= 4):
			#read in the header
			header = pd.read_csv(d+f, nrows=1)
	######################
	#NEED TO ACCOUNT FOR THE BINARY FRACTION when combining histograms
	#####################

			Nmult = header['clusterMass'][0]/mMean
			#Nmult = 1.

			RA.append(header['OpSimRA'])
			Dec.append(header['OpSimDec'])

			#read in rest of the file
			data = pd.read_csv(d+f, header = 2).fillna(-999)
			rF = 0.
			rN = 0.
			Nrec = 0.
			Nobs = 0.
			raN = 0.
			obN = 0.
			fiN = 0.
			fioN = 0.
			firN = 0.
			NallPrsa = 0.
			NobsPrsa = 0.
			NrecPrsa = 0.
			Nall = len(data.index)/intNorm ###is this correct? (and the only place I need to normalize?)
			prsa = data.loc[(data['p'] < 1000) & (data['p'] > 0.5)]

			# Selecting only WD candidates
			# print('foo: ', data.loc[(data['r1'] < wdMRrelation(data['m1'])) | (data['r2'] < wdMRrelation(data['m2']))]  )
			prsaWD = data.loc[(data['p'] < 1000) & (data['p'] > 0.5 ) & ((data['m1'] < 0.6) | (data['m2'] < 0.6))
			 & ((data['r1'] < wdMRrelation(data['m1'])) | (data['r2'] < wdMRrelation(data['m2'])))]
			all_WD.append(prsaWD)

			# Appending for Andrew
			eccAll.append(prsa['e'].values)
			pAll.append(prsa['p'].values)
			iAll.append(prsa['i'].values)
			m1All.append(prsa['m1'].values)
			m2All.append(prsa['m2'].values)
			r1All.append(prsa['r1'].values)
			r2All.append(prsa['r2'].values)
			magAll.append(prsa['appMagMean_r'].values)

			NallPrsa = len(prsa.index)
			if (Nall >= Nlim):
				#create histograms
				#All
				m1hAll0, m1b = np.histogram(data["m1"], bins=mbins)
				qhAll0, qb = np.histogram(data["m2"]/data["m1"], bins=qbins)
				ehAll0, eb = np.histogram(data["e"], bins=ebins)
				lphAll0, lpb = np.histogram(np.ma.log10(data["p"].values).filled(-999), bins=lpbins)
				dhAll0, db = np.histogram(data["d"], bins=dbins)
				maghAll0, magb = np.histogram(data["appMagMean_r"], bins=magbins)
				rhAll0, rb = np.histogram(data["r2"]/data["r1"], bins=rbins)

				if (doIndividualPlots):
					axmass.step(m1b[0:-1], m1hAll0/np.sum(m1hAll0), color='black', alpha=0.1)
					axqrat.step(qb[0:-1], qhAll0/np.sum(qhAll0), color='black', alpha=0.1)
					axecc.step(eb[0:-1], ehAll0/np.sum(ehAll0), color='black', alpha=0.1)
					axlper.step(lpb[0:-1], lphAll0/np.sum(lphAll0), color='black', alpha=0.1)
					axdist.step(db[0:-1], dhAll0/np.sum(dhAll0), color='black', alpha=0.1)
					axmag.step(magb[0:-1], maghAll0/np.sum(maghAll0), color='black', alpha=0.1)
					axrad.step(rb[0:-1], rhAll0/np.sum(rhAll0), color='black', alpha=0.1)

				#account for the binary fraction, as a function of mass
				dm1 = np.diff(m1b)
				m1val = m1b[:-1] + dm1/2.
				fb = np.sum(m1hAll0/len(data.index)*fbFit(m1val))
				#account for the hard-soft boundary
				Phs = getPhs(header['clusterVdisp'].iloc[0]*units.km/units.s).to(units.day).value
				fb *= RagNormal(np.log10(Phs), cdf = True)
				print("fb, Phs = ", fb, Phs)
				Nmult *= fb

							
				m1hAll += m1hAll0/Nall*Nmult
				qhAll += qhAll0/Nall*Nmult
				ehAll += ehAll0/Nall*Nmult
				lphAll += lphAll0/Nall*Nmult
				dhAll += dhAll0/Nall*Nmult
				maghAll += maghAll0/Nall*Nmult
				rhAll += rhAll0/Nall*Nmult

				#Obs
				obs = data.loc[data['LSM_PERIOD'] != -999]
				Nobs = len(obs.index)
				prsaObs = data.loc[(data['p'] < 1000) & (data['p'] >0.5) & (data['LSM_PERIOD'] != -999)]
				NobsPrsa = len(prsaObs.index)

				# White dwarf appending
				prsaObsWD = data.loc[(data['p'] < 1000) & (data['p'] > 0.5 ) & (data['LSM_PERIOD'] != -999)
					 & ((data['m1'] < 0.6) | (data['m2'] < 0.6))  & ((data['r1'] < wdMRrelation(data['m1'])) | (data['r2'] < wdMRrelation(data['m2'])))]
				obs_WD.append(prsaObsWD)

				# would like to see if there is a better way of doing this
				eccObs.append(prsaObs['e'].values)
				pObs.append(prsaObs['p'].values)
				iObs.append(prsaObs['i'].values)
				m1Obs.append(prsaObs['m1'].values)
				m2Obs.append(prsaObs['m2'].values)
				r1Obs.append(prsaObs['r1'].values)
				r2Obs.append(prsaObs['r2'].values)
				magObs.append(prsaObs['appMagMean_r'].values)

				if (Nobs >= Nlim):
					m1hObs0, m1b = np.histogram(obs["m1"], bins=mbins)
					qhObs0, qb = np.histogram(obs["m2"]/obs["m1"], bins=qbins)
					ehObs0, eb = np.histogram(obs["e"], bins=ebins)
					lphObs0, lpb = np.histogram(np.ma.log10(obs["p"].values).filled(-999), bins=lpbins)
					dhObs0, db = np.histogram(obs["d"], bins=dbins)
					maghObs0, magb = np.histogram(obs["appMagMean_r"], bins=magbins)
					rhObs0, rb = np.histogram(obs["r2"]/obs["r1"], bins=rbins)
					m1hObs += m1hObs0/Nall*Nmult
					qhObs += qhObs0/Nall*Nmult
					ehObs += ehObs0/Nall*Nmult
					lphObs += lphObs0/Nall*Nmult
					dhObs += dhObs0/Nall*Nmult
					maghObs += maghObs0/Nall*Nmult
					rhObs += rhObs0/Nall*Nmult

					#Rec
					recCombined = pd.DataFrame()
					prsaRecCombined = pd.DataFrame()
					for filt in filters:
						key = filt+'LSS_PERIOD'
						if (filt == 'all'):
							key = 'LSM_PERIOD'
						fullP = abs(data[key] - data['p'])/data['p']
						halfP = abs(data[key] - 0.5*data['p'])/(0.5*data['p'])
						twiceP = abs(data[key] - 2.*data['p'])/(2.*data['p'])
						rec = data.loc[(data[key] != -999) & ( (fullP < Pcut) | (halfP < Pcut) | (twiceP < Pcut))]
						prsaRec = data.loc[(data['p'] < 1000) & (data['p'] >0.5) & (data['LSM_PERIOD'] != -999) & ( (fullP < Pcut) | (halfP < Pcut) | (twiceP < Pcut))]
						Nrec = len(rec.index)

						#I'd like to account for all filters here to have more accurate numbers
						recCombined = recCombined.append(rec)
						prsaRecCombined = prsaRecCombined.append(prsaRec)

						# Going to use prsaRecCombined for ecc-p plots to account for all filters
						eccRec.append(prsaRec['e'].values)
						pRec.append(prsaRec['p'].values)
						iRec.append(prsaRec['i'].values)
						m1Rec.append(prsaRec['m1'].values)
						m2Rec.append(prsaRec['m2'].values)
						r1Rec.append(prsaRec['r1'].values)
						r2Rec.append(prsaRec['r2'].values)
						magRec.append(prsaRec['appMagMean_r'].values)

						# writeCornerFiles(prsaRec, ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r'], 'rec', 'M67','B','N')

						# White dwarf appending
						prsaRecWD = data.loc[(data['p'] < 1000) & (data['p'] > 0.5 ) & (data['LSM_PERIOD'] != -999)  & ( (fullP < Pcut) | (halfP < Pcut) | (twiceP < Pcut))
							 & ((data['m1'] < 0.6) | (data['m2'] < 0.6))  & ((data['r1'] < wdMRrelation(data['m1'])) | (data['r2'] < wdMRrelation(data['m2'])))]
						rec_WD.append(prsaRecWD)


						if (filt == 'all'):
							recCombined.drop_duplicates(inplace=True)
							prsaRecCombined.drop_duplicates(inplace=True)

						if (Nrec >= Nlim):
							m1hRec0, m1b = np.histogram(rec["m1"], bins=mbins)
							qhRec0, qb = np.histogram(rec["m2"]/rec["m1"], bins=qbins)
							ehRec0, eb = np.histogram(rec["e"], bins=ebins)
							lphRec0, lpb = np.histogram(np.ma.log10(rec["p"].values).filled(-999), bins=lpbins)
							dhRec0, db = np.histogram(rec["d"], bins=dbins)
							maghRec0, magb = np.histogram(rec["appMagMean_r"], bins=magbins)
							rhRec0, rb = np.histogram(rec["r2"]/rec["r1"], bins=rbins)
							m1hRec[filt] += m1hRec0/Nall*Nmult
							qhRec[filt] += qhRec0/Nall*Nmult
							ehRec[filt] += ehRec0/Nall*Nmult
							lphRec[filt] += lphRec0/Nall*Nmult
							dhRec[filt] += dhRec0/Nall*Nmult
							maghRec[filt] += maghRec0/Nall*Nmult
							rhRec[filt] += rhRec0/Nall*Nmult

							#for the mollweide
							if (filt == 'all'):
								Nrec = len(recCombined.index)
								rF = Nrec/Nall
								rN = Nrec/Nall*Nmult
								raN = Nmult
								obN = Nobs/Nall*Nmult
								fiN = Nall
								fioN = Nobs
								firN = Nrec

								NrecPrsa = len(prsaRecCombined.index)
								NrecPrsa = NrecPrsa/Nall*Nmult
								NobsPrsa = NobsPrsa/Nall*Nmult
								NallPrsa = NallPrsa/Nall*Nmult		



			recFrac.append(rF)
			recN.append(rN)
			rawN.append(raN)
			obsN.append(obN)
			fileN.append(fiN)
			fileObsN.append(fioN)
			fileRecN.append(firN)
			allNPrsa.append(NallPrsa)
			obsNPrsa.append(NobsPrsa)
			recNPrsa.append(NrecPrsa)
			#print(np.sum(lphRec), np.sum(recN), np.sum(lphRec)/np.sum(recN), np.sum(lphRec0), Nrec, np.sum(lphRec0)/Nrec, np.sum(lphObs), np.sum(obsN), np.sum(lphObs)/np.sum(obsN))

	# Concatenating param lists for 2D histograms -- eccentricity 1st (all 3 subpopulations)
	eccAll = np.concatenate(eccAll)
	eccObs = np.concatenate(eccObs)
	eccRec = np.concatenate(eccRec)

	# period
	pAll = np.concatenate(pAll)
	pObs = np.concatenate(pObs)
	pRec = np.concatenate(pRec)

	# Inclination
	iAll = np.concatenate(iAll)
	iObs = np.concatenate(iObs)
	iRec = np.concatenate(iRec)

	# Mass 1
	m1All = np.concatenate(m1All)
	m1Obs = np.concatenate(m1Obs)
	m1Rec = np.concatenate(m1Rec)

	# Mass 2
	m2All = np.concatenate(m2All)
	m2Obs = np.concatenate(m2Obs)
	m2Rec = np.concatenate(m2Rec)

	# Radius 1
	r1All = np.concatenate(r1All)
	r1Obs = np.concatenate(r1Obs)
	r1Rec = np.concatenate(r1Rec)

	# Radius 2
	r2All = np.concatenate(r2All)
	r2Obs = np.concatenate(r2Obs)
	r2Rec = np.concatenate(r2Rec)	

	# Apparent Magnitude
	magAll = np.concatenate(magAll)
	magObs = np.concatenate(magObs)
	magRec = np.concatenate(magRec)	


	# # print('Ecc lists:', eccAll, eccObs, eccRec)
	# # print('P lists:', pAll, pObs, pRec)
	# Appending lists with all the p/ecc values to our dataframes
	# All dataframe
	All['e'] = eccAll
	All['p'] = pAll
	All['i'] = iAll
	All['m1'] = m1All
	All['m2'] = m2All
	All['r1'] = r1All
	All['r2'] = r2All
	All['appMagMean_r'] = magAll

	# Observable dataframe
	Obs['e'] = eccObs
	Obs['p'] = pObs
	Obs['i'] = iObs
	Obs['m1'] = m1Obs
	Obs['m2'] = m2Obs
	Obs['r1'] = r1Obs
	Obs['r2'] = r2Obs
	Obs['appMagMean_r'] = magObs

	# Recovered dataframe
	Rec['e'] = eccRec
	Rec['p'] = pRec
	Rec['i'] = iRec
	Rec['m1'] = m1Rec
	Rec['m2'] = m2Rec
	Rec['r1'] = r1Rec
	Rec['r2'] = r2Rec
	Rec['appMagMean_r'] = magRec
	
	csv_cols = ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r']

	# 3 letter code corresponds to scenario (OC/GC, baseline/colossus, crowding/no crowding)
	All.to_csv('./data/all-GCN-histData.csv', header = ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r'])
	Obs.to_csv('./data/obs-GCN-histData.csv', header = ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r'])
	Rec.to_csv('./data/rec-GCN-histData.csv', header = ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r'])

	# Appending WD dataframes 
	WDall = pd.concat(all_WD)
	WDobs = pd.concat(obs_WD)
	WDrec = pd.concat(rec_WD)

	# Only want certain columns from df
	WDall = WDall[csv_cols]
	WDobs = WDobs[csv_cols]
	WDrec = WDrec[csv_cols]

	print('White Dwarf Candidates: ', WDall, WDobs, WDrec)

	WDall.to_csv('./data/wd/all-GCN-WD-histData.csv', header = ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r'])
	WDobs.to_csv('./data/wd/obs-GCN-WD-histData.csv', header = ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r'])
	WDrec.to_csv('./data/wd/rec-GCN-WD-histData.csv', header = ['p', 'm1', 'm2', 'r1', 'r2', 'e', 'i', 'appMagMean_r'])

	#plot and save the histograms
	saveHist(m1hAll, m1hObs, m1hRec, m1b, 'm1 (Msolar)', 'EBLSST_m1hist')
	saveHist(qhAll, qhObs, qhRec, qb, 'q (m2/m1)', 'EBLSST_qhist')
	saveHist(ehAll, ehObs, ehRec, eb, 'e', 'EBLSST_ehist')
	saveHist(lphAll, lphObs, lphRec, lpb, 'log(P [days])', 'EBLSST_lphist')
	saveHist(dhAll, dhObs, dhRec, db, 'd (kpc)', 'EBLSST_dhist')
	saveHist(maghAll, maghObs, maghRec, magb, 'mag', 'EBLSST_maghist')
	saveHist(rhAll, rhObs, rhRec, rb, 'r2/r1', 'EBLSST_rhist')


	#make the mollweide
	coords = SkyCoord(RA, Dec, unit=(units.degree, units.degree),frame='icrs')	
	lGal = coords.galactic.l.wrap_at(180.*units.degree).degree
	bGal = coords.galactic.b.wrap_at(180.*units.degree).degree
	RAwrap = coords.ra.wrap_at(180.*units.degree).degree
	Decwrap = coords.dec.wrap_at(180.*units.degree).degree

	f, ax = plt.subplots(subplot_kw={'projection': "mollweide"}, figsize=(8,5))
	ax.grid(True)
	#ax.set_xlabel(r"$l$",fontsize=16)
	#ax.set_ylabel(r"$b$",fontsize=16)
	#mlw = ax.scatter(lGal.ravel()*np.pi/180., bGal.ravel()*np.pi/180., c=np.log10(np.array(recFrac)*100.), cmap='viridis_r', s = 4)
	ax.set_xlabel("RA",fontsize=16)
	ax.set_ylabel("Dec",fontsize=16)
	mlw = ax.scatter(np.array(RAwrap).ravel()*np.pi/180., np.array(Decwrap).ravel()*np.pi/180., c=np.array(recFrac)*100., cmap='viridis_r', s = 4)
	cbar = f.colorbar(mlw, shrink=0.7)
	cbar.set_label(r'% recovered')
	f.savefig('./plots/analyse_plots/'  + 'mollweide_pct.pdf',format='pdf', bbox_inches = 'tight')

	f, ax = plt.subplots(subplot_kw={'projection': "mollweide"}, figsize=(8,5))
	ax.grid(True)
	#ax.set_xlabel(r"$l$",fontsize=16)
	#ax.set_ylabel(r"$b$",fontsize=16)
	#mlw = ax.scatter(lGal.ravel()*np.pi/180., bGal.ravel()*np.pi/180., c=np.log10(np.array(recN)), cmap='viridis_r', s = 4)
	ax.set_xlabel("RA",fontsize=16)
	ax.set_ylabel("Dec",fontsize=16)
	mlw = ax.scatter(np.array(RAwrap).ravel()*np.pi/180., np.array(Decwrap).ravel()*np.pi/180., c=np.log10(np.array(recN)), cmap='viridis_r', s = 4)
	cbar = f.colorbar(mlw, shrink=0.7)
	cbar.set_label(r'log10(N) recovered')
	f.savefig('./plots/analyse_plots/'  + 'mollweide_N.pdf',format='pdf', bbox_inches = 'tight')

	if (doIndividualPlots):
		fmass.savefig('./plots/analyse_plots/'  + 'massPDFall.pdf',format='pdf', bbox_inches = 'tight')
		fqrat.savefig('./plots/analyse_plots/'  + 'qPDFall.pdf',format='pdf', bbox_inches = 'tight')
		fecc.savefig('./plots/analyse_plots/'  + 'eccPDFall.pdf',format='pdf', bbox_inches = 'tight')
		flper.savefig('./plots/analyse_plots/'  + 'lperPDFall.pdf',format='pdf', bbox_inches = 'tight')
		fdist.savefig('./plots/analyse_plots/'  + 'distPDFall.pdf',format='pdf', bbox_inches = 'tight')
		fmag.savefig('./plots/analyse_plots/'  + 'magPDFall.pdf',format='pdf', bbox_inches = 'tight')
		frad.savefig('./plots/analyse_plots/'  + 'radPDFall.pdf',format='pdf', bbox_inches = 'tight')

	print("###################")
	print("number of binaries in input files (raw, log):",np.sum(fileN), np.log10(np.sum(fileN)))
	print("number of binaries in tested with gatspy (raw, log):",np.sum(fileObsN), np.log10(np.sum(fileObsN)))
	print("number of binaries in recovered with gatspy (raw, log):",np.sum(fileRecN), np.log10(np.sum(fileRecN)))
	print("recovered/observable*100 with gatspy:",np.sum(fileRecN)/np.sum(fileObsN)*100.)
	print("###################")
	print("total in sample (raw, log):",np.sum(rawN), np.log10(np.sum(rawN)))
	print("total observable (raw, log):",np.sum(obsN), np.log10(np.sum(obsN)))
	print("total recovered (raw, log):",np.sum(recN), np.log10(np.sum(recN)))
	print("recovered/observable*100:",np.sum(recN)/np.sum(obsN)*100.)
	print("###################")
	print("total in Prsa 15.8<r<19.5 P<1000d sample (raw, log):",np.sum(allNPrsa), np.log10(np.sum(allNPrsa)))
	print("total observable in Prsa 15.8<r<19.5 P<1000d sample (raw, log):",np.sum(obsNPrsa), np.log10(np.sum(obsNPrsa)))
	print("total recovered in Prsa 15.8<r<19.5 P<1000d sample (raw, log):",np.sum(recNPrsa), np.log10(np.sum(recNPrsa)))
	print("Prsa 15.8<r<19.5 P<1000d rec/obs*100:",np.sum(recNPrsa)/np.sum(obsNPrsa)*100.)


