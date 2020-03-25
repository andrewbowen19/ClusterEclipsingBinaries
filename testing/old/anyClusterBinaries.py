"""
Class derived from getClusterBinaries - want it to be able to read in any table from the web
Want it to take in same inputs as getClusterBinaries, plus a url (to read-in from web)

"""


import pandas as pd
import numpy as np
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.sample.sampler import independent
from cosmic.sample.sampler import multidim
from cosmic.evolve import Evolve


class anyClusterBinaries(object):
	
	def __init__(self, url):
		self.url = url
		self.age = 0
		self.Z = 0
		self.sigma = 0
		self.Nbin = 0

	# Method to read-in data from any cluster datafile: needs to pull correct info
	def cluster_readin(url):
		# Reading in from web

		clusters = pd.read_html(url, header = 0)[0]
		# self.clusters = clusters
test_url = 'https://webda.physics.muni.cz/tadross.html'

test_run = anyClusterBinaries.cluster_readin(url = test_url)

print(type(test_run), test_run)