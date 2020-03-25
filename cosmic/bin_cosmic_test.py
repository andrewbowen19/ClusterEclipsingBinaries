# test script for cosmic runs - let's generate some binaries!

import cosmic

from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.sample.sampler import independent
from cosmic.evolve import Evolve


final_kstar1 = [1]
final_kstar2 = [1]
InitialBinaries, sampled_mass, n_sampled = InitialBinaryTable.sampler('independent', final_kstar1, final_kstar2, \
	primary_model='kroupa93', ecc_model='thermal', SFH_model='const', \
	component_age=10000.0, met=0.02, size=10000)


print(InitialBinaries)