lr_Cl.out: lr.dal Cl.mol mc_Cl.tar.gz
	dalton -get 'AOPROPER AO2SOINT' -f mc_Cl lr Cl
	rm -f lr_Cl.AOPROPER lr_Cl.AO2SOINT

mc_Cl.tar.gz: mc.dal Cl.mol
	dalton mc Cl

