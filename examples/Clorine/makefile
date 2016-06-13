mol=Cl

lr_$(mol).tar.gz: lr.dal $(mol).mol mc_$(mol).tar.gz
	dalton -get "AOPROPER AO2SOINT" -f mc_$(mol) lr $(mol)
	rm -f lr_$(mol).AO*

mc_$(mol).tar.gz: mc.dal $(mol).mol
	dalton mc $(mol)

clean:
	rm -f lr_$(mol).hso* lr_$(mol).out lr_$(mol).tar.gz mc_$(mol).tar.gz mc_$(mol).out
