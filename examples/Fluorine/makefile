lr_F.tar.gz: lr.dal F.mol mc_F.tar.gz
	dalton -get "AOPROPER AO2SOINT" -f mc_F lr F
	rm -f lr_F.AO*

mc_F.tar.gz: mc.dal F.mol
	dalton mc F

clean:
	rm -f lr_F.hso* lr_F.out lr_F.tar.gz mc_F.tar.gz mc_F.out
