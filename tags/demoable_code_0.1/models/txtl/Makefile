# Makefile for MATLAB TXTL modeling library
# RMM, 29 Sep 2012

VERSION = 0.1a
LIBRARY = findspecies.m txtl_addspecies.m txtl_buffer.m txtl_combine.m \
  txtl_dna.m txtl_extract.m txtl_newtube.m txtl_prom_p70.m txtl_prom_ptet.m \
  txtl_protein_deGFP.m txtl_protein_tetR.m txtl_rnap_rnap70.m txtl_template.m \
  txtl_utr_rbs.m
EXAMPLES = gamS_plot.m geneexpr.m induction.m negautoreg.m

build:
	mkdir -p txtl-$(VERSION)
	cp -pf $(LIBRARY) $(EXAMPLES) txtl-$(VERSION)
	tar zcf txtl-$(VERSION).tgz txtl-$(VERSION)

