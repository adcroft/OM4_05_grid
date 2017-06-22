# Makefile to create supergrid.nc and interpolated_topog.nc
# To use:
#   module load python
#   setenv PYTHONPATH $cwd/MIDAS
#
# then
#   make all

SHELL=tcsh -f
NP=8

all: ocean_hgrid.nc
	md5sum -c md5sums.txt

showenv:
	env
	-set
	-module list
	which python
	-python --version

# Grids
supergrid.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python merge_grids.py

mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc: local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python create_grids.py

# Sets char tile='tile1'
ocean_hgrid.nc: supergrid.nc
	\cp $< $@
	./changeChar.py ocean_hgrid.nc tile tile1

MIDAS:
	git clone https://github.com/mjharriso/MIDAS.git
	(cd MIDAS; git checkout cb00d284c7b3ac652ba1dcf93ca42f02bc3c4a91)

local: MIDAS 
	-rm -rf $</build/*
	mkdir -p $@
	cd $<; make -f Makefile_GFDL INSTALL_PATH=../local
	touch $@


%.cdl: %.nc
	ncdump $< | egrep -v 'code_version|history' > $@

md5sums.txt: ocean_hgrid.nc antarctic_spherical_supergrid.nc mercator_supergrid.nc ncap_supergrid.nc supergrid.nc
	echo Grids > $@
	md5sum *supergrid.nc ocean_hgrid.nc >> $@

clean:
	rm -rf MIDAS local *.nc
