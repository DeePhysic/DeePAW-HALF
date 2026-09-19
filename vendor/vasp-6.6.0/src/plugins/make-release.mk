.PHONY: release force_build

FILES = interface_cpp.cpp interface_fortran.F ../plugins.F

UNIFDEF = unifdef -x2
CPP_OPTIONS ?=

release: $(FILES)
force_build: ;

$(FILES): force_build
	$(UNIFDEF) $(CPP_OPTIONS) $@ > $@.tmp
	mv $@.tmp $@

