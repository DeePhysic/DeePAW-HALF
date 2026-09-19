.PHONY = all \
	 html \
	 pdf \
	 clean-all \
	 clean \
	 sources \
	 clean-doxygen \
	 clean-doxygen-awesome-css

DOXYGEN_VERSION = "Release_1_10_0"
THEME_VERSION = "v2.3.2"

include ../makefile.include

PREDEFINED=$(subst -D,,$(CPP_OPTIONS))
SRC_LIB=../src/libvaspml
SRC_APP=../src/applications
SRC_TEST=../test/src/libvaspml

all: html

html: doxygen doxygen-awesome-css sources
	#( cat Doxyfile; echo "PREDEFINED += $(PREDEFINED)" ) | ./doxygen -
	( cat Doxyfile; echo "PREDEFINED += $(PREDEFINED)" ) > Doxyfile.mod
	# Maybe include output of all macros defined by the compiler:
	# gfortran -dM -E - < /dev/null
	./doxygen Doxyfile.mod
	#./doxygen

pdf: doxygen doxygen-awesome-css sources
	( cat Doxyfile; echo "PREDEFINED += $(PREDEFINED)" ) > Doxyfile.mod
	@sed -i.bak "s;^GENERATE_LATEX .*$$;GENERATE_LATEX         = YES;" Doxyfile.mod
	./doxygen Doxyfile.mod
	@sed -i.bak "s;^    top=2.5cm;    top=1.5cm;" build/latex/refman.tex
	@sed -i.bak "s;^    bottom=2.5cm;    bottom=1.5cm;" build/latex/refman.tex
	@sed -i.bak "s;^    left=2.5cm;    left=1.0cm;" build/latex/refman.tex
	@sed -i.bak "s;^    right=2.5cm;    right=1.0cm;" build/latex/refman.tex
	@sed -i.bak "s;courier;dejavu;" build/latex/refman.tex
	@sed -i.bak "s;scriptsize;small;" build/latex/doxygen.sty
	@sed -i.bak "s;pdflatex;pdflatex -interaction=nonstopmode;" build/latex/Makefile
	-cd build/latex && make pdf
	-cd build/latex && pdflatex -interaction=nonstopmode refman.tex
	cp build/latex/refman.pdf vaspml.pdf

clean-all: clean clean-doxygen clean-doxygen-awesome-css

clean:
	$(RM) -r build
	$(RM) -r src-vaspml
	$(RM) Doxyfile.mod
	$(RM) doxy.log
	$(RM) vaspml.pdf

sources:
	$(RM) -r src-vaspml && \
	mkdir src-vaspml && \
	rsync -u --ignore-missing-args $(SRC_LIB)/*.cpp $(SRC_LIB)/*.hpp src-vaspml && \
	rsync -u --ignore-missing-args $(SRC_APP)/*.cpp $(SRC_APP)/*.hpp src-vaspml && \
	rsync -u --ignore-missing-args $(SRC_TEST)/*.cpp $(SRC_TEST)/*.hpp src-vaspml

doxygen:
	git clone https://github.com/doxygen/doxygen.git src-doxygen
	cd src-doxygen && \
	git checkout ${DOXYGEN_VERSION} && \
	mkdir build && \
	cd build && \
	cmake -G "Unix Makefiles" .. && \
	$(MAKE) -j16 && \
	cp bin/doxygen ../..

clean-doxygen:
	$(RM) -r src-doxygen doxygen

doxygen-awesome-css:
	git clone https://github.com/jothepro/doxygen-awesome-css.git
	cd doxygen-awesome-css && \
	git checkout ${THEME_VERSION}

clean-doxygen-awesome-css:
	$(RM) -r doxygen-awesome-css
