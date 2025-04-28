PACKAGE := ConservationLawsParticles


.PHONY: default
default: makedocs


.PHONY: clean
clean:
	$(RM) -r build biblio/aux


.PHONY: makedocs
makedocs: setup-project
	julia --project make.jl


DOCTEST_CMD = \
using Documenter: DocMeta, doctest; \
using $(PACKAGE); \
DocMeta.setdocmeta!($(PACKAGE), :DocTestSetup, :(using $(PACKAGE)); recursive=true); \
doctest($(PACKAGE))

.PHONY: doctest
doctest: setup-project
	julia --project -e '$(DOCTEST_CMD)'


SETUP_PROJECT_CMD = \
using Pkg ;\
Pkg.develop(path="..") ;\
Pkg.instantiate()

.PHONY: setup-project
setup-project:
	julia --project -e '$(SETUP_PROJECT_CMD)'
