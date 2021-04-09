.PHONY: all
all: test


.PHONY: test
test:
	julia --project -e 'import Pkg; Pkg.test()'


.PHONY: makedocs
makedocs:
	$(MAKE) -C docs makedocs


.PHONY: doctest
doctest:
	$(MAKE) -C docs doctest
