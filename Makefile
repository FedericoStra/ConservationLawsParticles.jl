.PHONY: default
default: test


.PHONY: all
all: test makedocs


.PHONY: clean
clean:
	$(MAKE) -C docs clean


.PHONY: test
test:
	julia --project -e 'import Pkg; Pkg.test()'


.PHONY: makedocs
makedocs:
	$(MAKE) -C docs makedocs


.PHONY: doctest
doctest:
	$(MAKE) -C docs doctest
