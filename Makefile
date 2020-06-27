.PHONY: test

test:
	julia --project --color=yes -e 'import Pkg; Pkg.test()'
