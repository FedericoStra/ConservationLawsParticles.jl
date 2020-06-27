using Documenter
using Particles

makedocs(
    sitename = "Particles",
    authors = "Federico Stra",
    format = Documenter.HTML(prettyurls = false),
    modules = [Particles]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
