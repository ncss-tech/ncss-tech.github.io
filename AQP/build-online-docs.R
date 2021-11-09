# use pkgdown in each package folder
# output is stored in pkg/docs

# https://pkgdown.r-lib.org/articles/pkgdown.html
# configuration in inst/_pkgdown.yml

# interesting note: par('bg') is not set to 'white' when running build_site() ... why?

# install from CRAN
pkgdown::build_site()
unlink('docs/reference/*.pdf')



