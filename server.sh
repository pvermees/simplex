Rscript -e "remotes::install_github(c('tim-band/shinylight', 'pvermees/isoplotr'))"
Rscript -e "install.packages('.', repos=NULL)"
Rscript -e "simplex::simplex($1)"
