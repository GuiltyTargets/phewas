install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))

install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))

library(synapser)

# To access to Synapse, visit https://www.synapse.org/
synLogin(email="youremail", password="yourpassword")

## Download files for MSBB:
# syn8691099 - all counts
# syn6100548 - covariates
# syn6101474 - clinical

basedir = 'data'

files = c('syn8691099', 'syn6100548', 'syn6101474')

for (i in files) {
 synGet(i, downloadLocation=file.path(basedir, 'MSBB'))
}

## Download files for MayoRnaSeq
# CBE - cerebellum
# syn8690904 - all counts (CBE)
# syn5223705 - covariates (CBE) 
# syn6126119 - details (CBE)

files = c('syn8690904', 'syn5223705', 'syn6126119')

for (i in files) {
 synGet(i, downloadLocation=file.path(basedir, 'MayoCBE'))
}

# TCX - temporal cortex
# syn8690799 - all counts (TCX)
# syn3817650 - covariates (TCX)
# syn6126114 - details

files = c('syn8690799', 'syn3817650', 'syn6126114')

for (i in files) {
 synGet(i, downloadLocation=file.path(basedir, 'MayoTCX'))
}

## Download files for ROSMAP

# syn8691134 - all counts
# syn3382527 - id key
# syn3191087 - clinical

files = c('syn8691134', 'syn3382527', 'syn3191087')

for (i in files) {
 synGet(i, downloadLocation=file.path(basedir, 'ROS'))
}
