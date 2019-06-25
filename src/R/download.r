# To install, use:
# install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

library(synapser)

downloadSynapse = function(basedir, email, password) {
    # To access to Synapse, visit https://www.synapse.org/
    synLogin(email = email, password = password)

    # MSBB:
    # syn8691099 - all counts
    # syn6100548 - covariates
    # syn6101474 - clinical

    msbb_files = c('syn8691099', 'syn6100548', 'syn6101474')

    for (file in msbb_files) {
        synGet(file, downloadLocation = file.path(basedir, 'MSBB'))
    }

    # MayoRnaSeq
    # CBE - cerebellum
    # syn8690904 - all counts (CBE)
    # syn5223705 - covariates (CBE)
    # syn6126119 - details (CBE)

    mayo_rnaseq_files = c('syn8690904', 'syn5223705', 'syn6126119')

    for (file in mayo_rnaseq_files) {
        downloadLocation = file.path(basedir, 'MayoCBE')
        synGet(file, downloadLocation = downloadLocation)
    }

    # TCX - temporal cortex
    # syn8690799 - all counts (TCX)
    # syn3817650 - covariates (TCX)
    # syn6126114 - details

    mayo_tcx_files = c('syn8690799', 'syn3817650', 'syn6126114')

    for (file in mayo_tcx_files) {
        downloadLocation = file.path(basedir, 'MayoTCX')
        synGet(file, downloadLocation = downloadLocation)
    }

    # ROSMAP
    # syn8691134 - all counts
    # syn3382527 - id key
    # syn3191087 - clinical

    rosmap_files = c('syn8691134', 'syn3382527', 'syn3191087')

    for (file in rosmap_files) {
        downloadLocation = file.path(basedir, 'ROS')
        synGet(file, downloadLocation = downloadLocation)
    }
}

main = function() {
    args = commandArgs(trailingOnly = TRUE)
    basedir = args[0]
    synapseUser = args[1]
    synapsePassword = args[2]

    if (! file.exists(basedir)) {
        dir.create(basedir)
    }

    downloadSynapse(basedir, synapseUser, synapsePassword)
}

main()
