
# MAE container: 
library://timniem/rna_outliers/test:tmae

# OUTRIDER FRASER container:
library://timniem/rna_outliers/test:hg19hg38

# Add the Sylabscloud to the remotes on Apptainer
apptainer remote add --no-login SylabsCloud cloud.sylabs.io
apptainer remote use SylabsCloud

# if not already configured
export APPTAINER_CACHEDIR=/path/to/tmp
apptainer pull --dir 'path/to/cache/dir' container_name.sif library://timniem/rna_outliers/XXX:XXX
