
Images (on Sylabs cloud): 
```
library://timniem/rna_outliers/test:tmae #MAE container
library://timniem/rna_outliers/test:hg19hg38 #FRASER and OUTRIDER container
library://timniem/rna_outliers/test:report #HTML report container
```

Add the SylabsCloud to the remotes on Apptainer:
```
apptainer remote add --no-login SylabsCloud cloud.sylabs.io
apptainer remote use SylabsCloud
```

If not already configured:
```
export APPTAINER_CACHEDIR=/path/to/tmp
```
Get the Singularity images:
```
apptainer pull --dir 'path/to/cache/dir' container_name.sif library://timniem/rna_outliers/XXX:XXX
```
