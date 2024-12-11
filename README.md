
MAE container: 
<code>
library://timniem/rna_outliers/test:tmae
</code>
OUTRIDER FRASER container:
<code>
library://timniem/rna_outliers/test:hg19hg38
</code>

HTML report container:
<code>
library://timniem/rna_outliers/test:report
</code>

Add the Sylabscloud to the remotes on Apptainer
<code>
apptainer remote add --no-login SylabsCloud cloud.sylabs.io
apptainer remote use SylabsCloud
</code>

if not already configured
export APPTAINER_CACHEDIR=/path/to/tmp
apptainer pull --dir 'path/to/cache/dir' container_name.sif library://timniem/rna_outliers/XXX:XXX
