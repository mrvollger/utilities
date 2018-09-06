To use this it is probably best to load the anaconda enviorment on the cluster
```
#!/usr/bin/env bash
unset PYTHONPATH
module purge
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs/prod modules-eichler
module load anaconda/201710
```
or just 
```
source conda.cfg
```

