# scripts
Scripts for calculating observables from hydro model

## dN/dη
###### `dndeta.c`
This script reads three parameters:
- `directory` - directory of your root files
- `Npart_min` - lower limit of your cut on number of participants
- `Npart_max` - upper limit of your cut on number of participants

The script can be run either by
```root 'dndeta.c("directory",Npart_min,Npart_max)'```
or in batch mode
```root -q -b 'dndeta.c("directory",Npart_min,Npart_max)'```

Other parameters can be changed in the beginning of the script and are set to:
- Pseudorapidity interval = (-5.0; 5.0)
- Number of bins = 50

Output file contains the path to data and two columns:
```η dN/dη```
