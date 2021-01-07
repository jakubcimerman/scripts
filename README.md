# scripts
Scripts for calculating observables from hydro model

## dN/dη
###### `dndeta.c`
This script reads three parameters:
- `directory` - directory of your root files
- `Npart_min` - lower limit of your cut on number of participants
- `Npart_max` - upper limit of your cut on number of participants

The script can be run either by
```bash
root 'dndeta.c("directory",Npart_min,Npart_max)'
```
or in batch mode
```bash
root -q -b 'dndeta.c("directory",Npart_min,Npart_max)'
```

Other parameters can be changed in the beginning of the script and are set to:
- Pseudorapidity interval = (-5.0; 5.0)
- Number of bins = 50

Output file contains the path to data and two columns:

| η | dN/dη |
|---|---|

## v<sub>n</sub> {EP} (y)
###### `vn_EP_rapidity.c`
This script calculate the flow of any order using Event Plane method as a function of rapidity. The script reads four parameters:
- `directory` - directory of your root files
- `Npart_min` - lower limit of your cut on number of participants
- `Npart_max` - upper limit of your cut on number of participants
- `order` - order of the flow coefficient you want to calculate

The script can be run either by
```bash
root 'vn_EP_rapidity.c("directory",Npart_min,Npart_max,order)'
```
or in batch mode
```bash
root -q -b 'vn_EP_rapidity.c("directory",Npart_min,Npart_max,order)'
```

## v<sub>n</sub> {EP} (y) for identified hadrons
###### `vn_EP_rapidity_IH.c`
This script calculate the flow of any order using Event Plane method as a function of rapidity for identified hadrons. The script reads four parameters:
- `directory` - directory of your root files
- `Npart_min` - lower limit of your cut on number of participants
- `Npart_max` - upper limit of your cut on number of participants
- `order` - order of the flow coefficient you want to calculate
- `pid` - particle ID

The script can be run either by
```bash
root 'vn_EP_rapidity_IH.c("directory",Npart_min,Npart_max,order,pid)'
```
or in batch mode
```bash
root -q -b 'vn_EP_rapidity_IH.c("directory",Npart_min,Npart_max,order,pid)'
```

## v<sub>n</sub> {RP} (y) for identified hadrons
###### `vn_RP_rapidity_IH.c`
This script calculate the flow of any order using Reaction Plane method as a function of rapidity for identified hadrons. The script reads four parameters:
- `directory` - directory of your root files
- `Npart_min` - lower limit of your cut on number of participants
- `Npart_max` - upper limit of your cut on number of participants
- `order` - order of the flow coefficient you want to calculate
- `pid` - particle ID

The script can be run either by
```bash
root 'vn_RP_rapidity_IH.c("directory",Npart_min,Npart_max,order,pid)'
```
or in batch mode
```bash
root -q -b 'vn_RP_rapidity_IH.c("directory",Npart_min,Npart_max,order,pid)'
```
