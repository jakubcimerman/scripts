# scripts
Scripts for calculating observables from hydro model

## r<sub>n</sub> (η)
###### `decorrelation.c`
This script reads ten parameters:
- `directory` - directory of your root files
- `order` - order of the decorrelation coefficient you want to calculate
- `eventStep` - this calcultion uses oversampling of final events by creating so-called superevents. This parameter sets the number of merged events
- `etaRefMin` & `etaRefMax` - the reference region of pseudorapidity
- `etaTestMin` & `etaTestMax` - the test region of pseudorapidity
- `Npart_min` - lower limit of your cut on number of participants
- `Npart_max` - upper limit of your cut on number of participants
- `isSym` - set 1 for symmetric collision and 0 for asymmetric collision
- `option` - set 0 for calculating basic r<sub>n</sub>, 1 for calculating r<sub>n</sub><sup>v</sup>, 2 for calculating r<sub>n</sub><sup>ψ</sup>

The script can be run either by
```bash
root 'decorrelation.c("directory",order,eventStep,etaRefMin,etaRefMax,etaTestMin,etaTestMax,Npart_min,Npart_max,isSym,option)'
```
or in batch mode
```bash
root -q -b 'decorrelation.c("directory",order,eventStep,etaRefMin,etaRefMax,etaTestMin,etaTestMax,Npart_min,Npart_max,isSym,option)'
```

Other parameters can be changed in the beginning of the script and are set to:
- Transverse momentum cut = (0.4; 4.0)
- Number of bins = 5

Output file contains the path to data, order and three columns:

| η | r<sub>n</sub> | error |
|---|---|---|


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

## d<sup>2</sup>N/(2πp<sub>T</sub>dp<sub>T</sub>dy) - p<sub>T</sub> spectrum of π<sup>+</sup>, K<sup>+</sup> and protons
###### `spectrum.c`
This script reads three parameters:
- `directory` - directory of your root files
- `Npart_min` - lower limit of your cut on number of participants
- `Npart_max` - upper limit of your cut on number of participants

The script can be run either by
```bash
root 'spectrum.c("directory",Npart_min,Npart_max)'
```
or in batch mode
```bash
root -q -b 'spectrum.c("directory",Npart_min,Npart_max)'
```

Other parameters can be changed in the beginning of the script and are set to:
- Transverse momentum cut = (0.2; 2.0)
- Rapidity cut = (-0.1; 0.1)
- Number of bins = 36

Output file contains the path to data and four columns:

| p<sub>T</sub> | π<sup>+</sup> spectrum | K<sup>+</sup> spectrum | p spectrum |
|---|---|---|---|

## v<sub>n</sub> {EP} (η)
###### `vn_EP_pseudorapidity.c`
This script calculate the flow of any order using Event Plane method as a function of rapidity. The script reads four parameters:
- `directory` - directory of your root files
- `Npart_min` - lower limit of your cut on number of participants
- `Npart_max` - upper limit of your cut on number of participants
- `order` - order of the flow coefficient you want to calculate

The script can be run either by
```bash
root 'vn_EP_pseudorapidity.c("directory",Npart_min,Npart_max,order)'
```
or in batch mode
```bash
root -q -b 'vn_EP_pseudorapidity.c("directory",Npart_min,Npart_max,order)'
```

Output file contains the path to data and three columns:

| η | v<sub>n</sub> | error |
|---|---|---|

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
