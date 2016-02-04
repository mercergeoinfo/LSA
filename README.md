# LSA
Late season snow accumulation model  

## Background

Due to logistical contraints the timing of the 'end of season' accumulation survey rarely coincides with the true end of season.
In 2012 this produced a noticeable error in the summer balance, resulting in negative ablation on many parts of the glacier.
The net balance is not affected by this as the unobserved accumulation is exactly balanced by unobserved ablation.
This fact allows us to estimate this unobserved accumulation by firstly modelling expected ablation and differencing it against the observed ablation.

## Model

The basic model is a positive degree day model.
The temperature on glacier is adjusted from that measured at a weather station using an adiabatic lapse rate and a shading model.
The lapse rate is itself elevation dependent, empirically derived but justified in the literature.
The shading is a hillshade based adjusted summing the "daylight" hours shade at each pixel and normalising this against the maximum value found for the whole glacier on each day.

## Files

### ModifiedDegreeDay.py

__DEAD__The main model. Requires spatialfiles.py and standard.py from my other repo

### ParameterAnalysis.py

_DEAD__Used only for assessment of model parameters

### Shade.py

Creates a geotiff file containing multiple layers, each layer the shade cast on glacier for each julian day

### lsa.py

The melt model script. Primary usage is for deriving optimum and best scoring parameter vales

### lsaModel.py

Version of lsa.py to be main model with iterative recalculation of snow depth

### Other files

- ddf_model.py
- DiffTrigger.py
- KrigTrigger.py
- README.md
- smhiPrecipPlot.py
- SumTrigger.py


