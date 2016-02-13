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

1.	Comp2Sets.py
2.	DiffTrigger.py
3.	GetBest.py 
4.	KrigTrigger.py
5.	lsa.py
6.	lsaModel.py
7.	LSASO.py
    - Main Model file. Runs simplex algorithm to find best parameter combination unless parameters provided in settings and called with "p".  Iteratively calculates unmeasured snow.
8.	nao.py
9.	opparan.py
    - Analyse data from optimisation runs on all years, then runs using averages for parameters derived from all years, b1 years and b2 years
10.	README.md
11.	Shade.py
  -  Creates a geotiff file containing multiple layers, each layer the shade cast on glacier for each julian day
12.	smhiPrecipPlot.py
13.	smhisort.sh
14.	SumTrigger.py
