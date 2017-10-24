# vaisalaflux

Calculate chamber CO2 flux from Vaisala IRGA GMP343. Calculations do not account for changing chamber humidity.

* Check script for local adjustments (csv format, chamber size etc).
* Put script in a folder together with downloaded csv files (save as csv when you download GMP343 data)
* Make sure that working directory is the same as the folder ``` getwd()```. Run script
```
source("vaisala_flux_calc.R"")
```
```
flux.dat #contains the flux values. 
```
A [pdf](https://github.com/ggranath/vaisalaflux/blob/master/co2VStime.pdf) is created to inspect data and model fit.

# DISCLAIMER
The code has not been reviewed or verified and it may contain errors. Use at your own risk.