## Mälaren model setup
This folder contains a pyGETM-Selmaprotbas setup for Mälaren. The run.py script has been set up to run in chunks of 140 days. In the run_malaren.py script, you can control most of the model's behaviour in the upper part of the script, although more advanced changes would require more intricate changes to the script. The selmaprotbas model has been calibrated on [a multi-basin 1D GOTM-Selmaprotbas model](http://www.hydroshare.org/resource/fe8e0b2728324c7892da9a3312f358c2). This Github repository only has the ERA5 forcing for 2001-2003, although the model could run for a longer period if provided with meteorological forcing over a longer period.

The outflow from the lake is regulagulated at runtime, with QH-relationships for three different outflow locations (based on observed data and adjusted to reprodoce observed water level dynamics). 

1D output files are being generated for several locations in the lake, at a higher frequency that the all-lake 3D output.

### Additional files
Additional files related to this model setup can be found on a [Hydroshare repository](http://www.hydroshare.org/resource/608fe6b4eb9d46de8f700db96c9f6664)
