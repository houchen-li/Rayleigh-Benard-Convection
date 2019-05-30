#!/bin/zsh
 
./lhc_model --Prandtl_number=1.0 --Rayleigh_number=100000000.0 --wave_number=0.5 --start_time=0.0 --end_time=0.05 --output_h5file=data.h5
python plotNuEvolution.py
python RBCField.py
python plotTemperatureProfile.py
