# DCOPF
DC Optimal Power Flow (DCOPF): (i) one lossless DCOPF model and (ii) two lossy DCOPF models.

This set of AMPL codes implements a regular lossless DCOPF model and two lossy DCOPF model. The test case used here is a standard IEEE 73-bus system; but these codes can work on any other systems.

Model 1 (M1) for lossy DCOPF:
# implemented the model (2)-(14) of the following paper: O. W. Akinbode and K. W. Hedman "Fictitious losses in the DCOPF with a piecewise linear approximation of losses," IEEE PES General Meeting, Jul. 2013.

Model 2 (M2) for lossy DCOPF:
# Iterative Loss Approximation based model
# Each iteration: calc total losses and then the change in loss than previous iteration, and then assign the losses to all/selected buses with same/different participation factors.
# This method converges when the difference of total losses between two consecutive iterations is less than a threshold.
# Possibly, two iterations repeat themselves, which cause method to diverge.


License:
This work is licensed under the terms of the Creative Commons Attribution 4.0 (CC BY 4.0) license. 
https://creativecommons.org/licenses/by/4.0/


Disclaimer:
The author doesn’t make any warranty for the accuracy, completeness, or usefulness of any information disclosed and the author assumes no liability or responsibility for any errors or omissions for the information (data/code/results etc) disclosed.