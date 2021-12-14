## DCOPF
DC Optimal Power Flow (DCOPF): (i) one lossless DCOPF model and (ii) two lossy DCOPF models.

This set of AMPL codes implements a normal lossless DCOPF model and two lossy DCOPF model. The test case used here is a modified IEEE 73-bus system; but these codes can work on any other systems.

Model 1 (M1) for lossy DCOPF:
M1 implemented the model (2)-(14) of the following paper: O. W. Akinbode and K. W. Hedman "Fictitious losses in the DCOPF with a piecewise linear approximation of losses," IEEE PES General Meeting, Jul. 2013.

Model 2 (M2) for lossy DCOPF:
M2 is an iterative Loss approximation based model. For each iteration: it calculates total losses and then the change in losses than previous iteration; and then assign the losses to all/selected buses with same/different participation factors. This method converges when the difference of total losses between two consecutive iterations is less than a threshold. Possibly, two iterations repeat themselves, which cause method to diverge.

The following paper provides more details about Representation of Transmission Losses (how to assign losses to buses): 

[Xingpeng Li and Kory W. Hedman, “Enhanced Energy Management System with Corrective Transmission Switching Strategy— Part I: Methodology,” IEEE Transactions on Power Systems, vol. 34, no. 6, pp. 4490-4502, Nov. 2019. (DOI: 10.1109/TPWRS.2019.2922880)](https://ieeexplore.ieee.org/document/8736407)

## Contact:
Dr. Xingpeng Li

University of Houston

Email: xli83@central.uh.edu

Website: https://rpglab.github.io/


## License:
This work is licensed under the terms of the Creative Commons Attribution 4.0 (CC BY 4.0) license. 
https://creativecommons.org/licenses/by/4.0/


## Disclaimer:
The author doesn’t make any warranty for the accuracy, completeness, or usefulness of any information disclosed; and the author assumes no liability or responsibility for any errors or omissions for the information (data/code/results etc) disclosed.
