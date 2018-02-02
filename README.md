# Multiband RF pulse design for realistic gradient performance

## Content
This repository contains code to produce multiband VERSE (MBv) and VERSE multiband RF and gradient pulses.
The VERSE code used is an modified version version of https://github.com/mriphysics/reVERSE-GIRF by Shaihan Malik, which in turn
is an implementation of Time-optimal VERSE pulse design [(Lee et al, 2009)](http://doi.org/10.1002/mrm.21950). This has been implemented by modifying code released by Miki Lustig for designing time-optimal gradients, publically available on the  [authors' website](http://www.eecs.berkeley.edu/~mlustig/Software.html). 

Therefore, acknowledgements to the time-optimal gradient framework should be attributed to [Lustig IEEE-TMI 2008](https://www.ncbi.nlm.nih.gov/pubmed/18541493) whilst acknowledgement to time-optimal VERSE goes to Lee et al 2009.

## Purpose
Multiband RF pulses are an essential building block in multiband/simultaneous multi-slice imaging sequences. However conventional RF 
pulse design strategies lead to long RF pulse durations, which leads to long echo-times and thus lower SNR. One way of significantly reducing
pulse durations is to use time-variable selection gradients, like those used in the ISMRM pulse design challenge in 2017 (1). 

However, time-optimal design methods assume RF and gradient systems have similar output bandwidth, but gradient systems actually have much lower temporal bandwidth.
When ignored, this can lead to slice distortions and image artefacts.

## Theory
TBD

## Results
TBD

## Contact

02/02/2018

For contact:

Samy Abo Seada - Email:
samy.abo_seada@kcl.ac.uk

Shaihan Malik - Email:
shaihan.malik@kcl.ac.uk

King's College London, 2018.