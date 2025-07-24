# LinearInverseKsn
Use a linear inverse method to calculate the normalized steepness index (Ksn) at a grid of locations over a digital elevation model. Test multiple concavity values to find the one that produces the best fit to the stream profiles.

This is a heavily modified version of the INCA code by Adam Smith.
The original INCA code can be found here: https://github.com/adamsmith142/INCA
The method is described by Smith and Fox (2024, https://doi.org/10.1029/2023JF007584) and Smith et al. (2022, https://doi.org/10.1016/j.earscirev.2022.103970).
The Clearwater DEM used as an example can be found at DOI: 10.17632/pbppkd7gwn.1, and the outlet.kml file to use with it can be found in the INCA repository.

To run, modify the file RunModel.m as desired, and then run it.

This code requires TopoToolbox to run, which can be downloaded here: https://topotoolbox.wordpress.com/
It has only been tested with TopoToolbox 2, and it may or may not work with other versions.
