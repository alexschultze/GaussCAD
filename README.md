# GaussCAD / v0.3
A simulator for evaluating Gaussian optics and simulate projection onto a photodiode. Visual representation of intensity/phase fields as described by Gaussian optics.

DOI: 10.13140/RG.2.2.30260.27529

A. Schultze
30/09/2020

A simulator for evaluating Gaussian optics and simulate projection onto
a photodiode. Visual representation of intensity/phase 
fields as described by Gaussian optics.
It makes uses of Paraxial optics to propagate Gaussian Rays through optical elements. A simple Paraxial optics tracer functionality is included as well.

Features:
-Symmetrical Gaussian Beams of TEM00 mode
-Calculation of interference on a SED/QPD photodiode detector.
-Determination of Contrast (AC,DC), phase angle and DWS on detector.
-Detector shape (round, gap, per quadrant evaluation) supported.
-Propagation according to Matrix Optics (Paraxial approximation) of Gauss beam

It does not feature:
-Modes besides TEM00
-Raytracing
-Polarisation

Tested with Matlab v2018b.

Changelog:
----------
v3: - Fix Field Phase Calculation (fixes DWS)
