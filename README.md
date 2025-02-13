# GaussCAD / v1.1
A simulator for evaluating Gaussian optics and simulate projection onto a photodiode. Visual representation of intensity/phase fields as described by Gaussian optics.
https://www.doi.org/10.13140/RG.2.2.30260.27529/1

A. Schultze
14/02/2025

A simulator for evaluating Gaussian optics and simulate projection onto
a photodiode. Visual representation of intensity/phase 
fields as described by Gaussian optics.
It makes uses of Paraxial optics to propagate Gaussian Rays through optical elements. A simple Paraxial optics tracer functionality is included as well.

Features (Gaussian):
-Symmetrical Gaussian Beams of TEM00 mode
-Calculation of interference on a SED/QPD photodiode detector.
-Determination of Contrast (AC,DC), phase angle and DWS on detector.
-Detector shape (round, gap, per quadrant evaluation) supported.

Features (ABCD - Ray Transfer Matrix Analysis)
-Propagation according to Ray transfer matrix analysis (Paraxial approximation)
-Propagation of Gaussian beams with ABCD matrices
-Elementar optical elements (lenses, curved surfaces, propagation)

It does not feature:
-Modes besides TEM00
-Raytracing
-Polarisation

Tested with Matlab v2023.

Changelog:
----------
v0.3: - Fix Field Phase Calculation (fixes DWS)
v1.1: - Fix width bug for ray transfer matrix for gaussian beam
