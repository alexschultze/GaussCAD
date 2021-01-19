%% Matrix Optics (ABCD) example 01
% to calculate simple ray transfer optics example.
% A. Schultze 01/10/2020 (GaussCAD toolbox)

%Define a bench .
this_bench= paraxial.bench_abcd();
%Add (1) lens with f=0.05 @ 0.1m
this_bench.add(0.1, paraxial.element('lens',0.05));
%Add (2) curved mirror with r=0.02 @ 0.14m
this_bench.add(0.14, paraxial.element('mirror_curved',-0.02));
%Add (3) screen for terminating the bench
this_bench.add(0.28, paraxial.element('screen',0.02));
% Trace a  initial ray of 0.01m and alpha=0 rad
this_bench.plot([0.01; 0]);