%% Matrix Optics (ABCD) example 06
% to calculate new Gaussian complex beam parameter (q) with ray optics
% A. Schultze 01/10/2020 (GaussCAD toolbox)

% All paraxial approximation calculations (ABCD) are performed on axis.

this_bench= paraxial.bench_abcd(0.01,0);
this_bench.add(0.1, paraxial.element('lens',0.05));
this_bench.add(0.14, paraxial.element('mirror_curved',-0.02));
this_bench.add(0.28, paraxial.element('screen',0.02));
this_bench.plot();