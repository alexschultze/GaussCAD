%% Matrix Optics (ABCD) example 02
% to calculate new Gaussian complex beam parameter (q) with ray optics
% A. Schultze 01/10/2020 (GaussCAD toolbox)

q=-1+0.5i; % initial beam with complex beam parameter z=1 and z_r=0.4
% The beam has a rayleigh length of 0.4 m with a beam waist at 1m

this_bench= paraxial.bench_abcd(533e-9); %new bench with 533nm light

this_bench.add(0, paraxial.element('screen'));
this_bench.add(0.5, paraxial.element('screen'));
this_bench.add(1, paraxial.element('screen'));
this_bench.add(1.5, paraxial.element('screen'));
this_bench.add(10, paraxial.element('screen'));



this_bench.plot_gauss(q); % draw it
[q, pos, R, w,w0]=this_bench.plot_gauss(q); %data

table(pos', w0',w', q',R', 'VariableNames',{'z', 'w0','w', 'q','R'})