%% Matrix Optics (ABCD) example 02
% to calculate new Gaussian complex beam parameter (q) with ray optics
% A. Schultze 01/10/2020 (GaussCAD toolbox)


this_bench= paraxial.bench_abcd();
this_bench.add(0.1, paraxial.element('lens',0.05));
this_bench.add(0.14, paraxial.element('lens',-0.05));
this_bench.add(0.28, paraxial.element('screen',0.02));
q=1+0.4i; % initial beam with complex beam parameter z=1 and z_r=0.4
this_bench.plot_gauss(q); % draw it
[q, pos, R]=this_bench.plot_gauss(q); %data

w0=sqrt(imag(q)*1064e-9/pi); 
w = w0.*sqrt(1+(imag(q)./real(q)).^2);
table(pos', w0',w', q',R', 'VariableNames',{'z', 'w0','w', 'q','R'})