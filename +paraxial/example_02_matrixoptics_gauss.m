%% Matrix Optics (ABCD) example 02
% to calculate new Gaussian complex beam parameter (q) with ray optics
% A. Schultze 01/10/2020 (GaussCAD toolbox)



%Part 1. Trace a Gaussian Beam through a bench.
this_bench= paraxial.bench_abcd(0.01,0);
this_bench.add(0.1, paraxial.element('lens',0.05));
this_bench.add(0.14, paraxial.element('lens',-0.05));
this_bench.add(0.28, paraxial.element('screen',0.02));
this_bench.plot_gauss(-4.8+1i*0.45); % draw it
[q, pos, R]=this_bench.plot_gauss(-4.8+1i*0.45); %data

w0=sqrt(imag(q)*1064e-9/pi);
w = w0.*sqrt(1+(imag(q)./real(q)).^2);
table(pos', w0',w', q',R', 'VariableNames',{'z', 'w0','w', 'q','R'})

% Part 2.  Parameter study impact of original beam parameter distance

figure();
clear z_new this_qs;
%z0=wo^2*pi/lambda
z=-5:0.1:5;
for i_z=1:length(z)
        this_qs=this_bench.plot_gauss(z(i_z)+1i*0.30);
        z_new(i_z)=this_qs(end);
end
figure();
plot(z,imag(z_new));
xlabel('Distance z (m)');
ylabel('New zr');
title('Gauss Beam (zr orig=0.3)');
