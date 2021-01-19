%% Matrix Optics (ABCD) example 11
% to calculate new Gaussian complex beam parameter (q) with ray optics
% A. Schultze 01/10/2020 (GaussCAD toolbox)

% All paraxial approximation calculations (ABCD) are performed on axis.

% Setup: Lens , Spherical Surface, Lens (Double Pass)

f=0.10;
r=0.02;
d= 0.099;

% Without lenses, normal propagation
this_bench= paraxial.bench_abcd();
this_bench.add(0.5, paraxial.element('screen',0.02));
[q, pos, R]=this_bench.plot_gauss(-4.8+1i*0.45);

w0=sqrt(imag(q)*1064e-9/pi);
w = w0.*sqrt(1+(imag(q)./real(q)).^2);
table(pos', w0',w', q',R', 'VariableNames',{'z', 'w0','w', 'q','R'})

%This are beam parameters at screen by default without any optical
%components
w_orig=w(end);
R_orig=R(end);

this_bench= paraxial.bench_abcd();
this_bench.add(0.2-d, paraxial.element('lens',f));
this_bench.add(0.2, paraxial.element('mirror_curved',-r));
this_bench.add(0.2+d, paraxial.element('lens',f));
this_bench.add(0.5, paraxial.element('screen',0.02));
this_bench.plot([0.01;0]);

figure();
this_bench.plot_gauss(-4.8+1i*0.45);

[q, pos, R]=this_bench.plot_gauss(-4.8+1i*0.45);
w0=sqrt(imag(q)*1064e-9/pi);
table(pos', w0', q',R', 'VariableNames',{'z', 'w0', 'q','R'})


figure();
clear z_new this_qs;
%z0=wo^2*pi/lambda
d=0.07:0.0002:0.11;
q0=-4.8+1i*0.45;
for i_d=1:length(d)
        this_bench.elements_pos(1) = 0.2-d(i_d);
        this_bench.elements_pos(3) = 0.2+d(i_d);
        this_qs=this_bench.plot_gauss(q0);
        q_new(i_d)=this_qs(end); %final beam parameter
end
figure();
subplot(3,1,1)
plot(d,imag(q_new));hold on;
yline(imag(q0),'--b');
legend('zr (new)','zr (old)');
xlabel('Distance d (Lens-Mirror) (m)');
ylabel('New zr');
title('Gauss Beam (zr orig=0.45)');

subplot(3,1,2)
R=this_bench.gauss_curvature_q(q_new);
plot(d,R);hold on;
yline(R_orig);
xlabel('Distance d (Lens-Mirror) (m)');
ylabel('New R (m) at Screen');


subplot(3,1,3)
w0=sqrt(imag(q_new)*1064e-9/pi);
w = w0.*sqrt(1+(imag(q_new)./real(q_new)).^2);
plot(d,w);hold on;
yline(w_orig);
xlabel('Distance d (Lens-Mirror) (m)');
ylabel('Beam Width (m) at Screen');
