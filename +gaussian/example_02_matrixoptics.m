% GaussCad Example 02

% to calculate new Gaussian complex beam parameter (q) with ray optics
% A. Schultze 01/10/2020 (GaussCAD toolbox)

% All paraxial approximation calculations (ABCD) are performed on axis.

% This example mixes Gaussian Domain with Paraxial Domain for 
% Gaussian Beam Propagation

%% Part 1: Calculate new rayleigh parameters based on propagation
w0=320e-6;
lambda=1064e-9;
zd=-0.7; %zd = z-z0 distance from beam waist

beam1= gaussian.gauss_beam([0 0 0],[1 0 0],w0, 1, lambda );
lens= paraxial.element('lens',0.1);

q1=beam1.calc_q(zd); % q1=(zd)+1i*zr;
q2= lens.propagate_gauss(q1);
zd2=real(q2) %new (z-zd)
zr2=imag(q2) %new rayleigh


%% Part 2. Use real rays and a lens 100mm
%% Lens focusing gauss ray
% ><---------|--------------<o
% -----------|><---()-------<o
% |  - Screen
% >< - Beam Waist
% <o  - Beam origin

%LO1
w0= 400e-6;
zd= -150e-3; %@QPD
beam1= gaussian.gauss_beam([zd 0 0],[1 0 0],w0 ,1, lambda );
q_old= beam1.calc_q(zd);
q_new= lens.propagate_gauss(q_old);
beam11 =beam1.transform_beam(q_old, q_new);

%MEAS1 (Reflection from Center - Wavefront matching but not perfect)
w0= 400e-6;
zd= -800e-3; %@QPD 
beam2= gaussian.gauss_beam([zd 0 0],[1 0 0],w0, 1, lambda );
q_old= beam2.calc_q(zd);
q_new= lens.propagate_gauss(q_old);
beam22 =beam2.transform_beam(q_old, q_new);

gscreen = gaussian.field_screen([-0.1 0 0],[0 0 0],[1 1]*2e-3,[128 128]);
gscreen.set_mask_round();

gscreen.plot_mask();

gscreen.add_beam( beam11 );
gscreen.add_beam( beam22 );
gscreen.render();

[ac,dc]=gscreen.calc_contrast();
disp(['Max Contrast [%]' num2str(100*ac./(ac+dc))]);


figure();gscreen.plot_intensity();
figure();gscreen.plot_interference();
