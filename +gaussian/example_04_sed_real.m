% GaussCad Example 04

% This test superimposes two gaussian beams and observes
% the intereference on the photodiode (class field screen) element.
% The beams are taken from an actual interferometer setup.
% A. Schultze 01/10/2020 (GaussCAD toolbox)

close all;

lambda = 1064e-9;
w0     = 390e-6; % collimator output w0= [390/400]e-6, z=[-740,750 ] e-3

% collimator output w0= [390,400]e-6, z=[-740,750 ] e-3
% collimator qpd2   w0= [78, 150]e-6, z=[-43,-390 ] e-3

%LO1
w0= 400e-6;
zd= -150e-3;
beam1= gaussian.gauss_beam([zd 0 0],[1 0 0],w0, 1, lambda );
%MEAS1 (Reflection from Center - Wavefront matching)
w0= 400e-6;
zd= -800e-3;  %from optocad, max contrast is not great
beam2= gaussian.gauss_beam([zd 0 0],[1 0 0],w0, 1, lambda );

%MEAS2 (Reflection from Surface - Wavefront not matching)
w0= 80e-6;
zd= 43e-3;
beam3= gaussian.gauss_beam([zd 0 0],[1 0 0],w0, 1, lambda );



gscreen = gaussian.field_screen([0.0 0 0],[0 0 0],[1 1]*2e-3,[256 256]);
gscreen.rotang=4*pi/180;
gscreen.rotax =[0 1 0];

gscreen.set_mask_round();
gscreen.add_beam( beam1 );
gscreen.add_beam( beam2 );

gscreen.render();
figure(); gscreen.plot_intensity();
figure(); gscreen.plot_interference();
figure(); gscreen.plot_contrast();

[ac,dc]=gscreen.calc_contrast();
disp(['Max Contrast [%]' num2str(100*ac./(dc+ac))]);



% Parameter study. Beam1 has a offset
par = linspace(-50e-6, 50e-6,21);

p0=gscreen.beams(1).p;
for i=1:length(par)
    gscreen.beams(1).p=p0+[0 par(i) 0 ];
    gscreen.render();
    [~,this_ph]=gscreen.calc_int_phase();
    ph(i)=this_ph;
end
figure();
plot(par,(ph-max(ph))*lambda/2*pi);
xlabel('Beam Displacement Y (m)');
ylabel('Pathlength (m)');

suptitle('GaussCAD Example 04 - SED Offset-to-Length Coupling');


%% Make a video of interference on screen

figure();
v = VideoWriter('example04.avi','Motion JPEG AVI');
v.FrameRate = 10;
open(v);
par = linspace(0,4*pi,36);
for i = 1:length(par)
   fprintf('*');
   
   %Slow Variant
   %beam1.phase0=par(i);
   %gscreen.render();
   
   %Fast Variant
   gscreen.shift_phase_beam1(4*pi/36);
   
   gscreen.plot_interference();
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v)