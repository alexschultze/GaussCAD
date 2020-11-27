%% example 04
% - Demonstrate Tilt to Length Effects

% This test superimposes two gaussian beams and observes
% the intereference on the photodiode (class field screen) element.
% The beams are taken from an actual interferometer setup.
% A. Schultze 01/10/2020 (GaussCAD toolbox)

close all;

lambda = 1064e-9;
w0     = 390e-6;

%LO1
w0= 400e-6;
zd= -150e-3;

beam1= gaussian.gauss_beam([zd+5e-3 0 0],[1 0 0],w0 ,1, lambda );

%MEAS1 (Reflection from Center - Wavefront matching)
beam2= gaussian.gauss_beam([zd 0 0],[1 0 0],w0, 1, lambda );

gscreen = gaussian.field_screen([0.0 0 0],[0 0 0],[1 1]*2e-3,[512 512]);
gscreen.rotang=4*pi/180;
gscreen.rotax =[0 1 0];
gscreen.rotax =[0 0 1];
gscreen.rotax =[0 1 0.1];
gscreen.rotax = gscreen.rotax/norm(gscreen.rotax);

gscreen.set_mask_round();
gscreen.set_mask_gap(100e-6);

figure();
gscreen.plot_mask();




gscreen.add_beam( beam1 );
gscreen.add_beam( beam2 );





gscreen.render();
[ac,dc]=gscreen.calc_contrast();
disp(['Max Contrast [%]' num2str(100*ac./(ac+dc))]);
figure();
gscreen.plot_interference();


% Parameter study A. Beam1 has a angle. Coupling of Offset into Pathlength.
par = linspace(-1e-4,1e-4,10+1);
clear ph dws;
for i=1:length(par)
    theta=par(i);
    beam2.n=[cos(theta) sin(theta) 0];

    gscreen.render();
    [~,ph(i)]=gscreen.calc_int_phase;
    [~,ph_q(i,:)]=gscreen.calc_int_phase_quadrants();
    [dws(i,1),dws(i,2)]=gscreen.calc_dws();
    
    gscreen.plot_interference();
    fprintf('*');

end
figure();
subplot(2,1,1);
plot(par,(ph-max(ph))*lambda/2*pi);
xlabel('Tilt Angle Beam(rad)');ylabel('Path length (m)');
title('GaussCAD- Example 04 - TTL');
subplot(2,1,2);
plot(par,dws);
xlabel('Tilt Angle Beam(rad)');ylabel('DWS Signal (rad)');
legend('y','z');

%Calculate Tilt to Phase Difference (DWS Factor)
dws_fit = polyfit(dws(:,2),par',1);
disp([ 'DWS Factor: ' num2str(dws_fit(1)) ' (rad/rad) (Tilt/Phase Difference)']);

