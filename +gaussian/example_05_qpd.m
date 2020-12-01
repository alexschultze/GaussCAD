% example 05
% -Demonstrate QPD Masking and DWS evaluation
% -Demonstrate Beam Offset to Length Effects
% -Demonstrate Tilt to Length Effects

% This test superimposes two gaussian beams and observes
% the intereference on the photodiode (class field screen) element.
% The beams are taken from an actual interferometer setup.
% A. Schultze 01/10/2020 (GaussCAD toolbox)



close all;

lambda = 1064e-9;
w0     = 390e-6; % collimator output w0= [390/400]e-6, z=[-740,750 ] e-3



%Beam 1 (LO)
w0= 400e-6;
zd= -150e-3;
zd=-800e-3;
beam1= gaussian.gauss_beam([zd 0 0],[1 0 0],w0 ,1, lambda );


%Beam 2 (MEAS) (Reflection from Center - Wavefront matching)
w0= 400e-6;
zd= -800e-3;  %from optocad
beam2= gaussian.gauss_beam([zd 0 0],[1 0 0],w0, 1, lambda );


gscreen = gaussian.field_screen([0.0 0 0],[0 0 0],[1 1]*2e-3,[256 256]);
gscreen.rotang=0*4*pi/180;
gscreen.rotax =[0 0 1];

gscreen.set_mask_round();
gscreen.set_mask_gap(100e-6);

gscreen.plot_mask();

gscreen.add_beam( beam1 );
gscreen.add_beam( beam2 );

gscreen.render();

[ac,dc]=gscreen.calc_contrast();

disp(['Max Contrast [%]' num2str(100*ac./(ac+dc))]);



% Parameter study 1. Beam1 has a offset. Coupling of Offset into Pathlengths.
par = linspace(0,100e-6,20);
p0=beam1.p;

clear ph q_ph dws;
for i=1:length(par)
    beam1.p=p0+[0 par(i) 0 ];
    gscreen.render();
    [~,this_ph]=gscreen.calc_int_phase();
    ph(i)=this_ph;    
    [~,this_q_ph]=gscreen.calc_int_phase_quadrants();
    q_ph(:,i)=this_q_ph;
    [dws(i,1),dws(i,2)]=gscreen.calc_dws();
    fprintf('*');
end
%figure();
%gscreen.plot_interference();

figure();
subplot(3,1,1);
plot(par,(ph-max(ph))*lambda/2*pi);
xlabel('Offset Beam y (m)');
ylabel('Pathlength (m)');
suptitle('GaussCAD Example 05 - QPD Offset-to-Length Coupling');

subplot(3,1,2);
plot(par,(q_ph-max(q_ph))*lambda/2*pi);
xlabel('Offset Beam y (m)');
ylabel('Pathlength (m)');
legend('Q1','Q2','Q3','Q4');

subplot(3,1,3);
plot(par,dws);
xlabel('Offset Beam y (m)');ylabel('DWS Signal (rad)');
legend('y','z');


% Parameter Study 2
% Parameter study A. Beam1 has a angle. Coupling  into Pathlength/DWS.
par = linspace(-1e-3,1e-3,10+1);
clear ph dws;
%beam1.p = p0; % reset to original
for i=1:length(par)
    theta=par(i);
    beam2.n=[cos(theta) sin(theta) 0];

    gscreen.render();
    [~,ph(i)]=gscreen.calc_int_phase;
    [~,ph_q(i,:)]=gscreen.calc_int_phase_quadrants();
    [dws(i,1),dws(i,2)]=gscreen.calc_dws();
    
    %gscreen.plot_interference();
    fprintf('*');

end
%figure();
%gscreen.plot_interference();

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