%% example 05
% -Demonstrate Beam Offset to Length Effects

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
beam1= gauss_beam([zd 0 0],[1 0 0],w0 ,1, lambda );


%MEAS1 (Reflection from Center - Wavefront matching)
w0= 400e-6;
zd= -800e-3;  %from optocad
beam2= gauss_beam([zd 0 0],[1 0 0],w0, 1, lambda );

%MEAS2 (Reflection from Surface - Wavefront not matching)
w0= 80e-6;
zd= 43e-3;
beam3= gauss_beam([zd 0 0],[1 0 0],w0, 1, lambda );



gscreen = field_screen([0.0 0 0],[0 0 0],[1 1]*2e-3,[512 512]);
gscreen.rotang=0*4*pi/180;
gscreen.rotax =[0 1 0];

gscreen.set_mask_round();
gscreen.set_mask_gap(100e-6);

gscreen.add_beam( beam1 );
gscreen.add_beam( beam2 );

gscreen.render();

[ac,dc]=gscreen.calc_contrast();

disp(['Max Contrast [%]' num2str(100*ac./(ac+dc))]);


%% Section two, evaluate in X Y Z Coordinates (Cartesian)
x = linspace(-lambda*100000,lambda*100000,1000);
y = linspace(-lambda*1000,lambda*1000,1000);
z = 0;
[xx,yy,zz]=meshgrid(x,y,z);
points = [xx(:), yy(:) zz(:)];
E=beam1.calc_field_xyz(points)+beam2.calc_field_xyz(points);
EE = reshape(E,size(xx));

figure();
subplot(2,1,1);
surf(xx,yy,abs(EE),'EdgeColor','none');
xlabel('x');ylabel('y');zlabel('Magnitude');


% Parameter study. Beam1 has a offset. Coupling of Offset into Pathlength.
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
end
figure();
subplot(3,1,1);
plot(par,(ph-max(ph))*lambda/2*pi);
xlabel('Offset Beam y (m)');
ylabel('Pathlength (m)');
title('GaussCAD Example 05 - Offset-to-Length');

subplot(3,1,2);
plot(par,(q_ph-max(q_ph))*lambda/2*pi);
xlabel('Offset Beam y (m)');
ylabel('Pathlength (m)');
legend('Q1','Q2','Q3','Q4');

subplot(3,1,3);
plot(par,dws);
xlabel('Offset Beam y (m)');ylabel('DWS Signal (rad)');
legend('y','z');