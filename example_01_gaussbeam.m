%% Gauss Beam Example 01

% This example demonstrated the visualisation of a Gaussian field
% in its local (r,z) and cartesian (x,y,z) coordinate system.

% A. Schultze 01/10/2020 (GaussCAD toolbox)


close all;

lambda = 1064e-9;
w0     = 2*lambda;
beam1= gauss_beam([0.01 0 0],[1 1 0],w0, 1, lambda );

% Section one, evaluate in R-Z Coordinates (Beam reference)

r = linspace(-lambda*10,lambda*10,1000);
z = linspace(-lambda*100,lambda*100,100);
[rr,zz]=meshgrid(r,z);
E=beam1.calc_field_rz(zz,rr);
figure();
subplot(2,1,1);
surf(zz,rr,real(E),'EdgeColor','none');
xlabel('r');ylabel('z');zlabel('Magnitude');
subplot(2,1,2);
surf(zz,rr,angle(E),'EdgeColor','none');
xlabel('z');ylabel('r');zlabel('Phase');
%Phases Plot Represent as "Quiver"

figure();
quiver(zz,rr,real(E),imag(E));
xlabel('z');ylabel('r');
title('Complex Phase');

%Check Guoy Phase
E=beam1.calc_field_rz(z,0);
figure();
plot(z,angle(E));
xlabel('z');ylabel('Phase (rad)');


%Phases Plot Represent as "Wavefront"
r = linspace(-lambda*10,lambda*10,2000);
z = linspace(-lambda*10,lambda*10,2000);
[rr,zz]=meshgrid(r,z);
E=beam1.calc_field_rz(zz,rr);
indx_wf=find(abs(angle(E))<0.05);
scatter(zz(indx_wf),rr(indx_wf));
xlabel('z');ylabel('r');zlabel('Wavefront');

%Phases Plot Represent as "Wavefront with Intensity"
r = linspace(-lambda*10,lambda*10,2000);
z = linspace(-lambda*10,lambda*10,2000);
[rr,zz]=meshgrid(r,z);
E=beam1.calc_field_rz(zz,rr);
indx_wf=find(abs(angle(E))<0.05);
plot3(zz(indx_wf),rr(indx_wf),real(E(indx_wf)),'.');
xlabel('z');ylabel('r');zlabel('Wavefront Intensity');


%% Section two, evaluate in X Y Z Coordinates (Cartesian)
x = linspace(-lambda*100,lambda*100,1000)+0.01;
y = linspace(-lambda*100,lambda*100,1000);
z = 0;
[xx,yy,zz]=meshgrid(x,y,z);
points = [xx(:), yy(:) zz(:)];
E=beam1.calc_field_xyz(points);
EE = reshape(E,size(xx));

figure();
subplot(2,1,1);
surf(xx,yy,abs(EE),'EdgeColor','none');
xlabel('r');ylabel('z');zlabel('Mhagnitude');
subplot(2,1,2);
surf(xx,yy,angle(EE),'EdgeColor','none');
xlabel('z');ylabel('r');zlabel('Phase');



%% Check how it looks like when the phase is shifting for some point
beam1= gauss_beam([0.0 0 0],[1 1 0],w0, 1, lambda );
beam2= gauss_beam([0.0 0 0],[1 1 0],w0, 1, lambda );

point = [0*lambda 0 0];
ph=linspace(0,2*pi,20);
for i=1:length(ph)
    beam1.phase0=ph(i);
    E1=beam1.calc_field_xyz(point);
    E2=beam2.calc_field_xyz(point);
    Eph(i) =1/2*abs(E1+E2).^2;
end
figure();
plot(ph,Eph);
xlabel('Phase Offset (rad)');
ylabel('Intensity at [p]');

