% Gauss Beam Example 01

% This example demonstrated the visualisation of a Gaussian field
% in its local (r,z) and cartesian (x,y,z) coordinate system.

% A. Schultze 01/10/2020 (GaussCAD toolbox)


close all;

lambda = 1064e-9;
w0     = 2*lambda;
beam1= gaussian.gauss_beam([0.01 0 0],[1 1 0],w0, 1, lambda );

% Section one, evaluate in R-Z Coordinates (Beam reference)
r = linspace(-lambda*10,lambda*10,1000);
z = linspace(-lambda*100,lambda*100,100);
[rr,zz]=meshgrid(r,z);
E=beam1.calc_field_rz(zz,rr);
figure();
subplot(2,1,1);
surf(zz,rr,abs(E),'EdgeColor','none');
xlabel('r');ylabel('z');zlabel('Magnitude');
title('Intensity and Phase (Beam Coordinates)');
subplot(2,1,2);
surf(zz,rr,angle(E),'EdgeColor','none');
xlabel('z (m)');ylabel('r (m)');zlabel('Phase');


%Phases Plot Represent as "Quiver"
figure();
quiver(zz,rr,real(E),imag(E));
xlabel('z (m)');ylabel('r (m)');
title('Complex Phase (Beam Coordinates)');


%Vector Field Phase
E=beam1.calc_field_rz(z,0);
figure();
plot(z,angle(E));
xlabel('z');ylabel('Phase (rad)');
title('Wavefronts (Beam Coordinates)');

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


%% Section two, evaluate in X Y Z Coordinates (Cartesian) Cross Section
w0 = 100*lambda;
beam2= gaussian.gauss_beam([0.00 0 0],[0 0 1],w0, 1, lambda );

x = linspace(-lambda*500,lambda*500,1000);
y = linspace(-lambda*500,lambda*500,1000);
z = 0.1;
[xx,yy,zz]=meshgrid(x,y,z);
points = [xx(:), yy(:) zz(:)];
E=beam2.calc_field_xyz(points);
EE = reshape(E,size(xx));

EE_unwrap = lib.unwrap_2d_fast.unwrap_phase(angle(EE));

figure();
subplot(3,1,1);
surf(xx,yy,abs(EE),'EdgeColor','none');
xlabel('x (m)');ylabel('y (m)');zlabel('Magnitude');
title('Intensity and Phase (Cartesian Coordinates) Cross Section');
subplot(3,1,2);
surf(xx,yy,angle(EE),'EdgeColor','none');
xlabel('x (m)');ylabel('y (m)');zlabel('Phase');

subplot(3,1,3);
surf(xx,yy,EE_unwrap,'EdgeColor','none');
xlabel('x (m)');ylabel('y (m)');zlabel('Phase (Cont)');



%% Section three, evaluate in X Y Z Coordinates (Cartesian)
beam1= gaussian.gauss_beam([0.01 0 0],[1 1 0],w0, 1, lambda );
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
xlabel('x (m)');ylabel('y (m)');zlabel('Magnitude');
title('Intensity and Phase (Cartesian Coordinates)');
subplot(2,1,2);
surf(xx,yy,angle(EE),'EdgeColor','none');
xlabel('x (m)');ylabel('y (m)');zlabel('Phase');


