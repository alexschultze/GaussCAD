%% field screen example 03

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



gscreen = gaussian.field_screen([0.0 0 0],[0 0 0],[1 1]*2e-3,[512 512]);
gscreen.rotang=4*pi/180;
gscreen.rotax =[0 1 0];

gscreen.set_mask_round();
gscreen.add_beam( beam1 );
gscreen.add_beam( beam2 );

gscreen.render();
figure();
gscreen.plot_intensity();
figure();
gscreen.plot_interference();

[ac,dc,phi,int]=gscreen.calc_contrast();
figure();
plot(phi,int);
disp(['Max Contrast [%]' num2str(100*ac./(ac+dc))]);
xlabel('Intensity');ylabel('Phase (rad)');

%% Section two, evaluate in X Y Z Coordinates (Cartesian)
x = linspace(-lambda*100000,lambda*100000,1000);
y = linspace(-lambda*1000,lambda*1000,1000);
z = 0;
[xx,yy,zz]=meshgrid(x,y,z);
points = [xx(:), yy(:) zz(:)];
E=beam1.calc_field_xyz(points)+beam3.calc_field_xyz(points);
EE = reshape(E,size(xx));

figure();
subplot(2,1,1);
surf(xx,yy,abs(EE),'EdgeColor','none');
xlabel('x');ylabel('y');zlabel('Magnitude');


% Parameter study. Beam1 has a offset
if(0)
    p0=beam1.p;
    for i=1:20;
    beam1.p=p0+[0 (i-10)*5e-6 0 ];
    gscreen.render();
    [~,this_ph]=gscreen.calc_int_phase();
    ph(i)=this_ph;
    end
    figure();
    plot((ph-max(ph))*lambda/2*pi);
end


%% Make a video of interference on screen

figure();

v = VideoWriter('sed_interference_b.avi','Motion JPEG AVI');
v.FrameRate = 10;
open(v);
for i = 1:36
   %disp(i/);
   fprintf('*');
   ph = 4*pi*i/36;
   beam1.phase0=ph;
   gscreen.render();
   gscreen.plot_interference();
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v)