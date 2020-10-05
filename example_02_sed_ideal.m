% Example 02 - Field Screen Test
% This example uses two ideal, similar beams and creates a video
% of the interference on the SED.
% A. Schultze 01/10/2020 (GaussCAD toolbox)

close all;

lambda = 1064e-9;
w0     = 10*lambda;
beam1= gauss_beam([0 -0e-4 0],[1 0 0],w0, 1, lambda );
beam2= gauss_beam([0 0e-4 0],[1 0 0],w0, 1, lambda );

gscreen = field_screen([0.01 0 0],[0 0 0],[1 1]*1e-3,[512 512]);

gscreen.set_mask_round();
gscreen.add_beam( beam1 );
gscreen.add_beam( beam2 );

gscreen.render();
figure(); gscreen.plot_intensity();
figure(); gscreen.plot_interference();

return;
%Make a video of interference on screen
figure();
v = VideoWriter('sed_interference.avi','Motion JPEG AVI');
v.FrameRate = 10;
open(v);
par = linspace(0,4*pi,36);
for i = 1:length(par)
   fprintf('*');
   beam1.phase0=par(i);
   gscreen.plot_interference();
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v)