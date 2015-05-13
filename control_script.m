close all; clear all;clc;

%%%% To run this, from a console type: sudo comsol42 server matlab
%%%% then run it as if you were in real matlab

%%%%%%% Set Geometry %%%%%%%%
geom.n_eff = 3.5; %+.002*1i;
geom.wavelength = 1;

radius = 4;
geom.system_size = 10;

% radius = 20;
% geom.system_size = 44;

res=360;
phi = 0:2*pi/res:2*pi;

%epsilon = 0.11;
%r0 = radius/(1+epsilon);
%geom.x_coords = r0*(1+epsilon*cos(2*phi)).*cos(phi);
%geom.y_coords = r0*(1+epsilon*cos(2*phi)).*sin(phi);

flat_position = 0.5; %units of radius (1 is a circle, 0 a semi-circle)
geom.x_coords = min(radius.*cos(phi),radius*flat_position);
geom.y_coords = radius.*sin(phi);

gen_geometry(geom);
%%

num_modes=100;
%%
tic
solve_it(['scratch_file'],geom.wavelength,num_modes);
toc

%%
folder = '~/acerjan/comsol_results/Dcav_R4um_neff3p5/';

Qthresh=200;
[Q lambda] = extract_field('scratch_file_solved',Qthresh,num_modes,folder,[num2str(-radius-1),',.01,',num2str(radius+1)]);

xx=-(radius+1):0.01:(radius+1);yy=xx;
dlmwrite([folder 'grid_xy'],[xx; yy]);

return;
