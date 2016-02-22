%%%% To run this, from a console type: sudo comsol42 server matlab
%%%% then run it as if you were in real matlab

close all; clear all; clc;

%%%% Did you already run COMSOL?
DO_COMSOL = 1; % 1 = do comsol, 0 = skip comsol.

%%%% Here we will set the system parameters:

R = 5;        % Largest distance from the center of the cavity to
              % the edge in um.
n_inside = 3.5; % index of refraction inside the cavity. 
lambda_a = 1; % wavelength of the atomic resonance transition in
              % um.
gamma_perp_length = .01; % width of the gain curve in um.


directory = '~/acerjan/comsol_results/ellipse_R5um_test2/';
              % Directory to save results to. Make sure to include
              % the final '/'.
Q_thresh = 500; % minimum Q value for modes to save.
num_modes = 100; % number of modes to solve for from COMSOL.

angular_resolution = 360; % COMSOL angular resolution.


%%%% and the system geometry:

phi = 0:2*pi/angular_resolution:2*pi;

%% for D shaped cavity:
%geom_switch = 'D';
%flat_position = 0.5; % units of radius (1 is a circle, 0 a semi-circle)

%% for Quadrupole cavity:
%geom_switch = 'Quad';
%epsilon = 0.11; % deformation parameter.

%% for Elliptical cavity:
geom_switch = 'Ellipse';
aa = 5;
bb = 4;


%%%%%%% BEGIN COMSOL %%%%%%%%
%% Don't touch things in here.

switch geom_switch
  case 'D'
    assert((flat_position>=0)&&(flat_position<=1));

    geom.x_coords = min(R.*cos(phi),R*flat_position);
    geom.y_coords = R.*sin(phi);
    geom_element = flat_position;
    
  case 'Quad'
    assert(epsilon>=0);

    r0 = R/(1+epsilon);
    geom.x_coords = r0*(1+epsilon*cos(2*phi)).*cos(phi);
    geom.y_coords = r0*(1+epsilon*cos(2*phi)).*sin(phi);
    geom_element = epsilon;
    
  case 'Ellipse'
    geom.x_coords = aa*cos(phi);
    geom.y_coords = bb*sin(phi);
    geom_element(1) = aa;
    geom_element(2) = bb;
end    
    
geom.n_eff = n_inside;
geom.wavelength = lambda_a;

geom.system_size = 2*R + 2*lambda_a;

if (DO_COMSOL == 1)
    comsol_gen_geometry(geom,directory);

    %%
    tic
    comsol_solve_it(['scratch_file'],geom.wavelength,num_modes,directory);
    toc

    %%
    xx=-(R+1):0.01:(R+1);yy=xx;
    dlmwrite([directory, 'grid_xy'],[xx; yy]);

    [Q lambda] = comsol_extract_field('scratch_file_solved',Q_thresh,num_modes, ...
                                      directory,[num2str(-R-1),[',.01,'],num2str(R+1)],... 
                                      geom_switch, geom_element);

end
    
%%%%%%% BEGIN SPASALT %%%%%%%%

spasalt_setup(directory, R, n_inside, Q_thresh, geom_switch, geom_element);
spasalt_calc(directory, R, lambda_a, Q_thresh, gamma_perp_length, geom_switch, geom_element);

