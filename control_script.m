%%%% To run this, from a console type: sudo comsol42 server matlab
%%%% then run it as if you were in real matlab

close all; clear all; clc;

%%%% Did you already run COMSOL?
DO_COMSOL = 0; % 1 = do comsol, 0 = skip comsol.

%%%% Do you want to run adaptive pumping?
DO_ADAPT = 1; % 1 = do adaptive pumping, 0 = skip.

%%%% Do you want to find the condensed threshold?
DO_CONDEN = 0; % 1 = find condensed threshold, 0 = skip.

%%%% Here we will set the system parameters:

R = 4;        % Largest distance from the center of the cavity to
              % the edge in um.

n_inside = 3.5; % index of refraction inside the cavity. 
lambda_a = 1; % wavelength of the atomic resonance transition in
              % um.

gamma_perp_length = .01; % width of the gain curve in um.

directory = '~/Data/2d_salt/Dcav50_R4_rp5_dr20/comsol_new/';
              % Directory to save results to. Make sure to include
              % the final '/'.

Q_thresh = 800; % minimum Q value for modes to save.
num_modes = 100; % number of modes to solve for from COMSOL.

angular_resolution = 360; % COMSOL angular resolution.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% choose the system geometry:

phi = 0:2*pi/angular_resolution:2*pi;

%% for D shaped cavity:
geom_switch = 'D';
flat_position = 0.5; % units of radius (1 is a circle, 0 a semi-circle)

%% for Quadrupole cavity:
%geom_switch = 'Quad';
%epsilon = 0.11; % deformation parameter.

%% for Elliptical cavity:
%geom_switch = 'Ellipse';
%aa = 5; % length of semi-major axis in um.
%bb = 4; % length of semi-minor axis in um.
%assert((R >= aa) && (R >= bb));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for Adaptive and Condensed pumping methods:

numModes = 1; 
% FOR numModes >= 2: look for configuration that optimizes pump for
% this many modes.
% FOR numModes = 1: look for a pump configuration that promotes
% single-mode behavior while maintaining a constant intensity.

numR = 1; % number of 'breaks' in R.
numTH = 3; % number of 'breaks' in theta.

% we'll eventually construct the cavity as linspace(0,R,numR+2), so
% if you put in numR=0, everything will still work, and it means
% that you don't want to break the pump up in the radial
% direction.


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
                                      R, geom_switch, geom_element);

end
    
%%%%%%% BEGIN SPASALT %%%%%%%%

spasalt_setup(directory, R, n_inside, Q_thresh, geom_switch, geom_element);
spasalt_calc(directory, R, lambda_a, Q_thresh, gamma_perp_length, geom_switch, geom_element);

%%%%%%% BEGIN ADAPTIVE PUMPING %%%%%%%%

%assert(geom_switch=='D'); % adaptive pumping should work for all
% geometries now.

if (DO_ADAPT == 1)

    if ((numR==0)&&(numTH==0))
        error('no free pump variables.');
    end

    lb = zeros((numR+1)*(numTH+1)-1,1);
    %ub = ones((numR+1)*(numTH+1)-1,1);

    options = gaoptimset('Generations',50,'PopulationSize',(numR+1)*(numTH+1),'EliteCount',2);

    optPumpVec = ga(@(x)spasalt_adaptive(x, directory, R, lambda_a, Q_thresh, ...
                                         gamma_perp_length, ...
                                         numModes, numR, numTH, 0), ...
                    (numR+1)*(numTH+1)-1, [], [], [], [], lb, [], [], ...
                    [], options);


    %load([directory,'spasalt_adaptive.mat'],'optPumpVec');

    spasalt_adaptive(optPumpVec, directory, R, lambda_a, Q_thresh, ...
                     gamma_perp_length, numModes, numR, numTH, 1);

end

%%%%%%% BEGIN CONDENSED PUMPING %%%%%%%%

if (DO_CONDEN == 1)
    spasalt_condensed(directory, R, lambda_a, Q_thresh, gamma_perp_length, ...
                      geom_switch, geom_element, numModes);
end