README.txt
Alex Cerjan

Welcome to the readme file for using the single pole approximation of the
steady-state ab initio laser theory (SPA-SALT) with cavity wave functions
generated using COMSOL multiphysics. First and foremost, I'm here to help
you use this code. If something doesn't make sense, or the code isn't running,
please send me an email: alexcerjan@gmail.com, and I'd be happy to help
(without any need for any sort of authorship status).

If you use this code for academic work, please cite the original SALT
paper, the SPA-SALT paper, the papers demonstrating that SALT can be
correctly applied to realistic gain media, and the paper corresponding 
to this work demonstrating the correspondence between COMSOL solutions 
and SALT solutions. (All of the bibtex entries for these are given below.)



Version: v0.0 Alpha release, 5/13/15

Congratulations, you've downloaded my first release of this repository of
code. This is riddled with bugs, that I'll be ironing out in early June.


How to use:

For the moment, I have not yet written a full control script for this code,
this will be corrected very soon, but again, just email me for anything that
doesn't make sense.

1) 
control_script.m - defines the cavity for COMSOL, runs the COMSOL simulations
and then extracts the field profiles. Run using the comsol scripting command,
such as:
sudo comsol42 server matlab
then in the MATLAB program run control_script.

In this script, you'll need to set:

geom.n_eff = the effective index of refraction inside your structure.
geom.wavelength = the wavelength of your system, in um.
radius = maximum distance of your structure from the center point, in um.
geom.system_size = pick a number slightly larger than 2*R.

you'll then need to define the boundary of your cavity in Cartesian coordinates,
setting geom.x_coords and geom.y_coords. Script is currently set for
simulating a D-shaped cavity.

num_modes = number of modes you want COMSOL to solve for, in order of
highest to lowest Q.
Q_thresh = the minimum Q of a mode you want to save.
folder = directory for files generated in the simulation.

2) 
Once you've run COMSOL, you'll want to run comsol_spasalt.m from MATLAB.
In this code, you'll need to give as arguments folder, from above, and N,
which is the number of modes above threshold you want to solve for.

You'll also need to set:
lambda_a = geom.wavelength from before.
gammaPerpLambda = width of the gain curve in wavelength.
R = radius from before.
Q_thresh = same as above.

Finally, you'll again need to input the cavity boundary, near line 54.
cavityLocs(xii,yii) = 1 for points within the cavity, and 0 outside.
For a uniform pump across the entire cavity, pumpProfile = cavityLocs.
Otherwise, you should set pumpProfile(xii,yii) = 1 where you want 
to pump the cavity.

3)
Now that you have your SPA-SALT overlap integrals, calculated in step (2),
you're basically done. Either run lambdaG_gen.m or spasalt_inten_calc.m
to generate the generalize mode competition parameter lambda, or plots
of intensity as a function of pump as lasing modes are predicted to
turn on. You'll need to again specify data directory, folder from above,
and N, the number of modes above threshold you want to look at.


BIBTEX entries:
@article{tureci06,
author = {H. E. T\"{u}reci and A. D. Stone and B. Collier},
journal = {Phys. Rev. A},
pages = {043822},
title = {Self-consistent multimode lasing theory for complex or random lasing media},
volume = {74},
year = {2006},
}

@article{ge10,
author = {L. Ge and Y. D. Chong and A. D. Stone},
journal = {Phys. Rev. A},
pages = {063824},
title = {Steady-state ab initio laser theory: generalizations and analytic results},
volume = {82},
year = {2010},
}

@article{cerjan12,
author = {A. Cerjan and Y. D. Chong and L. Ge and A. D. Stone},
journal = {Opt. Express},
pages = {474-488},
title = {Steady-state ab initio laser theory for N-level lasers},
volume = {20},
year = {2012},
}

@article{cerjan_csalt_2015,
author = {A. Cerjan and Y. D. Chong and A. D. Stone},
title = {Steady-state ab initio laser theory for complex gain media},
volume = {23},
journal = {Opt. Express},
pages = {6455--6477},
year = {2015},
}

@unpublished{cerjan_comsol_salt_2015,
author = {A. Cerjan and B. Redding and H. Cao and A. Douglas Stone},
note = {in preparation},
}
