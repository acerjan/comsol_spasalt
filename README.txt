README.txt
Alex Cerjan

Welcome to the readme file for using the single pole approximation of the
steady-state ab initio laser theory with cavity wave functions
generated using COMSOL multiphysics (resonance SPA-SALT). First and foremost, I'm here to help
you use this code. If something doesn't make sense, or the code isn't running,
please send me an email: alexcerjan@gmail.com, and I'd be happy to help
(without any need for any sort of authorship status).

At the present, the manuscript describing this work is currently in submission.
So, if you currently find yourself using this code for academic work, please cite
the original SPA-SALT paper and the manuscript in submission:

@article{ge10,
author = {L. Ge and Y. D. Chong and A. D. Stone},
journal = {Phys. Rev. A},
pages = {063824},
title = {Steady-state ab initio laser theory: generalizations and analytic results},
volume = {82},
year = {2010},
}

@unpublished{cerjan_resonance_salt,
author = {A. Cerjan and B. Redding and L. Ge and S. F. Liew and H. Cao and A. Douglas Stone},
title = {Controlling mode competition by tailoring the spatial pump distribution in a laser: A resonance-based approach},
note = {in submission},
}

Version: v0.1, 8/25/16

Updated a few small bugs, and this README.

Version: v0.0 Beta release, 2/23/16

Congratulations, you've downloaded my second release of this repository of
code. A great deal of functionality has been added since the original
release, but this new version is also more sensitive to memory management,
and is sadly incompatible with the alpha release, as they store the data
in different ways.


How to use:

control_script.m - defines the cavity for COMSOL, runs the COMSOL simulations
and then extracts the field profiles. Run using the comsol scripting command,
such as:
sudo comsol42 server matlab
then in the MATLAB program run control_script.

This is the script where you will set all of the necessary parameters for the
cavity and the gain medium. If these variables are not clear in the comments
in the script, please let me know. If you need to add a new geometry, this currently
needs to be added in a few different ways in separate locations, please email me and
I'd be happy to assist.
