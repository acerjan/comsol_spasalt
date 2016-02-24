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
