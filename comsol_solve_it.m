function [ filename ] = comsol_solve_it( filename, wavelength,num_modes,directory )
%FCTNEIGENSOLVER Summary of this function goes here
%   Detailed explanation goes here
    format long;
    import com.comsol.model.*
    import com.comsol.model.util.*
    % filename = 'disclinationN1024lx_0.168ly_0.12k_1c_0.7854';
    model = mphload([filename,'.mph']);
    % a = 0.3;
    % wavelength = 0.9;
    %n_back = 2.83;
    model.param.set('wavelength',num2str(wavelength));
    %model.param.set('PMLorder',num2str(n_back));
    %model.param.set('n',num2str(n_back));
    model.study.create('std1');
    model.study('std1').feature.create('eig', 'Eigenfrequency');
    model.study('std1').feature('eig').set('geomselection', 'geom1');
    model.study('std1').feature('eig').set('physselection', 'emw');
    model.study('std1').feature('eig').set('neigs', num2str(num_modes));
    model.study('std1').feature('eig').set('shift', 'c_const/wavelength');
    model.study('std1').feature('eig').set('geomselection', 'geom1');
    model.study('std1').feature('eig').set('physselection', 'emw');
    
    model.sol.create('sol1');
    model.sol('sol1').study('std1');
    model.sol('sol1').feature.create('st1', 'StudyStep');
    model.sol('sol1').feature('st1').set('study', 'std1');
    model.sol('sol1').feature('st1').set('studystep', 'eig');
    model.sol('sol1').feature.create('v1', 'Variables');
    model.sol('sol1').feature.create('e1', 'Eigenvalue');
    model.sol('sol1').feature('e1').set('shift', 'c_const/wavelength');
    model.sol('sol1').feature('e1').set('neigs', num_modes);
    model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
    model.sol('sol1').feature('e1').set('control', 'eig');
    model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
    model.sol('sol1').attach('std1');

    model.sol('sol1').runAll;
    % filename = strcat(filename,'solvedwavelength',num2str(wavelength),'.mph');
    mphsave(model,[filename,'_solved']);
    %command = ['mv scratch_file_solved.mph ',directory,'scratch_file_solved.mph'];
    %system(command);

end