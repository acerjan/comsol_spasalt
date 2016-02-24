function [ filename ] = comsol_gen_geometry( geom, directory )
%FCTNDISCLINATIONGEN140523 Summary of this function goes here
%   Detailed explanation goes here
    format long;
    import com.comsol.model.*
    import com.comsol.model.util.*

    model = ModelUtil.create('Model');
    model.modelNode.create('mod1');
    model.geom.create('geom1', 2);
    model.mesh.create('mesh1', 'geom1');
    model.physics.create('emw', 'ElectromagneticWaves', 'geom1');

    % N = 1024;
    n_back = geom.n_eff;
    L = geom.system_size;
    wavelength = geom.wavelength;

    PMLt = wavelength;

    % startx = (-sqrt(N)/2+0.5)*a;
    % starty = (-sqrt(N)/2+0.5)*a;
    model.param.set('L', num2str(L));
    model.param.set('PMLt', 'wavelength');
    model.param.set('wavelength',num2str(wavelength));
    model.param.set('PMLorder',num2str(n_back));
    model.param.set('n',num2str(n_back));
    % model.param.set('a',num2str(a));

    model.geom('geom1').feature.create('r1', 'Rectangle');
    model.geom('geom1').feature('r1').setIndex('size', 'L', 0);
    model.geom('geom1').feature('r1').setIndex('size', 'L', 1);
    model.geom('geom1').feature('r1').setIndex('pos',  '-L/2', 0);
    model.geom('geom1').feature('r1').setIndex('pos', '-L/2', 1);
    model.geom('geom1').feature.create('r2', 'Rectangle');
    model.geom('geom1').feature('r2').setIndex('size', 'PMLt', 0);
    model.geom('geom1').feature('r2').setIndex('size', 'L+2*PMLt', 1);
    model.geom('geom1').feature('r2').setIndex('pos', '-L/2-PMLt', 0);
    model.geom('geom1').feature('r2').setIndex('pos', '-L/2-PMLt', 1);
    model.geom('geom1').feature.create('r3', 'Rectangle');
    model.geom('geom1').feature('r3').setIndex('size', 'PMLt', 0);
    model.geom('geom1').feature('r3').setIndex('size', 'L+2*PMLt', 1);
    model.geom('geom1').feature('r3').setIndex('pos', 'L/2', 0);
    model.geom('geom1').feature('r3').setIndex('pos', '-L/2-PMLt', 1);
    model.geom('geom1').feature.create('r4', 'Rectangle');
    model.geom('geom1').feature('r4').setIndex('size', 'L+2*PMLt', 0);
    model.geom('geom1').feature('r4').setIndex('size', 'PMLt', 1);
    model.geom('geom1').feature('r4').setIndex('pos',  '-L/2-PMLt', 0);
    model.geom('geom1').feature('r4').setIndex('pos',  '-L/2-PMLt', 1);
    model.geom('geom1').feature.create('r5', 'Rectangle');
    model.geom('geom1').feature('r5').setIndex('size', 'L+2*PMLt', 0);
    model.geom('geom1').feature('r5').setIndex('size', 'PMLt', 1);
    model.geom('geom1').feature('r5').setIndex('pos',  '-L/2-PMLt', 0);
    model.geom('geom1').feature('r5').setIndex('pos',  'L/2', 1);
    D = [];
    count = 0;

    x_coords = geom.x_coords;
    y_coords = geom.y_coords;

    count = count + 1;
    model.geom('geom1').feature.create('pol1', 'Polygon');
    model.geom('geom1').feature('pol1').set('x', num2str(x_coords));
    model.geom('geom1').feature('pol1').set('y', num2str(y_coords));

    D = [D 9+count];

    model.geom('geom1').runAll;
    model.geom('geom1').run;

    model.material.create('mat1');
    model.material.create('mat2');
    model.material('mat1').propertyGroup('def').set('electricconductivity', {'0'});
    model.material('mat1').propertyGroup('def').set('relpermeability', {'1'});
    model.material('mat1').propertyGroup('def').set('relpermittivity', {'n^2'});
    model.material('mat2').propertyGroup('def').set('electricconductivity', {'0'});
    model.material('mat2').propertyGroup('def').set('relpermittivity', {'1'});
    model.material('mat2').propertyGroup('def').set('relpermeability', {'1'});
    model.material('mat2').selection.set([1 2 3 4 5 6 8 9 10]);
    model.physics('emw').feature.create('pml1', 'PML', 2);
    model.physics('emw').feature('pml1').selection.set([1 2 3 4 6 8 9 10]);
    model.physics('emw').feature('pml1').set('PMLfactor', 1, 'PMLorder');

    % model.physics('emw').prop('components').set('components', 1, 'inplane');
    model.physics('emw').prop('components').set('components', 1, 'outofplane');

    model.mesh('mesh1').feature.create('ftri1', 'FreeTri');
    model.mesh('mesh1').feature('ftri1').feature.create('size1', 'Size');
    model.mesh('mesh1').feature('ftri1').set('method', 'del');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('custom', 'on');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', 'on');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hminactive', 'on');
    % model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', 'wavelength/10/n');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hmax', 'wavelength/10');
    model.mesh('mesh1').feature('ftri1').feature('size1').set('hmin', '0.0001');
    model.mesh('mesh1').feature('ftri1').selection.geom('geom1');
    model.mesh('mesh1').feature('ftri1').selection.set([1 2 3 4 5 6 7 8 9]);
    model.mesh('mesh1').feature('ftri1').feature('size1').selection.set([1 2 3 4 5 6 7 8 9]);

    model.mesh('mesh1').feature.create('ftri2', 'FreeTri');
    model.mesh('mesh1').feature('ftri2').feature.create('size1', 'Size');
    model.mesh('mesh1').feature('ftri2').set('method', 'del');
    model.mesh('mesh1').feature('ftri2').feature('size1').set('custom', 'on');
    model.mesh('mesh1').feature('ftri2').feature('size1').set('hmaxactive', 'on');
    model.mesh('mesh1').feature('ftri2').feature('size1').set('hminactive', 'on');

    %model.mesh('mesh1').feature('ftri2').feature('size1').set('hmax', 'wavelength/10');

    model.mesh('mesh1').feature('ftri2').feature('size1').set('hmax', 'wavelength/10');
    model.mesh('mesh1').feature('ftri2').feature('size1').set('hmin', '0.0001');
    model.mesh('mesh1').feature('ftri2').selection.geom('geom1');
    model.mesh('mesh1').feature('ftri2').selection.set([D]);
    model.mesh('mesh1').feature('ftri2').feature('size1').selection.set([D]);



    %filename = strcat('disclinationN',num2str(N),'lx_',num2str(lx),'ly_',num2str(ly),'k_',num2str(k),'c_',num2str(c),'v_',suffix,'.mph');
    mphsave(model,['scratch_file.mph']);
    %command = ['mv scratch_file.mph ',directory,'scratch_file.mph'];
    %system(command);
    
end

