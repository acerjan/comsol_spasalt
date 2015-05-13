function [Q lambda] = extract_field( filename,threshold,num_modes,field_fname,rangeString)
%FCTNEXTRACTFIELD001 Summary of this function goes here
%   Detailed explanation goes here

import com.comsol.model.*
import com.comsol.model.util.*


model = mphload(strcat(filename,'.mph'));


model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical('gev1').set('descr', ' Frequency Frequency');
model.result.numerical('gev1').set('expr', 'freq');
model.result.table.create('tbl1', 'Table');
model.result.table('tbl1').comments('Global Evaluation 1 (freq)');
model.result.numerical('gev1').set('table', 'tbl1');
model.result.numerical('gev1').setResult;
model.result.table('tbl1').save(['scratcheigenwavelength','.txt']);
model.result.table.remove('tbl1');
model.result.numerical.remove('gev1');

filename = strcat('scratcheigenwavelength','.txt');
A = dlmread(filename,'%',5,0);
C = 3e8;
Q = real(A(:,1))./imag(A(:,1))/2;
lambda = C./real(A(:,1));

figure(1);clf;
stem(lambda,Q);
dlmwrite([field_fname 'lambda_Q'],[lambda Q]);


for ii=1:num_modes
 
    timestart = clock
    

 if(Q(ii)>threshold)
    model.result.export.create('data1', 'Data');
    model.result.export('data1').set('solnum', {num2str(ii)});
    model.result.export('data1').set('expr', {'Ez'});
    model.result.export('data1').set('descr', {'Electric field, z component'});
    model.result.export('data1').set('filename','scratch_Ez');
    model.result.export('data1').set('location', 'grid');
    model.result.export('data1').set('gridstruct', 'grid');
    model.result.export('data1').set('gridx2', ['range(',rangeString,')']);
    model.result.export('data1').set('gridy2', ['range(',rangeString,')']);
    model.result.export('data1').run;
    model.result.export.remove('data1');
    
    
    
%     model.result.export.create('data1', 'Data');
%     model.result.export('data1').set('solnum', {num2str(ii)});
%     model.result.export('data1').set('expr', {'Ex'});
%     model.result.export('data1').set('descr',  {'Electric field, x component'});
%     model.result.export('data1').set('filename', 'scratch_Ex');
%     model.result.export('data1').set('location', 'grid');
%     model.result.export('data1').set('gridstruct', 'grid');
%     model.result.export('data1').set('gridx2', ['range(',rangeString,')']);
%     model.result.export('data1').set('gridy2', ['range(',rangeString,')']);
%     model.result.export('data1').run;
%     model.result.export.remove('data1');
%     
% 
%     model.result.export.create('data1', 'Data');
%     model.result.export('data1').set('solnum', {num2str(ii)});
%     model.result.export('data1').set('expr', {'Ey'});
%     model.result.export('data1').set('descr',  {'Electric field, y component'});
%     model.result.export('data1').set('filename', 'scratch_Ey');
%     model.result.export('data1').set('location', 'grid');
%     model.result.export('data1').set('gridstruct', 'grid');
%     model.result.export('data1').set('gridx2', ['range(',rangeString,')']);
%     model.result.export('data1').set('gridy2', ['range(',rangeString,')']);
%     model.result.export('data1').run;
%     model.result.export.remove('data1');
    
    
    Ez = dlmread('scratch_Ez','%',12,0);
%     Ex = dlmread('scratch_Ex','%',12,0);
%     Ey = dlmread('scratch_Ey','%',12,0);

    figure(1);clf;
    imagesc(abs(Ez));
    
    dlmwrite([field_fname 'Ez_sol' num2str(ii)],Ez);    
    
 end
    timeend = clock
end

    model.result.export.create('data1', 'Data');
    model.result.export('data1').set('solnum', '1');
    model.result.export('data1').set('expr', {'emw.epsrAv'});
    model.result.export('data1').set('descr', {'Relative permittivity, average'});
    model.result.export('data1').set('filename','scratch_structure');  
    model.result.export('data1').set('location', 'grid');
    model.result.export('data1').set('gridstruct', 'grid');
    model.result.export('data1').set('gridx2', ['range(',rangeString,')']);
    model.result.export('data1').set('gridy2', ['range(',rangeString,')']);
    model.result.export('data1').run;
    model.result.export.remove('data1');
    
    
    structure = dlmread('scratch_structure','%',12,0);
    dlmwrite([field_fname 'structure'],structure);
 
end

