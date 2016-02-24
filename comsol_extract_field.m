function [Q lambda] = comsol_extract_field( filename,threshold, ...
                                            num_modes,directory, ...
                                            rangeString, R, geom_switch, ...
                                            geom_element)

    %% SOME PARTS OF SPA-SALT MOVED HERE FOR SPEED:
    
    %% construct the cavity edge:
    
    tmp = dlmread([directory, 'grid_xy']);
    xpts = tmp(1,:);
    ypts = tmp(2,:);
    dx = abs(xpts(2) - xpts(1));
    dy = abs(ypts(2) - ypts(1));
    clear tmp;

    cavityLocs = zeros(length(xpts),length(ypts));
    switch geom_switch
      case 'D'
        r0 = geom_element * R;

        for xii=1:length(xpts)
            x = xpts(xii);
            for yii=1:length(ypts)
                y = ypts(yii);
                
                if ( (x^2 + y^2 <= R^2) && (y <= r0) )
                    cavityLocs(xii,yii) = 1;
                end
            end
        end
        
      case 'Quad'
        r0 = R/(1+geom_element);
        for xii=1:length(xpts)
            x = xpts(xii);
            for yii=1:length(ypts)
                y = ypts(yii);
                
                theta = atan(y/x)+pi/2;
                r = sqrt(x^2 + y^2);
                if (r <= r0*(1+geom_element*cos(2*theta)))
                    cavityLocs(xii,yii) = 1;
                end

            end
        end
       
      case 'Ellipse'
        aa = geom_element(1);
        bb = geom_element(2);
        for xii=1:length(xpts)
            x = xpts(xii);
            for yii=1:length(ypts)
                y = ypts(yii);
                
                if ( (x/bb)^2 + (y/aa)^2 <= 1)
                    cavityLocs(xii,yii) = 1;
                end
            end
        end
        
      otherwise
        error('I do not recognize your choice of geometry.');
    end    

    %% BEGIN COMSOL:
    %FCTNEXTRACTFIELD001 Summary of this function goes here
    %   Detailed explanation goes here

    import com.comsol.model.*
    import com.comsol.model.util.*


    model = mphload([filename,'.mph']);


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

    system(['mv scratcheigenwavelength.txt ',directory]);
    filename = [directory,'scratcheigenwavelength','.txt'];
    A = dlmread(filename,'%',5,0);
    C = 3e8;
    Q = real(A(:,1))./imag(A(:,1))/2;
    lambda = C./real(A(:,1));

    figure(1);clf;
    stem(lambda,Q);
    dlmwrite([directory,'lambda_Q'],[lambda Q]);


    parfor ii=1:num_modes
        
        if(Q(ii)>threshold)
            disp(Q(ii));
            
            model.result.export.create('data1', 'Data');
            model.result.export('data1').set('solnum', {num2str(ii)});
            model.result.export('data1').set('expr', {'Ez'});
            model.result.export('data1').set('descr', {'Electric field, z component'});
            model.result.export('data1').set('filename',['Ez_sol',num2str(ii),'.mat']);
            model.result.export('data1').set('location', 'grid');
            model.result.export('data1').set('gridstruct', 'grid');
            model.result.export('data1').set('gridx2', ['range(',rangeString,')']);
            model.result.export('data1').set('gridy2', ['range(',rangeString,')']);
            model.result.export('data1').run;
            model.result.export.remove('data1');
            
            system(['mv Ez_sol',num2str(ii),'.mat ',directory]);
            
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
            
            
            Ez = dlmread([directory,'Ez_sol',num2str(ii),'.mat'],'%',12,0);
            %     Ex = dlmread('scratch_Ex','%',12,0);
            %     Ey = dlmread('scratch_Ey','%',12,0);

            %figure(1);clf;
            %imagesc(abs(Ez));

            Ez = Ez .* cavityLocs;
            normFac = sum(reshape(Ez,[],1) .* reshape(Ez,[],1))*dx*dy;
            Ez = Ez/sqrt(normFac);

            parsave([directory,'Ez_sol',num2str(ii),'.mat'],Ez);
            %dlmwrite([directory,'Ez_sol',num2str(ii),'.mat'],Ez);    
            
        end
    end

    model.result.export.create('data1', 'Data');
    model.result.export('data1').set('solnum', '1');
    model.result.export('data1').set('expr', {'emw.epsrAv'});
    model.result.export('data1').set('descr', {'Relative permittivity, average'});
    model.result.export('data1').set('filename',['scratch_structure']);  
    model.result.export('data1').set('location', 'grid');
    model.result.export('data1').set('gridstruct', 'grid');
    model.result.export('data1').set('gridx2', ['range(',rangeString,')']);
    model.result.export('data1').set('gridy2', ['range(',rangeString,')']);
    model.result.export('data1').run;
    model.result.export.remove('data1');
    
    
    structure = dlmread(['scratch_structure'],'%',12,0);
    dlmwrite([directory,'structure'],structure);
    
end

function [] = parsave(fname,Ez)
    save(fname,'Ez');
end