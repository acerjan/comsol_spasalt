%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% comsol_spasalt.m
%% Alex Cerjan with some original code from Brandon Redding
%% 7/3/14
%%
%% Code for analyzing and performing the SPA-SALT calculation
%% on data acquired from COMSOL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = spasalt_setup(datadir, R, n_inside, Q_thresh, geom_switch, geom_element)

    
    %% init:
    n_eff = n_inside; 

    %% load comsol data:
    
    tmp = dlmread([datadir, 'lambda_Q']);
    lambda = tmp(:,1);
    Q = tmp(:,2); 
                    
    aboveZeroIdx = find(Q>Q_thresh);
    N = length(aboveZeroIdx);
    Q = Q(aboveZeroIdx);
    lambda = lambda(aboveZeroIdx);    
    
    %% construct the cavity edge and dielectric:
    
    tmp = dlmread([datadir, 'grid_xy']);
    xpts = tmp(1,:);
    ypts = tmp(2,:);
    dx = abs(xpts(2) - xpts(1));
    dy = abs(ypts(2) - ypts(1));

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
    epsCav = n_eff^2 * cavityLocs;
    
    %% calculate overlap integrals and epsvec:
    
    epsVec = zeros(length(aboveZeroIdx),1);
    chiMat = zeros(N);

    parfor ii=1:N
        idxI = aboveZeroIdx(ii);
        %EzI = dlmread([datadir,'Ez_sol',num2str(idxI)]);
        EzI = parload([datadir,'Ez_sol',num2str(idxI),'.mat']);
        EzI = reshape(EzI, [], 1);
        epsVec(ii) = real(sum(EzI .* EzI .* reshape(epsCav,[],1))*dx*dy);

        for jj=1:N
            idxJ = aboveZeroIdx(jj);
            %EzJ = dlmread([datadir,'Ez_sol',num2str(idxJ)]);
            EzJ = parload([datadir,'Ez_sol',num2str(idxJ),'.mat']);
            EzJ = reshape(EzJ, [], 1);
            
            chi = sum(EzI .* EzI .* EzJ .* conj(EzJ) .* reshape(epsCav,[],1))*dx*dy/epsVec(ii);            
            chiMat(ii,jj) = real(chi);
        end
    end
    
    save([datadir,'spasalt_setup.mat'],'chiMat','epsVec','Q','lambda','cavityLocs','aboveZeroIdx');    
    
end

function [Ez] = parload(fname)
    load(fname,'Ez');
end