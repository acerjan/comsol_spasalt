%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lambdaG_gen.m
%% Alex Cerjan
%% 5/21/14
%%
%% Calculates the generalized overlap coefficient lambdaG
%% given knowledge of the overlap elements chi_mu,nu
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lambdaVec] = spasalt_calc(datadir, R, lambda_a, Q_thresh,...
                                    gamma_perp_length, geom_switch, geom_element)

    %% INIT:
    %% load comsol data:
    
    tmp = dlmread([datadir, 'lambda_Q']);
    lambda = tmp(:,1);
    Q = tmp(:,2); 
                    
    aboveZeroIdx = find(Q>Q_thresh);
    clear tmp Q;
    
    %% load spasalt_setup:
    load([datadir,'spasalt_setup.mat'],'chiMat','epsVec','Q','lambda','cavityLocs');
    
    k_a = 2*pi/lambda_a;
    gamma_perp = (gamma_perp_length/lambda_a) * (2*pi/lambda_a);

    %% calculate pump profile and fvec:
    tmp = dlmread([datadir, 'grid_xy']);
    xpts = tmp(1,:);
    ypts = tmp(2,:);
    dx = abs(xpts(2) - xpts(1));
    dy = abs(ypts(2) - ypts(1));

    pumpProfile = zeros(length(xpts),length(ypts));
    switch geom_switch
      case 'D'
        r0 = geom_element * R;

        for xii=1:length(xpts)
            x = xpts(xii);
            for yii=1:length(ypts)
                y = ypts(yii);

                theta = atan(y/x)+pi/2;
                r = sqrt(x^2 + y^2);
                
                if ( (r^2 <= R^2) && (y <= r0) )
                    pumpProfile(xii,yii) = 1;
                end
                
                if (r <= 1.625)
                    pumpProfile(xii,yii) = 1;
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
                    pumpProfile(xii,yii) = 1;
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
                    pumpProfile(xii,yii) = 1;
                end
            end
        end

      otherwise
        error('I do not recognize your choice of geometry.');
    end
        
    %sum(reshape(cavityLocs,[],1)*dx*dy)
    %sum(reshape(pumpProfile,[],1)*dx*dy)    
    pumpProfile = pumpProfile*(sum(reshape(cavityLocs,[],1)*dx*dy)/sum(reshape(pumpProfile,[],1)*dx*dy));

    fVec = zeros(length(lambda),1);
    parfor ii=1:length(lambda)
        idxI = aboveZeroIdx(ii);
        %EzI = dlmread([datadir,'Ez_sol',num2str(idxI)]);
        EzI = parload([datadir,'Ez_sol',num2str(idxI),'.mat']);
        EzI = reshape(EzI, [], 1);
        
        fVec(ii) = real(sum(EzI .* EzI .* reshape(pumpProfile,[],1))*dx*dy);
    end
        
    %% calculate D from Q and lambda:

    k = 2*pi./lambda;
    k_n = k - 1i*k./(2*Q);
    et = ( (k_n.^2)./(k.^2) - 1 );
    D = et .* ((k - k_a + 1i*gamma_perp) / gamma_perp) .* (epsVec./fVec);
    D = real(D);
    [~,idx] = sort(D,'ascend');
    
    %% reorder everything:
    D = D(idx);
    k = k(idx);
    %EzMat = EzMat(:,idx);
    chiMat = chiMat(idx,idx);

    %% calculate generalized mode competition parameters:
    N = length(D);    
    Amat = zeros(N);
    for nii=1:N
        Amat(:,nii) = (gamma_perp^2)/( (k(nii)-k_a)^2 + gamma_perp^2 ) * chiMat(1:N,nii);
    end

    lambdaVec = zeros(N,1);
    for nii=2:N
        
        bvec = zeros(nii-1,1);
        cvec = zeros(nii-1,1);
        for mii=1:(nii-1)
            bvec(mii) = b_gen(mii, nii-1, Amat);
            cvec(mii) = c_gen(mii, nii-1, Amat, D);
        end
        
        %size(Amat(nii,1:(nii-1)))
        %size(bvec)
        
        Ab = Amat(nii,1:(nii-1)) * bvec;
        AcD = (Amat(nii,1:(nii-1)) * cvec) * D(nii);
        
        lambdaVec(nii) = (1/(1 - Ab))*(AcD - Ab);
        
    end     

    lambdaVec(lambdaVec >= .9999) = .9999;
    lambdaVec(lambdaVec < 0) = .9999;
    idx = find(lambdaVec<1);
    DinterUni = D(idx)./(1-lambdaVec(idx));
    [~,idx] = sort(DinterUni,'ascend');
    DinterUni = DinterUni(idx);

    save([datadir,'spasalt_results.mat'],'lambdaVec','DinterUni');
    
    %% calculate the intensities:
    %int_vec = zeros(N);
    %d0_inter = zeros(N,1);
    %d0_inter(1) = d0(1);
    %for ii=2:N
    %    if (lambdaVec(ii) < 1)
    %        d0_inter(ii) = d0(ii)*(1/(1-lambdaVec(ii)));
    %    else
    %        d0_inter(ii) = 999;
    %    end

    %    if (d0_inter(ii) > 0.03)
    %        d0_inter(ii) = 0.03;
    %    end
        
    %    for jj=1:(ii-1)
    %        b = b_gen(jj,(ii-1),Amat);
    %        c = c_gen(jj,(ii-1),Amat,d0);
    %        
    %        int_vec(ii,jj) = c*d0_inter(ii) - b;
    %    end
    %end
    
end



%%%%%%% helper functions:

function [c_mu] = c_gen(mu, n, Amat, d0vec)
    
    Amat = Amat(1:n,1:n);
    aInvMat = inv(Amat);
    
    c_mu = sum(aInvMat(mu,:).' ./ d0vec(1:n));
    
end

function [b_mu] = b_gen(mu, n, Amat)
    
    Amat = Amat(1:n,1:n);
    aInvMat = inv(Amat);
    
    b_mu = sum(aInvMat(mu,:));
    
end

function [Ez] = parload(fname)
    load(fname,'Ez');
end