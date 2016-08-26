%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lambdaG_gen.m
%% Alex Cerjan
%% 5/21/14
%%
%% Calculates the generalized overlap coefficient lambdaG
%% given knowledge of the overlap elements chi_mu,nu
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lambdaVec] = spasalt_condensed(datadir, R, lambda_a, Q_thresh, ...
                                         gamma_perp_length, geom_switch, ...
                                         geom_element, nModes)
    
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

      case 'Stadium'
        L = geom_element(1);
        r0 = geom_element(2);
        for xii=1:length(xpts)
            x = xpts(xii);
            for yii=1:length(ypts)
                y = ypts(yii);

                if ( (x<=r0)&&(x>=-r0)&&(y>=-L/2)&&(y<=L/2) )
                    pumpProfile(xii,yii) = 1;
                elseif ( (x<=r0)&&(x>=-r0)&&(y>L/2) )
                    r = sqrt((y-L/2)^2 + x^2);
                    if (r<=r0)                        
                        pumpProfile(xii,yii) = 1;
                    end
                elseif ( (x<=r0)&&(x>=-r0)&&(y<-L/2) )
                    r = sqrt((y+L/2)^2 + x^2);
                    if (r<=r0)                        
                        pumpProfile(xii,yii) = 1;
                    end                    
                end                    
            end
        end
      
      otherwise
        error('I do not recognize your choice of geometry.');
    end
        
    pumpProfile = pumpProfile*(sum(reshape(cavityLocs,[],1)*dx*dy)/sum(reshape(pumpProfile,[],1)*dx*dy));
    pumpUni = pumpProfile;
    
    fVec = zeros(length(lambda),1);
    parfor ii=1:length(lambda)
        idxI = aboveZeroIdx(ii);
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
    chiMat = chiMat(idx,idx);
    aboveZeroIdx = aboveZeroIdx(idx);

    %% calculate generalized mode competition parameters:
    N = length(D);    
    Amat = zeros(N);
    for nii=1:N
        Amat(:,nii) = (gamma_perp^2)/( (k(nii)-k_a)^2 + gamma_perp^2 ) * chiMat(1:N,nii);
    end

    lambdaVec = zeros(N,1);
    DinterUni = zeros(N,1);
    DinterUni(1) = D(1);
    
    trig=0;
    for nii=2:N
        if (trig==0)
            [D_uni_interacting, lambdaVec, trig] = gen_next_thresh(D,Amat,D_uni_interacting,lambdaVec);
        end
    end

    %% calculate the intensities when nMode turns on:

    [~,idx] = sort(DinterUni,'ascend');
    DinterUni = DinterUni(idx);
    D = D(idx);
    Amat = Amat(idx,idx);
    aboveZeroIdx = aboveZeroIdx(idx);
    Duni = D;
    
    int_vec = zeros(nModes,1);
    d0=DinterUni(nModes);
    for ii=1:nModes-1
        b = b_gen(ii,(nModes-1),Amat);
        c = c_gen(ii,(nModes-1),Amat,D);
        int_vec(ii) = c*d0-b;
    end
        
    %% calculate the holeburning term:
    
    holeTerm = zeros(length(xpts)*length(ypts),1);
    for ii=1:nModes-1
        idxI = aboveZeroIdx(ii);
        EzI = parload([datadir,'Ez_sol',num2str(idxI),'.mat']);
        EzI = reshape(EzI, [], 1);
        Gamma = (gamma_perp^2)/( (k(ii)-k_a)^2 + gamma_perp^2 );
        holeTerm = holeTerm + Gamma*int_vec(ii)*abs(EzI).^2;
    end
    holeTerm = reshape(cavityLocs,[],1)./(1+holeTerm);
    pumpHole = holeTerm*(sum(reshape(cavityLocs,[],1)*dx*dy)/sum(reshape(holeTerm,[],1)*dx*dy));
        
    % ok, at this point we have recovered the holeburning term when
    % the nMode has reached threshold, and defined a new pump
    % profile based upon this. We now need to redo the initial
    % calculation to recover the new set of thresholds.
    
    %% reload comsol data: (need to refresh aboveZeroIdx)
    
    tmp = dlmread([datadir, 'lambda_Q']);
    lambda = tmp(:,1);
    Q = tmp(:,2); 
                    
    aboveZeroIdx = find(Q>Q_thresh);
    clear tmp Q;

    %% calculate the thresholds with the new pump profile:
    load([datadir,'spasalt_setup.mat'],'chiMat','epsVec','Q','lambda','cavityLocs');

    fVec = zeros(length(lambda),1);
    parfor ii=1:length(lambda)
        idxI = aboveZeroIdx(ii);
        EzI = parload([datadir,'Ez_sol',num2str(idxI),'.mat']);
        EzI = reshape(EzI, [], 1);
        
        fVec(ii) = real(sum(EzI .* EzI .* reshape(pumpHole,[],1))*dx*dy);
    end
            
    D = et .* ((k - k_a + 1i*gamma_perp) / gamma_perp) .* (epsVec./fVec);
    D = real(D);
    [~,idx] = sort(D,'ascend');
    
    D = D(idx);
    Dhole = D;
    k = k(idx);
    chiMat = chiMat(idx,idx);
    aboveZeroIdx = aboveZeroIdx(idx);

    %% calculate generalized mode competition parameters:
    N = length(D);    
    Amat = zeros(N);
    for nii=1:N
        Amat(:,nii) = (gamma_perp^2)/( (k(nii)-k_a)^2 + gamma_perp^2 ) * chiMat(1:N,nii);
    end

    lambdaVec = zeros(N,1);
    DinterHole = zeros(N,1);
    DinterHole(1) = D(1);
    
    trig=0;
    for nii=2:N
        if (trig==0)
            [D_uni_interacting, lambdaVec, trig] = gen_next_thresh(D,Amat,D_uni_interacting,lambdaVec);
        end
    end

    [~,idx] = sort(DinterHole,'ascend');
    DinterHole = DinterHole(idx);
    aboveZeroIdxCond = aboveZeroIdx(idx);
    AmatCond = Amat(idx,idx);
    chiMatCond = chiMat(idx,idx);
    
    pumpConden = reshape(pumpHole,length(xpts),length(ypts));

    D_uni = Duni;
    D_uni_interacting = DinterUni;
    D_conden = Dhole(idx);
    D_conden_interacting = DinterHole;
    
    save([datadir,'spasalt_condensed.mat'],'D_uni','D_conden', ...
         'D_uni_interacting','D_conden_interacting','pumpUni', ...
         'pumpConden', 'aboveZeroIdxCond', 'AmatCond','chiMatCond');   

end



%%%%%%% helper functions:

function [Dint, lambdaSAVE, trig] = gen_next_thresh(DnonInt, Amat, Dint,lambdaSAVE)
    
    idxOn = find(Dint);
    idxOff = find(Dint==0);
    
    nCur = length(idxOn);
    Don = DnonInt(idxOn);
    AmatOn = Amat(idxOn,idxOn);
    bvec = zeros(nCur,1);
    cvec = zeros(nCur,1);
    for ii=1:nCur
        bvec(ii) = b_gen(ii,nCur,AmatOn);
        cvec(ii) = c_gen(ii,nCur,AmatOn,Don);
    end
    
    nTest = length(idxOff);    
    Dcheck = zeros(nTest,1);
    lambdaVec = zeros(nTest,1);
    DintCheck = zeros(nTest,1);
    for ii=1:nTest
        Ab = Amat(idxOff(ii),idxOn)*bvec;
        AcD = (Amat(idxOff(ii),idxOn)*cvec)*DnonInt(idxOff(ii));
        
        lambdaVec(ii) = (1/(1-Ab))*(AcD-Ab);
        if ((lambdaVec(ii) > 0.9999) || (lambdaVec(ii) < 0))
            lambdaVec(ii) = 0.9999;
        end
        DintCheck(ii) = DnonInt(idxOff(ii))/(1-lambdaVec(ii));
    end
        
    [~,idx] = sort(abs(DintCheck),'ascend');
    Dint(idxOff(idx(1))) = DintCheck(idx(1));
    lambdaSAVE(idxOff(idx(1))) = lambdaVec(idx(1));
    
    trig=0;
    if (lambdaVec(idx(1)) == 0.9999)
        trig=1;
        lambdaSAVE(idxOff) = 0.9999;
        Dint(idxOff) = DnonInt(idxOff)./(1-lambdaSAVE(idxOff));
    end
        
end

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
