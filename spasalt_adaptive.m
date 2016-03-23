%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lambdaG_gen.m
%% Alex Cerjan
%% 5/21/14
%%
%% Calculates the generalized overlap coefficient lambdaG
%% given knowledge of the overlap elements chi_mu,nu
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [minimizeThis] = spasalt_adaptive(pumpRatios, datadir, R, lambda_a, ...
                                           Q_thresh, gamma_perp_length, ...
                                           numModes, numR, numTH, ...
                                           FINAL_BOOL)

    pumpVec = [1, pumpRatios];
    
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

    [X,Y] = meshgrid(xpts,ypts);
    
    RR = sqrt(X.^2 + Y.^2);
    TH = atan2(Y,X);
    TH(isnan(TH)) = 0;
    
    rBreaks = linspace(0,R,numR+2);
    thBreaks = linspace(0,pi,numTH+2);
    
    for rii=2:length(rBreaks)
        for tii=2:length(thBreaks)
            
            pIdx = (length(thBreaks)-1)*(rii-2)+(tii-1);
            
            rPrev = rBreaks(rii-1);
            rNext = rBreaks(rii);
            
            thPrev = thBreaks(tii-1);
            thNext = thBreaks(tii);
            
            if (tii==length(thBreaks))
                pumpProfile( (rPrev<=RR) & (RR<rNext) & (thPrev<=TH) & ...
                             (TH<=thNext) ) = pumpVec(pIdx);
                pumpProfile( (rPrev<=RR) & (RR<rNext) & (-thNext<=TH) & ...
                             (TH<=-thPrev) ) = pumpVec(pIdx);
            else
                pumpProfile( (rPrev<=RR) & (RR<rNext) & (thPrev<=TH) & ...
                             (TH<thNext) ) = pumpVec(pIdx);
                pumpProfile( (rPrev<=RR) & (RR<rNext) & (-thNext<TH) & ...
                             (TH<=-thPrev) ) = pumpVec(pIdx);
            end
        end
    end
    pumpProfile = pumpProfile .* cavityLocs;
    pumpProfile = pumpProfile*(sum(reshape(cavityLocs,[],1)*dx*dy)/sum(reshape(pumpProfile,[],1)*dx*dy));
    
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
    
    Duni = et .* ((k - k_a + 1i*gamma_perp) / gamma_perp) .* (epsVec);
    Duni = real(Duni);
    
    %% reorder everything:
    D = D(idx);
    k = k(idx);
    chiMat = chiMat(idx,idx);
    aboveZeroIdx = aboveZeroIdx(idx);
    
    Duni = Duni(idx);
    DuniThr = min(Duni);

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
        
        Ab = Amat(nii,1:(nii-1)) * bvec;
        AcD = (Amat(nii,1:(nii-1)) * cvec) * D(nii);
        
        lambdaVec(nii) = (1/(1 - Ab))*(AcD - Ab);
        
    end     
    
    %% optimization location:

    lambdaVec(lambdaVec >= .9999) = .9999;
    lambdaVec(lambdaVec < 0) = .9999;
    D_interacting = D./(1-lambdaVec);

    [~,idx] = sort(D_interacting,'ascend');
    D_interacting = D_interacting(idx);

    if (numModes == 1)
        minimizeThis = (D_interacting(1)/D_interacting(2))+abs(1-(D_interacting(1)/DuniThr));
    else    
        D_ratio = D_interacting/DuniThr;
        minimizeThis = prod(Drat(1:numModes));    
    end
    disp(minimizeThis);
    
    if (FINAL_BOOL == 1)
        imagesc(pumpProfile);
        D_adaptive_nonInt = D(idx);
        D_adaptive_interacting = D_interacting;
        D_uniform_nonInt_sorted = Duni(idx);
        optPumpVec = pumpVec(2:end);
        pumpAdap = pumpProfile;
        aboveZeroIdxAdap = aboveZeroIdx(idx);
        AmatAdap = Amat(idx,idx);
        chiMatAdap = chiMat(idx,idx);
        
        save([datadir,'spasalt_adaptive.mat'],'optPumpVec', ...
             'D_adaptive_nonInt', 'D_adaptive_interacting', ...
             'D_uniform_nonInt_sorted', 'pumpAdap', 'aboveZeroIdxAdap', ...
             'AmatAdap','chiMatAdap');
        
        minimizeThis = [];
    end   
    
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
