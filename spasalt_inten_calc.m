%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lambdaG_gen.m
%% Alex Cerjan
%% 5/21/14
%%
%% Calculates the generalized overlap coefficient lambdaG
%% given knowledge of the overlap elements chi_mu,nu
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = spasalt_inten_calc(datadir, N)

    %load([datadir,'overlap_matrix.mat'], 'chiMat');
    load([datadir,'overlap_matrix_wFvec_partial_v2.mat'], 'chiMat', 'fvec');
    load([datadir,'output_iterations_partial_v2.mat'], 'k_save', 'd0_save');
    %load([datadir,'output_iterations.mat'], 'k_save', 'd0_save', 'a_save');
    k_cur = k_save(:,1);
    d0_cur = d0_save(:,1);
    %a_cur = a_save{1};
    %load([datadir,'init_complete.mat'],'k_cur','d0_cur','a_cur');
    load([datadir,'cavity_info.mat'], 'gammaPerpEffective','lambda_a');
    %load([datadir,'cavity_info.mat'], 'gammaPerpEffective','nkeep', ...
    %     'lambda_a','kmin','kmax','kLength','integrationMat','whereAreAtoms','R','Nr','Ntheta','dr','dtheta');
    %kvec = linspace(kmin, kmax, kLength);

    
    [val, idx] = sort(real(d0_cur));
    d0 = val;
    k = k_cur(idx);
    ka = (2*pi)/lambda_a;
    gp = gammaPerpEffective;

    Amat = zeros(N);
    
    % chiMat was ordered in its construction.
    for nii=1:N
        Amat(:,nii) = (gp^2)/( (k(nii)-ka)^2 + gp^2 ) * chiMat(1:N,nii);
    end

    lambdaVec = zeros(N,1);
    for nii=2:N
        
        bvec = zeros(nii-1,1);
        cvec = zeros(nii-1,1);
        for mii=1:(nii-1)
            bvec(mii) = b_gen(mii, nii-1, Amat);
            cvec(mii) = c_gen(mii, nii-1, Amat, d0);
        end
        
        %size(Amat(nii,1:(nii-1)))
        %size(bvec)
        
        Ab = Amat(nii,1:(nii-1)) * bvec;
        AcD = (Amat(nii,1:(nii-1)) * cvec) * d0(nii);
        
        lambdaVec(nii) = (1/(1 - Ab))*(AcD - Ab);
        
    end     

    lambdaVec(lambdaVec > 1) = 1;
    
    %uInt = zeros(N,1);
    %for mu=1:N
    %    [et, u] = tcf_lookup(k(mu),datadir,kvec);
    %    [~,uidx] = max(abs(a_cur(:,idx(mu))));
    %    figure(mu);
    %    SALTContourPlot(u(:,uidx),R,Nr,Ntheta,dr,dtheta);
    %    uInt(mu) = sum(u(:,uidx)' * integrationMat * u(:,uidx));
    %end
    uInt = ones(N,1);
    
    int_vec = zeros(N);
    d0_inter = zeros(N,1);
    d0_inter(1) = d0(1);
    for ii=2:N
        if (lambdaVec(ii) < 1)
            d0_inter(ii) = d0(ii)*(1/(1-lambdaVec(ii)));
        else
            d0_inter(ii) = 999;
        end
        d0_inter(3) = 0.0216; %0.0216; %0.017474;

        if (d0_inter(ii) > 0.03)
            d0_inter(ii) = 0.03;
        end
        
        for jj=1:(ii-1)
            b = b_gen(jj,(ii-1),Amat);
            c = c_gen(jj,(ii-1),Amat,d0);
            
            int_vec(ii,jj) = c*d0_inter(ii) - b;
        end
    end
        
    figure(N+1);
    cc = distinguishable_colors(N);
    for ii=1:N
        plot(d0_inter,int_vec(:,ii),'Color',cc(ii,:));
        hold on;
    end


end

%%%%%%%

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