%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lambdaG_gen.m
%% Alex Cerjan
%% 5/21/14
%%
%% Calculates the generalized overlap coefficient lambdaG
%% given knowledge of the overlap elements chi_mu,nu
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lambdaVec] = lambdaG_gen(datadir, N)

    load([datadir,'overlap_matrix_wFvec.mat'], 'chiMat');
    load([datadir,'output_iterations.mat'], 'k_save', 'd0_save');
    load([datadir,'cavity_info.mat'], 'gammaPerpEffective','lambda_a');
    [val, idx] = sort(d0_save(:,1));
    d0vec = val;
    k = k_save(idx,1);
    ka = (2*pi)/lambda_a;
    gp = gammaPerpEffective;

    Amat = zeros(N);
    for nii=1:N
        Amat(:,nii) = (gp^2)/( (k(nii)-ka)^2 + gp^2 ) * chiMat(1:N,nii);
    end

    lambdaVec = zeros(N,1);
    for nii=2:N
        
        bvec = zeros(nii-1,1);
        cvec = zeros(nii-1,1);
        for mii=1:(nii-1)
            bvec(mii) = b_gen(mii, nii-1, Amat);
            cvec(mii) = c_gen(mii, nii-1, Amat, d0vec);
        end
        
        %size(Amat(nii,1:(nii-1)))
        %size(bvec)
        
        Ab = Amat(nii,1:(nii-1)) * bvec;
        AcD = (Amat(nii,1:(nii-1)) * cvec) * d0vec(nii);
        
        lambdaVec(nii) = (1/(1 - Ab))*(AcD - Ab);
        
    end     

    lambdaVec(lambdaVec > 1) = 1;
    
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