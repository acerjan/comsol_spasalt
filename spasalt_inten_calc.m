%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lambdaG_gen.m
%% Alex Cerjan
%% 5/21/14
%%
%% Calculates the generalized overlap coefficient lambdaG
%% given knowledge of the overlap elements chi_mu,nu
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = spasalt_inten_calc(datadir, nModes, which_results)

    switch which_results
      case 1
        load([datadir,'spasalt_uniform_results.mat'],'AmatUni','D_uni', ...
             'D_uni_interacting');
        D = D_uni;
        Dint = D_uni_interacting;
        Amat = AmatUni;
        
      case 2
        load([datadir,'spasalt_adaptive.mat'],'AmatAdap','D_adaptive_nonInt', ...
             'D_adaptive_interacting');
        D = D_adaptive_nonInt;
        Dint = D_adaptive_interacting;
        Amat = AmatAdap;
                
      case 3
        load([datadir,'spasalt_condensed.mat'],'AmatCond','D_conden', ...
             'D_conden_interacting');
        D = D_conden;
        Dint = D_conden_interacting;
        Amat = AmatCond;
        
      otherwise
        error('I do not know which set of results to calculate the intensity for.');
    end
    
    int_vec = zeros(nModes);
    for ii=2:nModes
        %Dint(3) = 0.0216; %0.0216; %0.017474;

        for jj=1:(ii-1)
            b = b_gen(jj,(ii-1),Amat);
            c = c_gen(jj,(ii-1),Amat,D);
            
            int_vec(ii,jj) = c*Dint(ii) - b;
        end
    end
        
    figure(1);
    cc = distinguishable_colors(nModes);
    for ii=1:nModes
        plot(Dint(1:nModes),int_vec(:,ii),'Color',cc(ii,:));
        hold on;
    end
    hold off;

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