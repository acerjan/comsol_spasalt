%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% comsol_spasalt.m
%% Alex Cerjan with some original code from Brandon Redding
%% 7/3/14
%%
%% Code for analyzing and performing the SPA-SALT calculation
%% on data acquired from COMSOL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = comsol_spasalt(datadir, N)

    %% init:
    folder = datadir;

    lambda_a = 1; % um
    n_eff = 3.5;
    %n_eff = 3 + .002*1i;
    k_a = 2*pi/lambda_a;
    %gammaPerpLambda = .03;
    gammaPerpLambda = .01; % um
    gamma_perp = (gammaPerpLambda/lambda_a) * (2*pi/lambda_a);
    gammaPerpEffective = gamma_perp;
    
    R = 4;
    r0 = .5 * R;
    %R = 4;
    %epsilon = .16;
    %r0 = R/(1+epsilon);
    
    Q_thresh = 200;
    
    save([datadir,'cavity_info.mat'], 'gammaPerpEffective','lambda_a');

    %% load comsol data:
    
    tmp = dlmread([folder 'lambda_Q']);
    lambda = tmp(:,1); %conversion to um.
    Q = tmp(:,2); % need the 2 to convert to photon decay rate
                    % from field decay rate.
    aboveZeroIdx = find(Q>Q_thresh);
    Q = Q(aboveZeroIdx);
    lambda = lambda(aboveZeroIdx);
    
    %% construct the cavity edge and pump profile:
    
    tmp = dlmread([datadir, 'grid_xy']);
    xpts = tmp(1,:);
    ypts = tmp(2,:);
    dx = abs(xpts(2) - xpts(1));
    dy = abs(ypts(2) - ypts(1));

    cavityLocs = zeros(length(xpts),length(ypts));
    pumpProfile = zeros(length(xpts),length(ypts));
    for xii=1:length(xpts)
        x = xpts(xii);
        for yii=1:length(ypts)
            y = ypts(yii);
            
            if ( (x^2 + y^2 <= R^2) && (y <= r0) )
                cavityLocs(xii,yii) = 1;
                %pumpProfile(xii,yii) = 1;
            end
            
            theta = atan(y/x)+pi/2;
            r = sqrt(x^2 + y^2);
            %if (r <= r0*(1+epsilon*cos(2*theta)))
            %    cavityLocs(xii,yii) = 1;
            %    pumpProfile(xii,yii) = 1;
            %end

            if (r <= 1.625)
                pumpProfile(xii,yii) = 1;
            end
            
        end
    end
    
    %sum(reshape(cavityLocs,[],1)*dx*dy)
    %sum(reshape(pumpProfile,[],1)*dx*dy)
    
    pumpProfile = pumpProfile*(sum(reshape(cavityLocs,[],1)*dx*dy)/sum(reshape(pumpProfile,[],1)*dx*dy));
    
    %% normalize Ez fields and construct fvec:

    EzMat = zeros(length(reshape(cavityLocs,[],1)),length(aboveZeroIdx));
    fvec = zeros(length(aboveZeroIdx),1);
    
    for ii=1:length(aboveZeroIdx)
        idxI = aboveZeroIdx(ii);
       
        EzI = dlmread([datadir,'Ez_sol',num2str(idxI)]);
        EzI = EzI .* cavityLocs;
        EzI = reshape(EzI, [], 1);
        
        normFac = sum(EzI .* EzI)*dx*dy;
        EzMat(:,ii) = EzI/sqrt(normFac);
        
        fvec(ii) = real(sum(EzMat(:,ii) .* EzMat(:,ii) .* reshape(pumpProfile,[],1) *dx*dy));
    end
       
    save([datadir,'comsol_fvec_debug_partial_v2.mat'],'EzMat','fvec');
    %load([datadir,'comsol_fvec_debug_partial_v2.mat'],'EzMat','fvec');
    
    
    %% calculate D from Q and lambda:

    k = 2*pi./lambda;
    k_n = k - 1i*k./(2*Q);
    et = ( (k_n.^2)./(k.^2) - 1 )*n_eff*n_eff;
    D = et .* (k - k_a + 1i*gamma_perp) / gamma_perp ./ fvec;
    D = real(D);
    [~,idx] = sort(D,'ascend');
    [idx,D(idx)]
    EzMat = EzMat(:,idx);
    %idx = aboveZeroIdx(idx);
    
    d0_save = D;
    k_save = k;
    
    save([datadir,'output_iterations_partial_v2.mat'], 'k_save', 'd0_save');
    
    %% can plot here to see if data makes sense:
    %plot(lambda, Q, 'bo');
    %hold on;
    %plot(lambda(idx(1:10)), Q(idx(1:10)), 'rx');
    %hold off;
    
    %% generate overlap integral values:

    %aboveZeroIdx(idx(1))
    Ez = dlmread([folder 'Ez_sol' num2str(aboveZeroIdx(idx(1)))]);
    %Ez = Ez .* cavityLocs;
    %Ez = reshape(Ez, [], 1);
    %normFac = sum(Ez .* Ez) * dx*dy;
    %Ez = Ez/sqrt(normFac);
    %sum(Ez .* Ez) * dx*dy
    
    %assert(false);
    %size(Ez)
    %size(cavityLocs)
    
    %EzI = reshape(EzMat(:,idx(2)),length(xpts),length(ypts));
    %EzI = reshape(EzMat(:,1),length(xpts),length(ypts));
    figure(1); imagesc(abs(Ez)); axis square;
    figure(2); imagesc(abs(Ez.*cavityLocs)); axis square;
    figure(3); imagesc(abs(Ez.*pumpProfile)); axis square;
    %imagesc(abs(Ez .* cavityLocs));
    %assert(false);
    
    chiMat = zeros(N);

    for ii=1:N
        EzI = EzMat(:,ii);
        
        for jj=1:N
            EzJ = EzMat(:,jj);
            
            chi = sum(EzI .* EzI .* EzJ .* conj(EzJ)) *dx*dy;
            %[ii,jj,chi,real(chi)]
            
            chiMat(ii,jj) = real(chi);
        end
    end
    
    save([datadir,'overlap_matrix_wFvec_partial_v2.mat'],'chiMat','fvec');    
    
    %figure(1);clf;
    %stem(lambda,Q);

    %frame_count=1;
    %for ii=1:length(Q)
        
    %if(Q(ii)>800)
    %Ez = dlmread([folder 'Ez_sol' num2str(ii)]);
            
    %fig1=figure(1);clf;
    
    %subplot 121;
    %stem(lambda.*1e3,Q,'b'); 
    %hold on; 
    %stem(lambda(ii).*1e3,Q(ii),'r','LineWidth',2); 
    %title(['Q=' num2str(Q(ii))]);xlabel('Wavelength (nm)');
    %ylabel('Quality Factor');
    %xlim([min(lambda*1e3) max(lambda*1e3)]);
    
    %subplot 122;
    %imagesc(xx,yy,abs(Ez)); ylabel('\mum');xlabel('\mum');title('Ez');axis equal tight;
    
    %end
    %end

end
