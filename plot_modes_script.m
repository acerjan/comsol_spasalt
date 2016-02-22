clear;clc;



% folder = 'data/disk_R5um_D_at50_neff3.5/'
% xx = -6:0.01:6; yy=xx;

% folder = 'disk_R5um_D_at50_neff3.5/'
% xx = -6:0.01:6; yy=xx;

folder = '~/Data/2d_salt/Qcav40_R4_dr22_L002/comsol_results/';
xx = -5:0.01:5; yy=xx;

tmp = dlmread([folder 'lambda_Q']);
lambda = tmp(:,1);
Q = tmp(:,2);

figure(1);clf;
stem(lambda,Q);

frame_count=1;
for ii=1:length(Q)
    
    if(Q(ii)>0)
        Ez = dlmread([folder 'Ez_sol' num2str(ii)]);
        
        fig1=figure(1);clf;
        subplot 121;
        stem(lambda.*1e3,Q,'b'); hold on; stem(lambda(ii).*1e3,Q(ii),'r','LineWidth',2); title(['Q=' num2str(Q(ii))]);xlabel('Wavelength (nm)'); ylabel('Quality Factor');xlim([min(lambda*1e3) max(lambda*1e3)]);
        subplot 122;
        imagesc(xx,yy,abs(Ez)); ylabel('\mum');xlabel('\mum');title('Ez');axis equal tight;
       
    end
end

return;