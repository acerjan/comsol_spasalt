%%%%%%%%%%%%%%%%%%
% plot_mode.m
% Alex Cerjan
% 2.23.16
%%%%%%%%%%%%%%%%%%

function [] = plot_mode(datadir,modeIdx)
    load([datadir,'Ez_sol',num2str(modeIdx),'.mat'],'Ez');
    imagesc(abs(Ez));
    colormap(jet);
end
