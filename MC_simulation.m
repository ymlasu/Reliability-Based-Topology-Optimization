clear all;clc;
% % % nelx = 240;
% % % nely = 120;
% % % penal = 3;
% % % rmin = 2;
% % % u_threshold = [35 15];
% % % u_threshold = [42 19];
% % % u_threshold = [46 17];

nelx = 160;
nely = 100;
penal = 3;
rmin = 2;
u_threshold = 18;
%load('PMA_th18_300_R3_final.mat', 'xPhys_all_after');
load('IS_PMA51_300_th18_R3_sqp.mat', 'xPhys_all');
%xPhys_all_after = xPhys_all_after(:,151:300);
xPhys_all_after = xPhys_all(:,51:end);
parpool(24);
for i = 1:size(xPhys_all_after,2)
	%i = i + 50;
   [Reliability_MC(i),MC_glimit_sam(:,i)] = reliability_mc(nelx,nely,penal,rmin,xPhys_all_after(:,i),u_threshold);
end
   save('results_IS_MC_300_R3_sqp.mat','Reliability_MC','MC_glimit_sam');
   
