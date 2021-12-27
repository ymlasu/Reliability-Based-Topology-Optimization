
clear all;close all;clc;

nelx = 160;
nely = 100;
penal = 3;
rmin = 2;
ft = 1;
u_threshold = 18;
DIM = 187;


%% MATERIAL PROPERTIES
global eigenVec  eigenVal_GMP GMP_Marginal_Var  GMP_quantile GMP_ygrid;
load('microstructure_data_RBTO_small.mat');
global KE Emin  iK jK sK freedofs U_unit_virtual_force F;

beta = 3;
relia = normcdf(beta);
pf_constraint = normcdf(-beta);

% Emin = 1e3;
Emin = 1e-9;

nu = 0.3;
constr_index_1 = 2*(nelx+1)*(nely+1);
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% Define loads and supports (cantilever with two point loads)
% F = sparse(2*(nelx+1)*(nely+1)-1:2*(nelx+1)*(nely+1),1,[-0.5 -1],2*(nely+1)*(nelx+1),1);
% F = sparse(2*(nely+1)*(nelx+1),1,-1e6);
F = sparse(2*(nely+1)*(nelx+1),1,-1);


fixeddofs = 1:2*(nely+1);
displacement_constraint_dofs = 2*(nelx+1)*(nely+1);
U = zeros(2*(nely+1)*(nelx+1),1);
U_unit_virtual_force = zeros(2*(nely+1)*(nelx+1),1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION



%% IS_right

load('PMA_th18_300_R3_final.mat', 'xPhys_all_after');
%load('PMA_th18_300_R3.mat', 'xPhys_all_after');
%load('IS_PMA50_200_th18_R3.mat', 'xPhys_all');

xPhys = reshape(xPhys_all_after(:,50),nely,nelx);
x = xPhys;
loop = 50;
loop_ini = loop;
clear xPhys_all_after;


change = 1;
lam = 1e10;
F_unit_virtual_force =  sparse(constr_index_1,1,-1);

%% initial step reliability constraint
A = [];
B = [];
Aeq = [];
Beq = [];
LB = [];
UB = [];
maxloop = 200;    % Maximum number of iterations
tolx = 0.005;      % Terminarion criterion

% OPTIONS = optimset('TolX',tolx, 'MaxIter',maxloop, 'Algorithm','sqp',...
%     'GradObj','on', 'GradConstr','on', ...
%     'Display','none', 'OutputFcn',@(x,optimValues,state) myOutputFcn(optimValues,state), 'PlotFcns',@optimplotfval);

% OPTIONS = optimset('TolX',tolx, 'Algorithm','sqp',...
%     'GradObj','on', 'GradConstr','on', ...
%     'Display','off', 'OutputFcn',@(x,optimValues,state) myOutputFcn(optimValues,state));



OPTIONS = optimoptions('fmincon');
OPTIONS.SpecifyConstraintGradient=true;
OPTIONS.SpecifyObjectiveGradient =true;
OPTIONS.Algorithm = 'sqp';
% OPTIONS.OptimalityTolerance = 1e-10;
OPTIONS.FunctionTolerance = 5e-3;
OPTIONS.ConstraintTolerance = 2e-3;
OPTIONS.OptimalityTolerance = 2e-3;
OPTIONS.MaxIterations = 70;

OPTIONS.Display = 'iter';

u0 = 0.1*ones(length(eigenVal_GMP),1);

%     u0  = unifrnd(-1,1,100,1);


% [u_solu,displace_constraint,flag,output] = fmincon(@(u)myObjFcn(xPhys,u,u_threshold,nelx, nely,penal), u0, A, B, Aeq, Beq, LB, UB, @(u)myConstrFcn(u,relia), OPTIONS);
%
% E0 = compute_E(u_solu);

%%
%%%% git maker
% filename = 'topology_min_compliance.gif';
h = figure;
% % %    colormap(jet);
% %    imagesc((xPhys).*E1);
% %    colorbar;
% %     caxis([0 1]);
% %     axis equal; axis off; drawnow;



%% START ITERATION

N_Sam_IS = 2e3;
sigma_importance_func = 1e-3;
fval = -1;
global ConstraintsCount;

%%% IS samples

%parpool(24);

pool = parpool(str2num(getenv('SLURM_NTASKS')));
while (change > 1e-2 && loop < 300 )
    loop = loop + 1;
    
%Sam = lhsdesign(N_Sam_IS,DIM);
%Sam = norminv(Sam);
%Sam = Sam';
Sam = normrnd(0,1,DIM,N_Sam_IS);
    %% Displace constraint calculus / Relibility constraints solved via PMA.
   
    A = [];
    B = [];
    Aeq = [];
    Beq = [];
    LB = [];
    UB = [];
     
    %% FORM for MPP
    option_form = optimoptions('fmincon');
    option_form.SpecifyConstraintGradient=true;
    %option_form.Algorithm = 'active-set';
    option_form.Algorithm = 'sqp';
    option_form.OptimalityTolerance = 2e-3;
    option_form.FunctionTolerance = 5e-3;
    option_form.ConstraintTolerance = 2e-3;
    option_form.Display = 'iter';
    option_form.MaxIterations = 60;
    u0_form = 0.1*ones(DIM,1);
%     u0_form = unifrnd(0,1,DIM,1);

    
        ConstraintsCount = 0;
        [u_form_check,R_index,flag_form,output_form] = fmincon(@(u) gx(u), u0_form, A, B, Aeq, Beq, LB, UB, @(u)nonlcon(xPhys,u,u_threshold,nelx, nely,penal,constr_index_1,F_unit_virtual_force), option_form);
    u_form_check_xphy_all(:,loop) = u_form_check;
    R_index_all(:,loop) = R_index;
    ConstraintsCount_all(:,loop) = ConstraintsCount;
    
    %% Importance sampling
    Sam_Im = Sam + u_form_check;
    E0_IS = [];
    U_inter = [];
    U_unit_virtual_force_IS = [];
    glimit_IS = [];
    d_ro_IS = [];
    
    parfor i = 1: N_Sam_IS
        E0_IS(:,i) = compute_E(Sam_Im(:,i));
        %         E0 =  E0_IS(:,i) ;
        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal.*(E0_IS(:,i)'-Emin)),64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U_inter(:,i) = K(freedofs,freedofs)\F(freedofs);
        U_unit_virtual_force_IS(:,i) = K(freedofs,freedofs)\F_unit_virtual_force(freedofs);
    end
    
    
    
    for i = 1: N_Sam_IS
        U(freedofs) = U_inter(:,i);
        glimit_IS(1,i) = -U(end,1) - u_threshold;
        U_unit_virtual_force(freedofs) = U_unit_virtual_force_IS(:,i);
        ce = reshape(sum((U_unit_virtual_force(edofMat)*KE).*U(edofMat),2),nely,nelx);
        E0_matrix = reshape(E0_IS(:,i),nely,nelx);
        d_ro = -penal.*(E0_matrix-Emin).*xPhys.^(penal-1).*abs(ce);
        d_ro_IS(:,i) = d_ro(:);
    end
    
    
    loc = find(glimit_IS>=0);
    pf_IS_sample = exp(-0.5*(sum(Sam_Im(:,loc).^2)-sum((Sam_Im(:,loc)-u_form_check).^2)));
    pf_IS = sum(pf_IS_sample)/N_Sam_IS;
    glimit_IS_indicator = zeros(1,N_Sam_IS);
    glimit_IS_indicator(loc) = pf_IS_sample;
    pf_IS_cov = abs(sqrt(var(glimit_IS_indicator)/N_Sam_IS)./mean(glimit_IS_indicator,2));
    pf_IS_cov_all(loop,1) = pf_IS_cov;
    
    
    %     sigma_importance_func_matrix = linspace(1e-4,5,5e3);
    %     pf_IS_cov_matrix = [];
    % for sigma_importance_func= sigma_importance_func_matrix
    %     pf_IS_importance_func_sample=exp(-0.5*(sum(Sam_Im.^2)-sum((Sam_Im-u_form_check).^2))).*normcdf(glimit_IS/sigma_importance_func);
    %     pf_IS_cov_matrix(end+1)=abs(sqrt(var(pf_IS_importance_func_sample)/N_Sam_IS)./mean(pf_IS_importance_func_sample,2));
    %
    % end
    %
    % [~,min_cov] = min(pf_IS_cov_matrix);
    % pf_IS_importance_func_sample=exp(-0.5*(sum(Sam_Im.^2)-sum((Sam_Im-u_form_check).^2))).*normcdf(glimit_IS/sigma_importance_func_matrix(min_cov));
    
    pf_IS_importance_func_sample=exp(-0.5*(sum(Sam_Im.^2)-sum((Sam_Im-u_form_check).^2))).*normcdf(glimit_IS/sigma_importance_func);
    pf_IS_importance_func = sum(pf_IS_importance_func_sample)/N_Sam_IS;
    
    pf_IS_func_cov = abs(sqrt(var(pf_IS_importance_func_sample)/N_Sam_IS)./mean(pf_IS_importance_func_sample,2));
    pf_IS_func_cov_all(loop,1) = pf_IS_func_cov;
    pf_IS_all(loop,1) = pf_IS;
    pf_IS_importance_func_all(loop,1) = pf_IS_importance_func;
    
    
    
    % % %     [u_solu_form,displace_constraint_form,flag,output] = fmincon(@(u)myObjFcn(xPhys,u,u_threshold,nelx, nely,penal), u_form_check, A, B, Aeq, Beq, LB, UB, @(u)myConstrFcn(u,relia), OPTIONS);
    % % %     displace_constraint_form = -displace_constraint_form;
    % % %     displace_constraint = displace_constraint_form;
    % % %     u_solu = u_solu_form;
    % % %     u_solu_form_all(:,loop) = u_solu_form;
    % % %     displace_constraint_all(:,loop) = displace_constraint;
    
    %     [Reliability_MC,MC_glimit(:,loop),u_MC] = reliability_mc(nelx,nely,penal,rmin,xPhys,u_threshold);
    %     pf_MC_all(loop,1)  = 1-Reliability_MC;
    
    
    % % %     if displace_constraint >= 0
    % % %         disp('The relaiblity constraint is not satisfied using PMA ');
    % % %     else
    % % %     end
    
    
    
    
    dv = ones(nely,nelx)/nelx/nely;
    
    
    %% sensitivity of reliability to deterministic design variables (density)
    

    
    
     grad_IS_importance_func_sample = exp(-0.5*(sum(Sam_Im.^2)-sum((Sam_Im-u_form_check).^2)))./(sigma_importance_func).*normpdf(glimit_IS./(sigma_importance_func)).*d_ro_IS;
    
    dc_vector_IS = sum(grad_IS_importance_func_sample,2)/N_Sam_IS;
    cov_IS_optimal = abs(sqrt(var(grad_IS_importance_func_sample,0,2)/N_Sam_IS)./mean(grad_IS_importance_func_sample,2));
    cov_IS_optimal_all(:,loop) = cov_IS_optimal;
    dc = reshape(dc_vector_IS,nely,nelx);
    
    
    f0val = mean(xPhys(:));
    fval_IS_all(loop,1) = pf_IS-pf_constraint;
    fval = pf_IS_importance_func-pf_constraint;
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
	dc = -abs(dc);
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    %     l1 = 0; l2 = 2*lam; move = 0.05 ; i = 0;
    l1 = 0; l2 = 1e10; move = 0.025 ; i = 0;
    while (l2-l1)/(l1+l2) > 1e-6
        % Check if uniform reduction violates linearized constraint
       % if (pf_IS_importance_func-pf_constraint - sum(sum(dc))*move < 0)
       %     xnew = x - move;
       %     lam = 1e10;
       %     break;
       % end
        
        
        % Non-uniform design change
        i = i + 1;
        lam = 0.5*(l2+l1);
	xnew = max(1e-10,min(1,x.*((-lam*dc./dv).^0.5)));
        %xnew = max(1e-10,max(x-move,min(1,min(x+move,x.*((-lam*dc./dv).^0.5)))));
        
        if  pf_IS_importance_func-pf_constraint  + dc(:)'*((xnew(:)-x(:)).*(x(:)./xnew(:))) > 0
            l1 = lam;
        else l2 = lam;
        end
        
        
        
        
    end
    
        Lagrainge_multi_all(:,loop) = pf_IS_importance_func-pf_constraint   + dc(:)'*((xnew(:)-x(:)).*(x(:)./xnew(:)));
%     Lagrainge_multi_all(:,loop) = pf_IS_importance_func-pf_constraint   + dc(:)'*((xnew(:)-x(:)));
    
    
    i_all(loop,1) = i;
    lam_all(loop,1) = lam;
    %     ConstraintsCount_all(loop) = ConstraintsCount;
    % % %     %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    % % %   l1 = 0; l2 = 1e9; move = 0.5;i=0;
    % % %   while (l2-l1)/(l1+l2) > 1e-6
    % % %       i = i + 1;
    % % %     lmid = 0.5*(l2+l1);
    % % %     xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    % % %     if ft == 1
    % % %       xPhys = xnew;
    % % %     elseif ft == 2
    % % %       xPhys(:) = (H*xnew(:))./Hs;
    % % %     end
    % % %     if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
    % % %   end
    
    
    
    change = full(max(abs(xnew(:)-x(:))));
    x = xnew;
    if ft == 1
        xPhys = xnew;
    elseif ft == 2
        xPhys(:) = (H*xnew(:))./Hs;
    end
    %% PRINT RESULTS
    fprintf('********************Topology updating begin*******************\n');
    fprintf(' Topology It.:%5i Con.:%11.3e Vol.:%11.3e Lam.: %11.3e ch.:%7.3f innerit: %3i\n',...
        loop,fval,f0val,lam,change,i);
    fprintf('********************Topology updating end*******************\n');
% % %     %% PLOT DENSITIES
% % %     %         colormap(gray);
% % %     
% % %     imagesc((xPhys));
% % %     %     colorbar;
% % %     caxis([0 1]);
% % %     axis equal; axis off; drawnow;
% % %     %         colormap(jet); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
% % %     
    %%
    % %%%% git maker
    %   % Capture the plot as an image
    %       frame = getframe(h);
    %       im = frame2im(frame);
    %       [imind,cm] = rgb2ind(im,256);
    %       % Write to the GIF File
    %       if loop == 1
    %           imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf);
    %       else
    %           imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    %       end
    
    xPhys_all(:,loop) = xPhys(:);
    fval_all(loop,1) = fval;
    
end
save('IS_PMA_cluster_test.mat','xPhys_all','fval_all','ConstraintsCount_all','pf_IS_importance_func_all');

%% Save results
% save(filename,'f0val','xPhys');

% save('rslt_R_1_8.mat','xPhys','fval','f0val');
% title('R=1.8 dis=1.5e-4');







