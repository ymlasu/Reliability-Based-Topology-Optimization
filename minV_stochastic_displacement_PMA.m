
clear all;clc;

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
beta = 3.2;
relia = normcdf(beta);
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

%%% PMA_for_IS_right_1
%%load('PMA_th18_1_100.mat', 'xPhys_all_after');
% loop = 100;
% xPhys = reshape(xPhys_all_after(:,100),nely,nelx);
% x = xPhys;
% clear xPhys_all_after;





x = ones(nely,nelx);
xPhys = x;
loop = 0;


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
OPTIONS.Algorithm = 'active-set';
% OPTIONS.OptimalityTolerance = 1e-10;
OPTIONS.FunctionTolerance = 5e-3;
OPTIONS.ConstraintTolerance = 2e-3;
OPTIONS.OptimalityTolerance = 2e-3;
OPTIONS.MaxIterations = 60;
OPTIONS.Display = 'iter';

u0 = 0.1*ones(length(eigenVal_GMP),1);

%     u0  = unifrnd(-1,1,100,1);


% [u_solu,displace_constraint,flag,output] = fmincon(@(u)myObjFcn(xPhys,u,u_threshold,nelx, nely,penal), u0, A, B, Aeq, Beq, LB, UB, @(u)myConstrFcn(u,relia), OPTIONS);
%
% E0 = compute_E(u_solu);

%%
%%%% git maker
% filename = 'topology_min_compliance.gif';
% h = figure;
% % %    colormap(jet);
% %    imagesc((xPhys).*E1);
% %    colorbar;
% %     caxis([0 1]);
% %     axis equal; axis off; drawnow;

%% START ITERATION
fval = -1;
global ConstraintsCount;
parpool(24);
while (change > 1e-2 && loop < 100 )
    %     for loop = 37:58
    loop = loop + 1;
    
    
    %% Displace constraint calculus / Relibility constraints solved via PMA.
    % displace_constraint = -U(end,1);
    % dw_dksi = sensitivity_density_random(u);
    %
    % for i = 1:size(dw_dksi,2)
    % grad_sK = reshape(KE(:)*(xPhys(:)'.^penal).*dw_dksi(:,i),64*nelx*nely,1);
    %   grad_K_ksi = sparse(iK,jK,grad_sK); grad_K_ksi = (grad_K_ksi+grad_K_ksi')/2;
    %
    %   grad_u_ksi(i,1) = -U_unit_virtual_force(freedofs)'*grad_K_ksi(freedofs,freedofs)*U(freedofs);
    % end
    A = [];
    B = [];
    Aeq = [];
    Beq = [];
    LB = [];
    UB = [];
    
    
    displayflag = 0;  % Display structure flag
    u_form_check = u0;
    
    
    ConstraintsCount = 0;
    [u_solu_form,displace_constraint_form,flag,output] = fmincon(@(u)myObjFcn(xPhys,u,u_threshold,nelx, nely,penal,constr_index_1,F_unit_virtual_force), u_form_check, A, B, Aeq, Beq, LB, UB, @(u)myConstrFcn(u,relia), OPTIONS);
    displace_constraint_form = -displace_constraint_form;
    displace_constraint = displace_constraint_form;
    u_solu = u_solu_form;
    u_solu_form_all(:,loop) = u_solu_form;
    ConstraintsCount_all(:,loop) = ConstraintsCount;
    
    
    E0 = compute_E(u_solu);
    
    
    u_solu_all(:,loop) = u_solu;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal.*(E0(:)'-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    U_unit_virtual_force(freedofs) = K(freedofs,freedofs)\F_unit_virtual_force(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U_unit_virtual_force(edofMat)*KE).*U(edofMat),2),nely,nelx);
    %     displace_constraint = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    % displace_constraint = F_unit_virtual_force'*U;
    E0_matrix = reshape(E0,nely,nelx);
    dc = -penal.*(E0_matrix-Emin).*xPhys.^(penal-1).*abs(ce);
    dv = ones(nely,nelx)/nelx/nely;
    %%
    %     if (loop == 1 && nargin < 7)
    %         compconst = 2*displace_constraint;
    %     end
    f0val = mean(xPhys(:));
    fval = displace_constraint;
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    %     l1 = 0; l2 = 2*lam; move = 0.05 ; i = 0;
    l1 = 0; l2 = 1e10; move = 0.05 ; i = 0;
    while (l2-l1)/(l1+l2) > 1e-6
        % Check if uniform reduction violates linearized constraint
        if (displace_constraint - sum(sum(dc))*move < 0)
            xnew = x - move;
            lam = 1e10;
            break;
        end
        
        
        % Non-uniform design change
        i = i + 1;
        lam = 0.5*(l2+l1);
	xnew = max(1e-10,min(1,x.*((-lam*dc./dv).^0.5)));
        %xnew = max(1e-10,max(x-move,min(1,min(x+move,x.*((-lam*dc./dv).^0.5)))));

        
        
                    if displace_constraint  + dc(:)'*((xnew(:)-x(:)).*(x(:)./xnew(:))) > 0
       
        

                        l1 = lam;
                   else l2 = lam;
                   end
    
        
        
        
    end
    
    Lagrainge_multi_all(:,loop) = displace_constraint  + dc(:)'*((xnew(:)-x(:)).*(x(:)./xnew(:)));
    
    i_all(loop,1) = i;
    lam_all(loop,1) = lam;
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
    
    xPhys_all_before(:,loop) = xPhys(:);
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
% %     %% PLOT DENSITIES
% %     %         colormap(gray);
% %     
% %     imagesc((xPhys));
% %     %     colorbar;
% %     caxis([0 1]);
% %     axis equal; axis off; drawnow;
% %     %         colormap(jet); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    
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
    
    xPhys_all_after(:,loop) = xPhys(:);
    fval_all(loop,1) = fval;
end
save('PMA_th18_100_R32.mat','xPhys_all_after','xPhys_all_before','ConstraintsCount_all','fval_all','Lagrainge_multi_all');




