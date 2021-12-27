function [r,f] = reliability_mc(nelx,nely,penal,rmin,xPhys,u_threshold)
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

% % % % F = sparse(2*(nely+1)*(nelx+1),1);
% % % % % F(2*(60*121+1),1) = -1;
% % % % % F(2*(180*121+1),1) = -1;
% % % % 
% % % % 
% % % % F(2*(nelx+1)*(nely+1),1) = -1;
% % % % F(2*14641,1) = -1;


F = sparse(2*(nely+1)*(nelx+1),1);
F(constr_index_1,1) = -1;

fixeddofs = 1:2*(nely+1);

U = zeros(2*(nely+1)*(nelx+1),1);
U_unit_virtual_force = zeros(2*(nely+1)*(nelx+1),1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);


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
dim = 187;
N_sam = 1e4;
u = normrnd(0,1,dim,N_sam);
parfor i = 1: N_sam
%     u0 = u(:,i);
        E0 = compute_E(u(:,i));
%  E0 = compute_E(u0);

%% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal.*(E0'-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U_inter(:,i) = K(freedofs,freedofs)\F(freedofs);
%     U(freedofs) = K(freedofs,freedofs)\F(freedofs);
%     f(i,1) = -U(end,1) - u_threshold;
    
end




for i = 1: N_sam


    U(freedofs) = U_inter(:,i);
     ff(i,1) = -U(constr_index_1,1);
    f(i,1) = -U(constr_index_1,1)- u_threshold;
    
end

r = sum(f(:,1)<=0)/N_sam;
end
