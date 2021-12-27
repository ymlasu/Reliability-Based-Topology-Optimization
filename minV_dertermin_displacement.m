
% material property field volume fractio = 0.3;
function minV_dertermin_displacement(E0,nelx,nely,penal,rmin,ft,compconst)
%% MATERIAL PROPERTIES

Emin = 1e-9;
% Emin = 1e3;
nu = 0.3;
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
x = ones(nely,nelx);
xPhys = x;
loop = 0;
change = 1;
lam = 1e10;
F_unit_virtual_force =  sparse(2*(nelx+1)*(nely+1),1,-1);

%%
%%%% git maker
filename = 'topology_min_compliance_deter.gif';
h = figure;
% % %    colormap(jet); 
% %    imagesc((xPhys).*E1); 
% %    colorbar;
% %     caxis([0 1]); 
% %     axis equal; axis off; drawnow;
    
%% START ITERATION
while (change > 1e-2 && loop < 200)
    loop = loop + 1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal.*(E0(:)'-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
     U_unit_virtual_force(freedofs) = K(freedofs,freedofs)\F_unit_virtual_force(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U_unit_virtual_force(edofMat)*KE).*U(edofMat),2),nely,nelx);
%     displace_constraint = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
% displace_constraint = F_unit_virtual_force'*U;
displace_constraint = -U(end,1);
    dc = -penal.*(E0-Emin).*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx)/nelx/nely;
    if (loop == 1 && nargin < 7)
        compconst = 2*displace_constraint;
    end
    f0val = mean(xPhys(:));
    volume_all(loop,:) = f0val;
    fval = displace_constraint/compconst - 1;
    fval_all(loop,1) = fval;
     xPhys_all_before(:,loop) = xPhys(:);
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    l1 = 0; l2 = 1e10; move = 0.05 ; i = 0;
    while (l2-l1)/(l1+l2) > 1e-6
        % Check if uniform reduction violates linearized constraint
        if (displace_constraint - compconst - sum(sum(dc))*move < 0)
            xnew = x - move;
            lam = 1e10;
            break;
        end
        % Non-uniform design change
        i = i + 1;
        lam = 0.5*(l2+l1);
        xnew = max(1e-10,max(x-move,min(1,min(x+move,x.*((-lam*dc./dv).^0.5)))));
        if displace_constraint - compconst + dc(:)'*((xnew(:)-x(:)).*(x(:)./xnew(:))) > 0, l1 = lam; else l2 = lam; end
    end
    change = full(max(abs(xnew(:)-x(:))));
    x = xnew;
    if ft == 1
        xPhys = xnew;
    elseif ft == 2
        xPhys(:) = (H*xnew(:))./Hs;
    end
    %% PRINT RESULTS
    fprintf(' It.:%5i Con.:%11.3e Vol.:%11.3e Lam.: %11.3e ch.:%7.3f innerit: %3i\n',...
        loop,fval,f0val,lam,change,i);
    %% PLOT DENSITIES
%     colormap(jet); 
    imagesc((xPhys)); 
%     colorbar;
    caxis([0 1]); 
    axis equal; axis off; drawnow;
%         colormap(jet); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;

    %%
%%%% git maker
  % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if loop == 1 
          imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
      end 
      xPhys_all_after(:,loop) = xPhys(:);
end
%% Save results
save(filename,'f0val','xPhys');


%
