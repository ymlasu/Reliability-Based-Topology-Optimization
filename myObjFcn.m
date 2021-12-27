function [f, gradf] = myObjFcn(xPhys,u,u_threshold,nelx, nely,penal,index,F_unit_virtual_force)

        global KE Emin iK jK sK freedofs  U_unit_virtual_force F;
        global ConstraintsCount;
        ConstraintsCount  = ConstraintsCount +1;
        
        E0 = compute_E(u);
        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal.*(E0(:)'-Emin)),64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U(freedofs) = K(freedofs,freedofs)\F(freedofs);
        
         U_unit_virtual_force(freedofs) = K(freedofs,freedofs)\F_unit_virtual_force(freedofs);
        U = U(:);
        performance_constraint = -( -U(index,1) - u_threshold);
        
        
        dw_dksi = sensitivity_density_random(u);
        
        
        
        parfor k = 1:size(dw_dksi,2)
            
            grad_sK = reshape(KE(:)*(xPhys(:)'.^penal).*dw_dksi(:,k)',64*nelx*nely,1);
            
            grad_K_ksi = sparse(iK,jK,grad_sK); grad_K_ksi = (grad_K_ksi+grad_K_ksi')/2;
            
            
            grad_u_ksi(k,1) = -U_unit_virtual_force(freedofs)'*grad_K_ksi(freedofs,freedofs)*U(freedofs);
            
        end
        
        
        f = performance_constraint;
        gradf = -grad_u_ksi;
%         gradf = gradf*1e9;  %%%% 2e11 Pa
        gradf = gradf*1e-2;  %%%% 200 GPa
        
    end % myfun