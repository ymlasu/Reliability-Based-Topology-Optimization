function dw_dksi = sensitivity_density_random(u)
                global eigenVec eigenVal_GMP GMP_quantile GMP_ygrid;

        GMP_Marginal_icdf = @(y)interp1(GMP_quantile,GMP_ygrid,y,'pchip');
        % eigenVec = eigenVec(:,end-419:end);
        % eigenVal_GMP = eigenVal_GMP(end-419:end);
        %%% The stochastic mesh is the transpose of the classical finite element
        %%% mesh.
        
        dw_dksi = eigenVec.*(repmat(sqrt(eigenVal_GMP)',size(eigenVec,1),1));
        Num_ksi = length(eigenVal_GMP);
        
        
        % n = 21390;
        %
        % n1 = 10680;
        mean=zeros(Num_ksi,1);
        sigma=ones(Num_ksi,1);
        
        u=u.*sigma+mean;
        
        dw_dksi = Inv_GM(u(1:Num_ksi),GMP_Marginal_icdf,eigenVec,eigenVal_GMP)*dw_dksi;  %%%Inv_GM(u) is the Gaussian mixture inverse CDF at u.
        
    end
