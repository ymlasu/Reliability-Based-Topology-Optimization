 function E = compute_E(u)
%         global eigenVec eigenVal_GMP GMP_Marginal_Var Vomlume_Frac GMP_Marginal_icdf
%         load('data_calculated.mat');
%         load('/home/ygao115/Documents/MATLAB/Random filed/Translation/validation_case2.mat','eigenVal_GMP','eigenVec','GMP_Marginal_icdf','GMP_Marginal_Var','Vomlume_Frac');
load('microstructure_data_RBTO_small.mat');
        Num_grid = size(eigenVec,1)^0.5;
        Num_ksi = length(eigenVal_GMP);
        Num_R = Num_grid;
        mean=zeros(Num_ksi,1);
        sigma=ones(Num_ksi,1);
        
        u=u.*sigma+mean;
        E= eigenVec * (sqrt(eigenVal_GMP).*u(1:Num_ksi));
        GMP_Marginal_icdf = @(y)interp1(GMP_quantile,GMP_ygrid,y,'pchip');
        E = GMP_Marginal_icdf(normcdf(E));
        E = E.*1e9;
        E = E/1e11;
    end