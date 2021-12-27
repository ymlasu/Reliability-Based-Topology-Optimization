function f = Inv_GM(u,GMP_Marginal_icdf,eigenVec,eigenVal_GMP)
u = eigenVec * (sqrt(eigenVal_GMP).*u);
% step = u./1e3;
step = 1e-4;

f = (GMP_Marginal_icdf(normcdf(u+step))-GMP_Marginal_icdf(normcdf(u)))./(step);


f = diag(f);
f = sparse(f);
end
