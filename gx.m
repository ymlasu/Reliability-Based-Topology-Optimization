function [f, gf] = gx(x)

f = sqrt(sum(x.^2));
gf = abs(x)./sqrt(sum(x.^2));


end