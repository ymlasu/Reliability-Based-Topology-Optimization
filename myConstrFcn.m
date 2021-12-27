function [cneq, ceq, gradc, gradceq] = myConstrFcn(u,relia)
%         % nonLinear inequalities Constraints
%         cneq     =  [];
%         gradc = [];
%         
%         % nonLinear equalities Constraints
%         
%         ceq =  sqrt(sum(u.^2)) - norminv(relia) ;
%         gradceq = abs(u)./sqrt(sum(u.^2));



% % %     cneq     =  sqrt(sum(u.^2)) - norminv(relia) ;
% % %       gradc = abs(u)./sqrt(sum(u.^2));
% % %       ceq = [];
% % %       gradceq = [];     


cneq     =  sum(u.^2) - norminv(relia)^2 ;
      gradc = 2*u;
      ceq = [];
      gradceq = [];  


    end % mycon
