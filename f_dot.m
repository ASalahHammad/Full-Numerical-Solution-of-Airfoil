function [deriv] = f_dot(x, states, nu, U, U_d)

deriv = nan(2, 1);
H = HofH1(states(2)/states(1)/U);
if(H>=3)
% if(~(imag(H)==0))
    return
end

H1 = G(H);
Rtheta = U*states(1)/nu;
cf = 0.246*10^(-0.678*H)*Rtheta^(-0.268);
deriv(1) = -states(1)/U*U_d*(H+2) + cf/2;
deriv(2) = U * F(H1);

end % endfunction
