function [states] = RK4(FUNC, x, states_0, nu, U, U_d)

states = nan(length(states_0), length(x));
states(:, 1) = states_0;
for n = 1:length(x)-1
    h = x(n+1) - x(n);
    K1 = h*FUNC(x(n), states(:, n), nu, U(n), U_d(n));
    K2 = h*FUNC(x(n)+h/2, states(:, n)+0.5*K1, nu, (U(n)+U(n+1))/2, (U_d(n)+U_d(n+1))/2); % this needs fixing
    K3 = h*FUNC(x(n)+h/2, states(:, n)+0.5*K2, nu, (U(n)+U(n+1))/2, (U_d(n)+U_d(n+1))/2);
    K4 = h*FUNC(x(n)+h, states(:, n)+K3, nu, U(n+1), U_d(n+1));
    states(:, n+1) = states(:, n) + 1/6*(K1 + 2*K2 + 2*K3 + K4);
    if(isnan(states(:, n+1)))
        states(:, n+1:end) = nan(size(states(:, n+1:end)));
        break;
    end
end

end % endfunction
