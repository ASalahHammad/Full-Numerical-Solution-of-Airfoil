function [lambda] = find_lambda(K, lambda0, N)
if nargin < 2
    lambda0 = 7;
end

epsilon = 1;
cnt = 0;
x = lambda0;

% newton-raphson
while(epsilon>1e-10)
if(cnt > 1000)
    warning(strcat("N = ", num2str(N), ", too many iterations to calculate Lambda, this might indicate instability\n"));
    break;
end

%    f = (37/315)^2*x - 2*37/315/945*x^2 + (-2*37/315/9072 + 1/945^2)*x^3 + 2/945/9072*x^4 + 1/9072^2*x^5 - K;
%    fd = (37/315)^2 - 2*2*37/315/945*x + 3*(-2*37/315/9072 + 1/945^2)*x^2 + 4*2/945/9072*x^3 + 5/9072^2*x^4;
    cnt = cnt + 1;
    f = (37/315 - 1/945*x - 1/9072*x^2)^2 * x - K;
    fd = (37/315 - 1/945*x - 1/9072*x^2)^2 + 2*(37/315 - 1/945*x - 1/9072*x^2)*(-1/945 - 2/9072*x)*x;
    lambda = x - f/fd;
    epsilon = abs(lambda-x);
    x = lambda;
end
end



