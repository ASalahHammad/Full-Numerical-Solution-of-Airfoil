function [U, U_d, U_dd] = velocity(XC, V, x, method)

if nargin < 4
    method = "CD4";
end

N = length(x);

U = interp1(XC, V, x, "spline", "extrap");

if strcmp(method, "CD4")

  dx = x(2) - x(1);

  x_ghost = [x(1)-2*dx; x(1)-dx; x; x(end)+dx; x(end)+2*dx];

  U_ghost = zeros(length(U)+4, 1);
  U_ghost(3:length(U)+2) = U;
  U_ghost(1:2) = interp1(x(1:5), U(1:5), x_ghost(1:2), "spline", "extrap");
  U_ghost(end-1:end) = interp1(x(end-5:end), U(end-5:end), x_ghost(end-1:end), "spline", "extrap");

  U_d = (-U_ghost(5:end) + 8*U_ghost(4:end-1) - 8*U_ghost(2:end-3) + U_ghost(1:end-4)) / (12*dx);

  U_d_ghost = zeros(length(U_d)+4, 1);
  U_d_ghost(3:length(U_d)+2) = U_d;
  U_d_ghost(1:2) = interp1(x(1:5), U_d(1:5), x_ghost(1:2), "spline", "extrap");
  U_d_ghost(end-1:end) = interp1(x(end-5:end), U_d(end-5:end), x_ghost(end-1:end), "spline", "extrap");

  U_dd = (-U_d_ghost(5:end) + 8*U_d_ghost(4:end-1) - 8*U_d_ghost(2:end-3) + U_d_ghost(1:end-4)) / (12*dx);

%  h_m2 = x - x_ghost(1:N);
%  h_m1 = x - x_ghost(2:N+1);
%  h_p1 = x_ghost(4:N+1) - x;
%  h_p2 = x_ghost(5:N+2) - x;
%  term1 = (h_p2.^2 * U_ghost(4:N+1) - h_p1.^2 * U(5:N+2)) ./ (h_p1(3:N) .* h_p2(3:N) .* (h_p1(3:N) + h_p2(3:N)));
%  term2 = (h_m2.^2 * U(2:N-1) - h_m1.^2 * U(1:N)) ./ (h_m1(3:N) .* h_m2(3:N) .* (h_m1(3:N) + h_m2(3:N)));
%  U_d = term1 - term2;

else
error("method not implemented yet!\n");
end

end % endfunction

