function [x, y] =  LOAD_AIRFOIL(filename, AoA, N, closed_airfoil, method)

if nargin < 5
    method = "xfoil";
end
if nargin < 4
    closed_airfoil = 0;
end
if nargin < 3
    N = 100;
end

if strcmp(method, "xfoil")
  if(any([strfind(filename, "NACA"),  strfind(filename, "naca")]))
    fileID = fopen('xfoil_input','w');
    fprintf(fileID, strcat("naca ",strrep(filename, "naca", ''),"\n"));
    fprintf(fileID, strcat("PPAR\nN ",num2str(N),"\n\n\nPSAV ","./MY_AIRFOILS/",filename,".dat\n"));
    if (exist(strcat(filename, ".dat"), "file"))
      fprintf(fileID,'y\n'); % Overwrite existing file
    end
    fprintf(fileID, "\n\nQuit\n\n");
    fclose(fileID);
    system("xfoil < xfoil_input >/dev/null");
    delete("xfoil_input");
  end

  data = importdata(strcat("./MY_AIRFOILS/", filename, ".dat"));
  x = data(:, 1).';
  y = data(:, 2).';
  x_u = x(y>=0);
  x_l = x(y<=0);
  y_u = y(y>=0);
  y_l = y(y<=0); % this should be <= because I omit the repeated point later
elseif strcmp(method, "cosine")
  if mod(N, 2)
    warning("Please provide number of points as an even number\n");
  end
  m = str2double(filename(end-3)) / 100;
  p = str2double(filename(end-2)) / 10;
  t = str2double(filename(end-1:end)) / 100;

  NSIDE = ceil(N/2);
  theta = linspace(0, pi, NSIDE);
  x = (1 - cos(theta)) / 2;

  y_t = 5*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4);
  y_c = nan(1, NSIDE);
  y_c(x<p) = m/p^2 * (2*p*x(x<p) - x(x<p).^2);
  y_c(x>=p) = m/(1-p)^2 * ((1-2*p) + 2*p*x(x>=p) - x(x>=p).^2);
  dyc_dx = nan(1, NSIDE);
  dyc_dx(x<p) = 2*m/p^2 * (p-x(x<p));
  dyc_dx(x>=p) = 2*m/(1-p)^2 * (p-x(x>=p));
  theta = atan(dyc_dx);
  y_u = y_c + y_t.*cos(theta);
  y_l = y_c - y_t.*cos(theta);
end

if(y_l(end) == 0)
    warning("This might produce number of points N-1, this needs fixing");
    x_l = x_l(y_l(1:end-1));
    y_l = y_l(y_l(1:end-1));
end

x = [x_u, x_l].';
y = [y_u, y_l].';

%% Close the airfoil by adding final point
if(closed_airfoil)
  y(1) = 0;
  y(end) = 0;
end


end % endfunction

