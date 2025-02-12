clc;
clear;
clear all;
%Parameters
c=3e8;
eps0 = 8.85e-12;
wavelength =linspace(1600,1650,5)*1e-9;
frequency=c ./wavelength;
wavenumber=2*pi ./wavelength;


n_core = 1.44792;
n_clad = 1.4;
n_medium=1;
epsilon_core = eps0 * n_core^2;
epsilon_clad = eps0 * n_clad^2;

az_num=1;

radius_core= (5.25 / 2) * 1e-6;
radius_clad=62.5e-6;


%solve the matrices for beta to get the neff 

syms neff_symbolic
beta = wavenumber * neff_symbolic;clc;
clear;

%Parameters
c=3e8;
eps0 = 8.85e-12;

wavelength = linspace(1600,1650,50)*1e-9;
frequency=c ./wavelength;
wavenumber=2*pi ./wavelength;


n_core = 1.44792;
n_clad = 1.4;
n_medium=1;
epsilon_core = eps0 * n_core^2;
epsilon_clad = eps0 * n_clad^2;

az_num=1;

radius_core= (5.25 / 2) * 1e-6;
radius_clad=62.5e-6;


%solve the matrices for beta to get the neff 

syms neff_symbolic real

searchpoints=50;

%{
for i = 1:length(wavenumber)
    beta_symbolic = wavenumber(i) * neff_symbolic;

    %matrix__calculator(prop_const,wave_num,mode,radius,epsilon)
    Mcore=matrix_core_calculator(beta_symbolic,wavenumber(i),az_num,radius_core,epsilon_core);

    disp(Mcore);
    det_N=det(Mcore);
    n_eff_solution=vpasolve(det_N == 0, neff_symbolic);
    disp(n_eff_solution);
end
%}


for i = 1:length(wavenumber)
    % Define determinant function as an anonymous function
    det_eq = @(neff) det(matrix_core_calculator(wavenumber(i) * neff, ...
                                                wavenumber(i), az_num, ...
                                                radius_core, epsilon_core));

    % Check sign at boundaries
    sign_start = sign(det_eq(n_clad));
    sign_end = sign(det_eq(n_core));

    if sign_start ~= sign_end
        % If signs are different, use fzero
        try
            neff_solution = fzero(det_eq, [n_clad, n_core]); 
        catch
            neff_solution = NaN; % If fzero fails, assign NaN
        end
    else
        % If no sign change, try minimizing |det_eq| using fminbnd
        try
            neff_solution = fminbnd(@(neff) abs(det_eq(neff)), n_clad+1e-3, n_core);
        catch
            neff_solution = NaN;
        end
    end

    % Store the solution
    neff_solutions(i) = neff_solution;
end

%disp(neff_solutions);

r_values = linspace(0, 2e-6, 10);
figure; hold on;
colors = lines(length(wavenumber)); % Get distinct colors for each wavenumber

for i = 1:length(wavenumber)
    if isnan(neff_solutions(i))
        continue; % Skip if neff was not found
    end
    
    % Compute E_r for each radius
    Er_values = arrayfun(@(r) E_calculator(neff_solutions(i), n_core, n_clad, r, radius_core, wavenumber(i)), r_values);
    
    % Plot Er vs. r
    plot(r_values * 1e6, abs(Er_values), 'Color', colors(i, :), 'LineWidth', 2, ...
         'DisplayName', sprintf('\\lambda = %.0f nm', wavelength(i) * 1e9));
end

xlabel('Radius r (\mum)');
ylabel('|E_r|');
title('Radial Electric Field Distribution');
legend('show');
grid on;
hold off;


function M = matrix_core_calculator(prop_const,wave_num,az_num,radius,epsilon)
    function dJ = besselj_derivative(n, x)
            dJ = (n/x)*besselj(n,x) - besselj(n+1,x);
        end
        u=sqrt((wave_num^2)*epsilon-prop_const^2);
        phi=prop_const*az_num/(wave_num);
        m11=(-u^2/epsilon)*besselj(az_num,u*radius);
        m31=(phi/radius*epsilon)*besselj(az_num,u*radius);
        m41=-1*u*besselj_derivative(az_num,u*radius);
        m23=(u^2)*besselj(az_num,u*radius);
        m33=(-1*u)*besselj_derivative(az_num,u*radius);
        m34=(phi/radius)*besselj(az_num,u*radius);

        M=[m11,0,0,0;
            0,0,m23,0;
            m31,0,m33,0;
            m41,0,m34,0;];
end


function Er=E_calculator(neff,ncore,nclad,radius,rcore,wavenumber)
    b=(neff^2-nclad^2)*(ncore^2-nclad^2);
    V=wavenumber*rcore*sqrt(ncore^2-nclad^2);
    E01pt1=sqrt((377*b)/(pi*nclad*sqrt(1+2*b)));
    E01pt2=1/(rcore*besselj(1,V*sqrt(1-b)));
    E01=E01pt1*E01pt2;
    Er=1i*E01*besselj(0,V*sqrt(1-b)*(radius/rcore));
end
