function [p_s_ql,beta,xi_sqr,k_x] = quasilinearCoeff(k,k_y,k_x,waveSpectrum,Tv_k, beta)

dk=[0,diff(k)']; % resize to correct size

%xi_sqr = beta.^2.*cumtrapz(dk, waveSpectrum.*abs(Tv_k).^2);
%xi_sqr = func.resize(xi_sqr(:,2:end),waveSpectrum(1,:));
xi_sqr = beta.^2.*trapz(dk,waveSpectrum.*abs(Tv_k).^2);
xi = sqrt(xi_sqr);
test_k_neg = -(k_x).^2;
k_x_cutoff = xi.^(-1);
negative = 0;

% Cutoff k_x
if mean(k_x) < 0
    k_x_cutoff = -k_x_cutoff;
    negative = 1;
end    
% Initialize the index variable
index = 0;
% Iterate through the array to find the first index greater than the threshold
if negative
    k_x_cutoff_index = k_x > k_x_cutoff;
    k_x = k_x.*k_x_cutoff_index;
else    
    k_x_cutoff_index = k_x < k_x_cutoff;
    k_x = k_x.*k_x_cutoff_index;
end

p_s_ql = exp((k_x).^2.*xi_sqr);
end