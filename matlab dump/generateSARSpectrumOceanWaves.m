function generatedSARSpectrum = generateSARSpectrumOceanWaves(k,k_y,k_x,E_k,E_k_inv,incidenceAngle,r,nonLinOrder, beta,Tv_k,Tv_k_inv, Tr_k, Tr_k_inv)

% Co and autocovariance functions
[f_v_r] = orbitalVelocityCovariance(k,E_k,r,incidenceAngle,Tv_k);
[f_r_r] = rarImageIntensityAutocovariance(k,E_k,E_k_inv,r,Tr_k, Tr_k_inv);
[f_rv_r] = rarImageIntensityCovariance(k,E_k,E_k_inv,r,Tr_k,Tr_k_inv,Tv_k,Tv_k_inv);
[p_s_coeff,beta,xi_sqr] = quasilinearCoeff(k,k_y,k_x,E_k, Tv_k, beta);

% Spectral Expansions
[p_s_2n] = spectralExpansion2n(f_v_r,nonLinOrder);
f_rv_r_inv = rarImageIntensityCovariance(k,E_k,E_k_inv,-1.*r,Tr_k,Tr_k_inv,Tv_k,Tv_k_inv);
[p_s_2n_1] = spectralExpansion2n_1(f_rv_r,f_rv_r_inv,f_v_r,nonLinOrder);
f_rv_r_zero = rarImageIntensityCovariance(k,E_k,E_k_inv,0.*r,Tr_k,Tr_k_inv,Tv_k,Tv_k_inv);
[p_s_2n_2] = spectralExpansion2n_2(f_r_r,f_v_r,f_rv_r,f_rv_r_zero,f_rv_r_inv,1);

% Full spectral expansion
P_s = zeros(512,512);

for nonLinOrder = min(nonLinOrder):max(nonLinOrder)
    for m = 2*nonLinOrder-2:2*nonLinOrder
        fprintf('m = %d, ',m)
        fprintf('n = %d \n',nonLinOrder)
        coeff = ((k_x.*beta).^m);
        if (m == 2*nonLinOrder)
            fprintf('P_(n,2n) used \n')
            P_s = P_s + coeff.*p_s_2n;
        end
        if (m == 2*nonLinOrder-1)
            fprintf('P_(n,2n-1) used \n')
            P_s = P_s + coeff.*p_s_2n_1;
        end        
        if (m == 2*nonLinOrder-2)
            fprintf('P_(n,2n-2) used\n')
            P_s = P_s + coeff.*p_s_2n_2;
        end
        
    end
end
P_s = p_s_coeff.*P_s;
generatedSARSpectrum = P_s;
end