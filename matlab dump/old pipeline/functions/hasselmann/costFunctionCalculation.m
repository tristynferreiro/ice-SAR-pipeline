function [J,deltaFn,Fn,deltaPn] = costFunctionCalculation(generatedSARSpectrum,observedSARSpectrum,firstGuessWaveSpectrum,inverseFirstGuessWaveSpectrum,B,mu,Ts_k,Ts_k_inv,imageVarianceSpectrum,quasilinearCoeff,invQuasilinearCoeff,numOfIterations,dk)

deltaFn = 0;
deltaPn = 0;
Fn = firstGuessWaveSpectrum;
Pn = observedSARSpectrum;
J = zeros(numOfIterations,1);

for i=1:numOfIterations
    Fn = Fn + deltaFn;
    Pn = Pn + deltaPn;
    deltaPk = observedSARSpectrum - Pn;
    deltaFk = firstGuessWaveSpectrum - Fn;
    deltaFnk = inverseFirstGuessWaveSpectrum - Fn;
    Wk = abs(Ts_k).^2.*quasilinearCoeff;
    Wnk = abs(Ts_k_inv).^2.*invQuasilinearCoeff;
    Ak = Wk.^2 + 2.*mu;
    Ank = Wnk.^2 + 2.*mu;
    Bk = Wk.*Wnk;
    deltaFn = (Ank.*(Wk.*deltaPk + mu.*deltaFk) - Bk.*(Wnk.*deltaPk + mu.*deltaFnk))./(Ak.*Ank - Bk.^2);
    deltaPn = 0.5.*quasilinearCoeff.*imageVarianceSpectrum;
    J = cumtrapz(dk,(deltaPn-deltaPk).^2) + mu.*cumtrapz(dk,((deltaFn-deltaFk)./(B+deltaFk)).^2);
end

end