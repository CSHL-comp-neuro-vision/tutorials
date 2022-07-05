function PCorrect = OptimalCombDiscrim(FlashProjs, NoFlashProjs, SinglesMean, SinglesSD, NoiseMean, NoiseSD, FlashStrength, NumPhotons, Verbose)
%
%  function PCorrect = OptimalCombDiscrim(FlashProjs, NoFlashProjs, SinglesMean, SinglesSD, NoiseMean, NoiseSD, FlashStrength, NumPhotons, Verbose)
%
%   OptimalCombDiscrim uses a maximum likelihood approach to discriminat
%   pooled rod responses to flash/no flash trials.
%
%   Created:  GDF:  06/11/04
%
%

[NumRods, NumTrials] = size(FlashProjs);

% Check that you can assume rods are only getting < 1 photoisomerization;
MultiphotonProb = 1 - poisscdf(NumPhotons, FlashStrength);
if MultiphotonProb > 0.02
    fprintf('\n WARNING: you are not summing to high enough flash strengths; increase NumPhotons input! \n')
end

for trial = 1:NumTrials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  flash trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    FlashedRF = FlashProjs(:,trial);
    if Verbose
        if trial == 1
            plot(FlashedRF)
            pause(2)
        end
    end
    
    CombProb = 0;
    for photon = 0:NumPhotons
        CombProb = CombProb + (normpdf(FlashedRF, NoiseMean + (photon * SinglesMean), sqrt(NoiseSD^2 +  (photon * SinglesSD^2))) .* poisspdf(photon, FlashStrength));
    end
    FlashDist_flash = CombProb;
    
    % Calculate the probability that its from the no flash distribution
    NoFlashDist_flash = normpdf(FlashedRF, NoiseMean, sqrt(NoiseSD^2));
        
    if Verbose
        if trial == 1
            plot(log10(FlashDist_flash./NoFlashDist_flash))
            pause
            axis([0, length(FlashedRF), -0.001 0.001])
            pause
        end
    end
    
    TempLikelihood = sum(log10(FlashDist_flash./NoFlashDist_flash));
    FlashedLikelihood(trial) = TempLikelihood;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  no flash trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	NoFlashRF = NoFlashProjs(:,trial);
    if Verbose
        if trial == 1
            plot(NoFlashRF)
            pause(2)
        end
    end

    % Calculate the probabilitiy that the no flash responses are from the
    % flash distribution
    CombProb = 0;
    for photon = 0:NumPhotons
        CombProb = CombProb + (normpdf(NoFlashRF, NoiseMean + (photon * SinglesMean), sqrt(NoiseSD^2 +  (photon * SinglesSD^2))) .* poisspdf(photon, FlashStrength));
    end
    FlashDist_noflash = CombProb;

    % Calculte the probability that the no flash responses are from the no
    % flash distribution
    NoFlashDist_noflash = normpdf(NoFlashRF, NoiseMean, sqrt(NoiseSD^2));
      
    if Verbose
        if trial == 1
            plot(log10(FlashDist_noflash./NoFlashDist_noflash))
            pause
            axis([0, length(FlashedRF), -0.001, 0.001])
            pause
        end
    end
    
   % Calculate the log likelihood ratio for the `no-flash' receptive fields
   % for each trial
    TempLikelihood = sum(log10(FlashDist_noflash./NoFlashDist_noflash));
    NoFlashedLikelihood(trial) = TempLikelihood;
end

NumCorrectA = length(find(FlashedLikelihood > 0));
NumCorrectB = length(find(NoFlashedLikelihood < 0));

AChanceCorrect = floor(length(find(FlashedLikelihood == 0)) / 2);
BChanceCorrect = ceil(length(find(NoFlashedLikelihood == 0)) / 2);

NumCorrectA = NumCorrectA + AChanceCorrect;
NumCorrectB = NumCorrectB + BChanceCorrect;

PCorrectA = NumCorrectA ./ NumTrials;
PCorrectB = NumCorrectB ./ NumTrials;

PCorrect = (PCorrectA + PCorrectB) / 2;
    