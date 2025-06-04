%% Seismic Velocity Change Analysis Using Stretching Method (Coda Wave Interferometry) in Rt component
% This script calculates relative seismic velocity changes (dv/v) using the
% stretching method applied to coda waves. The method compares daily wave
% with a reference waveform by applying systematic time stretching/compression
% to find the optimal correlation, which corresponds to velocity changes.
%
% - Velocity decreases cause waveforms to arrive later (positive stretch)
% - Velocity increases cause waveforms to arrive earlier (negative stretch)
% - Method provides high precision (0.01% level) velocity change detection
%
% Main objectives:
% 1. Extract coda wave portions from daily transfer function and reference 
% 2. Apply systematic time stretching to daily coda
% 3. Find optimal stretch ratio that maximizes correlation with reference
% 4. Convert stretch ratios to relative velocity changes (dv/v)
% 5. Calculate correlation coefficients for quality assessment
%
% Input files required:
% - Rt_data.mat: Daily seismogram data for Rt component
% - Coda_Rt.mat: Coda wave end times for each station
%
% Output: dv_v_Rt.mat containing velocity changes and correlations
%
% Method reference: Sens-Schönfelder & Wegler (2006), Geophys. J. Int.
%
% Author: FENG CHEN

clc
clear
% Load coda wave timing information
load('Rt_data.mat');
% Expected variable: Coda_end - array of coda end times for each station
load('Coda_Rt.mat');
fsamp=100;
coda_start=5;
coda_end=Coda_end; 
seisdix=size(daily_data,3);
%% Create daily and reference coda waves 
for seis_ch=1:seisdix
    cut_daily_data{seis_ch}=squeeze(daily_data(coda_start*fsamp+1:coda_end(seis_ch)*fsamp,:,seis_ch)); % 5000*322*14
    cut_ref_data{seis_ch}=ref_data(coda_start*fsamp+1:coda_end(seis_ch)*fsamp,seis_ch); %5000*14
end
%% Calculated dv/v 
dayidx=size(daily_data,2);
seisdix=size(daily_data,3);
dv_v_Rt=nan(dayidx,seisdix);  % 322*14
CC_Rt=nan(dayidx,seisdix);  % 322*14
for seis_ch=1:seisdix %14
    st=tic;
    referenceTF=cut_ref_data{seis_ch};
    for day_num=1:length(t)
        currentTF = cut_daily_data{seis_ch}(:,day_num);
        if ~any(isnan(currentTF))
            maxCorrelation = -Inf;
            bestStretchRatio = 0;
            for stretchRatio = -0.05:0.0001:0.05 % Search range: -5% to +5% in 0.01% increments
                correlation = stretchTF(referenceTF, currentTF,stretchRatio,0:length(coda_start*fsamp+1:coda_end(seis_ch)*fsamp)-1);
                if correlation(1,2) > maxCorrelation
                    maxCorrelation = correlation(1,2);
                    bestStretchRatio = stretchRatio;
                end
            end
            % Convert stretch ratio to percentage for dv/v
            dv_v_Rt(day_num,seis_ch) = bestStretchRatio*100;
            CC_Rt(day_num,seis_ch) = maxCorrelation;
        end
    end
    fprintf('%05.0f of %05.0f Done. %05.2f [sec]\n',seis_ch,seisdix,toc(st))
end

save dv_v_Rt cut_daily_data cut_ref_data t dv_v_Rt CC_Rt



function correlation= stretchTF(ref_TF,org_TF,stretchRatio, tlag)
% STRETCHTF Applies time stretching to seismic waveform and calculates correlation
t_stretched = tlag * (1 + stretchRatio);
new_tlag_stretched=linspace(0, max(t_stretched), length(tlag));
stretched_TF = interp1(t_stretched, org_TF, new_tlag_stretched, 'linear');
if stretchRatio>=0
    stretched_TF=interp1(new_tlag_stretched,stretched_TF, tlag, 'linear');
else
    ref_TF=interp1(tlag, ref_TF, new_tlag_stretched, 'linear');
end
% if stretchRatio>=0
%     plot(tlag/100,ref_TF,'k');
%     hold on
%     plot(tlag/100,stretched_TF,'r');
%     hold off
%     xlim([0 3]);
%     text(1,0.5,sprintf('Ratio=%4.3f',stretchRatio*100),'FontSize',22);
% else
%     plot(new_tlag_stretched/100,ref_TF,'k');
%     hold on
%     plot(new_tlag_stretched/100,stretched_TF,'r');
%     hold off
%     xlim([0 3]);
%     text(1,0.5,sprintf('Ratio=%4.3f',stretchRatio*100),'FontSize',22);
% end
correlation = corrcoef(stretched_TF, ref_TF);
% text(1,10^-6,sprintf('COF=%4.3f',correlation(1,2)),'FontSize',22);
end