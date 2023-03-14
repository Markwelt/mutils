function [mnFt,fE_Inv,fE_Invs] = mafu_fEIratio(Signal, windows, DFA_Overlap)
% function [DFA_ys,fE_I_ind_w,fE_I_ind_ws] = mafu_fEIratio(Signal, windows [,DFA_Overlap])
% simplified implementation of critical functional Excitation/Inhibition ratio
%
% Signal : samples x channels, usually a band-passed amplitude
% windows : vector of scales (i.e. window lengths). Unit is number of samples
% DFA_Overlap (optional) : overlap between windows (to obtain higher SNR for the fluctuation estimates)
% mnFt (output) : windows x channels, mean(nF(t)) average fluctuations of detrended amplitude-normalized signal profile windows
% fE_Inv (output) : windows x channels, 1-fE/I correlations amplitude - nF(t)
% fE_Invs (output, optional) : windows x channels, 1-fE/I Spearman correlations of amplitude - nF(t)
%
% Originally a Neurophysiological Biomarker Toolbox (http://www.nbtwiki.net) function nbt_Scaling_DFA.m
% created by Klaus Linkenkaer-Hansen (2001), improved code by Simon-Shlomo Poil (2008,2009)
% Modified by Ehtasham Javed (2020) as calculate_fEI_ratio.m for calculating fE/I ratio

%%
if nargin<3
    DFA_Overlap=0;
end

%%
[mnFt,fE_Inv] = deal(zeros([length(windows),size(Signal,2)]));
if nargout > 2,fE_Invs = zeros([length(windows),size(Signal,2)]);end

for ChannelID = 1:(size(Signal,2)) % loop over channels 
        y = Signal(:,ChannelID);
        yi = y-mean(y); % demean
        yi = cumsum(yi); % Integrate the above signal (signal profile)
        for i = 1:size(windows,2) % Loop through each window size
            window = windows(i);
            [D,mean_w_Amp] = deal(zeros([floor(size(yi,1)/(window*(1-DFA_Overlap))),1])); % Initialize vectors
            % The smaller the window, the more windows are there (longer D)
            tt = 1;
            for nn = 1:round(window*(1-DFA_Overlap)):size(yi,1)-window % we are going to normalize all windows in steps of 'n*(1-DFA_Overlap)'. % Same number of steps as length of D
                winpart = nn:round(nn+window); % ML Integer operands required for colon operator warning
                mean_w_Amp(tt) = mean(y(winpart)); % mean of original amplitude (in the window)
                D(tt) = std(fastdetrend(yi(winpart)./mean_w_Amp(tt))); %divide each window by original amplitude, then detrend
                % std = fluctuation of detrended amplitude-normalized signal profile windows
                tt = tt + 1;
            end
            mnFt(i,ChannelID) = mean(D);	% nF(t), mean std fluctuations of the detrended normalized integrated signal.
            fE_Inv(i,ChannelID) = corr(D,mean_w_Amp);
            if nargout > 2,fE_Invs(i,ChannelID) = corr(D,mean_w_Amp,'Type','Spearman');end
        end % -- 
end

end

%% Supporting function
function signal = fastdetrend(signal)
% A simple and fast detrend, see also the supporting function fastdetrend
% in the supporting functions folder
persistent a
persistent N
if (isempty(a) || size(signal,1) ~= N)
    N = size(signal,1);
    a = [zeros(N,1) ones(N,1)];
    a(1:N) = (1:N)'/N;
end
signal = signal - a*(a\signal);

end
