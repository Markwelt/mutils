function Stats = WSubSE(Data, varargin)
%FUNCTION STATS = WSubSE(DATA[,Cfg])
%   This function takes in an nSubjects x nConditions matrix and computes 
%   the within-subjects standard error, as demonstrated in Morey (2008) 
%   - correction of Cousineau (2005) and alike Loftus and Masson (1994)-
%   in order to make confidence intervals that ignore between-subjects
%   variability, which is generally present in repeated measures studies.
%
%   The equation presented by Cousineau is as follows:
%   Y = Xij - X'i + X'
%   Where:
%   Xij is the observation for the ith subject and the jth condition
%   X'i is the ith subject's mean (for all his conditions) 
%   X' is the grand mean
%
%   The correction presented by Morey is as follows:
%   s' = s * (M/(M-1))
%   Where:
%   s is the sample variances in each condition
%   M is the number of conditions of within-S variables
%
%	The configuration can optionally contain:
%	Cfg.Morey  = Boolean. Whether to apply Morey's correction (default: true).
%  	Cfg.CI  = Double. Percentage of confidence intervals (default: 95).
%
%--------------------------------------------------------------------------
%
%   This function outputs a structure containing descriptive statistics of
%   the data matrix that was input. It has the following 6 fields:
%
%   conditionsMeans: a 1 x nConditions vector of means for each condition
%   subjectMeans: a nSubjects x 1 vector of means for each subject
%   normData: a normalized version of the data matrix that was input
%   wsSE: a 1 x nConditions vector of within-subject SE for each condition
%   tScore: degrees of freedom, a scalar (calculated as nSubjects - 1)
%   wsCI: a 1 x nConditions vector of 95% (default) confidence intervals
%
%--------------------------------------------------------------------------
%
%   Last Edit: WSubSE by Marco Fusca'  16 Sept 2014
%   Based on: computeWithinSubSE.m by Seth Levine
%   Based on: summarySEwithin.R  by Winston Chang
%
%   EXAMPLE CALL
%       Stats = WSubSE(DataMatrix)

%% Config
if isempty(varargin)
    Cfg = [];
else
    Cfg = varargin{1};
end
if  ~isfield(Cfg, 'Morey'), Cfg.Morey = true; else if ~islogical(Cfg.Morey) || ~isdouble(Cfg.Morey), Cfg.Morey = true; end, end
if  ~isfield(Cfg, 'CI'), Cfg.CI = 95; else if ~isdouble(Cfg.CI), Cfg.CI = 95; end, end, if Cfg.CI>100 || Cfg.CI<0, Cfg.CI = 95;  end,


%% Initialize a statistics structure

%All stats regarding the data matrix go into this structure to be output
%together
Stats = [];

%% Calculate the means

%Obtain the number of rows and columns that comprise the data matrix
nSubjects = size(Data,1);
nConditions = size(Data,2);

%Means of the conditions
Stats.conditionMeans = nanmean(Data,1);

%Means of each subject
Stats.subjectMeans = nanmean(Data,2);

%The grand mean is the mean of all the condition means
grandMean = mean(Stats.conditionMeans);
%Note that the following assignment would be equally sufficient
%   grandMean = mean(subjectMeans);

%% De-meaning the data

%Preallocate a matrix for the normalized data
Stats.normData = zeros(nSubjects, nConditions);

%Subtract each subject's mean from each observation in the conditions and
%add the grand mean
%Y = Xij - X'i + X'
for i = 1:nSubjects
    Stats.normData(i,:) = Data(i,:) - Stats.subjectMeans(i) + grandMean;
end

%% Prepare stats for output

%Compute the within-subjects standard error for each condition
Stats.wsSE = nanstd(Stats.normData,1)/sqrt(nSubjects);

%Apply correction from Morey (2008) to the ws standard error
if Cfg.Morey
correctionFactor = sqrt(nConditions/(nConditions-1));
Stats.wsSE = Stats.wsSE * correctionFactor;
end

%Calculate t-statistic for confidence interval
% e.g., if conf.interval is 95, use .975 (above/below), and use df=N-1
confinterv = ((Cfg.CI/100)/2) + .5;
Stats.tScore = tinv(confinterv,nSubjects-1);

%Compute the % confidence intervals using the within-subject SE
Stats.wsCI = Stats.wsSE*Stats.tScore;
