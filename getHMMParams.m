function [Output] = getHMM_Params_INatN170_Unique(varargin)
% Example (Run either locally or via Screen on Rondo):
% Written 4/28/17
%__________________________________________________________________________

%% -- Get Conversion Factor for Spike Times --
binsize = 200; %200 in units of 10 kHz = 20 ms width for binning
%(binsize is really the # of samples/time bin assuming 10 kHz sampling
%rate. The original spike times are in units of 10 kHz)

%% -- Initilize Specifications --
A.Nc          = 1;            %Defines Nc-fold cross-validation (for HMM only)
A.niter       = 100;          %# of iterations of the EM algorithm
A.Reg         = NaN;          %NaN = use default (0.002) 
A             = parseargs(A, varargin{:});

if ~isnan(A.Reg)
    addpath([pwd '/Code_JSP/Tree_HMM_cv/Regularization_' ...
             strrep(num2str(A.Reg),'.','pt') '/']);  
else %NaN (Use Default Regularization Parameter Value)
    addpath([pwd '/Code_JSP/Tree_HMM_cv/']);  
end 

%% -- Load the Optimal # Basins --
loaddir  = [pwd '/ParamFits_TreeHMM/NaturalMovie_N170_Repeated/' ...
            'Unique/']; 
loadname = 'BestnBasins_73Unique.mat';  
load([loaddir loadname]); %#ok - loads Elbow & Peak

%% -- Main Computations -- 
nbasins = Peak; 
% Load spike times:
st_loaddir  = [pwd '/DataSets/NaturalMovie_N170_Repeated_JSP/'];
st_loadname = 'Resps_to_73UniqueMovieSegments.mat';  
st = load([st_loaddir st_loadname],'SUn_times'); st = st.SUn_times; %units: 10 kHz

goodcells = 1:length(st); 

if A.Nc>1 %HMM syntax w/ Cross-Val
    % -- Compute Time Bin Indices for Nc-Fold Cross-Validation: -- 
    tmax = max(cell2mat(cellfun(@(x) max(double(x)), st, ...
               'UniformOutput', 0))); %tmax is in units of 10 kHz

    bins         = 0:binsize:tmax;
    s            = RandStream('mt19937ar','Seed',0);
    shuffle_bins = randperm(s,length(bins));
    ntest        = floor(length(bins)/A.Nc);

    %- Computations: -- 
    for k = 1:A.Nc

    testbins   = shuffle_bins((k-1)*ntest+1:k*ntest);
    train_bins = zeros(1,length(bins));
    train_bins(testbins) = 1;

    unobserved_low = bins(diff([0,train_bins]) == 1);
    unobserved_hi  = bins(diff([0,train_bins]) == -1);
    if (length(unobserved_hi) < length(unobserved_low))
        unobserved_hi = [unobserved_hi, tmax]; %#ok
    end

    [logli(:,k), trans, P_emiss, alpha, pred_prob, hist, params, sample] = ...
        EMBasins(st(goodcells), [unobserved_low', unobserved_hi'], ...
        binsize, nbasins, A.niter); %#ok
    end %for
    
    % -- Save: --
    Output.logli      = logli; %*cross-validated* logli (validated on test set)
    Output.trans      = trans;
    Output.P_emiss    = P_emiss;
    Output.alpha      = alpha;
    Output.pred_prob  = pred_prob;
    Output.hist       = hist;
    Output.params     = params;
    Output.sample     = sample; 

    savename = ['HMM_Params_RepeatNatN170_20ms' ...
                '_nb' num2str(nbasins) '_' num2str(A.Nc) 'cv.mat'];
    save([loaddir savename],'Output');
    
elseif A.Nc==1 %Use ALL data to determine params
    
    [logli, trans, P_emiss, alpha, pred_prob, hist, params, sample] = ...
        EMBasins(st(goodcells), [], binsize, nbasins, A.niter); 
    
    % -- Save: --
    Output.logli      = logli; %*cross-validated* logli (validated on test set)
    Output.trans      = trans;
    Output.P_emiss    = P_emiss;
    Output.alpha      = alpha;
    Output.pred_prob  = pred_prob;
    Output.hist       = hist;
    Output.params     = params;
    Output.sample     = sample; 

    savename = ['HMM_Params_INatN170_Unique_20ms' ...
                '_nb' num2str(nbasins) '.mat'];
    save([loaddir savename],'Output');
end

end %main fn