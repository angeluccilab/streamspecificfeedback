% This statistical analysis compares the distances at which the highest 30%
% of intensity is found relative to the blob-to-interblob distance
% This tests the significance between stripe types

clear all;
close all;

% runDate = '2019-07-22';
runDate = '2020-01-06';
labels = {'GFP' 'TDT' 'CO'};
layers = {'L1A' 'L1B' 'L2_3' 'L4B' 'L5_6'};
stripes = {'Thick' 'Thin' 'Pale'};
fileDate = '2019-07-20';

path( fullfile(pwd,'scriptFuncs'), path );
path( fullfile(pwd,'CombinedCases'), path );
% pathName = 'D:\Semi\Script_FeedbackPatternOnCO\FB Anatomy Scripts Project\CombinedCases\';
pathName = strcat(fullfile(pwd,'CombinedCases'),'\');

%extract the pixel conversions data
convs = xlsread('FBtoCO um to pixel conversions.xlsx');
mk356scale = convs(1,1);
mk365scale = convs(1,2);
mk373scale = convs(1,3);
mk374scale = convs(1,4);

% Load heatmaps, group data per stripe type
for lyr = 1:5
    load(strcat(pathName,fileDate,'_',layers{lyr},'_Heatmaps.mat'));      
%     load(strcat(pathName,fileDate,'_',layers{lyr},'_Heatmaps_Pale.mat'));
    for stripe = 1:size(stripes,2)
        eval(['matr' stripes{stripe} ' = ' layers{lyr} '_heatmaps(stripe).bestBars;']);
        eval(['cases' stripes{stripe} '{lyr} = ' layers{lyr} '_heatmaps(stripe).cases;']);
    end
%     eval(['matrThin = ' layers{lyr} '_heatmaps(2).bestBars;']);
%     eval(['matrPale = ' layers{lyr} '_heatmaps(3).bestBars;']);
    
%     eval(['casesThin{lyr} = ' layers{lyr} '_heatmaps(2).cases;']);
%     eval(['casesPale{lyr} = ' layers{lyr} '_heatmaps_pale(1).cases;']);


    thickValid = []; thinValid = []; paleValid = [];
    for mky = 1:3
        tk = 1; tn = 1; p = 1;
        %normalize
        matrThick(mky,:) = (matrThick(mky,:) - min(matrThick(mky,:)))/(max(matrThick(mky,:)) - min(matrThick(mky,:)));
        matrThin(mky,:) = (matrThin(mky,:) - min(matrThin(mky,:)))/(max(matrThin(mky,:)) - min(matrThin(mky,:)));
        try
            matrPale(mky,:) = (matrPale(mky,:) - min(matrPale(mky,:)))/(max(matrPale(mky,:))-min(matrPale(mky,:)));
        catch
        end
        %look for bins > 70% intensity
        for bin = 1:20
            if matrThick(mky,bin) > 0.7
                thickValid(mky,tk) = bin;
                tk = tk + 1;
            end
            if matrThin(mky,bin) > 0.7
                thinValid(mky,tn) = bin;
                tn = tn + 1;
            end
            try
                if matrPale(mky,bin) > 0.7
                    paleValid(mky,p) = bin;
                    p = p + 1;
                end
            catch
            end
        end
    end
    validBinsThick{lyr,:} = thickValid;
    validBinsThin{lyr,:} = thinValid;   
    validBinsPale{lyr,:} = paleValid;
    
end
clear tk tn p lyr mky bin thick thin pale thickValid thinValid paleValid            
 
% return

%% convert bins to ums
for stripe = 1:3
%     keyboard
    clear cases distances distancesUm;
    if stripe == 1
        cases = casesThick;
        distances = validBinsThick;
    elseif stripe == 2
        cases = casesThin;
        distances = validBinsThin;
    else
        cases = casesPale;
        distances = validBinsPale;
    end
    
%     keyboard
    for lyr = 1:size(cases,2)
        for cs = 1:size(cases{1,lyr},1)
            if cases{1,lyr}{cs,1} == 356
                distancesUm{lyr,1}(cs,:) = (distances{lyr,1}(cs,:)*(cases{1,lyr}{cs,3}/20))/mk356scale; %bins to pixels to microns
            elseif cases{1,lyr}{cs,1} == 365
                distancesUm{lyr,1}(cs,:) = (distances{lyr,1}(cs,:)*(cases{1,lyr}{cs,3}/20))/mk365scale;
            elseif cases{1,lyr}{cs,1} == 373
                distancesUm{lyr,1}(cs,:) = (distances{lyr,1}(cs,:)*(cases{1,lyr}{cs,3}/20))/mk373scale;
            elseif cases{1,lyr}{cs,1} == 374
                distancesUm{lyr,1}(cs,:) = (distances{lyr,1}(cs,:)*(cases{1,lyr}{cs,3}/20))/mk374scale;
            end
        end
    end
    
    if stripe == 1
        validBinsThickUm = distancesUm;
    elseif stripe == 2
        validBinsThinUm = distancesUm;
    else
        validBinsPaleUm = distancesUm;
    end
end


%% Perform statistical analysis comparing stripes (layer vs layer)
%group names for kruskalwallis 
groups = {'Thick' 'Thin' 'Pale'};
for lyr = 1:5
    %extract non-zero elements
    validThick = nonzeros(validBinsThick{lyr});
    validThin = nonzeros(validBinsThin{lyr});
    validPale = nonzeros(validBinsPale{lyr});
    validThickUm = nonzeros(validBinsThickUm{lyr});
    validThinUm = nonzeros(validBinsThinUm{lyr});
    validPaleUm = nonzeros(validBinsPaleUm{lyr});
    
    % create matrix with each column having different stripe data, will be
    % used for non-parametric analysis (kruskalwallis)
    matSize = max([size(validThick,1) size(validThin,1) size(validPale,1)]);
    validAllLayers = zeros(matSize,3);
    validAllLayers(:,1) = [validThick' zeros(1,matSize-size(validThick,1))]';
    validAllLayers(:,2) = [validThin' zeros(1,matSize-size(validThin,1))]';
    validAllLayers(:,3) = [validPale' zeros(1,matSize-size(validPale,1))]';
    validAllLayers(validAllLayers==0) = NaN;
    layersDataSet{lyr} = validAllLayers;
    
    matSize = max([size(validThickUm,1) size(validThinUm,1) size(validPaleUm,1)]);
    validAllLayersUm = zeros(matSize,3);
    validAllLayersUm(:,1) = [validThickUm' zeros(1,matSize-size(validThickUm,1))]';
    validAllLayersUm(:,2) = [validThinUm' zeros(1,matSize-size(validThinUm,1))]';
    validAllLayersUm(:,3) = [validPaleUm' zeros(1,matSize-size(validPaleUm,1))]';
    validAllLayersUm(validAllLayersUm==0) = NaN; %replace zeros with NaN, making them non-relevant to the statistical analysis
    layersDataSetUm{lyr} = validAllLayersUm; 
    
    % alternate data collection (for ANOVAN), put all data into one vector
    validAllLayersVector = [validThick' validThin' validPale']';
    validAllLayersVectorUm = [validThickUm' validThinUm' validPaleUm']';
    
    % create vector identifying the stripe type that each element of the
    % previous vector belong to
    groupLayers = {[repmat("Thick",1,size(validThick,1)) repmat("Thin",1, ...
        size(validThin,1)) repmat("Pale",1,size(validPale,1))]};

       
    % perform ANOVA n-way analysis
    [AnovaLayersPvaluesUm(lyr) AnovaLayersTablesUm{lyr} AnovaLayersStatsUm{lyr}] = anovan(validAllLayersVectorUm, ...
        groupLayers,'varnames',{char(layers{lyr})});
    figure;
    [AnovaLayersComparisonResultsUm{lyr}, AnovaLayersMeansAndErrorsUm{lyr}, fig1] = multcompare(AnovaLayersStatsUm{lyr});
    
    figName = strcat("ANOVA_micron_",layers{lyr});
    set(fig1,'name',figName);
    
    [AnovaLayersPvalues(lyr) AnovaLayersTables{lyr} AnovaLayersStats{lyr}] = anovan(validAllLayersVector, ...
        groupLayers,'varnames',{char(layers{lyr})});
    figure;
    [AnovaLayersComparisonResults{lyr}, AnovaLayersMeansAndErrors{lyr},fig2] = multcompare(AnovaLayersStats{lyr});
    
    figName = strcat("ANOVA_pixel_",layers{lyr});
    set(fig2,'name',figName);
    
    % perform kruskalwallis, non parametric analysis
    [KwLayersPvaluesUm(lyr) KwLayersTablesUm{lyr} KwLayersStatsUm{lyr}] = kruskalwallis(validAllLayersUm,groups);
    [KwLayersComparisonResultsUm{lyr}, KwLayersMeansAndErrorsUm{lyr},fig3] = multcompare(KwLayersStatsUm{lyr});
    
    figName = strcat("KW_micron_",layers{lyr});
    set(fig3,'name',figName);
    
    [KwLayersPvalues(lyr) KwLayersTables{lyr} KwLayersStats{lyr}] = kruskalwallis(validAllLayers,groups);
    [KwLayersComparisonResults{lyr}, KwLayersMeansAndErrors{lyr},fig4] = multcompare(KwLayersStats{lyr});
    
    figName = strcat("KW_pixel_",layers{lyr});
    set(fig4,'name',figName);
    
    clear groupLayers
end

% return  
% save(fullfile(pwd,strcat('/CombinedCases/StatisticalAnalysisTangantialCases.mat')));

%% Perform statistical analysis comparing stripes (all layers, cumulative)

%collect all data from all layers and group into stripe type
for stripe = 1:3
    clear valids
    if stripe == 1
        validBins = validBinsThick;
        validBinsUm = validBinsThickUm;
    elseif stripe ==2
        validBins = validBinsThin;
        validBinsUm = validBinsThinUm;
    else
        validBins = validBinsPale;
        validBinsUm = validBinsPaleUm;
    end
    
    valids = validBins{1,1}(find(validBins{1,1}));
    validsUm = validBinsUm{1,1}(find(validBinsUm{1,1}));
    validBinsStripe{stripe} = valids;
    validBinsStripeUm{stripe} = validsUm;
    for lyr = 2:5
        valids = validBins{lyr,1}(find(validBins{lyr,1}));
        validBinsStripe{stripe} = vertcat(validBinsStripe{stripe},valids);
        validsUm = validBinsUm{lyr,1}(find(validBinsUm{lyr,1}));
        validBinsStripeUm{stripe} = vertcat(validBinsStripeUm{stripe},validsUm);
    end
    stripeCount(stripe) = size(validBinsStripe{stripe},1);
end

%put collected into a matrix (non-cell structure) for kruskalwallis
%analysis
validAllStripes = zeros(max(stripeCount),3);
validAllStripesUm = zeros(max(stripeCount),3);
for stripe = 1:3
    validAllStripes(1:size(validBinsStripe{stripe},1),stripe) = validBinsStripe{stripe};
    validAllStripesUm(1:size(validBinsStripeUm{stripe},1),stripe) = validBinsStripeUm{stripe};
end
validAllStripes(validAllStripes==0) = NaN; %replace all zeros with NaN so that they are not included in the statistical analysis
validAllStripesUm(validAllStripesUm==0) = NaN;

%Put data into a single vector for ANOVAN analysis
validAllStripesVector = [validBinsStripe{1}' validBinsStripe{2}' validBinsStripe{3}']';
validAllStripesVectorUm = [validBinsStripeUm{1}' validBinsStripeUm{2}' validBinsStripeUm{3}']';

% create vector identifying the stripe type that each element of the
% previous vector belong to
for b = 1:size(validAllStripesVector,1)
    if b < stripeCount(1) + 1
        groupStripes{b} = 'T';
    elseif b < (stripeCount(1) + stripeCount(2) + 1)
        groupStripes{b} = 't';
    else
        groupStripes{b} = 'P';
    end
end
groupStripes = char(groupStripes);
groupStripes = {groupStripes};

% Perform ANOVAN
[AnovaStripesPvaluesUm AnovaStripesTablesUm AnovaStripesStatsUm] = anovan(validAllStripesVectorUm, ...
        groupStripes,'varnames',{'StripeType'});
figure;
[AnovaStripesComparisonResultsUm, AnovaStripesMeansAndErrorsUm, fig1] = multcompare(AnovaStripesStatsUm);

figName = strcat("ANOVA_micron_ALL Layers");
    set(fig1,'name',figName);

[AnovaStripesPvalues AnovaStripesTables AnovaStripesStats] = anovan(validAllStripesVector, ...
        groupStripes,'varnames',{'StripeType'});
figure;
[AnovaStripesComparisonResults, AnovaStripesMeansAndErrors, fig2] = multcompare(AnovaStripesStats);

figName = strcat("ANOVA_pixel_ALL Layers");
    set(fig2,'name',figName);

% perform kruskalwallis, non parametric analysis
[KwStripesPvaluesUm KwStripesTablesUm KwStripesStatsUm] = kruskalwallis(validAllStripesUm,groups);
[KwStripesComparisonResultsUm, KwStripesMeansAndErrorsUm, fig3] = multcompare(KwStripesStatsUm);

figName = strcat("KW_pixel_ALL Layers");
    set(fig3,'name',figName);

[KwStripesPvalues KwStripesTables KwStripesStats] = kruskalwallis(validAllStripes,groups);
[KwStripesComparisonResults, KwStripesMeansAndErrors, fig4] = multcompare(KwStripesStats);

figName = strcat("KW_pixel_ALL Layers");
    set(fig4,'name',figName);

clear L1A_heatmaps L1A_heatmaps_pale L1B_heatmaps L1B_heatmaps_pale ...
    L2_3_heatmaps L2_3_heatmaps_pale L4B_heatmaps L4B_heatmaps_pale ...
    L5_6_heatmaps L5_6_heatmaps_pale b lyr matrPale matrThick matrThin ...
    stripes validAllVector  ...
    validPale validThick validThin pathName matr labels fileDate ...
    validPaleUm validThickUm validThinUm  ...
    cases distances distancesUm mk356scale mk365scale ...
    mk373scale mk374scale stripe cs convs casesThick casesThin casesPale ...
%     matSize validAll validAllUm validAllVectorUm stripeCount validAllLayers ...
%     validAllLayersUm validAllLayersVector validAllLayersVectorUm validBins ...
%     validBinsPale validBinsPaleUm validBinsStripe validBinsStripeUm valids ...
%     validBinsThick validBinsThickUm validBinsThin validBinsThinUm validsUm ...
%     validBinsUm groups groupStripes validAllStripes validAllStripesUm ...
%     validAllStripesVector validAllStripesVectorUm fig1 fgi2 fig3 fig4

save(fullfile(pwd,strcat('/CombinedCases/StatsWorkspace_',runDate,'.mat')));

%% Gather comparison data
statsLayers = {'AnovaLayersComparisonResultsUm' 'AnovaLayersComparisonResults' ...
    'KwLayersComparisonResultsUm' 'KwLayersComparisonResults' ...
    'AnovaStripesComparisonResultsUm' 'AnovaStripesComparisonResults' ...
    'KwStripesComparisonResultsUm' 'KwStripesComparisonResults'};
for stat = 1:size(statsLayers,2)
    comparisonStats = [];
    if stat < 5
        for lyr = 1:size(layers,2)
            temp = [];
            eval(['temp =' statsLayers{stat} '{1,lyr};']);
            comparisonStats = [comparisonStats; temp];
            layerMat = [layers{lyr} layers{lyr} layers{lyr}];
        end
        heatmapStatistics(stat).StatisticType = statsLayers{stat};
        heatmapStatistics(stat).Data = comparisonStats;        
    else
        eval(['temp =' statsLayers{stat} ';']);
        heatmapStatistics(stat).StatisticType = statsLayers{stat};
        heatmapStatistics(stat).Data = temp;
    end
    
end
clear layers temp comparisonStats stat lyr statsLayers
save(fullfile(pwd,strcat('/CombinedCases/StatsResultsBlobsAnalysis_',runDate,'.mat')),'heatmapStatistics');
  