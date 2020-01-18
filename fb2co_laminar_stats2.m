% This statistical analysis compares the intensities between the pale and
% thin stripes on a layer by layer basis. Also, the intensities of layer
% 4a, 4b, and 4c is compared for both stripes, individually and combined

clear all;
close all;

runDate = '2020-01-06';
paleInjections = {'GFP Lat' 'GFP Med'};
thinInjections = {'TDT Lat' 'TDT Med'};
corticalLayers = {"L1" "L2_3" "L4A" "L4B" "L4C" "L5" "L6"};
stripes = {'Thin' 'Pale'};
fileDate = '2020-01-06';
isThresholded = 0;
isWithCO = 0;
midsOnly = 1;

performANOVA = 0;

path( fullfile(pwd,'scriptFuncs'), path );
path( fullfile(pwd,'MK359LH'), path );
if isThresholded
    if isWithCO
%         pathName = 'D:\Semi\Script_FeedbackPatternOnCO\FB Anatomy Scripts Project\MK359LH\ResultsAnalysisVer5a\Thresholded\CO-stained\';
        pathName = fullfile(pwd,strcat('\MK359LH\Results\Thresholded\CO-stained\'));
        load(strcat(pathName,'MK359LH_CombinedResults_Thresholded_withCO_',fileDate,'.mat'));
        filenameWS = strcat(pathName,'MK359LH_StatsWorkspace_Thresholded_withCO_',runDate,'.mat');
        filenameStats = strcat(pathName,'MK359LH_Statistics_Thresholded_withCO_',runDate,'.mat');
        
    else
%         pathName = 'D:\Semi\Script_FeedbackPatternOnCO\FB Anatomy Scripts Project\MK359LH\ResultsAnalysisVer5a\Thresholded\PreCO-stained\';
        pathName = fullfile(pwd,strcat('\MK359LH\Results\Thresholded\PreCO-stained\'));
        load(strcat(pathName,'MK359LH_CombinedResults_Thresholded_withoutCO_',fileDate,'.mat'));
        filenameWS = strcat(pathName,'MK359LH_StatsWorkspace_Thresholded_withoutCO_',runDate,'.mat');
        filenameStats = strcat(pathName,'MK359LH_Statistics_Thresholded_withoutCO_',runDate,'.mat');
    end
else
    if isWithCO
%         pathName = 'D:\Semi\Script_FeedbackPatternOnCO\FB Anatomy Scripts Project\MK359LH\ResultsAnalysisVer5a\Not Thresholded\CO-stained\';
        pathName = fullfile(pwd,strcat('\MK359LH\Results\Not Thresholded\CO-stained\'));
        load(strcat(pathName,'MK359LH_CombinedResults_NonThresholded_withCO_',fileDate,'.mat'));
        filenameWS = strcat(pathName,'MK359LH_StatsWorkspace_NonThresholded_withCO_',runDate,'.mat');
        filenameStats = strcat(pathName,'MK359LH_Statistics_NonThresholded_withCO_',runDate,'.mat');
    else
%         pathName = 'D:\Semi\Script_FeedbackPatternOnCO\FB Anatomy Scripts Project\MK359LH\ResultsAnalysisVer5a\Not Thresholded\PreCO-stained\';
        pathName = fullfile(pwd,strcat('\MK359LH\Results\Not Thresholded\PreCO-stained\'));
        load(strcat(pathName,'MK359LH_CombinedResults_NonThresholded_withoutCO_',fileDate,'.mat'));
        filenameWS = strcat(pathName,'MK359LH_StatsWorkspace_NonThresholded_withoutCO_',runDate,'.mat');
        filenameStats = strcat(pathName,'MK359LH_Statistics_NonThresholded_withoutCO_',runDate,'.mat');
    end
end
  
injInd = 1;
for inj = 1:4
    if midsOnly
        layerBrdrs(inj,:) = round(layers2(1).PopulationDataMiddleSections(inj).CombinedLayerBorderBins);
        intensities(inj,:) = layers2(1).PopulationDataMiddleSections(inj).NormalizedCombinedIntensities;
        intensitiesInjections(injInd:injInd+2,:) = layers2(1).PopulationDataMiddleSections(inj).AllIntensities;
        injInd = injInd + 3;
    else
        layerBrdrs(inj,:) = round(layers2(1).PopulationData(inj).CombinedLayerBorderBins);
        intensities(inj,:) = layers2(1).PopulationData(inj).NormalizedCombinedIntensities;
        intensitiesInjections(injInd:injInd+6,:) = layers2(1).PopulationData(inj).AllIntensities;
        injInd = injInd + 7;
    end
end

% if midsOnly
%     intensities = cat(1,layers2(1).PopulationDataMiddleSections(1).NormalizedCombinedIntensities, ...
%         layers2(1).PopulationDataMiddleSections(2).NormalizedCombinedIntensities, ...
%         layers2(1).PopulationDataMiddleSections(3).NormalizedCombinedIntensities, ...
%         layers2(1).PopulationDataMiddleSections(4).NormalizedCombinedIntensities);
%     layerBrdrs = cat(1,round(layers2(1).PopulationDataMiddleSections(1).CombinedLayerBorderBins), ...
%         round(layers2(1).PopulationDataMiddleSections(2).CombinedLayerBorderBins), ...
%         round(layers2(1).PopulationDataMiddleSections(3).CombinedLayerBorderBins), ...
%         round(layers2(1).PopulationDataMiddleSections(4).CombinedLayerBorderBins));
% else
%     intensities = cat(1,layers2(1).PopulationData(1).NormalizedCombinedIntensities, ...
%         layers2(1).PopulationData(2).NormalizedCombinedIntensities, ...
%         layers2(1).PopulationData(3).NormalizedCombinedIntensities, ...
%         layers2(1).PopulationData(4).NormalizedCombinedIntensities);
%     layerBrdrs = cat(1,round(layers2(1).PopulationData(1).CombinedLayerBorderBins), ...
%         round(layers2(1).PopulationData(2).CombinedLayerBorderBins), ...
%         round(layers2(1).PopulationData(3).CombinedLayerBorderBins), ...
%         round(layers2(1).PopulationData(4).CombinedLayerBorderBins));
% end
% 
% for b = 1:size(paleBorders,2)
   
    

brdr = 1; thin4 = []; pale4 = []; thin4group = []; pale4group = [];
for lyr = 1:size(corticalLayers,2)
    paleIntens = []; thinIntens = []; tempMat = [];
    paleIntensInjects = []; thinIntensInjects =[];
    cs2 = 1;
    for cs = 1:size(layerBrdrs,1)
        brdrs = layerBrdrs(cs,:);
        currIntens = [];  curr2 = [];
        if cs == 1 || cs == 2          
            currIntens = intensities(cs,brdrs(brdr):brdrs(brdr+1)-1);
            paleIntens = [paleIntens,currIntens]; % gather combined intensities
            if midsOnly
                curr2 = intensitiesInjections(cs2:cs2+2,brdrs(brdr):brdrs(brdr+1)-1);
                paleIntensInjects = [paleIntensInjects,curr2]; % gather individual section intensities
                cs2 = cs2 + 3;
            else
                curr2 = intensitiesInjections(cs2:cs2+6,brdrs(brdr):brdrs(brdr+1)-1);
                paleIntensInjects = [paleIntensInjects,curr2];
                cs2 = cs2 + 7;
            end
            
            
        else
            currIntens = intensities(cs,brdrs(brdr):brdrs(brdr+1)-1);
            thinIntens = [thinIntens,currIntens]; % gather combined intensities
            if midsOnly % only 3 sections as compared to 7 sections for the entire width of the injection site
                curr2 = intensitiesInjections(cs2:cs2+2,brdrs(brdr):brdrs(brdr+1)-1);
                thinIntensInjects = [thinIntensInjects,curr2]; %gather individual section intensities
                cs2 = cs2 + 3;
            else
                curr2 = intensitiesInjections(cs2:cs2+6,brdrs(brdr):brdrs(brdr+1)-1);
                thinIntensInjects = [thinIntensInjects,curr2];
                cs2 = cs2 + 7;
            end
        end
    end
    
    %create matrix for kw analysis from combined data (smallest N)
    paleBrdrs = mean(layerBrdrs(1:2,:));
    thinBrdrs = mean(layerBrdrs(3:4,:));
    paleAvgIntens = mean(intensities(1:2,paleBrdrs(brdr):paleBrdrs(brdr+1)-1));
    thinAvgIntens = mean(intensities(3:4,thinBrdrs(brdr):thinBrdrs(brdr+1)-1));
    matSize3 = max([size(paleAvgIntens,2) size(thinAvgIntens,2)]);
    tempMat = zeros(matSize3,2);
    tempMat(:,1) = [thinAvgIntens ones(1,matSize3-size(thinAvgIntens,2))*-1]';
    tempMat(:,2) = [paleAvgIntens ones(1,matSize3-size(paleAvgIntens,2))*-1]';
    tempMat(tempMat==-1) = NaN;
    kwDataAvg{lyr,1} = tempMat;
            
%     tempPale = [paleIntensities(1,paleBorders(1):paleBorders(2)),paleIntensities(2,paleBorders(1):paleBorders(2))];
%     tempThin = thinIntensities(inj,thinBorders(1):thinBorders(2));
    
    % create matrix for kw analysis from individual case/stripe data(mid N)
    thinCount = size(thinIntens,2);
    paleCount = size(paleIntens,2);
    matSize = max([thinCount paleCount]);
    tempMat = zeros(matSize,2);
    tempMat(:,1) = [thinIntens ones(1,matSize-size(thinIntens,2))*-1]';
    tempMat(:,2) = [paleIntens ones(1,matSize-size(paleIntens,2))*-1]';
    tempMat(tempMat==-1) = NaN;
    kwData{lyr,1} = tempMat;
    
    %create matrix for kw analysis from all sections in each case (max N)
    thinCount2 = size(thinIntensInjects,2);
    thinCount3 = size(thinIntensInjects,1);
    paleCount2 = size(paleIntensInjects,2);
    paleCount3 = size(paleIntensInjects,1);
    matSize2 = max([thinCount2 paleCount2]);
    tempMat = zeros(matSize2,thinCount3+paleCount3);
    if midsOnly
        tempMat(:,1:3) = [thinIntensInjects ones(thinCount3,matSize2-size(thinIntensInjects,2))*-1]';
        tempMat(:,4:6) = [paleIntensInjects ones(paleCount3,matSize2-size(paleIntensInjects,2))*-1]';
    else
        tempMat(:,1:7) = [thinIntensInjects ones(thinCount3,matSize2-size(thinIntensInjects,2))*-1]';
        tempMat(:,8:14) = [paleIntensInjects ones(paleCount3,matSize2-size(paleIntensInjects,2))*-1]';
    end
    tempMat(tempMat==-1) = NaN;
    kwDataInjects{lyr,1} = tempMat;
    
    if performANOVA
        anovaData{lyr,1} = [thinIntens paleIntens]';
        anovaData{lyr,2} = {[repmat("Thin",1,thinCount) repmat("Pale",1,paleCount)]'};
    end
    
    if ismember(lyr, [3 4 5])
        thin4 = [thin4 thinIntens];
        pale4 = [pale4 paleIntens];
        thin4group = [thin4group repmat(corticalLayers{lyr},1,thinCount)];
        pale4group = [pale4group repmat(corticalLayers{lyr},1,paleCount)];
    end
    brdr = brdr + 1;     
end

if performANOVA
    anovaDataL4{1,1} = thin4';
    anovaDataL4{1,2} = {thin4group'};
    anovaDataL4{2,1} = pale4';
    anovaDataL4{2,2} = {pale4group'};
end

% create matrix for kw analysis of Layers 4A,b and C
tempSize = max([size(kwData{3,1},1) size(kwData{4,1},1) size(kwData{5,1},1)]);
ind2 = 1;
for stripe = 1:2
    ind = 1;
    
    kwDataL4Avg{stripe} = zeros(tempSize,3);
    kwDataL4{stripe} = zeros(tempSize,3);
    if midsOnly
        kwDataL4Injects{stripe} = zeros(tempSize,3*3);
    else
        kwDataL4Injects{stripe} = zeros(tempSize,3*7);
    end
    for lyr = 3:5  
        kwDataL4Avg{stripe}(:,lyr-2) = [kwDataAvg{lyr,1}(:,stripe)' ones(1,tempSize-size(kwDataAvg{lyr,1},1))*-1]';
        kwDataL4{stripe}(:,lyr-2) = [kwData{lyr,1}(:,stripe)' ones(1,tempSize-size(kwData{lyr,1},1))*-1]';
        if midsOnly
            kwDataL4Injects{stripe}(:,ind:ind+2) = [kwDataInjects{lyr,1}(:,ind2:ind2+2)' ones(3,tempSize-size(kwDataInjects{lyr,1},1))*-1]';
            ind = ind + 3;
        else
            kwDataL4Injects{stripe}(:,ind:ind+6) = [kwDataInjects{lyr,1}(:,ind2:ind2+6)' ones(7,tempSize-size(kwDataInjects{lyr,1},1))*-1]';
            ind = ind + 7;
        end
    end
    kwDataL4Avg{stripe}(kwDataL4Avg{stripe} == -1) = NaN;
    kwDataL4{stripe}(kwDataL4{stripe} == -1) = NaN;
    kwDataL4Injects{stripe}(kwDataL4Injects{stripe} == -1) = NaN;
    if midsOnly
        ind2 = ind2 + 3;
    else
        ind2 = ind2 + 7;
    end
end

%% Compare individual layers betweem pale and thin stripes
for lyr = 1:size(corticalLayers,2) 
    if performANOVA
        % perform ANOVA n-way analysis
        [AnovaLayersPvalues(lyr) AnovaLayersTables{lyr} AnovaLayersStats{lyr}] = anovan(anovaData{lyr,1}, ...
            anovaData{lyr,2},'varnames',{char(corticalLayers{lyr})});
        figure;
        [AnovaLayersComparisonResults{lyr}, AnovaLayersMeansAndErrors{lyr}, fig1] = multcompare(AnovaLayersStats{lyr});
        figName = strcat("ANOVA_",corticalLayers{lyr});
        set(fig1,'name',figName);
    end
    
    [KwLayersPvalues(lyr) KwLayersTables{lyr} KwLayersStats{lyr}] = kruskalwallis(kwData{lyr,1},{'Thin' 'Pale'});
    [KwLayersComparisonResults{lyr}, KwLayersMeansAndErrors{lyr}, fig2] = multcompare(KwLayersStats{lyr});
    figName = strcat("KW_",corticalLayers{lyr});
    set(fig2,'name',figName);
    
    [KwLayersAvgPvalues(lyr) KwLayersAvgTables{lyr} KwLayersAvgStats{lyr}] = kruskalwallis(kwDataAvg{lyr,1},{'Thin' 'Pale'});
    [KwLayersAvgComparisonResults{lyr}, KwLayersAvgMeansAndErrors{lyr}, fig2] = multcompare(KwLayersAvgStats{lyr});
    figName = strcat("KW_Avg",corticalLayers{lyr});
    set(fig2,'name',figName);
    
    if midsOnly
        kwGroups = {'Thin' 'Thin' 'Thin' 'Pale' 'Pale' 'Pale'};
        [KwLayersAllPvalues(lyr) KwLayersAllTables{lyr} KwLayersAllStats{lyr}] = kruskalwallis(kwDataInjects{lyr,1},kwGroups);
        [KwLayersAllComparisonResults{lyr}, KwLayersAllMeansAndErrors{lyr}, fig3] = multcompare(KwLayersAllStats{lyr});
        figName = strcat("KW_ALLData_",corticalLayers{lyr});
        set(fig2,'name',figName);
    else
        kwGroups = {'Thin' 'Thin' 'Thin' 'Thin' 'Thin' 'Thin' 'Thin' 'Pale' 'Pale' 'Pale' 'Pale' 'Pale' 'Pale' 'Pale'};
        [KwLayersAllPvalues(lyr) KwLayersAllTables{lyr} KwLayersAllStats{lyr}] = kruskalwallis(kwDataInjects{lyr,1},kwGroups);
        [KwLayersAllComparisonResults{lyr}, KwLayersAllMeansAndErrors{lyr}, fig3] = multcompare(KwLayersAllStats{lyr});
        figName = strcat("KW_ALLData_",corticalLayers{lyr});
        set(fig2,'name',figName);
    end
    
end

%% Compare layers 4a, 4b, and 4c for each stripe type
for stripe = 1:2
    if performANOVA
        % perform ANOVA n-way analysis
        [AnovaL4Pvalues(stripe) AnovaL4Tables{stripe} AnovaL4Stats{stripe}] = anovan(anovaDataL4{stripe,1}, ...
            anovaDataL4{stripe,2},'varnames',{char(stripes{stripe})});
        figure;
        [AnovaL4ComparisonResults{stripe}, AnovaL4MeansAndErrors{stripe}, fig1] = multcompare(AnovaL4Stats{stripe});
        figName = strcat("ANOVA_",stripes{stripe},'_Layer 4');
        set(fig1,'name',figName);
    end
    
    [KwL4Pvalues(stripe) KwL4Tables{stripe} KwL4Stats{stripe}] = kruskalwallis(kwDataL4{stripe},{'L4A' 'L4B' 'L4C'});
    [KwL4ComparisonResults{stripe}, KwL4MeansAndErrors{stripe}, fig2] = multcompare(KwL4Stats{stripe});
    figName = strcat("KW_",stripes{stripe},"_Layer 4");
    set(fig2,'name',figName);
    
    [KwL4AvgPvalues(stripe) KwL4AvgTables{stripe} KwL4AvgStats{stripe}] = kruskalwallis(kwDataL4Avg{stripe},{'L4A' 'L4B' 'L4C'});
    [KwL4AvgComparisonResults{stripe}, KwL4AvgMeansAndErrors{stripe}, fig2] = multcompare(KwL4AvgStats{stripe});
    figName = strcat("KW_",stripes{stripe},"_Layer 4_Avg");
    set(fig2,'name',figName);
    
    if ~midsOnly
        kwGroupsL4 = {'L4A' 'L4A' 'L4A' 'L4A' 'L4A' 'L4A' 'L4A' 'L4B' 'L4B' 'L4B' 'L4B' ...
            'L4B' 'L4B' 'L4B' 'L4C' 'L4C' 'L4C' 'L4C' 'L4C' 'L4C' 'L4C'};
        [KwL4AllPvalues(stripe) KwL4AllTables{stripe} KwL4AllStats{stripe}] = kruskalwallis(kwDataL4Injects{stripe},kwGroupsL4);
        [KwL4AllComparisonResults{stripe}, KwL4AllMeansAndErrors{stripe}, fig3] = multcompare(KwL4AllStats{stripe});
        figName = strcat("KW_",stripes{stripe},"_Layer 4_ALL");
        set(fig3,'name',figName);
    else
        kwGroupsL4 = {'L4A' 'L4A' 'L4A' 'L4B' 'L4B' 'L4B' 'L4C' 'L4C' 'L4C'};
        [KwL4AllPvalues(stripe) KwL4AllTables{stripe} KwL4AllStats{stripe}] = kruskalwallis(kwDataL4Injects{stripe},kwGroupsL4);
        [KwL4AllComparisonResults{stripe}, KwL4AllMeansAndErrors{stripe}, fig3] = multcompare(KwL4AllStats{stripe});
        figName = strcat("KW_",stripes{stripe},"_Layer 4_ALL");
        set(fig3,'name',figName);
    end
end
       
% clear brdr brdrs cs currIntens fig1 fig2 figName fileDate inj ...
%     lyr ...
%     matSize pale4 pale4group paleCount paleInjections paleIntens sectionImages ...
%     stripe stripes tempMat tempSize thin4 thin4group thinCount thinInjections ...
%     thinIntens pathName anovaData anovaDataL4 cs2 curr2 fig3 ind ind2 injInd ...
%     matSize2 paleCount2 paleCount3 thinCount2 thinCount3 paleIntensInjects ...
%     thinIntensInjects 

save(filenameWS);

%% Gather comparison data
if performANOVA
    statsLayers = {'AnovaLayersComparisonResults' 'KwLayersComparisonResults' ...
        'KwLayersAvgComparisonResults' 'KwLayersAllComparisonResults' ...
        'AnovaL4ComparisonResults' 'KwL4ComparisonResults' ...
        'KwL4AvgComparisonResults' 'KwL4AllComparisonResults'};
    randVar = 5;
else
    statsLayers = {'KwLayersComparisonResults' 'KwLayersAvgComparisonResults' ...
        'KwLayersAllComparisonResults' 'KwL4ComparisonResults' ...
        'KwL4AvgComparisonResults' 'KwL4AllComparisonResults'};  
    randVar = 4;
end
for stat = 1:size(statsLayers,2)
    comparisonStats = [];
    if stat < randVar
        for lyr = 1:size(corticalLayers,2)
            temp = [];
            eval(['temp =' statsLayers{stat} '{1,lyr};']);
            comparisonStats = [comparisonStats; temp];
            layerMat = [corticalLayers{lyr} corticalLayers{lyr} corticalLayers{lyr}];
        end
        if midsOnly
            laminarStatisticsMiddleSections(stat).StatisticType = statsLayers{stat};
            laminarStatisticsMiddleSections(stat).Data = comparisonStats; 
        else
            laminarStatistics(stat).StatisticType = statsLayers{stat};
            laminarStatistics(stat).Data = comparisonStats; 
        end
    else
        for lyr = 1:2
            temp = [];
            eval(['temp =' statsLayers{stat} '{1,lyr};']);
            comparisonStats = [comparisonStats; temp];
            layerMat = [corticalLayers{lyr} corticalLayers{lyr} corticalLayers{lyr}];
        end
        if midsOnly
            laminarStatisticsMiddleSections(stat).StatisticType = statsLayers{stat};
            laminarStatisticsMiddleSections(stat).Data = comparisonStats;
        else
            laminarStatistics(stat).StatisticType = statsLayers{stat};
            laminarStatistics(stat).Data = comparisonStats;
        end
    end
    
end
% clear layers temp comparisonStats stat lyr statsLayers
if midsOnly
    fN = strcat(filenameStats,'_MIDSonly');
    save(fN,'laminarStatisticsMiddleSections'); clear fN;
else
    save(filenameStats,'laminarStatistics');
end

% %% Save data to excel
% for lyr = 1:size(kwData,1)
%     xlswrite('StripeVsStripeStats',kwData
%     