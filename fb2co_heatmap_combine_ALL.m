% Combine all cases for the thick, thin, and pale stripes
% Saved File: Layer_Heatmap

clear all;
close all;

wantSaveFigs = 1; %save figures?

runDate = '2019-10-19';
monkeys = [356 365 373 374]; %can take out or add cases from the 4 cases used
hemi = 'RH';
rotAngle = 20;
barWidth = 200;
bins = 20;
labels = {'GFP' 'TDT' 'CO'};
layers = {'L1A' 'L1B' 'L2_3' 'L4B' 'L5_6' 'P'};
fileDate = '2019-07-19';
stripes = {'Thick' 'Thin' 'Pale'};
casesThick = {'MK356RH','GFP';'MK365RH','GFP';'MK374RHmedial','TDT'};
casesThin = {'MK356RH','TDT';'MK365RH','TDT';'MK374RHlateral','GFP'};
casesPale = {'MK374RHmedial','GFP';'MK374RHlateral','TDT';'MK373RHmedial','GFP'};
casesProcessed = {'MK356RH','CO'};

path( fullfile(pwd,'scriptFuncs'), path );

% scales
convs = xlsread('FBtoCO um to pixel conversions.xlsx');

%load all the heatmaps for all cases
for monkey = monkeys
    if monkey == 374 || monkey == 373
        locations = {'lateral' 'medial'};
        for loc = 1:size(locations,2)
            if monkey == 373 && strcmp(locations{loc},'lateral')
                continue
            end
            path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/',locations{loc},'/Mat Files'), path );
            path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/',locations{loc}), path );
            for lyr = 1:size(layers,2)
                if monkey == 373 && strcmp(layers{lyr},'L1A')
                    continue
                end
                load(strcat('MK',num2str(monkey),hemi,locations{loc},'_',fileDate,'_',layers{lyr},'_Heatmap.mat'));
            end
        end
    else
        path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/Mat Files'), path );
        path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/'), path );
        for lyr = 1:size(layers,2)
            load(strcat('MK',num2str(monkey),hemi,'_',fileDate,'_',layers{lyr},'_Heatmap.mat'));
        end
    end
end

%create a matrix for all the cases and their heatmaps (except processed
%images)
for lyr = 1:size(layers,2)
    if strcmp(layers{lyr},'P')
        numLabels = 1;
    else
        numLabels = 2;
    end
    ctr = 1; % accounts for all the monkeys in the layer
    for monkey = monkeys
        if monkey == 356
            px2um = convs(1,1);
        elseif monkey == 365
            px2um = convs(1,2);
        elseif monkey == 373
            px2um = convs(1,3);
        elseif monkey == 374
            px2um = convs(1,4);
        else
            px2um = 1;
        end
        
        if monkey == 374 || monkey == 373
            if monkey == 373 && strcmp(layers{lyr},'L1A')
                continue
            end
            for loc = 1:size(locations,2)
                if monkey == 373 && strcmp(locations{loc},'lateral')
                    continue
                end
                for label = 1:numLabels
                    eval([layers{lyr} '_heatmaps(ctr).monkey = monkey;']);
                    eval([layers{lyr} '_heatmaps(ctr).case = ''MK' num2str(monkey) hemi locations{loc} ''';']);
                    eval([layers{lyr} '_heatmaps(ctr).label = mk' num2str(monkey) locations{loc} '_' layers{lyr} '_heatmap(label).label;']);
                    eval([layers{lyr} '_heatmaps(ctr).heatMap = mk' num2str(monkey) locations{loc} '_' layers{lyr} '_heatmap(label).heatMap;']);
                    eval([layers{lyr} '_heatmaps(ctr).avgBlobRadiusPx = mk' num2str(monkey) locations{loc} '_' layers{lyr} '_heatmap(label).avgBlobRadius;']);                    
                    eval([layers{lyr} '_heatmaps(ctr).avgInterblobRadiusPx = mk' num2str(monkey) locations{loc} '_' layers{lyr} '_heatmap(label).avgInterblobRadius;']);
                    eval([layers{lyr} '_heatmaps(ctr).avgBlobRadiusMicron = (mk' num2str(monkey) locations{loc} '_' layers{lyr} '_heatmap(label).avgBlobRadius)/px2um;']);
                    eval([layers{lyr} '_heatmaps(ctr).avgInterblobRadiusMicron = (mk' num2str(monkey) locations{loc} '_' layers{lyr} '_heatmap(label).avgInterblobRadius)/px2um;']);
                    eval([layers{lyr} '_heatmaps(ctr).blobCenterPx = mk' num2str(monkey) locations{loc} '_' layers{lyr} '_heatmap(label).blobCenter;']);
                    eval(['bb = mk' num2str(monkey) locations{loc} '_' layers{lyr} '_heatmap(label).bestBarMean;']);
                    eval([layers{lyr} '_heatmaps(ctr).bestBarIntens = mk' num2str(monkey) locations{loc} '_' layers{lyr}... 
                        '_heatmap(label).barlineIntensities{1,bb};']);
                    ctr = ctr + 1;
                end                
            end   
        else
            for label = 1:numLabels
                eval([layers{lyr} '_heatmaps(ctr).monkey = monkey;']);
                eval([layers{lyr} '_heatmaps(ctr).case = ''MK' num2str(monkey) hemi ''';']);
                eval([layers{lyr} '_heatmaps(ctr).label = mk' num2str(monkey) '_' layers{lyr} '_heatmap(label).label;']);
                eval([layers{lyr} '_heatmaps(ctr).heatMap = mk' num2str(monkey) '_' layers{lyr} '_heatmap(label).heatMap;']);
                eval([layers{lyr} '_heatmaps(ctr).avgBlobRadiusPx = mk' num2str(monkey) '_' layers{lyr} '_heatmap(label).avgBlobRadius;']);
                eval([layers{lyr} '_heatmaps(ctr).avgInterblobRadiusPx = mk' num2str(monkey) '_' layers{lyr} '_heatmap(label).avgInterblobRadius;']);
                eval([layers{lyr} '_heatmaps(ctr).avgBlobRadiusMicron = (mk' num2str(monkey) '_' layers{lyr} '_heatmap(label).avgBlobRadius)/px2um;']);
                eval([layers{lyr} '_heatmaps(ctr).avgInterblobRadiusMicron = (mk' num2str(monkey) '_' layers{lyr} '_heatmap(label).avgInterblobRadius)/px2um;']);
                eval([layers{lyr} '_heatmaps(ctr).blobCenterPx = mk' num2str(monkey) '_' layers{lyr} '_heatmap(label).blobCenter;']);
                eval(['bb = mk' num2str(monkey) '_' layers{lyr} '_heatmap(label).bestBarMean;']);
                eval([layers{lyr} '_heatmaps(ctr).bestBarIntens = mk' num2str(monkey) '_' layers{lyr}... 
                    '_heatmap(label).barlineIntensities{1,bb};']);
                ctr = ctr + 1;
            end
        end
    end
end

clear label lyr monkey loc;

%find biggest heatmap and record its size (one side)
for lyr = 1:size(layers,2)
    eval(['numCases = size(' layers{lyr} '_heatmaps,2);']);
    eval(['temp = ' layers{lyr} '_heatmaps;']);

    maxWindow = 0;
    for cs = 1:numCases
        if size(temp(cs).heatMap,1) > maxWindow
            maxWindow = size(temp(cs).heatMap,1);
            maxInd = cs;
        end
    end

    eval([layers{lyr} '_heatmaps(1).maxSize = maxWindow;'])
    eval([layers{lyr} '_heatmaps(1).stripeBlobCenterPx =' layers{lyr} '_heatmaps(maxInd).blobCenterPx;'])
end

%rescale all the heatmaps to the size of the biggest heatmap case
for lyr = 1:size(layers,2)
    eval(['numCases = size(' layers{lyr} '_heatmaps,2);']);
    eval(['temp = ' layers{lyr} '_heatmaps;']);

    maxWindow = temp(1).maxSize;
    for cs = 1:numCases
        newHM = imresize(temp(cs).heatMap,[maxWindow,maxWindow]);
        eval([layers{lyr} '_heatmaps(cs).scaledHeatMap = newHM;'])
    end
end

%combine all heatmaps per label, get barlines, and select best bar
for lyr = 1:size(layers,2)
    eval(['temp = ' layers{lyr} '_heatmaps;']);
    maxMky = temp(maxInd).monkey;
    
    for stp = 1:size(stripes,2)
        if strcmp(stripes{stp},'Thick')
            if strcmp(layers{lyr},'P')
                cases = casesProcessed;
            else
                cases = casesThick;
            end
        elseif strcmp(stripes{stp},'Thin')
            if strcmp(layers{lyr},'P')
                cases = casesProcessed;
            else
                cases = casesThin;
            end
        elseif strcmp(stripes{stp},'Pale')
            if strcmp(layers{lyr},'P')
                cases = casesProcessed;
            else
                cases = casesPale;
            end
        end
        
        tempHM = zeros(temp(1).maxSize,temp(1).maxSize);
        tempBB = [];
        tempBL = [];
        tempBLmicron = [];
        tempIBL = [];
        tempIBLmicron = [];
        ind = 1;
        for cs = 1:size(cases,1)
            for k = 1:size(temp,2)
                if strcmp(cases{cs,1},temp(k).case) && strcmp(cases{cs,2},temp(k).label)
                    currHMsize = size(temp(k).heatMap,1);
                    tempHM = tempHM + temp(k).scaledHeatMap;
                    tempBB(ind,1:20) = temp(k).bestBarIntens(1:20);
                    tempBB(ind,1:20) = (tempBB(ind,1:20) - min(tempBB(ind,1:20)))/ ...
                        (max(tempBB(ind,1:20))-min(tempBB(ind,1:20)));
                    tempCases{ind,1} = temp(k).monkey;
                    tempCases{ind,2} = strcat(temp(k).case,'_',temp(k).label);
                    tempCases{ind,3} = temp(k).blobCenterPx(1);
                    if temp(1).maxSize > currHMsize
                        scl = temp(1).maxSize/size(temp(k).heatMap,1);
                        tempBL(ind) = scl*temp(k).avgBlobRadiusPx;
                        tempBLmicron(ind) = temp(k).avgBlobRadiusMicron;
                        tempIBL(ind) = scl*temp(k).avgInterblobRadiusPx;
                        tempIBLmicron(ind) = temp(k).avgInterblobRadiusMicron;
                    elseif temp(1).maxSize == currHMsize
                        tempBL(ind) = temp(k).avgBlobRadiusPx;
                        tempBLmicron(ind) = temp(k).avgBlobRadiusMicron;
                        tempIBL(ind) = temp(k).avgInterblobRadiusPx;
                        tempIBLmicron(ind) = temp(k).avgInterblobRadiusMicron;
                    end
                    ind = ind + 1;
                end
            end
        end
        eval([layers{lyr} '_heatmaps(stp).blank = [];']);
        eval([layers{lyr} '_heatmaps(stp).stripe = ''' stripes{stp} ''';']);
        eval([layers{lyr} '_heatmaps(stp).stripeHM = tempHM;']);
        tempAvgIntens = mean(tempBB);
        tempSEM = std(tempBB)./ sqrt(size(tempBB,1));
        eval([layers{lyr} '_heatmaps(stp).bestBars = tempBB;']);
        eval([layers{lyr} '_heatmaps(stp).cases = tempCases;']);
        eval([layers{lyr} '_heatmaps(stp).bestBarAvg = tempAvgIntens;']);
        eval([layers{lyr} '_heatmaps(stp).bestBarSEM = tempSEM;']);
        eval([layers{lyr} '_heatmaps(stp).stripeBlobRadiusPx = mean(tempBL);']);
        eval([layers{lyr} '_heatmaps(stp).stripeBlobRadiusMicron = mean(tempBLmicron);']);
        clear tempBB tempCases temAvgIntens tempSEM tempBL tempBLmicron tempHM
    end
    
    
%     eval([layers{lyr} '_heatmaps = getBarlinesStripeHeatmap(maxMky,' layers{lyr} '_heatmaps,rotAngle,barWidth,bins);']);
%     eval([layers{lyr} '_heatmaps = getBarlineIntensitiesStripeHeatmap(' layers{lyr} '_heatmaps);']);
%     eval([layers{lyr} '_heatmaps = getBestBarSingleHeatmap(' layers{lyr} '_heatmaps);']);

    fN = fullfile(pwd,strcat('/CombinedCases/',runDate,'_',layers{lyr},'_Heatmaps.mat'));
    fprintf(strcat('Saving heatmap for layer :',layers{lyr},'...'));
    save(fN,strcat(layers{lyr},'_heatmaps'));
    fprintf('DONE! \n');
end

% clear temp numCases newHM maxWindow lyr cs tempBB tempHM tempSEM tempBL;
%% Plot
for lyr = 1:size(layers,2)
    eval(['numCases = size(' layers{lyr} '_heatmaps,2);']);
    eval(['temp = ' layers{lyr} '_heatmaps;']);

    

    for fig = 1:2
        if fig == 2
            stripes = [3 3];
        else
            stripes = [1 2];
        end
        figName = strcat(layers{lyr});
        currFig = figure('Name',figName,'units','normalized','outerposition',[0 0 1 1]);

        % plot heatmaps
        for sp = 1:2
            if fig == 1
                stripe = sp;
            else
                stripe = 3;
            end
            
            subplot(2,2,sp);            
            imagesc(temp(stripe).stripeHM); colormap(jet(numel(temp(stripe).stripeHM))); colorbar;
            viscircles(temp(1).stripeBlobCenterPx, temp(stripe).stripeBlobRadiusPx,'Color','m');
            axis('square');
            title(temp(stripe).stripe);
%             yticks([1 230 459]);
%             xticks([1 230 459]);
%             yticklabels({'Interblob Center' 'Blob Center' 'Interblob Center'});
%             xticklabels({'Interblob Center' 'Blob Center' 'Interblob Center'});
            hmTicks = linspace(0,size(temp(stripe).stripeHM,1),41);
            xticks(hmTicks); set(gca,'xticklabel',[]);
            yticks(hmTicks); set(gca,'yticklabel',[]);
            set(gca,'fontsize',15);
        end
        % plot intensity curves
        lineColor = {'-b' '-r' '-g'};
        lineType = {'b.' 'r.' '.g'};
        for sp = 3:4
            if fig == 1
                stripe = sp-2;
            else
                stripe = 3;
            end
            
            subplot(2,2,sp)
            if strcmp(layers{lyr},'P')
                for cs = 1:size(P_heatmaps,2)
                    tTemp(cs,1:20) = temp(cs).bestBarIntens(1:20);
                end
                intens = mean(tTemp);
            else
                intens = temp(stripe).bestBarAvg;
                sem = temp(stripe).bestBarSEM;
            end
            intens = (intens - min(intens))/(max(intens)-min(intens));

            xval = linspace(1,size(intens,2),size(intens,2));
            tl = polyfit(xval,intens,3);
            intensTL = polyval(tl,xval);

            plot(intensTL,lineColor{stripe},'LineWidth',2); hold on;
            if ~strcmp(layers{lyr},'P')
                errorbar((1:20),intens,sem,'vertical',lineType{stripe},'linewidth',1)
            end
            xlim([1 20]);
            ylim([-0.1 1.2]);
            yticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 1.1]);
    %         xticks([1 20]);
    %         xticklabels({'Blob Center' 'Interblob Center'});
            set(gca,'fontsize',15);
        end
        if wantSaveFigs
            if fig == 1
                fN = fullfile(pwd,strcat('/Plots/PopulationThickThin_',figName,'.tiff'));
                saveas(currFig,fN,'tiff');
                fN = fullfile(pwd,strcat('/Plots/PopulationThickThin_',figName,'.eps'));
                saveas(fig,fN,'epsc');
            elseif fig == 2
                fN = fullfile(pwd,strcat('/Plots/PopulationPale_',figName,'.tiff'));
                saveas(currFig,fN,'tiff');
                fN = fullfile(pwd,strcat('/Plots/PopulationPale_',figName,'.eps'));
                saveas(fig,fN,'epsc');
            end
        end
    end
end
%     subplot(2,2,4)
%         if strcmp(layers{lyr},'P')
%             for cs = 1:size(P_heatmaps,2)
%                 tTemp(cs,1:20) = temp(cs).bestBarIntens(1:20);
%             end
%             intens = mean(tTemp);
%         else
%             sem = temp(2).bestBarSEM;
%             intens = temp(2).bestBarAvg;
%         end
%         intens = (intens - min(intens))/(max(intens)-min(intens));
% 
%         xval = linspace(1,size(intens,2),size(intens,2));
%         tl = polyfit(xval,intens,3);
%         intensTL = polyval(tl,xval);
%         
%         plot(intensTL,'-r','LineWidth',2); hold on;
%         if ~strcmp(layers{lyr},'P')
%             errorbar((1:20),intens,sem,'vertical','r.','linewidth',1)
%         end
%         xlim([1 20]);
%         ylim([-0.1 1.2]);
%         yticks([0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 1.1]);
% 
% %         xticks([1 20]);
% %         xticklabels({'Blob Center' 'Interblob Center'});
%         set(gca,'fontsize',15);
% 
% end

clear temp numCases maxWindow lyr cs figName fN ans barWidth bb ...
        bins cases convs ctr currHMsize hemi ind intens intensTL k ...
        labels layers loc locations maxInd maxMky mk356_L1A_heatmap ...
        mk356_L1B_heatmap mk356_L2_3_heatmap mk356_L4B_heatmap ...
        mk356_L5_6_heatmap mk356_P_heatmap mk365_L1B_heatmap ...
        mk365_L2_3_heatmap mk365_L4B_heatmap mk365_L5_6_heatmap ...
        mk365_P_heatmap mk365_L1A_heatmap mk374lateral_L1B_heatmap ...
        mk374lateral_L2_3_heatmap mk374lateral_L4B_heatmap mk374lateral_L5_6_heatmap ...
        mk374lateral_P_heatmap mk374lateral_L1A_heatmap mk374medial_L1B_heatmap ...
        mk374medial_L2_3_heatmap mk374medial_L4B_heatmap mk374medial_L5_6_heatmap ...
        mk374medial_P_heatmap mk374medial_L1A_heatmap monkey monkeys ...
        newHM px2um rotAngle fileDate scl sem stp stripes tempAvgIntens ...
        tempBB tempBL tempBLmicron tempCases tempHM tempIBL tempIBLmicron ...
        tempSEM tl xval tTemp numLabels casesProcessed casesThick casesThin ...
        casesPale fig lineColor lineType mk373medial_L1B_heatmap stripe sp ...
        mk373medial_L2_3_heatmap mk373medial_L4B_heatmap mk373medial_P_heatmap ...
        mk373medial_L5_6_heatmap

