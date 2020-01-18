% Plot heatmaps and intensities for each layer for each case
% Plots include both heatmaps and quantifications
% Use script "fb2co_heatmap_3" to plot just the heatmaps
% Saved figures are in the "Plots" folder

clear all; 
close all;
path( path, fullfile(pwd,'scriptFuncs') )

wantSaveFigs = 1; %save figs?

monkey = 374;               % monkey number
loc = 'lateral';             % only for mk373 and mk374
fileDate = '2019-07-19';    % date for saved layer heatmap.

rotAngle = 20;              % angle to rotate the lines
bins = 20;                  % number of bins for normalization of line lengths
barWidth = 200;             % microns
hemi = 'RH';                % left or right hemisphere

method = 'pick_mean';       %picked bar based on highest mean
% method = 'all_mean';        %averaged all bars
% method = 'adjust';          %minimized effects of bars with no label 
% method = 'pick_skew';       %picked bar based on skew
% method = 'pick_sum';       %picked bar based on skew

layers = {'L1A' 'L1B' 'L2_3' 'L4B' 'L5_6' 'P'};
% layers = {'L1B'};
if monkey == 373
    layers = {'L1B' 'L2_3' 'L4B' 'L5_6' 'P'};
end

if monkey == 373 || monkey == 374
    path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/',loc,'/Mat Files'), path );
    path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/',loc), path );
else
    path( fullfile(pwd,strcat('MK',num2str(monkey),hemi,'/Mat Files')), path );
    path( fullfile(pwd,strcat('MK',num2str(monkey),hemi)), path );
end

numLayers = size(layers,2);
for lyr = 1:numLayers
    if monkey(1) == 373 || monkey(1) == 374
        load(strcat('MK',num2str(monkey),hemi,loc,'_',fileDate,'_',layers{lyr},'_Heatmap.mat'));
        figName = strcat('mk',num2str(monkey),loc,'_',layers{lyr});
        eval(['hm = mk' num2str(monkey) loc '_' layers{lyr} '_heatmap;']);

    else
        load(strcat('MK',num2str(monkey),hemi,'_',fileDate,'_',layers{lyr},'_Heatmap.mat'));
        figName = strcat('mk',num2str(monkey),'_',layers{lyr});
        eval(['hm = mk' num2str(monkey) '_' layers{lyr} '_heatmap;']);
    end
    
    hmTicks = linspace(0,size(hm(1).heatMap,1),41);
    
    if strcmp(layers{lyr},'P')
        fig = figure('name',figName,'units','normalized','outerposition',[0 0 1 1]); 
        subplot(2,2,1);
            imagesc(hm(1).heatMap); colormap(jet(numel(hm(1).heatMap))); colorbar;
            viscircles(hm(1).blobCenter, hm(1).avgBlobRadius,'Color','m');
            axis square;
            xticks(hmTicks); set(gca,'xticklabel',[]);
            yticks(hmTicks); set(gca,'yticklabel',[]);
            title('Processed');
        subplot(2,2,3);
            if strcmp(method,'pick_mean')
                bestbar = hm(1).bestBarMean;
                intens = hm(1).barlineIntensities{1,bestbar};
            elseif strcmp(method,'pick_skew')
                bestbar = hm(1).bestBarSkew;
                intens = hm(1).barlineIntensities{1,bestbar};
            elseif strcmp(method,'all_mean')
                intens = avgIntens(hm(1).barlineIntensities);
            elseif strcmp(method,'adjust')
                intens = avgIntens(hm(1).barlineIntensitiesAdj);
            end
            intens = intens(1:20);
            intens = (intens - min(intens))/(max(intens)-min(intens));
            plot(intens,'-k','LineWidth',2);
            xlim([1 20]);
            ylim([0 1]);
        subplot(2,2,2);
%             imagesc(hm(3).heatMap); colormap(jet(numel(hm(3).heatMap))); colorbar;
            imagesc(hm(1).heatMap); colormap('gray'); colorbar;
            viscircles(hm(1).blobCenter, hm(1).avgBlobRadius,'Color','m');
            axis square;
            xticks(hmTicks); set(gca,'xticklabel',[]);
            yticks(hmTicks); set(gca,'yticklabel',[]);
            title('Processed');
        subplot(2,2,4);
            if strcmp(method,'pick_mean')
                bestbar = hm(1).bestBarMean;
                intens = hm(1).barlineIntensities{1,bestbar};
            elseif strcmp(method,'pick_skew')
                bestbar = hm(1).bestBarSkew;
                intens = hm(1).barlineIntensities{1,bestbar};
            elseif strcmp(method,'all_mean')
                intens = avgIntens(hm(1).barlineIntensities);
            elseif strcmp(method,'adjust')
                intens = avgIntens(hm(1).barlineIntensitiesAdj);
            end
            intens = intens(1:20);
            intens = (intens - min(intens))/(max(intens)-min(intens));
            plot(intens,'-k','LineWidth',2);
            xlim([1 20]);
            ylim([0 1]);
    else
        fig = figure('name',figName,'units','normalized','outerposition',[0 0 1 1]); 
        subplot(2,2,1);
            if strcmp(method,'pick_mean11')
                tempImg = hm(1).heatMap;
                bbSkew = hm(1).bestBarSkew;
                bbAvg = hm(1).bestBarMean;
%                 intensSkew = hm(1).barlineIntensities{1,bbSkew};
                intensAvg = hm(1).barlineIntensities{1,bbAvg};
                bars = [bbAvg];
                for br = bars
                    numBins = size(hm(1).barlines{1,br},2);       
                    for bn = 1:2:numBins
                        numPix = size(hm(1).barlines{1,br}{1,bn},1);
                        for px = 1:numPix
                        r = hm(1).barlines{1,br}{1,bn}(px,1);
                        c = hm(1).barlines{1,br}{1,bn}(px,2);
%                             if br == bbSkew       
%                                 tempImg(r,c)=1;
%                             else
                                tempImg(r,c) = 0;
%                             end
                        end
                    end
                end
                imagesc(tempImg); colormap(jet(numel(tempImg))); colorbar;
                axis('square');
                xticks(hmTicks); set(gca,'xticklabel',[]);
                yticks(hmTicks); set(gca,'yticklabel',[]);
                title('GFP');
            else
                imagesc(hm(1).heatMap); colormap(jet(numel(hm(1).heatMap))); colorbar;
                viscircles(hm(1).blobCenter, hm(1).avgBlobRadius,'Color','k','linewidth',2);
                axis square;
                xticks(hmTicks); set(gca,'xticklabel',[]);
                yticks(hmTicks); set(gca,'yticklabel',[]);
                title('GFP');
            end
        subplot(2,2,2);
            if strcmp(method,'pick_mean11')
                tempImg = hm(2).heatMap;
                bbSkew = hm(2).bestBarSkew;
                bbAvg = hm(2).bestBarMean;
%                 intensSkew = hm(2).barlineIntensities{1,bbSkew};
                intensAvg = hm(2).barlineIntensities{1,bbAvg};
                bars = [bbAvg];
                for br = bars
                    numBins = size(hm(1).barlines{1,br},2);       
                    for bn = 1:2:numBins
                        numPix = size(hm(1).barlines{1,br}{1,bn},1);
                        for px = 1:numPix
                        r = hm(1).barlines{1,br}{1,bn}(px,1);
                        c = hm(1).barlines{1,br}{1,bn}(px,2);
%                             if br == bbSkew       
%                                 tempImg(r,c)=1;
%                             else
                                tempImg(r,c) = 0;
%                             end
                        end
                    end
                end
                imagesc(tempImg); colormap(jet(numel(tempImg))); colorbar;
                axis('square');
                xticks(hmTicks); set(gca,'xticklabel',[]);
                yticks(hmTicks); set(gca,'yticklabel',[]);
                title('TDT');
            else
                imagesc(hm(2).heatMap); colormap(jet(numel(hm(2).heatMap))); colorbar;
                viscircles(hm(2).blobCenter, hm(2).avgBlobRadius,'Color','k','linewidth',2);
                axis square;
                xticks(hmTicks); set(gca,'xticklabel',[]);
                yticks(hmTicks); set(gca,'yticklabel',[]);
                title('TDT');
            end
        subplot(2,2,3);
            if strcmp(method,'pick_mean')
                bestbar = hm(1).bestBarMean;
                intens = hm(1).barlineIntensities{1,bestbar};
            elseif strcmp(method,'pick_skew')
                bestbar = hm(1).bestBarSkew;
                intens = hm(1).barlineIntensities{1,bestbar};
            elseif strcmp(method,'all_mean')
                intens = avgIntens(hm(1).barlineIntensities);
            elseif strcmp(method,'pick_sum')
                bestbar = hm(1).bestBarSum;
                intens = hm(1).barlineIntensities{1,bestbar};
            elseif strcmp(method,'adjust')
                intens = avgIntens(hm(1).barlineIntensitiesAdj);
            end 
            
            intens = (intens(1:20) - min(intens(1:20)))/(max(intens(1:20))-min(intens(1:20)));
            xval = linspace(1,size(intens(1:20),2),size(intens(1:20),2));
            tl = polyfit(xval,intens(1:20),3);
            intensTL = polyval(tl,xval);
            plot(intens(1:20),'*b','LineWidth',2); hold on;
            plot(intensTL,'-b','LineWidth',2);
            xlim([1 20]);
            ylim([0 1]);

        subplot(2,2,4);
            if strcmp(method,'pick_mean')
                bestbar = hm(2).bestBarMean;
                intens = hm(2).barlineIntensities{1,bestbar};
            elseif strcmp(method,'pick_skew')
                bestbar = hm(2).bestBarSkew;
                intens = hm(2).barlineIntensities{1,bestbar};
            elseif strcmp(method,'pick_sum')
                bestbar = hm(2).bestBarSum;
                intens = hm(2).barlineIntensities{1,bestbar};
            elseif strcmp(method,'all_mean')
                intens = avgIntens(hm(2).barlineIntensities);
            elseif strcmp(method,'adjust')
                intens = avgIntens(hm(2).barlineIntensitiesAdj);
            end
            
            intens = (intens(1:20) - min(intens(1:20)))/(max(intens(1:20))-min(intens(1:20)));
            xval = linspace(1,size(intens(1:20),2),size(intens(1:20),2));
            tl = polyfit(xval,intens(1:20),3);
            intensTL = polyval(tl,xval);
            plot(intens(1:20),'*r','LineWidth',2); hold on;
            plot(intensTL,'-r','LineWidth',2);
            xlim([1 20]);
            ylim([0 1]);

    end
    if wantSaveFigs
        fN = fullfile(pwd,strcat('/Plots/HeatmapsWithQuantification_',figName,'.tiff'));
        saveas(fig,fN,'tiff');
        fN = fullfile(pwd,strcat('/Plots/HeatmapsWithQuantification_',figName,'.eps'));
        saveas(fig,fN,'epsc');
    end
end

function output = avgIntens(matr)
    avgInten = [];
    for bin = 1:20
        numBars = size(matr,2);
        intens = [];
        for bar = 1:numBars
            intens(bar) = matr{1,bar}(bin);
        end
        avgInten(bin) = mean(intens);
    end
    output = avgInten;
end
    