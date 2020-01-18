% Get and plot single heatmap for each layer of the given monkey case
% Saved File: Mky_Layer_Heatmap

clear all; 
close all;
path( path, fullfile(pwd,'scriptFuncs') )

%Enter variables
rotAngle = 20;              % angle to rotate the lines
bins = 20;                  % number of bins for normalization of line lengths
barWidth = 200;             % microns
hemi = 'RH';                % left or right hemisphere

runDate = '2019-07-19';
monkeys = [356 365 373 374];               % monkey number
injections373 = {'medial'};             % only for mk373 and mk374
injections = {'lateral' 'medial'};
fileDate = '2019-07-18';    % date for blob heatmaps file

wantSaveFigs = 1; %save figures?

layersAll = {'L1A' 'L1B' 'L2_3' 'L4B' 'L5_6' 'P'};
layers373 = {'L1B' 'L2_3' 'L4B' 'L5_6' 'P'};

for monkey = monkeys
    if ismember(monkey,[356 365])
        locations = {'NA'};
    elseif monkey == 373
        locations = injections373;
    else
        locations = injections;
    end
    
    for loc = locations
        loc = char(loc);

        % assign the path to the storage data folder
        if monkey == 373 || monkey == 374
            fileName = strcat('/MK',num2str(monkey),hemi,'/',loc,'/Mat Files/MK',num2str(monkey),hemi,loc,'_',runDate,'_');
            path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/',loc,'/Mat Files'), path );
            path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/',loc), path );
        else
            fileName = strcat('/MK',num2str(monkey),hemi,'/Mat Files/MK',num2str(monkey),hemi,'_',runDate,'_');
            path( fullfile(pwd,strcat('MK',num2str(monkey),hemi,'/Mat Files')), path );
            path( fullfile(pwd,strcat('MK',num2str(monkey),hemi)), path );
        end

        % import the excel spreadsheet that contains the relevant section and blobs
        % that will be used for the analysis
        if monkey == 356
            blobData = importdata('MK356 label quality and locations.xlsx');
            layers = layersAll;
        elseif monkey == 365
            blobData = importdata('MK365 label quality and locations.xlsx');
            layers = layersAll;
        elseif monkey == 373
            if strcmp(loc,'medial')
                blobData = importdata('MK373 medial label quality and locations.xlsm');
            elseif strcmp(loc,'lateral')
                blobData = importdata('MK373 lateral label quality and locations.xlsm');
            end
            layers = layers373;

        elseif monkey == 374
            if strcmp(loc,'medial')
                blobData = importdata('MK374 medial label quality and locations.xlsm');
            elseif strcmp(loc,'lateral')
                blobData = importdata('MK374 lateral label quality and locations.xlsm');
            end
            layers = layersAll;
        end

        % combine heatmaps from the sections and plot; draw barlines around the blobs and determine
        % the barline with the highest average intensity
        numLayers = size(layers,2);
        for lyr = 1:numLayers
            if monkey == 373 || monkey == 374
                eval(['mk' num2str(monkey) loc '_' layers{lyr} '_heatmap = plotHeatMap(blobData,layers{lyr},fileDate,hemi,wantSaveFigs,loc);' ]);
                eval(['mk' num2str(monkey) loc '_' layers{lyr} '_heatmap = getBarlinesSingleHeatmap(monkey,mk'...
                    num2str(monkey) loc '_' layers{lyr} '_heatmap,rotAngle,barWidth,bins,loc);']);
                eval(['mk' num2str(monkey) loc '_' layers{lyr} '_heatmap = getBarlineIntensitiesSingleHeatmap(mk'...
                    num2str(monkey) loc '_' layers{lyr} '_heatmap);']);
                eval(['mk' num2str(monkey) loc '_' layers{lyr} '_heatmap = getBestBarSingleHeatmap(mk'...
                    num2str(monkey) loc '_' layers{lyr} '_heatmap);']);
        %         eval(['mk' num2str(monkey) loc '_' layers{lyr} '_heatmap = adjustBarlinesSingleHeatmap(mk'...
        %             num2str(monkey) loc '_' layers{lyr} '_heatmap);']);

                fN = fullfile(pwd,strcat(fileName,layers{lyr},'_Heatmap.mat'));
                fprintf(strcat('Saving heatmap for mk',num2str(monkey),loc,'_',layers{lyr},'...'));
                save(fN,strcat('mk',num2str(monkey),loc,'_',layers{lyr},'_heatmap'));
                fprintf('DONE! \n');

            else
                eval(['mk' num2str(monkey) '_' layers{lyr} '_heatmap = plotHeatMap(blobData,layers{lyr},fileDate,hemi,wantSaveFigs);' ]);
                eval(['mk' num2str(monkey) '_' layers{lyr} '_heatmap = getBarlinesSingleHeatmap(monkey,mk'...
                    num2str(monkey) '_' layers{lyr} '_heatmap,rotAngle,barWidth,bins);']);
                eval(['mk' num2str(monkey) '_' layers{lyr} '_heatmap = getBarlineIntensitiesSingleHeatmap(mk'...
                    num2str(monkey) '_' layers{lyr} '_heatmap);']);
                eval(['mk' num2str(monkey) '_' layers{lyr} '_heatmap = getBestBarSingleHeatmap(mk'...
                    num2str(monkey) '_' layers{lyr} '_heatmap);']);
        %         eval(['mk' num2str(monkey) '_' layers{lyr} '_heatmap = adjustBarlinesSingleHeatmap(mk'...
        %             num2str(monkey) '_' layers{lyr} '_heatmap);']);

                fN = fullfile(pwd,strcat(fileName,layers{lyr},'_Heatmap.mat'));
                fprintf(strcat('Saving heatmap for mk',num2str(monkey),'_',layers{lyr},'...'));
                save(fN,strcat('mk',num2str(monkey),'_',layers{lyr},'_heatmap'));
                fprintf('DONE! \n');

            end
        end
    end
end
clear barWidth bins blobData fileDate fileName fN hemi injections layers ...
    injections373 layers373 layersAll loc locations lyr monkey monkeys ...
    numLayers rotAngle 

