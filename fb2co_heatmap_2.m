% Get the heatmaps for each blob in each section image of the given case
% File: BlobHeatmaps
% saved in the monkey folder under Mat Files

% Extract data from file
clear all; 
close all;

% msgbox('Enter date and other pertinent info');

%Enter variables
rotAngle = 20;              % angle to rotate the lines
bins = 20;                  % number of bins for normalization of line lengths
barWidth = 200;             % microns
hemi = 'RH';                % left or right hemisphere

runDate = '2019-07-18';
monkeys = [356 365 373 374];  % monkey number
injections373 = {'medial'};               % only for mk373
injections = {'lateral' 'medial'};   % for mk374
fileDate = '2019-02-19';      % date for Lines file

% sections used for each case
sections356 = [4 5 6 9 11 19 31];
sections365 = [4 5 9 13 17 21 25 27 30];
sections373Med = [4 6 10 13 20 21];
sections373Lat = [4 6 8 10 13 25];
sections374Med = [9 12 15 16 17 18 21 24 25 26 27 30];
sections374Lat = [5 6 9 12 15 21 24];
 
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
        
        % make path to data storage folder and load the lines_V4 mat file
        if monkey == 373 || monkey == 374
            fileName = strcat('/MK',num2str(monkey),hemi,'/',loc,'/Mat Files/MK',num2str(monkey),hemi,loc,'_',runDate,'_');
            path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/',loc,'/Mat Files'), path );
            path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/',loc), path );
            load(strcat('MK',num2str(monkey),hemi,loc,'_',fileDate,'_Lines_v4.mat'));
            eval(['M = mk' num2str(monkey) loc '_lines;']);
        else
            fileName = strcat('/MK',num2str(monkey),hemi,'/Mat Files/MK',num2str(monkey),hemi,'_',runDate,'_');
            path( fullfile(pwd,strcat('MK',num2str(monkey),hemi,'/Mat Files')), path );
            path( fullfile(pwd,strcat('MK',num2str(monkey),hemi)), path );
            load(strcat('MK',num2str(monkey),hemi,'_',fileDate,'_Lines_v4.mat'));
            eval(['M = mk' num2str(monkey) '_lines;']);
        end

        path( path, fullfile(pwd,'scriptFuncs') )

        % read in images of blob centers, interblob centerlines, and blob outlines
        if monkey == 356
            BCImage = imread('MK356RHcrop ALL Blob Centers Blue Single Pixel ver2.tif');
            BC = color2bw(BCImage,'blue');
            sections = sections356;

        elseif monkey == 365
            BCImage = imread('MK365RH blob centers blue.tif');    
            BC = color2bw(BCImage,'blue');
            sections = sections365;

        elseif monkey == 373

            if strcmp(loc,'medial')
                BCImage = imread('_0013_Blob centers blue.tif');
                sections = sections373Med;

            elseif strcmp(loc, 'lateral')
                BCImage = imread('_0013_Blob centers blue.tif'); 
                sections = sections373Lat;

            end

            BC = color2bw(BCImage,'blue');

        elseif monkey == 374

            if strcmp(loc,'medial')
                BCImage = imread('MK374RH V1 MEDIAL CROP LARGE Blob Centers.tif');
                sections = sections374Med;

            elseif strcmp(loc, 'lateral')
                BCImage = imread('MK374RH V1 stack LATERAL CROP LARGE Blob Centers.tif'); 
                sections = sections374Lat;

            end

            BC = color2bw(BCImage,'blue');

        end

        % Get the coordinates of the blob centers
        [rBC, cBC] = find(BC);

        % Get the average distance from blob centers to interblob centers
        numBlobs = size(M,2);
        ibDist = [];
        for blob = 1:numBlobs
            ibDist(blob) = M(blob).avgIbLengths;
        end
        avgIbDist = ceil(mean(ibDist));

        % Get heat maps for each section
        if monkey == 356
            for sec = sections
                if sec == 9
                    eval(['labelImage = imread(''MK356RHcrop sec' num2str(sec) ' label-250umRollBallSubtract_ver2.tif'');']);

                elseif sec == 11        
                    eval(['labelImage = imread(''MK356RHcrop sec' num2str(sec) ' processed blobs raw.tif'');']);

                else
                    eval(['labelImage = imread(''MK356RHcrop sec' num2str(sec) ' label-250umRollBallSubtract.tif'');']);
                end

                eval(['mk' num2str(monkey) '_sec' num2str(sec) ' = getHeatMaps(mk356_lines,labelImage,avgIbDist,rBC,cBC);']);

                numBlobs = size(mk356_lines,2);
                for blob = 1:numBlobs
                    eval([ 'mk356_sec' num2str(sec) '(blob).heatMapNorm = mk356_sec' num2str(sec) ...
                        '(blob).heatMap/max(max(mk356_sec' num2str(sec) '(blob).heatMap));']);
                end

                fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), '_BlobHeatmaps.mat'));
                fprintf(strcat('Saving heatmaps for mk356sec',num2str(sec),'...'));
                save(fN,strcat('mk',num2str(monkey),'_sec',num2str(sec)));
                fprintf('DONE! \n');

            end

        elseif monkey == 365
            labels = {' GFP' ' TDT'};
            for sec = sections 
                for la = 1:2
                    try
                        label = char(labels(la));
                        if sec == 9
                            labelImage = imread(strcat('MK365RH ',label, ' sec', num2str(sec),' 12um background subtract.tif'));
                        else
                            labelImage = imread(strcat('MK365RH ',label, ' sec', num2str(sec),' 250um background subtract.tif'));
                        end
                        label = label(~isspace(label));
                        eval(['mk' num2str(monkey) '_sec' num2str(sec) label ' = getHeatMaps(mk365_lines,labelImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk365_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) '_sec' num2str(sec) label '(blob).heatMapNorm = mk' num2str(monkey) '_sec' num2str(sec) label ...
                                '(blob).heatMap/max(max(mk' num2str(monkey) '_sec' num2str(sec) label '(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), label,'_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk365sec',num2str(sec),label,'...'));
                        save(fN,strcat('mk',num2str(monkey),'_sec',num2str(sec),label));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('No ',label,' image for mk365sec',num2str(sec),'...skipped! \n'));
                    end
                end
                
                try
                    COImage = imread(strcat('MK365RH sec', num2str(sec),' processed.tif'));
                    eval(['mk' num2str(monkey) '_sec' num2str(sec) 'CO = getHeatMaps(mk365_lines,COImage,avgIbDist,rBC,cBC);']);

                    numBlobs = size(mk365_lines,2);
                    for blob = 1:numBlobs
                        eval([ 'mk' num2str(monkey) '_sec' num2str(sec) 'CO(blob).heatMapNorm = mk' num2str(monkey) '_sec' num2str(sec) ...
                            'CO(blob).heatMap/max(max(mk' num2str(monkey) '_sec' num2str(sec) 'CO(blob).heatMap));']);
                    end

                    fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'CO_BlobHeatmaps.mat'));
                    fprintf(strcat('Saving heatmaps for mk365sec',num2str(sec),'CO...'));
                    save(fN,strcat('mk',num2str(monkey),'_sec',num2str(sec),'CO'));
                    fprintf('DONE! \n');
                catch
                    fprintf(strcat('No CO image for mk365sec',num2str(sec),'...skipped! \n'));
                end
            end

        elseif monkey == 374
            if strcmp(loc,'medial')
                for sec = sections            
                    try
                        eval(['labelImage = imread(''MK374 RH Sec' num2str(sec) ' V1 GFP medial to 18-registered.tif'');']);
                        eval(['mk' num2str(monkey) 'medial_sec' num2str(sec) 'GFP = getHeatMaps(mk374medial_lines,labelImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk374medial_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'medial_sec' num2str(sec) 'GFP(blob).heatMapNorm = mk' num2str(monkey) 'medial_sec' num2str(sec) ...
                                'GFP(blob).heatMap/max(max(mk' num2str(monkey) 'medial_sec' num2str(sec) 'GFP(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'GFP_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk374medial_sec',num2str(sec),'GFP...'));
                        save(fN,strcat('mk',num2str(monkey),'medial_sec',num2str(sec),'GFP'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO GFP label for mk374medial sec',num2str(sec), ' \n'));
                    end

                    try
                        eval(['labelImage = imread(''MK374 RH Sec' num2str(sec) ' V1 TDT medial to 18-registered.tif'');']);
                        eval(['mk' num2str(monkey) 'medial_sec' num2str(sec) 'TDT = getHeatMaps(mk374medial_lines,labelImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk374medial_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'medial_sec' num2str(sec) 'TDT(blob).heatMapNorm = mk' num2str(monkey) 'medial_sec' num2str(sec) ...
                                'TDT(blob).heatMap/max(max(mk' num2str(monkey) 'medial_sec' num2str(sec) 'TDT(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'TDT_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk374medial_sec',num2str(sec),'TDT...'));
                        save(fN,strcat('mk',num2str(monkey),'medial_sec',num2str(sec),'TDT'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO TDT label for mk374medial sec',num2str(sec), ' \n'));
                    end
                    try
                        eval(['COImage = imread(''MK374RH V1 MEDIAL CO ' num2str(sec) ' crop processed.tif'');']);  
                        eval(['mk' num2str(monkey) 'medial_sec' num2str(sec) 'CO = getHeatMaps(mk374medial_lines,COImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk374medial_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'medial_sec' num2str(sec) 'CO(blob).heatMapNorm = mk' num2str(monkey) 'medial_sec' num2str(sec) ...
                                'CO(blob).heatMap/max(max(mk' num2str(monkey) 'medial_sec' num2str(sec) 'CO(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'CO_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk374medial_sec',num2str(sec),'CO...'));
                        save(fN,strcat('mk',num2str(monkey),'medial_sec',num2str(sec),'CO'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO CO image for mk374medial sec',num2str(sec), ' \n'));
                    end
                end
            elseif strcmp(loc,'lateral')
                for sec = sections
                    try
                        eval(['labelImage = imread(''MK374 RH Sec' num2str(sec) ' V1 GFP lateral to 12R-registered.tiff'');']);
                        eval(['mk' num2str(monkey) 'lateral_sec' num2str(sec) 'GFP = getHeatMaps(mk374lateral_lines,labelImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk374lateral_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'lateral_sec' num2str(sec) 'GFP(blob).heatMapNorm = mk' num2str(monkey) 'lateral_sec' num2str(sec) ...
                                'GFP(blob).heatMap/max(max(mk' num2str(monkey) 'lateral_sec' num2str(sec) 'GFP(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'GFP_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk374lateral_sec',num2str(sec),'GFP...'));
                        save(fN,strcat('mk',num2str(monkey),'lateral_sec',num2str(sec),'GFP'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO GFP label for mk374lateral sec',num2str(sec), ' \n'));
                    end

                    try
                        eval(['labelImage = imread(''MK374 RH Sec' num2str(sec) ' V1 TDT lateral to 12R-registered.tiff'');']);
                        eval(['mk' num2str(monkey) 'lateral_sec' num2str(sec) 'TDT = getHeatMaps(mk374lateral_lines,labelImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk374lateral_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'lateral_sec' num2str(sec) 'TDT(blob).heatMapNorm = mk' num2str(monkey) 'lateral_sec' num2str(sec) ...
                                'TDT(blob).heatMap/max(max(mk' num2str(monkey) 'lateral_sec' num2str(sec) 'TDT(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'TDT_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk374lateral_sec',num2str(sec),'TDT...'));
                        save(fN,strcat('mk',num2str(monkey),'lateral_sec',num2str(sec),'TDT'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO TDT label for mk374lateral sec',num2str(sec), ' \n'));
                    end

                    try
                        eval(['COImage = imread(''MK374RH V1 LATERAL CROP LARGE sec' num2str(sec) ' processed blobs ver2.tif'');']);  
                        eval(['mk' num2str(monkey) 'lateral_sec' num2str(sec) 'CO = getHeatMaps(mk374lateral_lines,COImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk374lateral_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'lateral_sec' num2str(sec) 'CO(blob).heatMapNorm = mk' num2str(monkey) 'lateral_sec' num2str(sec) ...
                                'CO(blob).heatMap/max(max(mk' num2str(monkey) 'lateral_sec' num2str(sec) 'CO(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'CO_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk374lateral_sec',num2str(sec),'CO...'));
                        save(fN,strcat('mk',num2str(monkey),'lateral_sec',num2str(sec),'CO'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO CO image for mk374lateral sec',num2str(sec), ' \n'));
                    end
                end
            end

        elseif monkey == 373
            if strcmp(loc,'medial')
                for sec = sections
                    try
                        labelImage = imread(strcat('MK373 med GFP sec', num2str(sec), ' to 10-registered.tif.tif'));         
                        eval(['mk' num2str(monkey) 'medial_sec' num2str(sec) 'GFP = getHeatMaps(mk373medial_lines,labelImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk373medial_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'medial_sec' num2str(sec) 'GFP(blob).heatMapNorm = mk' num2str(monkey) 'medial_sec' num2str(sec) ...
                                'GFP(blob).heatMap/max(max(mk' num2str(monkey) 'medial_sec' num2str(sec) 'GFP(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'GFP_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk373medial_sec',num2str(sec),'GFP...'));
                        save(fN,strcat('mk',num2str(monkey),'medial_sec',num2str(sec),'GFP'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO GFP label image for mk373medial sec',num2str(sec), ' \n'));
                    end

                    try
                        labelImage = imread(strcat('MK373 med TDT sec', num2str(sec), ' to 10-registered.tif.tif')); 
                        eval(['mk' num2str(monkey) 'medial_sec' num2str(sec) 'TDT = getHeatMaps(mk373medial_lines,labelImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk373medial_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'medial_sec' num2str(sec) 'TDT(blob).heatMapNorm = mk' num2str(monkey) 'medial_sec' num2str(sec) ...
                                'TDT(blob).heatMap/max(max(mk' num2str(monkey) 'medial_sec' num2str(sec) 'TDT(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'TDT_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk373medial_sec',num2str(sec),'TDT...'));
                        save(fN,strcat('mk',num2str(monkey),'medial_sec',num2str(sec),'TDT'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO TDT label image for mk373medial sec',num2str(sec), ' \n'));
                    end

                    try
                        COImage = imread(strcat('_0017_MK373 vGlut 1.25x sec', num2str(sec), ' ver3 processed blob.tif.tif')); 
                        eval(['mk' num2str(monkey) 'medial_sec' num2str(sec) 'CO = getHeatMaps(mk373medial_lines,COImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk373medial_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'medial_sec' num2str(sec) 'CO(blob).heatMapNorm = mk' num2str(monkey) 'medial_sec' num2str(sec) ...
                                'CO(blob).heatMap/max(max(mk' num2str(monkey) 'medial_sec' num2str(sec) 'CO(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'CO_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk373medial_sec',num2str(sec),'CO...'));
                        save(fN,strcat('mk',num2str(monkey),'medial_sec',num2str(sec),'CO'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO CO image for mk373medial sec',num2str(sec),' \n'));
                    end
                end
            elseif strcmp(loc,'lateral')
                for sec = sections
                    try
                        labelImage = imread(strcat('MK373 lat GFP sec', num2str(sec), ' to 10-registered.tif.tif'));
                        eval(['mk' num2str(monkey) 'lateral_sec' num2str(sec) 'GFP = getHeatMaps(mk373lateral_lines,labelImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk373lateral_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'lateral_sec' num2str(sec) 'GFP(blob).heatMapNorm = mk' num2str(monkey) 'lateral_sec' num2str(sec) ...
                                'GFP(blob).heatMap/max(max(mk' num2str(monkey) 'lateral_sec' num2str(sec) 'GFP(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'GFP_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk373lateral_sec',num2str(sec),'GFP...'));
                        save(fN,strcat('mk',num2str(monkey),'lateral_sec',num2str(sec),'GFP'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO GFP image for mk373lateral sec',num2str(sec),' \n'));
                    end

                    try
                        labelImage = imread(strcat('MK373 lat TDT sec', num2str(sec), ' to 10-registered.tif.tif'));
                        eval(['mk' num2str(monkey) 'lateral_sec' num2str(sec) 'TDT = getHeatMaps(mk373lateral_lines,labelImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk373lateral_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'lateral_sec' num2str(sec) 'TDT(blob).heatMapNorm = mk' num2str(monkey) 'lateral_sec' num2str(sec) ...
                                'TDT(blob).heatMap/max(max(mk' num2str(monkey) 'lateral_sec' num2str(sec) 'TDT(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'TDT_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk373lateral_sec',num2str(sec),'TDT...'));
                        save(fN,strcat('mk',num2str(monkey),'lateral_sec',num2str(sec),'TDT'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO TDT image for mk373lateral sec',num2str(sec),' \n'));
                    end

                    try
                        labelImage = imread(strcat('_0017_MK373 vGlut 1.25x sec', num2str(sec), ' ver3 processed blob.tif.tif')); 
                        eval(['mk' num2str(monkey) 'lateral_sec' num2str(sec) 'CO = getHeatMaps(mk373lateral_lines,labelImage,avgIbDist,rBC,cBC);']);

                        numBlobs = size(mk373lateral_lines,2);
                        for blob = 1:numBlobs
                            eval([ 'mk' num2str(monkey) 'lateral_sec' num2str(sec) 'CO(blob).heatMapNorm = mk' num2str(monkey) 'lateral_sec' num2str(sec) ...
                                'CO(blob).heatMap/max(max(mk' num2str(monkey) 'lateral_sec' num2str(sec) 'CO(blob).heatMap));']);
                        end

                        fN = fullfile(pwd,strcat(fileName,'sec', num2str(sec), 'CO_BlobHeatmaps.mat'));
                        fprintf(strcat('Saving heatmaps for mk373lateral_sec',num2str(sec),'CO...'));
                        save(fN,strcat('mk',num2str(monkey),'lateral_sec',num2str(sec),'CO'));
                        fprintf('DONE! \n');
                    catch
                        fprintf(strcat('NO CO image for mk373lateral sec',num2str(sec),' \n'));
                    end
                end
            end
        end    

        if monkey == 374 || monkey == 373
            fprintf(strcat('Done with monkey',num2str(monkey),loc,'!\n'));
        else
            fprintf(strcat('Done with monkey',num2str(monkey),'!\n'));
        end
    
    end
end
clear avgIbDist barWidth BC BCImage bins blob cBC COImage fileDate fileName ...
fN hemi ibDist injections injections373 la label labelImage labels loc M ...
locations monkey monkeys numBlobs rBC rotAngle sec sections sections356 ...
sections365 sections373Lat sections373Med sections374Lat sections374Med ...
mk356_lines mk373medial_lines mk365_lines mk373lateral_lines mk374lateral_lines ...
mk374medial_lines
%% Plot heat maps based on chosen blobs

% % sec4img = imread(strcat('MK356RHcrop sec',num2str(sec),' label-250umRollBallSubtract.tif'));
% heatMapAll = zeros(size(mk374medial_sec27TDT(1).heatMap,1),size(mk374medial_sec27TDT(1).heatMap,2));
% for blob = [3,9,17,19,20,23,25,26,29,32,33,35,36,37,38,42,44,45,46,50,54,55,56,57,59,62,63,69,73,77,78,79,82,83,84,87,90,92,94,97,101,102,104,106,108,110,111,118,119,125,127,133,136,137,144,151] %1:53 %numBlobs
%     mk374medial_sec27TDT(blob).heatMapNorm = mk374medial_sec27TDT(blob).heatMap/max(max(mk374medial_sec27TDT(blob).heatMap));
% %     mk356heatmap_sec4(blob).heatMapFiltered10 = medfilt2(mk356heatmap_sec4(blob).heatMap,[10 10]);
%     heatMapAll = heatMapAll + mk374medial_sec27TDT(blob).heatMap;
% end
% 
% heatMapAll2 = heatMapAll/max(max(heatMapAll));


% figure();
% imshow(heatMap);
