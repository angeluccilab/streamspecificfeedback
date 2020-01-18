%% Extract data from file
clear all; close all;

monkeys = [356 365 373 374 374];
flag = 0;

% for mk = 1:length(monkeys)
% %variables
% runDate = '2019-02-15';
% rotAngle = 20; %angle to rotate the lines
% bins = 20;  %number of bins for normalization of line lengths
% barWidth = 200; %microns
% monkey = monkeys(mk); %monkey number
% RHLH = 'RH'; %left or right hemisphere
% 
% if monkey > 372
%     if flag < 2
%         loc = 'medial';
%         flag = flag + 1;
%     else
%         loc = 'lateral';
%     end
%     fileName = strcat('/MK',num2str(monkey),RHLH,'/',loc,'/Mat Files/MK',num2str(monkey),loc,'_',runDate,'_');
% 
%     path( path, fullfile(pwd,'scriptFuncs') );
%     path( fullfile(pwd,strcat('MK',num2str(monkey),RHLH),loc) , path)
% else
%     fileName = strcat('/MK',num2str(monkey),RHLH,'/','/Mat Files/MK',num2str(monkey),RHLH,'_',runDate,'_');
% 
%     path( path, fullfile(pwd,'scriptFuncs') );
%     path( path, fullfile(pwd,strcat('MK',num2str(monkey),RHLH)) );
% end
% 
% %read in images of blob centers, interblob centerlines, and blob outlines
% if monkey == 356
%     BCImage = imread('MK356RHcrop ALL Blob Centers Blue Single Pixel ver2.tif');
%     BOImage = imread('MK356RHcrop Blob Outlines Yellow.tif');
%     IBCOImage = imread('MK356RHcrop Interblob Centers Outlines Green VER2.tif');
% elseif monkey == 365
%     BCImage = imread('MK365RH blob centers blue.tif');
%     BOImage = imread('MK365RH blob outlines yellow.tif');
%     IBCOImage = imread('MK365RH interblob outlines green VER2.tif');
% elseif monkey == 373
%     BCImage = imread('_0013_Blob centers blue.tif');
%     BOImage = imread('_0014_Blob outlines yellow.tif');
%     IBCOImage = imread('_0012_Interblob centers green VER2.tif');
% elseif monkey == 374
%     if strcmp(loc,'medial')
%         BCImage = imread('MK374RH V1 MEDIAL CROP LARGE Blob Centers.tif');
%         BOImage = imread('MK374RH V1 MEDIAL CROP LARGE Blob Outlines.tif');
%         IBCOImage = imread('MK374RH V1 MEDIAL CROP LARGE Interblob Outlines VER2.tif');
%     elseif strcmp(loc, 'lateral')
%         BCImage = imread('MK374RH V1 stack LATERAL CROP LARGE Blob Centers.tif');
%         BOImage = imread('MK374RH V1 stack LATERAL CROP LARGE Blob Outlines.tif');
%         IBCOImage = imread('MK374RH V1 stack LATERAL CROP LARGE Interblob Outlines VER2.tif');
%     end
% end
% %extract the pixel conversions data
% convs = xlsread('FBtoCO um to pixel conversions.xlsx');
% 
% %get pixel to micron conversions for each monkey
% if monkey == 356
%     px2um = convs(1,1);
% elseif monkey == 365
%     px2um = convs(1,2);
% elseif monkey == 373
%     px2um = convs(1,3);
% elseif monkey == 374
%     px2um = convs(1,4);
% else
%     px2um = 1;
% end
% 
% %turn color images into b/w images for individual data sets
% tic
% fprintf('find blobs, borders, and centers:   ');
% IBCO = color2bw(IBCOImage,'green');
% BO = color2bw(BOImage,'yellow');
% BC = color2bw(BCImage,'blue');
% toc
% 
% figure();imshow(BO); title('Blob Outlines');
% figure();imshow(IBCO); title('Interblob borders');
% 
% 
% 
% %% Get the coordinates of the blob centers
% [rBC, cBC] = find(BC);
% %% Get lines to blob borders at every 'rotAngle' degrees
% % This is done by rotating the image and extending a line vertically from 
% % the blob center until it hits an edge or line
% if monkey > 372
%     tic
%     fprintf('Get Lines to Blob borders: ');
%     eval(['Lines2Blob_mk' num2str(monkey) loc ' = getLines(BO, rBC, cBC, rotAngle);']);
%     toc
%     
%     fN = fullfile(pwd,strcat(fileName,'Lines2BlobBorder.mat'));
%     save(fN,strcat('Lines2Blob_mk',num2str(monkey),loc));
% else
%     tic
%     fprintf('Get Lines to Blob borders: ');
%     eval(['Lines2Blob_mk' num2str(monkey) ' = getLines(BO, rBC, cBC, rotAngle);']);
%     toc
%     
%     fN = fullfile(pwd,strcat(fileName,'Lines2BlobBorder.mat'));
%     save(fN,strcat('Lines2Blob_mk',num2str(monkey)));
% end
% 
% 
% end

%%

for mk = 1:length(monkeys)
    runDate = '2019-02-19';
    fileDate = '2019-01-16'; %for interblob lines
    fileDate2 = '2019-02-15'; %for blob lines
    rotAngle = 20; %angle to rotate the lines
    bins = 20;  %number of bins for normalization of line lengths
    barWidth = 200; %microns
    monkey = monkeys(mk); %monkey number
    RHLH = 'RH'; %left or right hemisphere
    
    if monkey > 372
        if flag < 2
            loc = 'medial';
            flag = flag + 1;
        else
            loc = 'lateral';
        end
        fileName = strcat('/MK',num2str(monkey),RHLH,'/',loc,'/Mat Files/MK',num2str(monkey),loc,'_',runDate,'_');
        
        fileLoc = strcat('/MK',num2str(monkey),RHLH,'/',loc,'/Mat Files/MK',num2str(monkey),loc,'_',fileDate,'_');
        load(fullfile(pwd,strcat(fileLoc,'GoodLines.mat')),strcat('Data5_GoodLines2Interblob_mk',num2str(monkey),loc));
        eval(['M2 = Data5_GoodLines2Interblob_mk' num2str(monkey) loc ';']);
        
        fileLoc = strcat('/MK',num2str(monkey),RHLH,'/',loc,'/Mat Files/MK',num2str(monkey),loc,'_',fileDate2,'_');
        load(fullfile(pwd,strcat(fileLoc,'Lines2BlobBorder.mat')));
        eval(['M1 = Lines2Blob_mk' num2str(monkey) loc ';']);
        
        path( path, fullfile(pwd,'scriptFuncs') );
        path( fullfile(pwd,strcat('MK',num2str(monkey),RHLH),loc) , path)        
        
    else
        fileName = strcat('/MK',num2str(monkey),RHLH,'/','/Mat Files/MK',num2str(monkey),RHLH,'_',runDate,'_');
        
        fileLoc = strcat('/MK',num2str(monkey),RHLH,'/','/Mat Files/MK',num2str(monkey),RHLH,'_',fileDate,'_');
        load(fullfile(pwd,strcat(fileLoc,'GoodLines.mat')),strcat('GoodLines2Interblob_mk',num2str(monkey)));
        eval(['M2 = GoodLines2Interblob_mk' num2str(monkey) ';']);
        
        fileLoc = strcat('/MK',num2str(monkey),RHLH,'/','/Mat Files/MK',num2str(monkey),RHLH,'_',fileDate2,'_');
        load(fullfile(pwd,strcat(fileLoc,'Lines2BlobBorder.mat')));
        eval(['M1 = Lines2Blob_mk' num2str(monkey) ';']);
        
        path( path, fullfile(pwd,'scriptFuncs') );
        path( path, fullfile(pwd,strcat('MK',num2str(monkey),RHLH)) );
    end
    
    convs = xlsread('FBtoCO um to pixel conversions.xlsx');

    %get pixel to micron conversions for each monkey
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
    
    if monkey > 372
        tic
        fprintf(strcat('Get Good Lines to Blob borders (mk', num2str(monkey),'): '));
        eval(['GoodLines2Blob_mk' num2str(monkey) loc ' = selectGoodLines(M1,M2,px2um,1);']);
        toc

        fN = fullfile(pwd,strcat(fileName,'GoodLines2BlobBorder.mat'));
        save(fN,strcat('GoodLines2Blob_mk',num2str(monkey),loc));
        
        fN = fullfile(pwd,strcat(fileName,'GoodLines2InterblobCenter.mat'));
        eval(['GoodLines2Interblob_mk' num2str(monkey) loc ' = Data5_GoodLines2Interblob_mk' num2str(monkey) loc ';']);
        save(fN,strcat('GoodLines2Interblob_mk',num2str(monkey),loc));
    else
        tic
        fprintf(strcat('Get Good Lines to Blob borders (mk', num2str(monkey),'): '));
        eval(['GoodLines2Blob_mk' num2str(monkey) ' = selectGoodLines(M1,M2,px2um,1);']);
        toc

        fN = fullfile(pwd,strcat(fileName,'GoodLines2BlobBorder.mat'));
        save(fN,strcat('GoodLines2Blob_mk',num2str(monkey)));
        
        fN = fullfile(pwd,strcat(fileName,'GoodLines2InterblobCenter.mat'));
%         eval(['GoodLines2Interlob_mk' num2str(monkey) ' = Data5_GoodLines2Interblob_mk' num2str(monkey) ';']);
        save(fN,strcat('GoodLines2Interblob_mk',num2str(monkey)));
    end
    
end
