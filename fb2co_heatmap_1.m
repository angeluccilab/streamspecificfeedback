% Get the line lengths to both blob borders and interblob centerlines
% This code does not need to be run again. Just use the Lines_v4.mat file
% that is already saved for each case
% File: Line_v4
clear all; close all;

path( path, fullfile(pwd,'scriptFuncs') )
monkeys = [356 365 373 374 374];
flag = 0;

for mk = 1:length(monkeys)
    monkey = monkeys(mk);
    runDate = '2019-02-19';
    fileDate = '2019-02-19'; %for good interblob and blob lines
    rotAngle = 20; %angle to rotate the lines
    bins = 20;  %number of bins for normalization of line lengths
    barWidth = 200; %microns
    hemi = 'RH'; %left or right hemisphere
    
    if monkey > 372
        if flag < 2
            loc = 'medial';
            flag = flag + 1;
        else
            loc = 'lateral';
        end
        fileLoc = strcat('/MK',num2str(monkey),hemi,'/',loc,'/Mat Files/MK',num2str(monkey),loc,'_',fileDate,'_');
        load(fullfile(pwd,strcat(fileLoc,'GoodLines2BlobBorder.mat')),strcat('GoodLines2Blob_mk',num2str(monkey),loc));
        load(fullfile(pwd,strcat(fileLoc,'GoodLines2InterblobCenter.mat')),strcat('GoodLines2Interblob_mk',num2str(monkey),loc));
    else
        fileLoc = strcat('/MK',num2str(monkey),hemi,'/','/Mat Files/MK',num2str(monkey),hemi,'_',fileDate,'_');
        load(fullfile(pwd,strcat(fileLoc,'GoodLines2BlobBorder.mat')),strcat('GoodLines2Blob_mk',num2str(monkey)));
        load(fullfile(pwd,strcat(fileLoc,'GoodLines2InterblobCenter.mat')),strcat('GoodLines2Interblob_mk',num2str(monkey)));
    end

    if monkey == 356

        numBlobs = size(GoodLines2Interblob_mk356,2);
        for blob = 1:numBlobs
            mk356_lines(blob).interblobLineLengths = GoodLines2Interblob_mk356(blob).lineLengths;
            mk356_lines(blob).blobLineLengths = GoodLines2Blob_mk356(blob).lineLengths;
        end

        for blob = 1:numBlobs
            mk356_lines(blob).avgIbLengths = mean(mk356_lines(blob).interblobLineLengths);
            mk356_lines(blob).avgBlobLengths = mean(mk356_lines(blob).blobLineLengths);
        end

    elseif monkey == 365

        numBlobs = size(GoodLines2Interblob_mk365,2);    
        for blob = 1:numBlobs
            mk365_lines(blob).interblobLineLengths = GoodLines2Interblob_mk365(blob).lineLengths;
            mk365_lines(blob).blobLineLengths = GoodLines2Blob_mk365(blob).lineLengths;
        end

        for blob = 1:numBlobs
            mk365_lines(blob).avgIbLengths = mean(mk365_lines(blob).interblobLineLengths);
            mk365_lines(blob).avgBlobLengths = mean(mk365_lines(blob).blobLineLengths);
        end

    elseif monkey == 373
        if strcmp(loc,'lateral')
            numBlobs = size(GoodLines2Interblob_mk373lateral,2);
            for blob = 1:numBlobs
                mk373lateral_lines(blob).interblobLineLengths = GoodLines2Interblob_mk373lateral(blob).lineLengths;
                mk373lateral_lines(blob).blobLineLengths = GoodLines2Blob_mk373lateral(blob).lineLengths;
            end

            for blob = 1:numBlobs
                mk373lateral_lines(blob).avgIbLengths = mean(mk373lateral_lines(blob).interblobLineLengths);
                mk373lateral_lines(blob).avgBlobLengths = mean(mk373lateral_lines(blob).blobLineLengths);
            end

        elseif strcmp(loc,'medial')
            numBlobs = size(GoodLines2Interblob_mk373medial,2);
            for blob = 1:numBlobs
                mk373medial_lines(blob).interblobLineLengths = GoodLines2Interblob_mk373medial(blob).lineLengths;
                mk373medial_lines(blob).blobLineLengths = GoodLines2Blob_mk373medial(blob).lineLengths;
            end

            for blob = 1:numBlobs
                mk373medial_lines(blob).avgIbLengths = mean(mk373medial_lines(blob).interblobLineLengths);
                mk373medial_lines(blob).avgBlobLengths = mean(mk373medial_lines(blob).blobLineLengths);
            end

        end

    elseif monkey == 374
        if strcmp(loc,'lateral')
            numBlobs = size(GoodLines2Interblob_mk374lateral,2);
            for blob = 1:numBlobs
                mk374lateral_lines(blob).interblobLineLengths = GoodLines2Interblob_mk374lateral(blob).lineLengths;
                mk374lateral_lines(blob).blobLineLengths = GoodLines2Blob_mk374lateral(blob).lineLengths;
            end

            for blob = 1:numBlobs
                mk374lateral_lines(blob).avgIbLengths = mean(mk374lateral_lines(blob).interblobLineLengths);
                mk374lateral_lines(blob).avgBlobLengths = mean(mk374lateral_lines(blob).blobLineLengths);
            end

        elseif strcmp(loc,'medial')
            numBlobs = size(GoodLines2Interblob_mk374medial,2);
            for blob = 1:numBlobs
                mk374medial_lines(blob).interblobLineLengths = GoodLines2Interblob_mk374medial(blob).lineLengths;
                mk374medial_lines(blob).blobLineLengths = GoodLines2Blob_mk374medial(blob).lineLengths;
            end

            for blob = 1:numBlobs
                mk374medial_lines(blob).avgIbLengths = mean(mk374medial_lines(blob).interblobLineLengths);
                mk374medial_lines(blob).avgBlobLengths = mean(mk374medial_lines(blob).blobLineLengths);
            end

        end
    end


    if monkey == 373 || monkey == 374
        fileName = strcat('/MK',num2str(monkey),hemi,'/',loc,'/Mat Files/MK',num2str(monkey),hemi,loc,'_',runDate,'_');
        path( fullfile(pwd,strcat('MK',num2str(monkey),hemi),'/',loc,'/Mat Files'), path );
        fN = fullfile(pwd,strcat(fileName,'Lines_v4.mat'));
        save(fN,strcat('mk',num2str(monkey),loc,'_lines'));
    else
        fileName = strcat('/MK',num2str(monkey),hemi,'/Mat Files/MK',num2str(monkey),hemi,'_',runDate,'_');
        path( fullfile(pwd,strcat('MK',num2str(monkey),hemi,'/Mat Files')), path );
        fN = fullfile(pwd,strcat(fileName,'Lines_v4.mat'));
        save(fN,strcat('mk',num2str(monkey),'_lines'));
    end

end