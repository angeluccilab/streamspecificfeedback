% Laminar case with bins instead of layers
% Bin number is defined arbitrarily, makes it easier to normalize
% This method uses a centerline of the label instead of a box
% This version (ver6) is the exact same as ver5c except it stores information 
% differently...more informative and detailed

clear all; close all;
runDate = '2020-01-06';
monkey = 359;
hemi = 'LH';
sections = [5 15 25 35 43 45 53 64 65 74 94 106 118];
locations = {'LAT' 'MED'};
labels = {'GFP' 'TDT'};
labelMidlines = {'GFP Lat' 'GFP Med' 'TDT Lat' 'TDT Med'};
isThresholded = 0;
isWithCO = 0;
wantSave = 1;
wantEPS = 0;
wantSeeIndividualPlots = 0;
isAllCases = 1;

path( fullfile(pwd,'scriptFuncs'), path );
path( fullfile(pwd,strcat('MK',num2str(monkey),hemi)), path );

% read images and rotate them 
for sec = 1:size(sections,2)
    sectionImages(sec).section = sections(sec);   
    layers2(sec).section = sections(sec);
    ind = 1;
    for lm = 1:size(labelMidlines,2)
        try
            eval(['temp = imread(''MK359LH Sec ' num2str(sections(sec)) ' ' labelMidlines{lm} ' CL.tif'');']);
            eval(['temp2 = imread(''MK359LH Sec ' num2str(sections(sec)) ' ' labelMidlines{lm} ' layers.tif'');']);
        catch
            continue
        end
        
        tempLM = imcomplement(im2bw(temp,0.8));       
        centerlines(ind).InjectionSite = labelMidlines{lm};
        centerlines(ind).LabelImage = tempLM;
        layers2(sec).LabelFieldCenterlines = centerlines;
        
        lyrs = logical(imcomplement(rgb2gray(temp2)));
        layeroutlines(ind).InjectionSite = labelMidlines{lm};
        layeroutlines(ind).LayerOutlineImage = lyrs;
        layers2(sec).LayerOutlines = layeroutlines;

        [midlineX midlineY] = find(tempLM);
        rotAngle = abs(atand((midlineY(1) - midlineY(end))/(midlineX(1) - midlineX(end))));
        midInd = round(size(midlineY,1)/2);
        rotAnchor = [midlineY(midInd) midlineX(midInd)];
        rotatedData(ind).InjectionSite = labelMidlines{lm};
        rotatedData(ind).RotationAngle = rotAngle;
        rotatedData(ind).RotationAnchorPoint = rotAnchor;
        layers2(sec).RotationInfo = rotatedData;
        
        tempLM2 = logical(tempLM);
        tempLM2 = rotateAround2(tempLM2,rotAnchor(2),rotAnchor(1),-rotAngle);
        [lmR lmC] = find(tempLM2);
        rotatedcenterlines(ind).InjectionSite = labelMidlines{lm};
        rotatedcenterlines(ind).RotatedImage = tempLM2;
        rotatedcenterlines(ind).RotatedCenterlineCoords = [lmR lmC];
        layers2(sec).RotatedLabelFieldCenterlines = rotatedcenterlines;
        
        lyrs2 = rotateAround2(lyrs,rotAnchor(2),rotAnchor(1),-rotAngle);
        [lyrR lyrC] = find(lyrs2);
        rotatedoutlines(ind).InjectionSite = labelMidlines{lm};
        rotatedoutlines(ind).RotatedImage = lyrs2;
        rotatedoutlines(ind).RotatedLayersCoords = [lyrR lyrC];
        layers2(sec).RotatedLayerOutlines = rotatedoutlines;
        
        ind = ind + 1;
    end
    clear centerlines layeroutlines rotatedData rotatedcenterlines rotatedoutlines
    
    ind = 1;
    for label = 1:size(labels,2)  
        try
            if isWithCO 
                eval(['temp = imread(''MK359LH Sec ' num2str(sections(sec)) ' ' labels{label} '.tif'');']);                
            else
                eval(['temp = imread(''MK359LH Sec' num2str(sections(sec)) ' ' labels{label} '-registered.tif'');']);
            end
            sectionImages(sec).labelImages{ind,1} = labels{label};
            sectionImages(sec).labelImages{ind,2} = temp;
            ind = ind + 1;
        catch
            continue
        end
        
    end

end
if isThresholded
    sectionImages = thresholdImages(sectionImages);
    if isWithCO
        fprintf('Analysis: CO images, thresholded\n');
    else
        fprintf('Analysis: Non-CO images, thresholded\n');
    end
else
    if isWithCO
        fprintf('Analysis: CO images, NOT thresholded\n');
    else
        fprintf('Analysis: Non-CO images, NOT thresholded\n');
    end
end

% get scales for the images
[scales,labelFields] = xlsread('MK359LH section scales.xlsx');
for sec = 1:size(sectionImages,2)
    for sc = 1:size(scales,1)
        if scales(sc,1) == sectionImages(sec).section
            sectionImages(sec).scale = scales(sc,3);
        end
    end
end

clear tempLM2 tempLM temp2 sec sc rotAngle rotAnchor midlineX midlineY ...
    midInd lyrs2 lyrs lyrR lyrC lmR lmC lm label ind temp labelFields

%% Find the layer boundaries
binWidth = 500; %um

for sec = 1:size(sections,2)
    for img = 1:size(sectionImages,2)
        if sectionImages(img).section == layers2(sec).section
            um2px = sectionImages(img).scale;
            break
        end
    end
    binWidthHalf = round(um2px*binWidth/2); %px
    
    numInjs = size(layers2(sec).RotatedLabelFieldCenterlines,2);
    for inj = 1:numInjs
        tempLB(inj).InjectionSite = layers2(sec).RotatedLayerOutlines(inj).InjectionSite;
        
        currLyrs = layers2(sec).RotatedLayerOutlines(inj).RotatedImage;
        currLyrs = medfilt2(currLyrs);       
        x1 = round(min(layers2(sec).RotatedLabelFieldCenterlines(inj).RotatedCenterlineCoords(:,2)));
        x2 = round(max(layers2(sec).RotatedLabelFieldCenterlines(inj).RotatedCenterlineCoords(:,2)));
        xStart = round(mean(x1,x2)) - binWidthHalf;
        xEnd = round(mean(x1,x2)) + binWidthHalf;
        yStart = round(min(layers2(sec).RotatedLabelFieldCenterlines(inj).RotatedCenterlineCoords(:,1)));
        yEnd = round(max(layers2(sec).RotatedLabelFieldCenterlines(inj).RotatedCenterlineCoords(:,1)));
        
        col = 1;
        for x = xStart:xEnd
            row = 1; flag = 0; tempBrdr = 0;
            for y = yStart:yEnd
                if currLyrs(y,x) == 1
                    if y > (tempBrdr + 11)
                        tempBrdr = y;
                        lyrBrdr(row,col) = y;
                        flag = 1;
                    end
                    if flag == 1
                        row = row + 1;
                        flag = 0;
                    end          
                end
            end
            col = col + 1;
        end
        clear x y flag row tempBrdr;

        [r c] = find(lyrBrdr == 0);
        for cI = 1:size(r,1)
            for rI = 1:size(c,1)
                rr = r(rI); cc = c(cI);
                newNum = lyrBrdr(rr,cc-1);
                lyrBrdr(rr,cc) = newNum;
            end
        end
        clear newNum rr cc rI cI;
           
        % Find layer 1 and 6 borders, average out layers in between
        for lyr = 1:size(lyrBrdr,1) 
            if lyr == 1
                brdr(lyr) = max(lyrBrdr(lyr,:));
            elseif lyr == size(lyrBrdr,1)
                brdr(lyr) = min(lyrBrdr(lyr,:));
            else
                brdr(lyr) = round(mean(lyrBrdr(lyr,:)));
            end
        end
        clear lyr;

        tempLB(inj).BorderRawYCoords = lyrBrdr;
        tempLB(inj).BorderDefinedYCoords = brdr;
        tempLB(inj).BorderDefinedXCoords = round(mean(x1,x2));
        layers2(sec).LayerBorders = tempLB;
        brdr = []; lyrBrdr = [];       
    end
    clear tempLB;
end      
clear yStart yEnd xStart xEnd x1 x2 um2px sec r lyrBrdr inj img c col ...
    currLyrs brdr binWidthHalf

%% Get coordinates of bin fills
numBins = 100; %unit: pixels
% binWidth = 500; %unit: pixels

for sec = 1:size(sections,2)
    
    for img = 1:size(sectionImages,2)
        if sectionImages(img).section == layers2(sec).section
            um2px = sectionImages(img).scale;
            break
        end
    end
    binWidthHalf = um2px*binWidth/2; %px
    clear img um2px;
    
    numInjs = size(layers2(sec).LabelFieldCenterlines,2);
    for inj = 1:numInjs

        tempBinCoords(inj).InjectionSite = layers2(sec).RotatedLabelFieldCenterlines(inj).InjectionSite;
        xStart = round(layers2(sec).LayerBorders(inj).BorderDefinedXCoords - binWidthHalf);
        xEnd = round(layers2(sec).LayerBorders(inj).BorderDefinedXCoords + binWidthHalf);
        yStart = round(layers2(sec).LayerBorders(inj).BorderDefinedYCoords(1));
        yEnd = round(layers2(sec).LayerBorders(inj).BorderDefinedYCoords(end));
        
        binHt = floor((yEnd - yStart)/numBins);
        binHtExtra = mod((yEnd - yStart),numBins);
        binHt2 = binHt + 1;
        firstBins = numBins - binHtExtra;

        bin = 1; ind2 = 2;
        for yPix = yStart:binHt:yEnd
            ind = 1; rows = []; cols = [];
            for yP = yPix:yPix+binHt   
                if yP > yEnd
                    break
                else
                    for x = xStart:xEnd
                        rows(ind) = yP;
                        cols(ind) = x;
                        ind = ind + 1;
                    end
                end
            end
            eval(['tempBinCoords(inj).Bin' num2str(bin) 'Coords = [rows'' cols''];']);
            bin = bin + 1;
            if bin > firstBins
                break
            end
        end
        yStart = yP+1; binHt = binHt2;
        for yPix = yStart:binHt:yEnd
            ind = 1; rows = []; cols = [];
            for yP = yPix:yPix+binHt2  
                if yP > yEnd
                    break
                else
                    for x = xStart:xEnd
                        rows(ind) = yP;
                        cols(ind) = x;
                        ind = ind + 1;
                    end
                end
            end
            eval(['tempBinCoords(inj).Bin' num2str(bin) 'Coords = [rows'' cols''];']);
            bin = bin + 1;
        end
    end
    layers2(sec).BinCoordinates = tempBinCoords;
    clear tempBinCoords
end
clear yStart yPix yP yEnd xStart xEnd x sec rows inj ind2 ind firstBins ...
    cols binWidthHalf binHtExtra binHt2 binHt bin


%% Determine layer borders
for sec = 1:size(sections,2)

    numInjs = size(layers2(sec).BinCoordinates,2);
    for inj = 1:numInjs
        % Borders by bin
        ind = 1; lyrBrdrBin = []; 
        for bin = 1:numBins
            eval(['binYCoords = layers2(sec).BinCoordinates(inj).Bin' num2str(bin) 'Coords(:,1);']); 
            if ismember(layers2(sec).LayerBorders(inj).BorderDefinedYCoords(ind),binYCoords)
                lyrBrdrBin(ind) = bin;
                ind = ind + 1;
                if ind > size(layers2(sec).LayerBorders(inj).BorderDefinedYCoords,2)
                    break
                end
            end
        end
        layers2(sec).LayerBorders(inj).LayerBorderBins = lyrBrdrBin;
        
        %Border by pixel        
        lyrBorders = layers2(sec).LayerBorders(inj).BorderDefinedYCoords;
        brdrFraction = [];
        for brdr = 1:size(lyrBorders,2)-1
            brdrFraction(brdr) = (lyrBorders(brdr+1)-lyrBorders(1))/ ...
                (lyrBorders(end)-lyrBorders(1));
        end
        layers2(sec).LayerBorders(inj).LayerBorderAsFractionOfDepth = brdrFraction;
    end
    
end
clear sec lyrBrdrBin lyrBorders inj ind brdrFraction brdr bin


%% get intensities
for sec = 1:size(sections,2)
    for label = 1:size(sectionImages(sec).labelImages,1)
        if strcmp(sectionImages(sec).labelImages{label,1},'GFP')
            gfpImg = sectionImages(sec).labelImages{label,2};
        elseif strcmp(sectionImages(sec).labelImages{label,1},'TDT')
            tdtImg = sectionImages(sec).labelImages{label,2};
        end
    end
    for inj = 1:size(layers2(sec).BinCoordinates,2) 
        if strcmp(layers2(sec).BinCoordinates(inj).InjectionSite,'GFP Med') || strcmp(layers2(sec).BinCoordinates(inj).InjectionSite,'GFP Lat')                               
            tempImg = rotateAround2(gfpImg,layers2(sec).RotationInfo(inj).RotationAnchorPoint(2),...
                layers2(sec).RotationInfo(inj).RotationAnchorPoint(1),-layers2(sec).RotationInfo(inj).RotationAngle);
            sectionImages(sec).rotatedImage{inj,1} = tempImg;
            sectionImages(sec).rotatedImage{inj,2} = tempImg;
        elseif strcmp(layers2(sec).BinCoordinates(inj).InjectionSite,'TDT Med') || strcmp(layers2(sec).BinCoordinates(inj).InjectionSite,'TDT Lat')
            tempImg = rotateAround2(tdtImg,layers2(sec).RotationInfo(inj).RotationAnchorPoint(2),...
                layers2(sec).RotationInfo(inj).RotationAnchorPoint(1),-layers2(sec).RotationInfo(inj).RotationAngle);
            sectionImages(sec).rotatedImage{inj,3} = tempImg;
            sectionImages(sec).rotatedImage{inj,4} = tempImg;
        end
        intensities(inj).InjectionSite = layers2(sec).BinCoordinates(inj).InjectionSite;
        for bin = 1:numBins 
            r = []; c = [];
            eval(['binCheck = layers2(sec).BinCoordinates(inj).Bin' num2str(bin) 'Coords;']);
            if ~isempty(binCheck)
                eval(['r = layers2(sec).BinCoordinates(inj).Bin' num2str(bin) 'Coords(:,1);']);
                eval(['c = layers2(sec).BinCoordinates(inj).Bin' num2str(bin) 'Coords(:,2);']);
            else
                continue
            end
            intens = [];
            for j = 1:size(r)
                intens(j) = tempImg(r(j),c(j));
            end               
            tempIntens(bin) = mean(intens);
        end
        intensities(inj).RawBinIntensities = tempIntens;
        clear tempIntens
    end
    layers2(sec).Intensities = intensities;
    clear intensities
end

% keyboard

%% Plot individual results
lyrsText = {'L1' 'L2/3' 'L4A' 'L4B' 'L4C' 'L5' 'L6'};

for sec = 1:size(sections,2)
    for inj = 1:size(layers2(sec).Intensities,2)
       intensity = layers2(sec).Intensities(inj).RawBinIntensities;
       intensity = (intensity - min(intensity))/(max(intensity)-min(intensity));
       layers2(sec).Intensities(inj).NormalizedBinIntensities = intensity;
       
       
       yPt = max(intensity);
       verts = linspace(0,yPt);
       yPt = yPt - 0.1;
       textY = [yPt yPt yPt yPt yPt yPt yPt];
       fig = figure; plot(intensity); hold on;
       
       for k = 2:(size(layers2(sec).LayerBorders(inj).LayerBorderBins,2)-1)
%            border = ones(size(verts,2),1)*(layers2(sec).LayerBorders(inj).LayerBorderBins(k));
%            plot(border,verts,'--c'); hold on;

           currlyrbrdr = layers2(sec).LayerBorders(inj).LayerBorderAsFractionOfDepth(k);
           prevlyrbrdr = layers2(sec).LayerBorders(inj).LayerBorderAsFractionOfDepth(k-1);
           brdrFrac = prevlyrbrdr*size(intensity,2);
           border2 = ones(size(verts,2),1)*brdrFrac;
           brdr2 = (((currlyrbrdr - prevlyrbrdr)/2) + prevlyrbrdr)*size(intensity,2);
           plot(border2,verts,'--m'); hold on;
           
           text(brdr2,textY(k),lyrsText{k},'HorizontalAlignment','center');
           
       end
       brdr2 = (((layers2(sec).LayerBorders(inj).LayerBorderAsFractionOfDepth(1))/2))*size(intensity,2);
       text(brdr2,textY(1),lyrsText{1},'HorizontalAlignment','center');

       title(strcat('Sec ',num2str(layers2(sec).section),...
        '...',layers2(sec).Intensities(inj).InjectionSite));
       xlim([1 numBins]); set(gca,'Xticklabel',[]);
       if ~wantSeeIndividualPlots
           close(fig);
       end
       clear currlyrbrdr prevlyrbrdr
    end
end
clear brdr intens bin border2 brdr2 brdrFrac c gfpImg inj intensity j k ...
    label r sec tdtImg tempImg tempIntens yPt fig
% keyboard
if ~isAllCases
    fprintf('Not enough data for combined results\n');
    return
end
%% combine results (per label)
for inj = 1:4
    sections = [];
    if strcmp(labelMidlines{inj},'GFP Lat')
        sections = [15,25,35,45,65,74,94];
        sectionsMidsOnly = [35,45,65];
    elseif strcmp(labelMidlines{inj},'GFP Med')
        sections = [43,53,64,74,94,106,118];
        sectionsMidsOnly  = [64,74,94];
    elseif strcmp(labelMidlines{inj},'TDT Lat')
        sections = [5,15,25,35,65,74,94];
        sectionsMidsOnly  = [25,35,65];
    elseif strcmp(labelMidlines{inj},'TDT Med')
        sections = [25,35,43,53,74,94,106];
        sectionsMidsOnly  = [43,53,74];
    end  
    
    intens = [];
    k = 1;
    for sec = 1:size(layers2,2)
        if max(ismember(layers2(sec).section,sections))
            for inj2 = 1:size(layers2(sec).Intensities,2)
                if strcmp(layers2(sec).Intensities(inj2).InjectionSite,labelMidlines{inj})
                    intens(k,:) = layers2(sec).Intensities(inj2).NormalizedBinIntensities;
                    brdr(k,:) = layers2(sec).LayerBorders(inj2).LayerBorderBins;
                    brdrF(k,:) = layers2(sec).LayerBorders(inj2).LayerBorderAsFractionOfDepth;
                    k = k + 1;
                end
            end
        end
    end
    
    intensMids = [];
    k = 1;
    for sec = 1:size(layers2,2)
        if max(ismember(layers2(sec).section,sectionsMidsOnly))
            for inj2 = 1:size(layers2(sec).Intensities,2)
                if strcmp(layers2(sec).Intensities(inj2).InjectionSite,labelMidlines{inj})
                    intensMids(k,:) = layers2(sec).Intensities(inj2).NormalizedBinIntensities;
                    brdrMids(k,:) = layers2(sec).LayerBorders(inj2).LayerBorderBins;
                    brdrFMids(k,:) = layers2(sec).LayerBorders(inj2).LayerBorderAsFractionOfDepth;
                    k = k + 1;
                end
            end
        end
    end       

    layers2(1).PopulationData(inj).InjectionSite = labelMidlines{inj};
    layers2(1).PopulationData(inj).Sections = sections;
    layers2(1).PopulationData(inj).AllIntensities = intens;
    layers2(1).PopulationData(inj).CombinedIntensities = mean(intens);
    normIntens = (mean(intens) - min(mean(intens)))/(max(mean(intens))-min(mean(intens)));
    layers2(1).PopulationData(inj).NormalizedCombinedIntensities = normIntens;  %normalized intensities
    layers2(1).PopulationData(inj).StandardError = std(intens)./sqrt(size(intens,1)); %standard error
    layers2(1).PopulationData(inj).CombinedLayerBorderBins = mean(brdr);
    layers2(1).PopulationData(inj).CombinedLayerBorderAsFractionOfDepth = mean(brdrF);
    
    layers2(1).PopulationDataMiddleSections(inj).InjectionSite = labelMidlines{inj};
    layers2(1).PopulationDataMiddleSections(inj).Sections = sectionsMidsOnly;
    layers2(1).PopulationDataMiddleSections(inj).AllIntensities = intensMids;
    layers2(1).PopulationDataMiddleSections(inj).CombinedIntensities = mean(intensMids);
    normIntens = [];
    normIntens = (mean(intensMids) - min(mean(intensMids)))/(max(mean(intensMids))-min(mean(intensMids)));
    layers2(1).PopulationDataMiddleSections(inj).NormalizedCombinedIntensities = normIntens; %normalized intensities
    layers2(1).PopulationDataMiddleSections(inj).StandardError = std(intensMids)./sqrt(size(intensMids,1)); %standard error    
    layers2(1).PopulationDataMiddleSections(inj).CombinedLayerBorderBins = mean(brdrMids);
    layers2(1).PopulationDataMiddleSections(inj).CombinedLayerBorderAsFractionOfDepth = mean(brdrFMids);
end


%% Plot combined per label results
% for figWin = 1:2
    for inj = 1:4
        fig = figure; hold on;

        intensities = layers2(1).PopulationDataMiddleSections(inj).NormalizedCombinedIntensities;
        stdError = layers2(1).PopulationDataMiddleSections(inj).StandardError;
        gX= 1:size(intensities,2);
        fill([gX';flipud(gX')],[intensities' - stdError'; flipud(intensities' + stdError')],[0.8 1 1],'linestyle','none');
        plot(intensities,'-c');  

        intensities = layers2(1).PopulationData(inj).NormalizedCombinedIntensities;
        stdError = layers2(1).PopulationData(inj).StandardError;

        gX= 1:size(intensities,2);
        fill([gX';flipud(gX')],[intensities' - stdError';flipud(intensities' + stdError')],[0.8 0.8 1],'linestyle','none');
        plot(intensities,'-b');

    %     plot(layers2(1).PopulationData(inj).NormalizedCombinedIntensities,'-b','LineWidth',2); hold on;
    %     plot(layers2(1).PopulationDataMiddleSections(inj).NormalizedCombinedIntensities,'-c','LineWidth',2); hold on;
    % 
    %     errorbar((1:size(layers2(1).PopulationData(inj).NormalizedCombinedIntensities,2)),...
    %         layers2(1).PopulationData(inj).NormalizedCombinedIntensities,...
    %         layers2(1).PopulationData(inj).StandardError,'vertical','.b','linewidth',2);
    %     errorbar((1:size(layers2(1).PopulationDataMiddleSections(inj).NormalizedCombinedIntensities,2)),...
    %         layers2(1).PopulationDataMiddleSections(inj).NormalizedCombinedIntensities,...
    %         layers2(1).PopulationDataMiddleSections(inj).StandardError,'vertical','.c','linewidth',1)
    %     set(b,{'FaceColor'},{'c'});
        brdr = layers2(1).PopulationData(inj).CombinedLayerBorderBins;
        brdrF = layers2(1).PopulationData(inj).CombinedLayerBorderAsFractionOfDepth;
        brdrMids = layers2(1).PopulationDataMiddleSections(inj).CombinedLayerBorderBins;
        brdrFMids = layers2(1).PopulationDataMiddleSections(inj).CombinedLayerBorderAsFractionOfDepth;
        yPt = 1.1;   
        textY = [yPt yPt yPt-0.1 yPt yPt yPt yPt];
    %     keyboard
        for k = 2:(size(brdr,2)-1)
            border = ones(size(verts,2),1)*(brdr(k));
            brdrFrac = brdrF(k-1)*size(layers2(1).PopulationData(inj).NormalizedCombinedIntensities,2);
            border2 = ones(size(verts,2),1)*(brdrFrac);
            brdr2 = (((brdrF(k) - brdrF(k-1))/2)+brdrF(k-1))*size(layers2(1).PopulationData(inj).NormalizedCombinedIntensities,2);
            plot(border,verts,'--m'); hold on;
            plot(border2,verts,'-m'); hold on;
            
            border = ones(size(verts,2),1)*(brdrMids(k));
            brdrFrac = brdrFMids(k-1)*size(layers2(1).PopulationDataMiddleSections(inj).NormalizedCombinedIntensities,2);
            border2 = ones(size(verts,2),1)*(brdrFrac);
%             brdr2 = (((brdrFMids(k) - brdrFMids(k-1))/2)+brdrFMids(k-1))*size(layers2(1).PopulationDataMiddleSections(inj).NormalizedCombinedIntensities,2);
            plot(border,verts,'--g'); hold on;
            plot(border2,verts,'-g'); hold on;
            
            text(brdr2,textY(k),lyrsText{k},'HorizontalAlignment','center', ...
                'FontSize',17,'FontWeight','bold');  
        end
        brdr2 = (brdrF(1)/2)*size(layers2(1).PopulationData(inj).NormalizedCombinedIntensities,2);
        text(brdr2,textY(1),lyrsText{1},'HorizontalAlignment','center', ...
                'FontSize',17,'FontWeight','bold');

        title({'MK359LH Combined Results',labelMidlines{inj}});

        xlim([1 numBins]); ylim([0 1.2]);
        ylabel('Intensity');
        set(gca,'Xticklabel',[]) 

        if isThresholded
            if isWithCO
                filename = fullfile(pwd,strcat('/MK359LH/Results/Thresholded/CO-stained/MK359LH_CombinedResults_',...
                    labelMidlines{inj},'_Thresholded_withCO'));
            else
                filename = fullfile(pwd,strcat('/MK359LH/Results/Thresholded/PreCO-stained/MK359LH_CombinedResults_',...
                labelMidlines{inj},'_Thresholded_withoutCO'));
            end
        else
            if isWithCO
                filename = fullfile(pwd,strcat('/MK359LH/Results/Not Thresholded/CO-stained/MK359LH_CombinedResults_',...
                    labelMidlines{inj},'_NonThresholded_withCO'));
            else
                filename = fullfile(pwd,strcat('/MK359LH/Results/Not Thresholded/PreCO-stained/MK359LH_CombinedResults_',...
                labelMidlines{inj},'_Thresholded_withoutCO'));
            end
        end
        if wantSave
            print( fig, '-dtiff', filename ); 
            if wantEPS
               saveas(fig,filename,'epsc');
            end
        end
    end
%     keyboard
% end
clear alt altInd bin border2 brdr brdr2 brdrF brdrFrac c fig flag inj ...
    intens intensity j k label inj numBins plotStyle r sec  tdtImg tempImg ...
    tempIntens yPt
% keyboard
%% Combine results per stripe (all sections
 
gInd = 1; tInd = 1; gfpIntens = []; tdtIntens = []; gfpBrdr = []; tdtBrdr = [];
gfpIntensMids = []; tdtIntensMids = [];
for label = 1:4
%     keyboard
    if strcmp(layers2(1).PopulationData(label).InjectionSite,'GFP Lat') || ...
            strcmp(layers2(1).PopulationData(label).InjectionSite,'GFP Med')
        gfpIntens(gInd,:) = layers2(1).PopulationData(label).NormalizedCombinedIntensities;
        gfpIntensMids(gInd,:) = layers2(1).PopulationDataMiddleSections(label).NormalizedCombinedIntensities;
        gfpBrdr(gInd,:) = layers2(1).PopulationData(label).CombinedLayerBorderBins;
        gfpBrdrF(gInd,:) = layers2(1).PopulationData(label).CombinedLayerBorderAsFractionOfDepth;
        gfpBrdrMids(gInd,:) = layers2(1).PopulationDataMiddleSections(label).CombinedLayerBorderBins;
        gfpBrdrFMids(gInd,:) = layers2(1).PopulationDataMiddleSections(label).CombinedLayerBorderAsFractionOfDepth;
        gInd = gInd + 1;
    else
        tdtIntens(tInd,:) = layers2(1).PopulationData(label).NormalizedCombinedIntensities;
        tdtIntensMids(tInd,:) = layers2(1).PopulationDataMiddleSections(label).NormalizedCombinedIntensities;
        tdtBrdr(tInd,:) = layers2(1).PopulationData(label).CombinedLayerBorderBins;
        tdtBrdrF(tInd,:) = layers2(1).PopulationData(label).CombinedLayerBorderAsFractionOfDepth;
        tdtBrdrMids(tInd,:) = layers2(1).PopulationDataMiddleSections(label).CombinedLayerBorderBins;
        tdtBrdrFMids(tInd,:) = layers2(1).PopulationDataMiddleSections(label).CombinedLayerBorderAsFractionOfDepth;
        tInd = tInd + 1;
    end
end

gfpIntensity = mean(gfpIntens); gfpIntensity = (gfpIntensity - min(gfpIntensity))/(max(gfpIntensity) - min(gfpIntensity));
gfpIntensityMids = mean(gfpIntensMids); gfpIntensityMids = (gfpIntensityMids - min(gfpIntensityMids))/(max(gfpIntensityMids) - min(gfpIntensityMids));
tdtIntensity = mean(tdtIntens); tdtIntensity = (tdtIntensity - min(tdtIntensity))/(max(tdtIntensity) - min(tdtIntensity));
tdtIntensityMids = mean(tdtIntensMids); tdtIntensityMids = (tdtIntensityMids - min(tdtIntensityMids))/(max(tdtIntensityMids) - min(tdtIntensityMids));
gfpSEM = std(gfpIntens)./sqrt(size(gfpIntens,1)); gfpMidsSEM = std(gfpIntensMids)./sqrt(size(gfpIntensMids,1));
tdtSEM = std(tdtIntens)./sqrt(size(tdtIntens,1)); tdtMidsSEM = std(tdtIntensMids)./sqrt(size(tdtIntensMids,1));
gfpBorder = round(mean(gfpBrdr)); gfpBorder2 = (mean(gfpBrdrF));
tdtBorder = round(mean(tdtBrdr)); tdtBorder2 = (mean(tdtBrdrF));
gfpBorderMids = round(mean(gfpBrdrMids)); gfpBorder2Mids = (mean(gfpBrdrFMids));
tdtBorderMids = round(mean(tdtBrdrMids)); tdtBorder2Mids = (mean(tdtBrdrFMids));

yPt = 1.1;   
textY = [yPt yPt yPt-0.1 yPt yPt yPt yPt];

fig=figure; hold on;
gX= 1:size(gfpIntensityMids,2);
fill([gX';flipud(gX')],[gfpIntensityMids'-gfpMidsSEM';flipud(gfpIntensityMids'+gfpMidsSEM')], ...
        [0.8 1 1],'linestyle','none');

gX= 1:size(gfpIntensity,2);
fill([gX';flipud(gX')],[gfpIntensity'-gfpSEM';flipud(gfpIntensity'+gfpSEM')], ...
        [0.8 0.8 1],'linestyle','none');
plot(gfpIntensityMids,'-c');    
plot(gfpIntensity,'-b');

for k = 2:(size(gfpBrdr,2)-1)
    border = ones(size(verts,2),1)*(gfpBorder(k));
    brdrFrac = gfpBorder2(k-1)*size(gfpIntensity,2);
    border2 = round(ones(size(verts,2),1)*brdrFrac);
    plot(border,verts,'-m'); hold on;
%     plot(border2,verts,'-m'); hold on;
    
    border = ones(size(verts,2),1)*(gfpBorderMids(k));
    brdrFrac = gfpBorder2Mids(k-1)*size(gfpIntensityMids,2);
    border2 = round(ones(size(verts,2),1)*brdrFrac);
    plot(border,verts,'-g'); hold on;
%     plot(border2,verts,'-g'); hold on;
    
    brdr2 = ((gfpBorder2(k) - gfpBorder2(k-1))/2+gfpBorder2(k-1))* ...
        size(gfpIntensityMids,2);
    text(brdr2,textY(k),lyrsText{k},'HorizontalAlignment','center', ...
            'FontSize',17,'FontWeight','bold');
end
brdr2 = (gfpBorder2(1)/2)*size(gfpIntensityMids,2);
text(brdr2,textY(1),lyrsText{1},'HorizontalAlignment','center', ...
            'FontSize',17,'FontWeight','bold');
title({'MK359LH Pale Stripe Projections (GFP)'});
ylabel('Intensity'); ylim([0 1.2]);
set(gca,'Xticklabel',[]);

clear brdrFrac border2 brdr2 gfpBorder2Mids gfpBorderMids gfpBrdrFMids gfpBrdrMids

if isThresholded
    if isWithCO
        filename = fullfile(pwd,strcat('/MK359LH/Results/Thresholded/CO-stained/MK359LH_CombinedResults_Stripe PALE(GFP)_Thresholded'));
    else
        filename = fullfile(pwd,strcat('/MK359LH/Results/Thresholded/PreCO-stained/MK359LH_CombinedResults_Stripe PALE(GFP)_Thresholded'));
    end
else
    if isWithCO
        filename = fullfile(pwd,strcat('/MK359LH/Results/Not Thresholded/CO-stained/MK359LH_CombinedResults_Stripe PALE(GFP)_NonThresholded'));
    else
        filename = fullfile(pwd,strcat('/MK359LH/Results/Not Thresholded/PreCO-stained/MK359LH_CombinedResults_Stripe PALE(GFP)_NonThresholded'));
    end
end
if wantSave
    print(fig, '-dtiff', filename );
     if wantEPS
        saveas(fig,filename,'epsc');
     end
end

fig=figure; hold on;
tX= 1:size(tdtIntensityMids,2);
fill([tX';flipud(tX')],[tdtIntensityMids'-tdtMidsSEM';flipud(tdtIntensityMids'+tdtMidsSEM')], ...
        [0.8 1 1],'linestyle','none');
tX= 1:size(tdtIntensity,2);
fill([tX';flipud(tX')],[tdtIntensity'-tdtSEM';flipud(tdtIntensity'+tdtSEM')], ...
        [0.8 0.8 1],'linestyle','none');
plot(tdtIntensityMids,'-c');
plot(tdtIntensity,'-b');
for k = 2:(size(tdtBrdr,2)-1)
%         keyboard
    border = ones(size(verts,2),1)*(tdtBorder(k));
    brdrFrac = tdtBorder2(k-1)*size(tdtIntensity,2);
    border2 = round(ones(size(verts,2),1)*brdrFrac);
    plot(border,verts,'-m'); hold on;
%     plot(border2,verts,'-m'); hold on;
    
    border = ones(size(verts,2),1)*(tdtBorderMids(k));
    brdrFrac = tdtBorder2Mids(k-1)*size(tdtIntensityMids,2);
    border2 = round(ones(size(verts,2),1)*brdrFrac);
    plot(border,verts,'-g'); hold on;
%     plot(border2,verts,'-g'); hold on;
    
    brdr2 = ((tdtBorder2(k) - tdtBorder2(k-1))/2+tdtBorder2(k-1))* ...
        size(tdtIntensityMids,2);
    text(brdr2,textY(k),lyrsText{k},'HorizontalAlignment','center', ...
        'FontSize',17,'FontWeight','bold');
end

brdr2 = (tdtBorder2(1)/2)*size(tdtIntensity,2);
text(brdr2,textY(1),lyrsText{1},'HorizontalAlignment','center', ...
        'FontSize',17,'FontWeight','bold');
title({'MK359LH Thin Stripe Projections (TDT)'});
ylabel('Intensity'); ylim([0 1.2]);
set(gca,'Xticklabel',[]) 

% %create a new figure for legend
% cloneFig = figure;
% compCopy(fig,cloneFig);
% set(cloneFig,'DefaultFigureVisible','off');
% set(binsAll(k-1),'visible','on');
% set(binsMids(k-1),'visible','on');

clear brdrFrac border2 brdr2 tdtBorder2Mids tdtBorderMids tdtBrdrFMids tdtBrdrMids stdError

if isThresholded
    if isWithCO
        filename = fullfile(pwd,strcat('/MK359LH/Results/Thresholded/CO-stained/MK359LH_CombineResults_Stripe THIN(TDT)_Thresholded_withCO'));
        fN = fullfile(pwd,strcat('/MK359LH/Results/Thresholded/CO-stained/MK359LH_CombinedResults_Thresholded_withCO_',runDate,'.mat'));
    else
        filename = fullfile(pwd,strcat('/MK359LH/Results/Thresholded/preCO-stained/MK359LH_CombinedResults_Stripe THIN(TDT)_Thresholded_withoutCO'));
        fN = fullfile(pwd,strcat('/MK359LH/Results/Thresholded/preCO-stained/MK359LH_CombinedResults_Thresholded_withoutCO_',runDate,'.mat'));
    end
else
    if isWithCO
        filename = fullfile(pwd,strcat('/MK359LH/Results/Not Thresholded/CO-stained/MK359LH_CombinedResults_Stripe THIN(TDT)_NonThresholded_withCO'));
        fN = fullfile(pwd,strcat('/MK359LH/Results/Not Thresholded/CO-stained/MK359LH_CombinedResults_NonThresholded_withCO_',runDate,'.mat'));
    else
        filename = fullfile(pwd,strcat('/MK359LH/Results/Not Thresholded/preCO-stained/MK359LH_CombinedResults_Stripe THIN(TDT)_NonThresholded_withoutCO'));
        fN = fullfile(pwd,strcat('/MK359LH/Results/Not Thresholded/preCO-stained/MK359LH_CombinedResults_NonThresholded_withoutCO_',runDate,'.mat'));
    end
end
if wantSave
    save(fN,'layers2','sectionImages');
    print( fig, '-dtiff', filename );
    if wantEPS
        saveas(fig,filename,'epsc');
    end
end

clear gInd tInd label lyr yEnd yStart x1 x2 verts sec sc numBins ...
    monkey midlineY midlineX midInd lyrBrdr locations loc ...
    k j intensity intens inds ind hemi gfpIntensity gfpIntens fN currLyrs ...
    brdr2 brdr border binWidthHalf binWidth binHt tdtIntensity tdtIntens ...
    sections runDate rotAngle rotAnchor bin binHt2 binHtExtra c col cols ...
    fig filename firstBins gfpBorder gfpImg ind2 inj labelFields labels lm lmC lmR ...
    lyrBrdrBin lyrC lyrR lyrs lyrs2 normIntens r rows tdtBorder tdtImg temp ...
    temp2 tempImg tempIntens tempLM tempLM2 x xEnd xStart yP yPix yPt gfpBorder2 ...
    tdtBorder2 brdrF gfpBrdrF tdtBrdrF gfpBrdr tdtBrdr plotStyle sec sections sectionsMidsOnly ...
    binCheck binYCoords gfpIntensityMids gfpIntensMids gfpSEM gfpMidsSEM gX inj2 ...
    intensMids isAllCases labelMidlines lyrsText numInjs scales tX ...
    tdtSEM tdtMidsSEM tdtIntensMids tdtIntensityMids isAllCases brdrFMids ...
    brdrMids textY wantEPS wantSave wantSeeIndividualPlots intensities

function compCopy(op, np)
%COMPCOPY copies a figure object represented by "op" and its % descendants to another figure "np" preserving the same hierarchy.
     ch = get(op, 'children');
     if ~isempty(ch)   
         nh = copyobj(ch,np);
         for k = 1:length(ch)
             compCopy(ch(k),nh(k));
         end
     end
     return
end
