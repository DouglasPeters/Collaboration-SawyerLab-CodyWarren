clear; clc; tic; cd(userpath);
%% Editables %%

Folder = 'F:\DATA\01.27.21 Phenix for Cody\'; 
%Where are the images located? IMPORTANT: Use format 'FILEPATH\'. The apostrophes and ending slash are important.

FigShow = 1; %Do you want to show the Figure with segmentation overlay? (1=yes, 0=no)
FigAnnotate = 1; %Do you want to show the number designation for nuclei and cell objects on the Figure? (1=yes, 0=no)

FigSave = 1; %Do you want to save the Figure that is generated during analysis? (1=yes, 0=no)
MIP_ImageSave = 1; %Do you want to save the max intensity projection image hyperstack (XYC)? (1=yes, 0=no)

Channels = 3; %How many fluorescent channels are in the image?

CH_DAPI = 1; %Which channel corresponds to DAPI signal?

CH_CellMask = 2; %Which channel corresponds to CellMask (CM) signal?
CM_LowSizeFilt = 1000; %What is the minimum area of cell objects that you want to include in your analysis? (pixels^2)



%% Analysis Pre-Reqs and Metadata %%

cd(Folder);
srcFiles_Folder = dir(Folder);
srcFiles_Folder = srcFiles_Folder(3:end);

for b = 5:length(srcFiles_Folder)
clearvars -except Folder FigShow FigAnnotate FigSave MIP_ImageSave Channels CH_DAPI CH_CellMask CM_LowSizeFilt b srcFiles_Folder;
cd(Folder); cd(srcFiles_Folder(b).name); cd Images; cd ImageStacks;
srcFiles = dir('*.tiff');
Image_Path = strcat(Folder,srcFiles_Folder(b).name,'\Images\ImageStacks\');

START = 1;
FINISH = length(srcFiles);

for u = 1:length(srcFiles)
    srcFiles_FOV(u,1) = str2num(srcFiles(u).name(8:9));
    srcFiles_WELL(u,1) = string(srcFiles(u).name(1:6));
end
Fields = max(srcFiles_FOV);
Wells = numel(srcFiles)/Fields;

if FigShow == 1, figure,
else end

%% Analysis %%
for f = START:FINISH
    time(f,1).ElapsedSeconds = toc;
    
    %try
    
clc
filename = strcat(Image_Path,srcFiles(f).name);
progress = (((FINISH-START+1)-(FINISH-f))/FINISH)*100;
progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,FINISH,progress,'%');
disp(progress2);

if f == START
        cd(Image_Path); mkdir('Analysis'); cd(Image_Path);
else end

if progress < 10
    disp('Estimated time remaining will display after 10% of images are analyzed...');
else
    time(f).AverageSecondsPerLoop = time(f).ElapsedSeconds/((FINISH-START+1)-(FINISH-f));
    time(f).EstimatedTotalSeconds = time(f).AverageSecondsPerLoop*(FINISH-START+1);
    time(f).EstimatedSecondsLeft = time(f).EstimatedTotalSeconds-time(f).ElapsedSeconds;
    time(f).EstimatedMinutesLeft = time(f).EstimatedSecondsLeft/60;
    time(f).EstimatedMinutesElapsed = time(f).ElapsedSeconds/60;
    estimate = sprintf('Run time: %0.2f minutes.',time(f).EstimatedMinutesElapsed);
    estimate2 = sprintf('Estimated time remaining: %0.2f minutes.',time(f).EstimatedMinutesLeft);
    disp(estimate);
    disp(estimate2);
end
    
Results(f).FileName = srcFiles(f).name;
Results(f).Well = srcFiles_WELL{f};
Results(f).FieldofView = srcFiles_FOV(f);

I = bfopen(filename);

ResY = size(I{1,1}{1,1},1);
ResX = size(I{1,1}{1,1},2);
ZPlanes = (length(I{1,1})/Channels);
Blank3D = zeros(ResY,ResX,ZPlanes);
Slices = Channels*ZPlanes; 

%% Parse Channels and Create Max Intensity Projection Images %%

disp('Generating in-focus images...');
[Ch1,Ch2,Ch3,Ch4,Ch1_MIP,Ch2_MIP,Ch3_MIP,Ch4_MIP] = MaxIntProj(I,ZPlanes,Channels);

if MIP_ImageSave == 1
    disp('Saving max intensity projection images...');
    if f == START
        cd(Image_Path); cd Analysis; mkdir('MIPImages'); cd MIPImages;
    else cd(Image_Path); cd Analysis; cd MIPImages; 
    end
        
    MIP_IMAGE = zeros(ResY,ResX,Channels, 'uint16');  
    MIP_IMAGE(:,:,1) = Ch1_MIP;
    MIP_IMAGE(:,:,2) = Ch2_MIP;
    if Channels > 2, MIP_IMAGE(:,:,3) = Ch3_MIP; else end
    if Channels > 3, MIP_IMAGE(:,:,4) = Ch4_MIP; else end
    MIP_IMAGENAME = strcat(srcFiles(f).name(1:9),' MIP Image.tiff');
    bfsave(MIP_IMAGE(:,:,:),MIP_IMAGENAME);
else
end

%% Cell Mask Segmentation %%

disp('Segmenting Cell Mask signal...');
if CH_CellMask == 1, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Ch1_MIP,CM_LowSizeFilt,ResY,ResX);
elseif CH_CellMask == 2, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Ch2_MIP,CM_LowSizeFilt,ResY,ResX);
elseif CH_CellMask == 3, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Ch3_MIP,CM_LowSizeFilt,ResY,ResX);
elseif CH_CellMask == 4, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Ch4_MIP,CM_LowSizeFilt,ResY,ResX); 
else end

%% Nuclear Segmentation %%

disp('Segmenting DAPI signal...');
if CH_DAPI == 1, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Ch1_MIP,ResY,CM_Watershed_BW2);
elseif CH_DAPI == 2, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Ch2_MIP,ResY,CM_Watershed_BW2);
elseif CH_DAPI == 3, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Ch3_MIP,ResY,CM_Watershed_BW2);
elseif CH_DAPI == 4, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Ch4_MIP,ResY,CM_Watershed_BW2); 
else end

%% Cellular Analysis %%

disp('Quantifying image properties for each cell/nucleus...');
if CH_DAPI == 1, DAPI = Ch1_MIP;
elseif CH_DAPI == 2, DAPI = Ch2_MIP;
elseif CH_DAPI == 3, DAPI = Ch3_MIP;
elseif CH_DAPI == 4, DAPI = Ch4_MIP;
else end

if CH_CellMask == 1, CellMask = Ch1_MIP;
elseif CH_CellMask == 2, CellMask = Ch2_MIP;
elseif CH_CellMask == 3, CellMask = Ch3_MIP;
elseif CH_CellMask == 4, CellMask = Ch4_MIP;
else end

clearvars NearestNuc;
[Results_CellAnalysis,NearestNucDistanceFiltered] = CellularAnalysis(DAPI_75Percentile,DAPI_Watershed_BW2,Channels,CMseg_props,CM_IndCells,DAPI,Ch2_MIP,Ch3_MIP,CellMask);

clearvars NucleiNumber AdjustedNucleiNumber;
for q = 1:size(Results_CellAnalysis,2)
    NucleiNumber(q,1) = Results_CellAnalysis(q).NucleiNumber; 
    AdjustedNucleiNumber(q,1) = Results_CellAnalysis(q).AdjustedNucleiNumber;
end
   
%% Figure %%

if FigShow > 0,
    disp('Generating Figure...');
totalfilteredcellperims = false(ResY,ResX);
for v = 1:size(Results_CellAnalysis,2)
    if v == 1, totalfilteredcellperims = Results_CellAnalysis(v).LogPerim;
    else totalfilteredcellperims = or(totalfilteredcellperims,Results_CellAnalysis(v).LogPerim);
    end
end

figimage1 = imoverlay(imadjust(DAPI),imdilate(DAPI_Watershed_Perim,strel('disk',1)),[0 0 1]);
figimage1 = imoverlay(figimage1,imdilate(totalfilteredcellperims,strel('disk',1)),[0.5 0 0]);
figimage2 = imoverlay(imadjust(CellMask),imdilate(totalfilteredcellperims,strel('disk',1)),[0.5 0 0]);
figimage2 = imoverlay(figimage2,imdilate(DAPI_Watershed_Perim,strel('disk',1)),[0 0 0.5]);
figimage3 = Blank3D(:,:,1:3);

C = [figimage1 figimage2 figimage3];
imshow(C); title('Segmentation and Analysis Summary');

dtclock = fix(clock);
scriptname = mfilename('fullpath');
scriptname_slash = strfind(scriptname,'\');
scriptname_str = scriptname(scriptname_slash(end)+1:end);
scriptname_underscore = strfind(scriptname_str,'_');
scriptname_str(scriptname_underscore) = ' ';
title('Segmentation and Analysis Summary');
figtext1 = ['Script: ' scriptname_str];
if dtclock(5)>10, figtext2 = ['Date and Time of Analysis: ' num2str(dtclock(2)) '/' num2str(dtclock(3)) '/' num2str(dtclock(1)) ' at ' num2str(dtclock(4)) ':' num2str(dtclock(5))];
else figtext2 = ['Date and Time of Analysis: ' num2str(dtclock(2)) '/' num2str(dtclock(3)) '/' num2str(dtclock(1)) ' at ' num2str(dtclock(4)) ':0' num2str(dtclock(5))];
end
figtext3 = ['Image Filename: ' srcFiles(f).name];
figtext4 = [' '];
figtext6 = ['Nuclei Detected: ' num2str(sum(NucleiNumber,'all'))];
figtext7 = ['Nuclei Detected (Adjusted): ' num2str(sum(AdjustedNucleiNumber,'all'))];
figtext8 = ['Cell Objects Detected: ' num2str(size(Results_CellAnalysis,2))];
figtext_all = {figtext1,figtext2,figtext3,figtext4,figtext6,figtext7,figtext8};
t = text((ResX*2)+100,500,figtext_all);
t.Color = [1 1 1];
t.FontSize = 12;

if FigAnnotate >0,
    for b = 1:size(Results_CellAnalysis,2)
        cmtext = text((ResX+Results_CellAnalysis(b).CMCentroid(1)),(Results_CellAnalysis(b).CMCentroid(2)),num2str(b));
        cmtext.Color = [1 1 0];
        cmtext.FontSize = 10;
        cmtext.FontWeight = 'bold';
        for n = 1:numel(Results_CellAnalysis(b).NuclearProps)
            nuctext = text(Results_CellAnalysis(b).NuclearProps(n).Centroid(1),Results_CellAnalysis(b).NuclearProps(n).Centroid(2),num2str(n));
            nuctext.Color = [1 1 0];
            nuctext.FontSize = 10;
            nuctext.FontWeight = 'bold';
        end
    end       
else
end

drawnow; hold off;
else
end

if FigSave == 1
    disp('Saving Figure...');
    if f == START,
        cd(Image_Path); cd Analysis; mkdir('FigureImages'); cd FigureImages;
    else cd(Image_Path); cd Analysis; cd FigureImages;
    end
    ax = gca;
    Fig_Name = append(filename(end-17:end-9),' Segmentation.tif');
    exportgraphics(ax,Fig_Name,'ContentType','image','Resolution','400');
else
% elseif f > START
%     cd(Folder); cd Analysis; cd FigureImages;
end

%% Results %%

clearvars CellularMeans CellularSums CMAreas;
disp('Collating results...')

for r = 1:size(Results_CellAnalysis,2)
    CellularMeans(r,1) = Results_CellAnalysis(r).MeanCh1;
    CellularMeans(r,2) = Results_CellAnalysis(r).MeanCh2;
    if Channels>2, CellularMeans(r,3) = Results_CellAnalysis(r).MeanCh3; else end;
    if Channels>3, CellularMeans(r,4) = Results_CellAnalysis(r).MeanCh4; else end;
    CellularSums(r,1) = Results_CellAnalysis(r).SumCh1;
    CellularSums(r,2) = Results_CellAnalysis(r).SumCh2;
    if Channels >2, CellularSums(r,3) = Results_CellAnalysis(r).SumCh3; else end;
    if Channels >3, CellularSums(r,4) = Results_CellAnalysis(r).SumCh4; else end;
    CMAreas(r,1) = Results_CellAnalysis(r).CellArea;    
end

Results(f).NucNearestNeighborDistance = NearestNucDistanceFiltered;
Results(f).TotalNuclei = sum(NucleiNumber,'all');
Results(f).NucleiPerGroup = NucleiNumber;
Results(f).AdjustedTotalNuclei = sum(AdjustedNucleiNumber,'all');
Results(f).AdjustedNucGroupSizes = AdjustedNucleiNumber;
Results(f).MeanIntensities = CellularMeans; %Each column is a channel.
Results(f).SumIntensities = CellularSums; %Each column is a channel.
Results(f).CMNumber = size(Results_CellAnalysis,2);
Results(f).CMAreas = CMAreas;
Results(f).CMTotalArea = sum(CMAreas);

%     catch
%          warning('An error occurred during analysis. Saving Results and skipping to next image.'); pause(2);
%          Results = Results(START:f);
%          cd(Image_Path); cd Analysis;
%          save('AnalysisResults.mat','Results', '-v7.3');
%     end
end %End of Analysis Loop

cd(Image_Path); cd Analysis;
save('AnalysisResults.mat','Results', '-v7.3');

    %% Result Sorting By Well%%

disp('Finished analyzing images, now sorting results by well...');
[Results_Sorted,WellNumber] = WellSort(Channels,Results);
cd(Image_Path); cd Analysis;
save('AnalysisResultsSortedByWell.mat','Results_Sorted','-v7.3');

end
toc