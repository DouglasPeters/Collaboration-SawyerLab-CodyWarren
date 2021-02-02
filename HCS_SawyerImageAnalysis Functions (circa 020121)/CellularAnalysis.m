
function [Results_CellAnalysis,NearestNucDistanceFiltered] = CellularAnalysis(DAPI_75Percentile,DAPI_Watershed_BW2,Channels,CMseg_props,CM_IndCells,DAPI,Ch2_MIP,Ch3_MIP,CellMask)

for i = 1:size(CM_IndCells,2)    
    Results_CellAnalysis(i).LogImage = CM_IndCells(i).LogImage;
    Results_CellAnalysis(i).LogPerim = CM_IndCells(i).LogPerim;
    Results_CellAnalysis(i).CellArea = CMseg_props(i).Area;
    Results_CellAnalysis(i).CMCentroid = CMseg_props(i).Centroid;
    Results_CellAnalysis(i).NuclearImage = DAPI_Watershed_BW2.*CM_IndCells(i).LogFilledHoles;
    Results_CellAnalysis(i).NucleiCC = bwconncomp(Results_CellAnalysis(i).NuclearImage,8);
    Results_CellAnalysis(i).NuclearProps = regionprops(Results_CellAnalysis(i).NucleiCC,DAPI,'Area','Circularity','Centroid');
  
    NuclearArea = [];
    
    for p = 1:Results_CellAnalysis(i).NucleiCC.NumObjects
        NuclearArea(p,1) = Results_CellAnalysis(i).NuclearProps(p).Area;
    end
    Results_CellAnalysis(i).TotalNuclearArea = sum(NuclearArea,'all');
    Results_CellAnalysis(i).NucleiNumber = length(Results_CellAnalysis(i).NuclearProps);
    Results_CellAnalysis(i).NucAreaThreshForAdjustment = DAPI_75Percentile*3;
    
    if size(Results_CellAnalysis(i).NuclearProps,1) > 0,
        for a = 1:size(Results_CellAnalysis(i).NuclearProps,1)
            RawNucAreas(a,1) = Results_CellAnalysis(i).NuclearProps(a).Area;
        end
    else RawNucAreas = 0;
    end
    
    NucAdjustmentFilter = RawNucAreas > Results_CellAnalysis(i).NucAreaThreshForAdjustment;
    
    Results_CellAnalysis(i).NuclearAreaToBeAdjusted = sum(RawNucAreas(NucAdjustmentFilter),'all');
    Results_CellAnalysis(i).AdjustedNucleiNumber = Results_CellAnalysis(i).NucleiNumber + round(Results_CellAnalysis(i).NuclearAreaToBeAdjusted/(DAPI_75Percentile));
    
    Results_CellAnalysis(i).MeanCh1 = mean(DAPI(CM_IndCells(i).LogFilledHoles),'all');
    Results_CellAnalysis(i).SumCh1 = sum(DAPI(CM_IndCells(i).LogFilledHoles),'all');
    Results_CellAnalysis(i).MeanCh2 = mean(Ch2_MIP(CM_IndCells(i).LogFilledHoles),'all');
    Results_CellAnalysis(i).SumCh2 = sum(Ch2_MIP(CM_IndCells(i).LogFilledHoles),'all');
    if Channels>2, 
        Results_CellAnalysis(i).MeanCh3 = mean(Ch3_MIP(CM_IndCells(i).LogFilledHoles),'all');
        Results_CellAnalysis(i).SumCh3 = sum(Ch3_MIP(CM_IndCells(i).LogFilledHoles),'all');
    else
    end;
    if Channels>3, 
        Results_CellAnalysis(i).MeanCh4 = mean(CellMask(CM_IndCells(i).LogFilledHoles),'all');
        Results_CellAnalysis(i).SumCh4 = sum(CellMask(CM_IndCells(i).LogFilledHoles),'all');
    else
    end;
end

QCFilter = false(0);
for m = 1:size(Results_CellAnalysis,2)
QCFilter(m,1) = size(Results_CellAnalysis(m).NuclearProps,1)>0 && Results_CellAnalysis(m).CellArea<1000000;
end

Results_CellAnalysis = Results_CellAnalysis(QCFilter); %Only CellMask objects containing 1 or more DAPI objects will pass this filter. This removes cells on the very edges of images, as well as fragments of CellMask signal within images.

%% Nearest Neighbor Analysis %%

for j = 1:size(Results_CellAnalysis,2)
    clearvars knnCoord;
    for n = 1:Results_CellAnalysis(j).NucleiNumber
        knnCoord(n,1) = Results_CellAnalysis(j).NuclearProps(n).Centroid(1,1);
        knnCoord(n,2) = Results_CellAnalysis(j).NuclearProps(n).Centroid(1,2);
    end

Results_CellAnalysis(j).NearestNuc = zeros(size(Results_CellAnalysis(j).NuclearProps,1),3);

    for t = 1:Results_CellAnalysis(j).NucleiNumber
    knnCoord2 = knnCoord(knnCoord~=knnCoord(t,:));
    
        if numel(knnCoord2)>0,
            knnCoord3 = reshape(knnCoord2,size(knnCoord2,1)/2,2);
            [knnIdx,knnDist] = knnsearch(knnCoord3,knnCoord(t,:));
            Results_CellAnalysis(j).NearestNuc(t,1:2) = [knnIdx,knnDist];
            Results_CellAnalysis(j).NearestNuc(t,3) = Results_CellAnalysis(j).NucleiNumber;
        else
            knnCoord3 = 0;
            knnIdx = 0; knnDist = 0;
            Results_CellAnalysis(j).NearestNuc(t,1:2) = 0;
            Results_CellAnalysis(j).NearestNuc(t,3) = Results_CellAnalysis(j).NucleiNumber;
        end

    end
    
    for s = 1:size(Results_CellAnalysis(j).NearestNuc,1)
        if s==1 && j==1,
            NearestNucDistance(1,1) = Results_CellAnalysis(j).NearestNuc(s,2);
            NearestNucDistance(1,2) = Results_CellAnalysis(j).NearestNuc(s,3);
        else
            NearestNucDistance(end+1,1) = Results_CellAnalysis(j).NearestNuc(s,2);
            NearestNucDistance(end,2) = Results_CellAnalysis(j).NearestNuc(s,3);
        end
    end

Results_CellAnalysis(j).NearestNucFiltered(:,1) = unique(Results_CellAnalysis(j).NearestNuc(:,2));
NNFiltered(j,1) = Results_CellAnalysis(j).NearestNucFiltered(1,1);
end
NNFilteredLog = NNFiltered>0;
if sum(NNFilteredLog) > 0,
    [NearestNucDistanceFiltered1, ia, ic] = unique(NearestNucDistance(:,1),'stable');
    NNZeroFilt = NearestNucDistanceFiltered1>0;
    NearestNucDistanceFiltered(:,1) = NearestNucDistanceFiltered1(NNZeroFilt);
    NearestNucDistanceFiltered(:,2) = NearestNucDistance(ia(NNZeroFilt),2);
else
    NearestNucDistanceFiltered(1,1) = 0;
    NearestNucDistanceFiltered(1,2) = 0;
end


