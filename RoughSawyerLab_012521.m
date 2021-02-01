clear all; clc; tic

Folder = 'F:\DATA\01.25.21 SDC Imaging for Cody\';

cd(Folder);
srcFiles = dir('*.nd2');

Channels = 3;
ShowFig = 0;

START = 1;
FINISH = length(srcFiles);

if ShowFig ==1, figure; else end
for f = START:FINISH;
        clc
    filename = strcat(Folder,srcFiles(f).name);
    progress = (((FINISH-START+1)-(FINISH-f))/FINISH)*100;
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,FINISH,progress,'%');
    disp(progress2);
    
    clearvars Raw Ch1 Ch2 Ch3;
    Raw = bfopen(filename);
    
    ResY = size(Raw{1,1}{1,1},1);
    ResX = size(Raw{1,1}{1,1},2);
    Planes = (length(Raw{1,1})/Channels);
    
    
    for i = 1:Planes
    Ch1_planes(i,1) = 1+(Channels*i-Channels);
    Ch1(:,:,i) = Raw{1,1}{Ch1_planes(i,1),1};
    
    Ch2_planes(i,1) = 2+(Channels*i-Channels);
    Ch2(:,:,i) = Raw{1,1}{Ch2_planes(i,1),1};
    
    if Channels>2,
        Ch3_planes(i,1) = 3+(Channels*i-Channels);
        Ch3(:,:,i) = Raw{1,1}{Ch3_planes(i,1),1};
    else
    end
    if Channels>3, 
        Ch4_planes(i,1) = 4+(Channels*i-Channels); 
        Ch4(:,:,i) = Raw{1,1}{Ch4_planes(i,1),1};
    else
    end
    
    end
    
    Ch1_MIP = uint16(sum(Ch1,3));
    Ch2_MIP = uint16(sum(Ch2,3));
    if Channels >2, Ch3_MIP = uint16(sum(Ch3,3)); else end
    if Channels >3, Ch4_MIP = uint16(sum(Ch4,3)); else end
    
    disp('Making Binary Mask...');
    BW = imbinarize(imopen(Ch2_MIP,strel('disk',10)));
    BW2 = imfill(BW,'holes');
    BW3 = imerode(imdilate(BW2,strel('disk',5)),strel('disk',2));
    BW4 = bwareaopen(BW3,10000);
    BW4perim = imdilate(bwperim(BW4),strel('disk',2));
    
    BW4_cc = bwconncomp(BW4,8);
    
    FISH_props = regionprops('table',BW4_cc,Ch1_MIP,'Area','MeanIntensity','MaxIntensity','PixelValues');
    Results(f).Filename = filename;
    Results(f).MeanFISH = FISH_props.MeanIntensity;
    Results(f).MaxFISH = FISH_props.MaxIntensity;
    
    for c = 1:size(FISH_props,1)
    Results(f).PrctileVals(c,1) = prctile(FISH_props.PixelValues{c,1},99.9,'all');
    Results(f).PrctileLogical(c).Logical = FISH_props.PixelValues{c,1}(:)>Results(f).PrctileVals(c,1);
    Results(f).PrctileMean(c,1) = mean(FISH_props.PixelValues{c,1}(Results(f).PrctileLogical(c).Logical));
    
    end
    
    if ShowFig == 1,
    disp('Generating Figure...');
    figimage1 = imoverlay(Ch2_MIP,BW4perim);
    figimage2 = imoverlay(Ch1_MIP,BW4perim);
    figimage3 = imoverlay(Ch3_MIP,BW4perim);
    C = [figimage1 figimage2 figimage3];
    imshow(C);
    else end;
    
end

toc