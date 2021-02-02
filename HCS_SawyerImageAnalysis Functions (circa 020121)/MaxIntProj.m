

function [Ch1,Ch2,Ch3,Ch4,Ch1_MIP,Ch2_MIP,Ch3_MIP,Ch4_MIP] = MaxIntProj(I,ZPlanes,Channels)

for m = 1:ZPlanes
    Ch1(:,:,m) = I{1,1}{m,1}(:,:);
    Ch2(:,:,m) = I{1,1}{ZPlanes+m,1}(:,:);
    if Channels>2, Ch3(:,:,m) = I{1,1}{(2*ZPlanes)+m,1}(:,:); 
    else Ch3 = zeros(1); end
    if Channels>3, Ch4(:,:,m) = I{1,1}{(3*ZPlanes)+m,1}(:,:); 
    else Ch4 = zeros(1); end
end

Ch1_MIP = uint16(sum(Ch1,3));
Ch2_MIP = uint16(sum(Ch2,3));
if Channels >2, Ch3_MIP = uint16(sum(Ch3,3)); else end
if Channels >3, Ch4_MIP = uint16(sum(Ch4,3)); else end

if exist('Ch3_MIP') == 0, Ch3_MIP = Ch3; else end
if exist('Ch4_MIP') == 0, Ch4_MIP = Ch4; else end

end