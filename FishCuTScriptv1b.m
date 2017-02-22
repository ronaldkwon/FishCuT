function FishCuTScriptv9(inPath,inFileName,threshold)

for i = 1:9
    figure(i)
    close(i)
end

Resolution = 21;
f = filesep;
Version = 'beta';
inPath = [inPath filesep];
Algorithm = 'Default';
thresholdMultiplier = threshold;
rotAngle = 0;
%start at first verterbra following webberian apparatus, and count until
%two vertebrae from tail fins 

%Read in the image
im = squeeze(dicomread([inPath inFileName]));

%Turn the matrix into grayscale values.
grayim = mat2gray(im,[-2^15 2^15]);

%Plot the resulting intensity projection.
sipDummy = max(grayim,[],2);

sip = reshape(sipDummy,size(sipDummy,1),size(sipDummy,3));

figure(1);
imshow(sip,[0 1],'Border','tight');

% standardLength = (sqrt(diff(z).^2 + diff(y).^2)*21)/1000;
% thresholdMultiplier = -50.67*standardLength^(-1.767)+1;

sipEdges = edge(sip,'Canny');
se90 = strel('line', 8.5, 90);
se0 = strel('line', 8.5, 0);
BWsdil = imdilate(sipEdges, [se90 se0]);
BWdfill = imfill(BWsdil, 'holes');
CC = bwconncomp(BWdfill);
[~,fishIdx] = max(cellfun(@numel,CC.PixelIdxList));
filteredFishIdxs = CC.PixelIdxList{fishIdx};
filteredFish = sip;
filteredFish(filteredFishIdxs) = 1;
filteredFish(filteredFish~=1) = 0;
seD = strel('diamond',3);
BWfinal = imerode(filteredFish,seD);
BWfinal = imerode(BWfinal,seD);
BWoutline = bwperim(BWfinal);
outlinedSip = sip; outlinedSip(BWoutline) = 1;

figure(1);
siphandle1 = imshow(outlinedSip);

disp('Displayed is automatic fish detection')
userInput = input('Enter "m" to manually outline, or "e" to proceed: ','s');
keepAsking = 1;
while keepAsking
    if strcmp(userInput,'m')
        h = impoly(gca);
        BW = createMask(h);
        fishim = sip;
        fishim(~BW) = 0;
        keepAsking = 0;
    elseif strcmp(userInput,'e')
        fishim = sip;
        fishim(~logical(BWfinal)) = 0;
        keepAsking = 0;
    else
        disp('Please enter "m" or "e"')
        userInput = input('Enter "m" to manually outline, or "e" to proceed: ','s');
    end
end
% 
% bone = 5540;
% [a,b] = imhist(fishim(:),2^16);
% a = a(2:end); b = b(2:end);
% cumDensity = cumsum(a)./sum(a(:));
% [~,percIndx] = min(abs(b-(bone/2^16+0.5)));
% bonePerc = cumDensity(percIndx);

% fishindxs = find(fishim(:));
% GMModel1 = fitgmdist(sqrt(fishim(fishindxs)),2);
% 
% bins = linspace(-1,0,2^16)';
% wholeDist = pdf(GMModel1,bins);
% means = GMModel1.mu; stds = sqrt(squeeze(GMModel1.Sigma));
% dist1 = normpdf(bins,means(1),stds(1));
% dist2 = normpdf(bins,means(2),stds(2));
% 
% figure
% bar(bins,dist1,'r')
% hold on
% bar(bins,dist2);

mijirun = 1;
classpath = javaclasspath('-all');
i=1;
while i <= size(classpath,1)
    [~,name,~] = fileparts(classpath{i});
    if strcmp(name,'z_spacing-0.0.1-SNAPSHOT')
        mijirun = 0;
        break
    end
    i = i+1;
end

if mijirun
    Miji(false);
end

fishim(fishim==0.5) = 0;
MIJ.createImage('SIP',uint16(fishim*2^16),true);
MIJ.run('Auto Threshold',['method=' Algorithm ' ignore_black ignore_white white']);
threshindxs = MIJ.getCurrentImage;
MIJ.run('Close All')
threshindxs(threshindxs~=0) = 1;
fishim(~logical(threshindxs)) = 0;
threshold = (((2^16*min(fishim(fishim~=0))-2^15)*thresholdMultiplier)+2^15)/2^16;

cutspine = 1;

newim = grayim;
newim(newim < threshold) = 0;
newim(newim ~= 0) = 1;

bwImg = newim;

manualMarks = zeros(size(newim,1),size(newim,2),size(newim,3));

vk = 1;
VertebralxCoordsAll = [];
VertebralzCoordsAll = [];

AllCoords = cell(0);

vloc = 0;

while cutspine
    
    %compute connectivity
    tic;
    disp('Computing connectivity...');
    regions = bwconncomp(bwImg,6); %regions is a structure containing connected regions
    pixelIDs = regions.PixelIdxList; %compute pixel index list
    numPixels = cellfun(@numel,pixelIDs);    
    regProps = regionprops(regions,'Centroid'); %regProps is an array of structures with region properties
    centroids = cat(1,regProps.Centroid);
    
    sortedpixelIDs = sortrows([pixelIDs.' num2cell(numPixels.') num2cell(centroids)],2);
    centroids = cell2mat(sortedpixelIDs(:,3:end));
    sortedpixelIDs = sortedpixelIDs(:,1:2);
    
    newregions = regions;
    newregions.PixelIdxList = sortedpixelIDs(:,1).';
    labelImg = labelmatrix(newregions); %compute labelled matrix
    numRegions = regions.NumObjects;
    toc;

    %compute MIP for label matrix, visualize using custom colormap
    tic;
    disp('Making labelled projections...'); 
    mipDummy = max(labelImg,[],2); %compute maxima along dimension 2 
    mip = reshape(mipDummy,size(mipDummy,1),size(mipDummy,3)); %reshape mip
    omip = mip;
    figure(3)
    close(3)
    figure(3)
    mipHandle = imshow(mip,'DisplayRange',[0 numRegions],'Border','tight');
    map = colorcube;
    map = map(1:end-64,:);
    map = flipud(repmat(map,round(numRegions/length(map)),1));
    map(1,:) = [0 0 0];
    colormap(map);
    
    %compute SIP for segmented image, visualize
    graynewim = grayim;
    graynewim(~bwImg)=0;
    sipDummy = max(graynewim,[],2);
    sip = reshape(sipDummy,size(sipDummy,1),size(sipDummy,3));
    figure(4)
    close(4)
    figure(4)
    sipHandle = imshow(sip,'DisplayRange',[0 1],'Border','tight');
    colormap(bone);
    toc;
    
    askinput = 1;
    while askinput
        remainInLoop1 = input('Enter "s" to segment, or "e" to exit: ','s');
        askinput = ~ismember(remainInLoop1,{'s','e'});
        if askinput
            disp('');
            disp('Make sure to enter either "s" or "e"!')
        end
    end
    
    AutoSegment = 0;
    ManualSegment = 0;
    if remainInLoop1 == 'e';
        break
    else
        askinput = 1;
        while askinput
            lineType = input('Enter "v" for vertebral segmentation, "m" for manual segmentation or "e" to exit: ','s');
            askinput = ~ismember(lineType,{'v','m','e'});
            if askinput
                disp('')
                disp('Make sure to enter either "v", "m" or "e"!')
            end
        end
        
        if strcmp(lineType,'v')
            AutoSegment = 1;
        elseif strcmp(lineType,'m') 
            ManualSegment = 1;
        end
    end
    
    if AutoSegment
        bwImg = newim;
        bwImg(logical(manualMarks)) = 0;
        
        if vloc ~= 0
            AllCoords(vloc) = [];
        end
        
        vloc = size(AllCoords,2)+1;
        
        vk = 0;
        VertebralxCoordsAll = [];
        VertebralzCoordsAll = [];
        %obtain line coordinates from user
        disp('Enter line segments, press space to finish');        
        remainInLoop2 = 1;
        lineData = [];
        numPoints = 0;
        while remainInLoop2 == 1
            figure(4)
            hold on
            [zCoords,xCoords,buttonPress] = ginputc(1,'Color',[1 1 1]); %get x and z line coordinates
            if buttonPress == 1
                numPoints = numPoints+1;
                plot(zCoords,xCoords,'*','Color','w','Tag',['Point' num2str(numPoints)])
                lineData = [lineData;xCoords zCoords];
                
            elseif buttonPress == 98 
                if size(lineData,1) == 1
                    lineData = [];
                    tempPoint1 = findobj('Tag',['Point' num2str(vk*2+1)]);
                    delete(tempPoint1)
                    numPoints = numPoints - 1;
                elseif size(lineData,1) == 0 && size(VertebralxCoordsAll,1) > 0
                    VertebralxCoordsAll(vk,:) = [];
                    VertebralzCoordsAll(vk,:) = [];
                    tempLine = findobj('Tag',['VertebralLine' num2str(vk)]);
                    delete(tempLine)
                    tempPoint1 = findobj('Tag',['Point' num2str(vk*2-1)]);
                    delete(tempPoint1)
                    tempPoint2 = findobj('Tag',['Point' num2str(vk*2)]);
                    delete(tempPoint2)
                    vk = vk-1;
                    numPoints = numPoints - 2;
                end
                
            elseif buttonPress == 32
                remainInLoop2 = 0;
                
            elseif ~ishandle(4)
                remainInLoop2 = 0;
            
            else
                disp('Invalid point section because of extraneous key or mouse press.')
                disp('Please re-enter point or hit "b" to backspace.')
            end

            if size(lineData,1) == 2
                vk = vk + 1;
                line(lineData(:,2),lineData(:,1),'Tag',['VertebralLine' num2str(vk)],'Color','w');
                VertebralxCoordsAll(vk,:) = lineData(:,1); %record
                VertebralzCoordsAll(vk,:) = lineData(:,2);
                lineData = [];
            end
        end
        
        if ~isempty(VertebralxCoordsAll)
            ids = VertebralxCoordsAll(:,1) > VertebralxCoordsAll(:,2);
            VertebralxCoordsAll(ids,:) = fliplr(VertebralxCoordsAll(ids,:));
            VertebralzCoordsAll(ids,:) = fliplr(VertebralzCoordsAll(ids,:));

            AllCoords = [AllCoords {[VertebralxCoordsAll VertebralzCoordsAll]}];

            fullxs = VertebralxCoordsAll;
            fullzs = VertebralzCoordsAll;

            avexCoords = mean(VertebralxCoordsAll,2);
            avezCoords = mean(VertebralzCoordsAll,2);

            ylist = zeros(size(VertebralxCoordsAll,1),1);
            for m = 1:size(VertebralxCoordsAll,1)
                [x1,xmin] = min(VertebralxCoordsAll(m,:));
                [x2,xmax] = max(VertebralxCoordsAll(m,:));

                x1 = floor(x1);
                if x1 <= 0
                    x1 = 1;
                end
                x2 = ceil(x2);
                if x2 > size(sipDummy,1)
                    x2 = size(sipDummy,1);
                end
                xBar = floor(avexCoords(m));
                zBar = floor(avezCoords(m));
                if m == 1
                    dz = avezCoords(2)-avezCoords(1);
                elseif m == size(VertebralxCoordsAll,1)
                    dz = avezCoords(end)-avezCoords(end-1);
                else
                    dz = 0.5*(avezCoords(m+1)-avezCoords(m-1));
                end
                dx = dz;
                dy = dz;
                mdz = VertebralzCoordsAll(m,xmax)-VertebralzCoordsAll(m,xmin);
                mdx = x2-x1;
                xcoords = round(x1):round(x2);
                zarray1 = round(mdz/mdx*(xcoords-x1)+VertebralzCoordsAll(m,xmin));
                zarray2 = zarray1+1;
                grow = 1;

                while grow
                    [s1,~,s3] = size(bwImg);
                    xel1 = round(xBar-dx);
                    xel2 = round(xBar+dx);
                    zel1 = round(zBar-dz);
                    zel2 = round(zBar+dz);
                    if xBar-dx < 0.5
                        xel1 = 1;
                    end
                    if xBar+dx > s1
                        xel2 = s1;
                    end
                    if zBar+dz > s3
                        zel2 = s3;
                    end
                    if zBar-dz < 0.5
                        zel1 = 1;
                    end
                    bwSubImage = bwImg(xel1:xel2,:,zel1:zel2);
                    [xDummy,yDummy,zDummy] = ind2sub(size(bwSubImage),find(bwSubImage));
                    y = floor(mean(yDummy));
                    ylist(m) = mean(yDummy);
                    %make line mask
                    if ~isnan(y)
                        for n = 1:size(bwImg,2);
                            sliceToMask = squeeze(bwImg(:,n,:));
                            sliceToMask([sub2ind(size(sliceToMask),xcoords,zarray1)' sub2ind(size(sliceToMask),xcoords,zarray2)']) = 0;
                            bwImg(:,n,:) = sliceToMask;
                        end
                    end

                    tempseparation = bwImg;

                    if m == 1
                        lowz = round(VertebralzCoordsAll(1)-(VertebralzCoordsAll(2)-VertebralzCoordsAll(1)));
                        if lowz <= 0
                            lowz = 1;
                        end
                        highz = round(VertebralzCoordsAll(m+1));
                    elseif m == size(VertebralxCoordsAll,1)
                        highz = round(VertebralzCoordsAll(end)+(VertebralzCoordsAll(end)-VertebralzCoordsAll(end-1)));
                        if highz > size(sipDummy,3)
                            highz = size(sipDummy,3);
                        end
                        lowz = round(VertebralzCoordsAll(m-1));
                    else
                        highz = round(VertebralzCoordsAll(m+1));
                        lowz = round(VertebralzCoordsAll(m-1));
                    end
                    region = tempseparation(:,:,lowz:highz);
                    region(region~=0) = 1;
                    regions = bwconncomp(region,6);
                    labelImg = labelmatrix(regions);
                    split1 = labelImg(xel1:xel2,round(y-dy):round(y+dy)...
                        ,1:round(zBar-lowz));
                    split1(split1 == 0) = [];
                    tallied1 = tabulate(split1(:));
                    [~,group1] = max(tallied1(:,2));
                    group1 = tallied1(group1,1);

                    split2 = labelImg(round(xBar-dx):round(xBar+dx),round(y-dy):round(y+dy)...
                        ,round(zBar-lowz):end);
                    split2(split2 == 0) = [];
                    tallied2 = tabulate(split2);
                    [~,group2] = max(tallied2(:,2));
                    group2 = tallied2(group2,1);
                    if group1 ~= group2
                        grow = 0;
                    else
                        x1 = x1-2;
                        x2 = x2+2;
                        if x1 <= 0
                            x1 = 1;
                            if x2 > size(sipDummy,1)
                                grow = 0;
                            end
                        end
                        if x2 > size(sipDummy,1)
                            x2 = size(sipDummy,1);
                        end
                        xcoords = x1:x2;
                        zarray1 = round(mdz/mdx*(xcoords-x1)+VertebralzCoordsAll(m,xmin));
                        zarray2 = zarray1+1;
                    end
                end
            end
        end
        
    elseif ManualSegment
        k = 0;
        xCoordsAll = [];
        zCoordsAll = [];
        disp('Enter line segments, press space to finish');        
        remainInLoop2 = 1;
        lineData = [];
        numPoints = 0;
        while remainInLoop2 == 1
            figure(4)
            hold on
            [zCoords,xCoords,buttonPress] = ginputc(1,'Color',[1 1 1]); %get x and z line coordinates
            if buttonPress == 1
                numPoints = numPoints+1;
                plot(zCoords,xCoords,'*','Color','w','Tag',['Point' num2str(numPoints)])
                lineData = [lineData;xCoords zCoords];
                
            elseif buttonPress == 98 
                if size(lineData,1) == 1
                    lineData = [];
                    tempPoint1 = findobj('Tag',['Point' num2str(k*2+1)]);
                    delete(tempPoint1)
                    numPoints = numPoints - 1;
                elseif size(lineData,1) == 0 && size(xCoordsAll,1) > 0
                    xCoordsAll(k,:) = [];
                    zCoordsAll(k,:) = [];
                    tempLine = findobj('Tag',['Line' num2str(k)]);
                    delete(tempLine)
                    tempPoint1 = findobj('Tag',['Point' num2str(k*2-1)]);
                    delete(tempPoint1)
                    tempPoint2 = findobj('Tag',['Point' num2str(k*2)]);
                    delete(tempPoint2)
                    k = k-1;
                    numPoints = numPoints - 2;
                end
                
            elseif buttonPress == 32
                remainInLoop2 = 0;
                
            elseif ~ishandle(4)
                remainInLoop2 = 0;
            
            else
                disp('Invalid point section because of extraneous key or mouse press.')
                disp('Please re-enter point or hit "b" to backspace.')
            end

            if size(lineData,1) == 2
                k = k + 1;
                line(lineData(:,2),lineData(:,1),'Tag',['Line' num2str(k)],'Color','w');
                xCoordsAll(k,:) = lineData(:,1); %record
                zCoordsAll(k,:) = lineData(:,2);
                lineData = [];
            end
        end
        
        if ~isempty(xCoordsAll)
            AllCoords = [AllCoords {[xCoordsAll zCoordsAll]}];
            for m = 1:size(xCoordsAll,1)
                [x1,xmin] = min(xCoordsAll(m,:));
                [x2,xmax] = max(xCoordsAll(m,:));

                if xmin == xmax
                    xmax = 2;
                end

                x1 = floor(x1);
                if x1 <= 0
                    x1 = 1;
                end
                x2 = ceil(x2);
                if x2 > size(sipDummy,1)
                    x2 = size(sipDummy,1);
                end

                mdz = zCoordsAll(m,xmax)-zCoordsAll(m,xmin);
                mdx = x2-x1;
                xcoords = round(x1):round(x2);
                if mdx == 0
                    zarray1 = round(min(zCoordsAll)):round(min(zCoordsAll))+abs(mdz);
                    zarray2 = zarray1+1;
                else
                    zarray1 = round(mdz/mdx*(xcoords-x1)+zCoordsAll(m,xmin));
                    zarray2 = zarray1+1;
                end
                %make line mask
                for n = 1:size(manualMarks,2);
                    sliceToMask = squeeze(manualMarks(:,n,:));
                    sliceToMask([sub2ind(size(sliceToMask),xcoords,zarray1)' sub2ind(size(sliceToMask),xcoords,zarray2)']) = 1;
                    manualMarks(:,n,:) = sliceToMask;
                end
                bwImg(logical(manualMarks)) = 0;
            end
        end
    end
end

n=1;
vertIDs = cell(size(VertebralxCoordsAll,1)-1,1);

while n <= size(vertIDs,1)
    vertIDs = cell(size(VertebralxCoordsAll,1)-1,1);
    rollbackmips = cell(size(vertIDs,2)+1,1);
    rollbackmips{1} = reshape(mipDummy,size(mipDummy,1),size(mipDummy,3)); %reshape mip
    mip = rollbackmips{1};
    figure(3)
    n = 1;
    h = cell(size(vertIDs,1)+1,1);
    keepSelecting = 1;
    
    while n <= size(vertIDs,1)
        figure(3);
        [x,y,button] = ginputc(1,'Color',[1 1 1]);

        if button == 1
            r = mip(round(y),round(x));
            if r > 0
                disp(['Specifying vertebral body ',num2str(n)]);
                vertIDs{n} = [vertIDs{n} r];
                mip(mip == r) = 0;
                hold on;
                currenttxt = h{n};
                currenttxt(length(vertIDs{n})) = text(centroids(r,3),centroids(r,2),num2str(n),'FontSize',12,'Color','white');
                h{n} = currenttxt;
                hold off;
            end
        elseif button == 98
            if ~isempty(vertIDs{n})
                vertIDs{n} = [];
            else
                n = n - 1;
                if n < 1
                    n = 1;
                end
                vertIDs{n} = [];
            end
            mip = rollbackmips{n};
            txtarray = h{n};
            for i = 1:length(txtarray)
                delete(txtarray(i))
            end
            h{n} = [];
        elseif button == 32
            n = n + 1;
            rollbackmips{n} = mip;
        elseif button == 109
            k = 0;
            xCoordsAll = [];
            zCoordsAll = [];
            disp('Enter line segments, press space to finish');        
            remainInLoop2 = 1;
            lineData = [];
            numPoints = 0;
            while remainInLoop2 == 1
                figure(4)
                hold on
                [zCoords,xCoords,buttonPress] = ginputc(1,'Color',[1 1 1]); %get x and z line coordinates
                if buttonPress == 1
                    numPoints = numPoints+1;
                    plot(zCoords,xCoords,'*','Color','w','Tag',['Point' num2str(numPoints)])
                    lineData = [lineData;xCoords zCoords];

                elseif buttonPress == 98 
                    if size(lineData,1) == 1
                        lineData = [];
                        tempPoint1 = findobj('Tag',['Point' num2str(k*2+1)]);
                        delete(tempPoint1)
                        numPoints = numPoints - 1;
                    elseif size(lineData,1) == 0 && size(xCoordsAll,1) > 0
                        xCoordsAll(k,:) = [];
                        zCoordsAll(k,:) = [];
                        tempLine = findobj('Tag',['Line' num2str(k)]);
                        delete(tempLine)
                        tempPoint1 = findobj('Tag',['Point' num2str(k*2-1)]);
                        delete(tempPoint1)
                        tempPoint2 = findobj('Tag',['Point' num2str(k*2)]);
                        delete(tempPoint2)
                        k = k-1;
                        numPoints = numPoints - 2;
                    end

                elseif buttonPress == 32
                    remainInLoop2 = 0;

                elseif ~ishandle(4)
                    remainInLoop2 = 0;

                else
                    disp('Invalid point section because of extraneous key or mouse press.')
                    disp('Please re-enter point or hit "b" to backspace.')
                end

                if size(lineData,1) == 2
                    k = k + 1;
                    line(lineData(:,2),lineData(:,1),'Tag',['Line' num2str(k)],'Color','w');
                    xCoordsAll(k,:) = lineData(:,1); %record
                    zCoordsAll(k,:) = lineData(:,2);
                    lineData = [];
                end
            end

            if ~isempty(xCoordsAll)
                AllCoords = [AllCoords {[xCoordsAll zCoordsAll]}];
                for m = 1:size(xCoordsAll,1)
                    [x1,xmin] = min(xCoordsAll(m,:));
                    [x2,xmax] = max(xCoordsAll(m,:));

                    if xmin == xmax
                        xmax = 2;
                    end

                    x1 = floor(x1);
                    if x1 <= 0
                        x1 = 1;
                    end
                    x2 = ceil(x2);
                    if x2 > size(sipDummy,1)
                        x2 = size(sipDummy,1);
                    end

                    mdz = zCoordsAll(m,xmax)-zCoordsAll(m,xmin);
                    mdx = x2-x1;
                    xcoords = round(x1):round(x2);
                    if mdx == 0
                        zarray1 = round(min(zCoordsAll)):round(min(zCoordsAll))+abs(mdz);
                        zarray2 = zarray1+1;
                    else
                        zarray1 = round(mdz/mdx*(xcoords-x1)+zCoordsAll(m,xmin));
                        zarray2 = zarray1+1;
                    end
                    %make line mask
                    for sn = 1:size(manualMarks,2);
                        sliceToMask = squeeze(manualMarks(:,sn,:));
                        sliceToMask([sub2ind(size(sliceToMask),xcoords,zarray1)' sub2ind(size(sliceToMask),xcoords,zarray2)']) = 1;
                        manualMarks(:,sn,:) = sliceToMask;
                    end
                    bwImg(logical(manualMarks)) = 0;
                end

                regions = bwconncomp(bwImg,6); %regions is a structure containing connected regions
                pixelIDs = regions.PixelIdxList; %compute pixel index list
                numPixels = cellfun(@numel,pixelIDs);    
                regProps = regionprops(regions,'Centroid'); %regProps is an array of structures with region properties
                centroids = cat(1,regProps.Centroid);

                sortedpixelIDs = sortrows([pixelIDs.' num2cell(numPixels.') num2cell(centroids)],2);
                centroids = cell2mat(sortedpixelIDs(:,3:end));
                sortedpixelIDs = sortedpixelIDs(:,1:2);

                newregions = regions;
                newregions.PixelIdxList = sortedpixelIDs(:,1).';
                labelImg = labelmatrix(newregions); %compute labelled matrix
                numRegions = regions.NumObjects;
                toc;

                %compute MIP for label matrix, visualize using custom colormap
                tic;
                disp('Making labelled projections...'); 
                mipDummy = max(labelImg,[],2); %compute maxima along dimension 2 
                mip = reshape(mipDummy,size(mipDummy,1),size(mipDummy,3)); %reshape mip
                omip = mip;
                figure(3)
                close(3)
                figure(3)
                mipHandle = imshow(mip,'DisplayRange',[0 numRegions],'Border','tight');
                map = colorcube;
                map = map(1:end-64,:);
                map = flipud(repmat(map,round(numRegions/length(map)),1));
                map(1,:) = [0 0 0];
                colormap(map);
                set(mipHandle,'cdata',mip)
                break
            end
        end
        set(mipHandle,'cdata',mip)
    end
end

%Tissue Properties
grayim2 = mat2gray(im,[-2^15 2^15]);

FishSegFishImage = zeros(size(im));
FishSegVertImage = zeros(size(im));

%Isolate Vertebral parts
numBodies = size(vertIDs,1);
VertebralImages{numBodies} = [];
CentrumImages{numBodies} = [];
HaemalArchImages{numBodies} = [];
NeuralArchImages{numBodies} = [];
MuscleImages{numBodies} = [];
VertebralVolumes = zeros(numBodies,1);
CentrumVolumes = zeros(numBodies,1);
HaemalVolumes = zeros(numBodies,1);
NeuralVolumes = zeros(numBodies,1);
VertebralSAs = zeros(numBodies,1);
CentrumSAs = zeros(numBodies,1);
HaemalSAs = zeros(numBodies,1);
NeuralSAs = zeros(numBodies,1);
VertebralMIs = zeros(numBodies,1);
CentrumMIs = zeros(numBodies,1);
HaemalMIs = zeros(numBodies,1);
NeuralMIs = zeros(numBodies,1);
VertebralISs = zeros(numBodies,1);
CentrumISs = zeros(numBodies,1);
HaemalISs = zeros(numBodies,1);
NeuralISs = zeros(numBodies,1);
VertebralMTs = zeros(numBodies,1);
CentrumMTs = zeros(numBodies,1);
HaemalMTs = zeros(numBodies,1);
NeuralMTs = zeros(numBodies,1);
VertebralTSs = zeros(numBodies,1);
CentrumTSs = zeros(numBodies,1);
HaemalTSs = zeros(numBodies,1);
NeuralTSs = zeros(numBodies,1);
% MuscleEccentricity = zeros(numBodies,1);
% MuscleMinAxLength = zeros(numBodies,1);
% MuscleMajAxLength = zeros(numBodies,1);
% MuscleArea = zeros(numBodies,1);
% MusclePerimiter = zeros(numBodies,1);
% MuscleMIs = zeros(numBodies,1);
% MuscleISs = zeros(numBodies,1);
CentrumLength = zeros(numBodies,1);

for n = 1:size(VertebralImages,2)

    CentrumLength(n) = sqrt((mean(fullxs(n,:)) - mean(fullxs(n+1,:)))^2+...
        (mean(fullzs(n,:)) - mean(fullzs(n+1,:)))^2+(ylist(n)-ylist(n+1))^2);
    wholeImg = ismember(labelImg, vertIDs{n});
    
    vertbounds = zeros(1,6);
    slicevert = zeros(1,size(wholeImg,1));
    for r = 1:size(wholeImg,1)
        slicevert(r) = sum(any(squeeze(wholeImg(r,:,:))));
    end 
    ia = find(slicevert);
    vertbounds(1) = min(ia);
    vertbounds(4) = max(ia);
    
    slicevert = zeros(1,size(wholeImg,2));
    for r = 1:size(wholeImg,2)
        slicevert(r) = sum(any(squeeze(wholeImg(:,r,:)==1)));
    end 
    ia = find(slicevert);
    vertbounds(2) = min(ia);
    vertbounds(5) = max(ia);
    
    slicevert = zeros(1,size(wholeImg,3));
    for r = 1:size(wholeImg,3)
        slicevert(r) = sum(any(squeeze(wholeImg(:,:,r))));
    end 
    ia = find(slicevert);
    vertbounds(3) = min(ia);
    vertbounds(6) = max(ia);
  
    box = wholeImg(vertbounds(1):vertbounds(4),vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6));
    subImg = graynewim(vertbounds(1):vertbounds(4),vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6));
    subImg(~box) = 0; 
    bwSubImg = subImg;
    bwSubImg(bwSubImg ~= 0) = 1;
    VertebralImages{n} = subImg;
    VertebralVolumes(n) = nnz(subImg);
    VertebralMIs(n) = mean(nonzeros(subImg(:)));
    VertebralISs(n) = std(nonzeros(subImg(:)));
    VertebralSAs(n) = nnz(bwperim(bwSubImg, 26));
    
    MIJ.createImage('Vertebrae',uint16(VertebralImages{n}*2^16),true);
    MIJ.run('8-bit');
    MIJ.run('Local Thickness (complete process)','threshold=1');
    thicknessIm = MIJ.getCurrentImage;
    MIJ.run('Close All');
    
    VertebralMTs(n) = mean(nonzeros(thicknessIm));
    VertebralTSs(n) = std(nonzeros(thicknessIm));
    
    centrumxlimits1 = round(fullxs(n,:)-vertbounds(1));
    lowercentrumx1 = centrumxlimits1(1)-2;
    if lowercentrumx1 <= 0
        lowercentrumx1 = 0.5;
    end
    uppercentrumx1 = centrumxlimits1(2)+2;
    if uppercentrumx1 > size(subImg,1)
        uppercentrumx1 = size(subImg,1)+0.5;
    elseif uppercentrumx1 <= 0
        uppercentrumx1 = 0.5;
    end
    
    centrumxlimits2 = round(fullxs(n+1,:)-vertbounds(1));
    lowercentrumx2 = centrumxlimits2(1)-2;
    if lowercentrumx2 <= 0
        lowercentrumx2 = 0.5;
    end
    uppercentrumx2 = centrumxlimits2(2)+2;
    if uppercentrumx2 > size(subImg,1)
        uppercentrumx2 = size(subImg,1)+0.5;
    elseif uppercentrumx2 <= 0
        uppercentrumx2 = 0.5;
    end
    
    centrumzlimits1 = round(fullzs(n,:)-vertbounds(3));
    lowercentrumz1 = centrumzlimits1(1)-2;
    if lowercentrumz1 <= 0
        lowercentrumz1 = 0.5;
    end
    uppercentrumz1 = centrumzlimits1(2)-2;
    if uppercentrumz1 > size(subImg,3)
        uppercentrumz1 = size(subImg,3)+0.5;
    elseif uppercentrumz1 <= 0
        uppercentrumz1 = 0.5;
    end
    
    centrumzlimits2 = round(fullzs(n+1,:)-vertbounds(3));
    lowercentrumz2 = centrumzlimits2(1)+2;
    if lowercentrumz2 <= 0
        lowercentrumz2 = 0.5;
    end
    uppercentrumz2 = centrumzlimits2(2)+2;
    if uppercentrumz2 > size(subImg,3)
        uppercentrumz2 = size(subImg,3)+0.5;
    elseif uppercentrumz2 <= 0
        uppercentrumz2 = 0.5;
    end
    
    centrummask = repmat(poly2mask([lowercentrumz1 uppercentrumz1 uppercentrumz2 lowercentrumz2 lowercentrumz1]...
        ,[lowercentrumx1 uppercentrumx1 uppercentrumx2 lowercentrumx2 lowercentrumx1],size(subImg,1),size(subImg,3)),...
        [1,1,size(subImg,2)]);
    centrummask = flip(permute(centrummask, [1 3 2]),2);
    subcentrum = subImg;
    subcentrum(~centrummask) = 0;
    bwSubCentrum = subcentrum;
    bwSubCentrum(bwSubCentrum ~= 0) = 1;
    props = regionprops(bwSubCentrum, 'Centroid');
    centrumcentroid = props.Centroid;
    radius = mean([abs(diff(fullxs(n+1,:))) abs(diff(fullxs(n,:)))])/2+5;
    yel1 = round(centrumcentroid(1)-radius);
    yel2 = round(centrumcentroid(1)+radius);
    if yel1 <= 0
        yel1 = 1;
    end
    if yel2 > size(subImg,2)
        yel2 = size(subImg,2);
    end
    
    centrumImage = subcentrum(:,yel1:yel2,:);
    centrumThicknessImage = thicknessIm;
    centrumThicknessImage(~centrummask) = 0;
    centrumThicknessImage = centrumThicknessImage(:,yel1:yel2,:);
    CentrumMTs(n) = mean(nonzeros(centrumThicknessImage));
    CentrumTSs(n) = std(nonzeros(centrumThicknessImage));
    
    bwCentrumImage = centrumImage;
    bwCentrumImage(bwCentrumImage ~= 0) = 1;
    
    centrum2whole = zeros(size(subcentrum));
    centrum2whole(:,yel1:yel2,:) = bwCentrumImage;
    centrum2whole(centrum2whole == 1) = n;
    
    FishSegFishImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6)) = FishSegFishImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6))+ centrum2whole;
    
    centrum2whole(centrum2whole == n) = 1;
    FishSegVertImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6)) = FishSegVertImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6))+ centrum2whole;
    
    bounds = regionprops(bwCentrumImage, 'BoundingBox');
    bounds = bounds.BoundingBox;
    bounds = [floor(bounds(1:3)) round(bounds(4:6))+1];
    maxbounds1 = bounds(2)+bounds(5);
    maxbounds2 = bounds(1)+bounds(4);
    maxbounds3 = bounds(3)+bounds(6);
    
    if bounds(1) <= 0
        bounds(1) = 1;
    end
    if bounds(2) <= 0
        bounds(2) = 1;
    end
    if bounds(3) <= 0
        bounds(3) = 1;
    end
    if maxbounds1 > size(centrumImage,1);
        maxbounds1 = size(centrumImage,1);
    end
    if maxbounds2 > size(centrumImage,2);
        maxbounds2 = size(centrumImage,2);
    end
    if maxbounds3 > size(centrumImage,3);
        maxbounds3 = size(centrumImage,3);
    end
    CentrumImages{n} = centrumImage(bounds(2):maxbounds1,bounds(1):maxbounds2,bounds(3):maxbounds3);
    CentrumVolumes(n) = nnz(centrumImage);
    CentrumSAs(n) = nnz(bwperim(bwCentrumImage, 26));
    CentrumMIs(n) = mean(nonzeros(centrumImage(:)));
    CentrumISs(n) = std(nonzeros(centrumImage(:)));
    
    neuralmask = repmat(poly2mask([0.5 0.5 size(subImg,3)+0.5 size(subImg,3)+0.5 0.5]...
        ,[0.5 lowercentrumx1+0.5 lowercentrumx2+0.5 0.5 0.5],size(subImg,1),size(subImg,3)),...
        [1,1,size(subImg,2)]);
    neuralmask = flip(permute(neuralmask, [1 3 2]),2);
    
    neuralImage = subImg;
    neuralImage(~neuralmask) = 0;
    
    neuralThicknessImage = thicknessIm;
    neuralThicknessImage(~neuralmask) = 0;
    NeuralMTs(n) = mean(nonzeros(neuralThicknessImage));
    NeuralTSs(n) = std(nonzeros(neuralThicknessImage));
    
    bwNeuralImage = neuralImage;
    bwNeuralImage(bwNeuralImage ~= 0) = 1;
    
    vertebralRegion = FishSegFishImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6));
    vertebralRegion(vertebralRegion~=0) = 1;
    
    overlaps = find(bwNeuralImage + vertebralRegion == 2);
    
    neural2whole = bwNeuralImage;
    neural2whole(neural2whole == 1) = n;
    
    neural2whole(overlaps) = 0;
    
    FishSegFishImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6)) = ...
            FishSegFishImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6))+neural2whole;
        
    neural2whole(neural2whole == n) = 2;
    
    neural2whole(overlaps) = 0;
    
    FishSegVertImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6)) = ...
            FishSegVertImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6)) + neural2whole;
        
    bounds = regionprops(bwNeuralImage, 'BoundingBox');
    
    if isempty(bounds)
        NeuralArchImages{n} = [];
        NeuralVolumes(n) = 0;
        NeuralSAs(n) = 0;
        NeuralMIs(n) = 0;
        NeuralISs(n) = 0;
    else
        bounds = bounds.BoundingBox;
        bounds = [floor(bounds(1:3)) round(bounds(4:6))+1];
        maxbounds1 = bounds(2)+bounds(5);
        maxbounds2 = bounds(1)+bounds(4);
        maxbounds3 = bounds(3)+bounds(6);

        if bounds(1) <= 0
            bounds(1) = 1;
        end
        if bounds(2) <= 0
            bounds(2) = 1;
        end
        if bounds(3) <= 0
            bounds(3) = 1;
        end
        if maxbounds1 > size(neuralImage,1);
            maxbounds1 = size(neuralImage,1);
        end
        if maxbounds2 > size(neuralImage,2);
            maxbounds2 = size(neuralImage,2);
        end
        if maxbounds3 > size(neuralImage,3);
            maxbounds3 = size(neuralImage,3);
        end
        NeuralArchImages{n} = neuralImage(bounds(2):maxbounds1,bounds(1):maxbounds2,bounds(3):maxbounds3);
        NeuralVolumes(n) = nnz(neuralImage);
        NeuralSAs(n) = nnz(bwperim(bwNeuralImage, 26));
        NeuralMIs(n) = mean(nonzeros(neuralImage(:)));
        NeuralISs(n) = std(nonzeros(neuralImage(:)));
    end
    
    haemalImage = subImg;
    croppedcentrummask = centrummask; croppedcentrummask(:,1:yel1-1,:) = 0;
    croppedcentrummask(:,yel2+1:end,:) = 0;
    haemalImage(logical(croppedcentrummask+neuralmask)) = 0;
    haemalImage(1:maxbounds1+1,:,:) = 0;
    haemalThicknessImage = thicknessIm;
    haemalThicknessImage(logical(croppedcentrummask+neuralmask)) = 0;
    haemalThicknessImage(1:maxbounds1+1,:,:) = 0;    
    HaemalMTs(n) = mean(nonzeros(haemalThicknessImage));
    HaemalTSs(n) = std(nonzeros(haemalThicknessImage));
    bwHaemalImage = haemalImage;
    bwHaemalImage(bwHaemalImage ~= 0) = 1;
    
    vertebralRegion = FishSegFishImage(vertbounds(1):vertbounds(4),...
    vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6));
    vertebralRegion(vertebralRegion~=0) = 1;

    overlaps = find(bwHaemalImage + vertebralRegion == 2);

    haemal2whole = bwHaemalImage;
    haemal2whole(haemal2whole == 1) = n;
    
    haemal2whole(overlaps) = 0;
    
    FishSegFishImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6)) = ...
            FishSegFishImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6))+haemal2whole;
            
    haemal2whole(haemal2whole == n) = 3;
    
    haemal2whole(overlaps) = 0;
    
    FishSegVertImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6)) = ...
            FishSegVertImage(vertbounds(1):vertbounds(4),...
        vertbounds(2):vertbounds(5),vertbounds(3):vertbounds(6))+haemal2whole;
        
    bounds = regionprops(bwHaemalImage, 'BoundingBox');
    if isempty(bounds)
        HaemalArchImages{n} = [];
        HaemalVolumes(n) = 0;
        HaemalSAs(n) = 0;
        HaemalMIs(n) = 0;
        HaemalISs(n) = 0;
    else
        bounds = bounds.BoundingBox;

        bounds = [floor(bounds(1:3)) round(bounds(4:6))+1];
        maxbounds1 = bounds(2)+bounds(5);
        maxbounds2 = bounds(1)+bounds(4);
        maxbounds3 = bounds(3)+bounds(6);

        if bounds(1) <= 0
            bounds(1) = 1;
        end
        if bounds(2) <= 0
            bounds(2) = 1;
        end
        if bounds(3) <= 0
            bounds(3) = 1;
        end
        if maxbounds1 > size(haemalImage,1);
            maxbounds1 = size(haemalImage,1);
        end
        if maxbounds2 > size(haemalImage,2);
            maxbounds2 = size(haemalImage,2);
        end
        if maxbounds3 > size(haemalImage,3);
            maxbounds3 = size(haemalImage,3);
        end
        HaemalArchImages{n} = haemalImage(bounds(2):maxbounds1,bounds(1):maxbounds2,bounds(3):maxbounds3);
        HaemalVolumes(n) = nnz(haemalImage);
        HaemalSAs(n) = nnz(bwperim(bwHaemalImage, 26));
        HaemalMIs(n) = mean(nonzeros(haemalImage(:)));
        HaemalISs(n) = std(nonzeros(haemalImage(:)));
    end
%     
%     regionsToAnalyze = vertIDs{n};
%     kIndex = round(mean([mean(fullzs(n,:)) mean(fullzs(n+1,:))]));
%     fishSlice = squeeze(grayim2(:,:,kIndex));
%     fishSlice(fishSlice == 0.5) = 0;
% 
%     MIJ.createImage('Muscle',uint16(fishSlice*2^16),true);
%     MIJ.run('Auto Threshold',['method=' Algorithm ' ignore_black ignore_white white']);
%     threshindxs = MIJ.getCurrentImage;
%     MIJ.run('Close All')
%     threshindxs(threshindxs~=0) = 1;
%     fishSlice(~logical(threshindxs)) = 0;
%     
%     muscleThreshImg = fishSlice;
%     
%     muscleImgSliceRegions = bwconncomp(muscleThreshImg); %regions is a structure containing connected regions
%     muscleImgSlicePixelIDs = muscleImgSliceRegions.PixelIdxList; %compute pixel index list
%     muscleImgSliceNumPixels = cellfun(@numel,muscleImgSlicePixelIDs);
%     [~,muscleImgSliceMaxRegion] = max(muscleImgSliceNumPixels);
%     
%     if isempty(muscleImgSliceMaxRegion)
%         MuscleEccentricity(n) = 0;
%         MuscleMinAxLength(n) = 0;
%         MuscleMajAxLength(n) = 0;
%         MuscleArea(n) = 0;
%         MusclePerimiter(n) = 0;
%         MuscleMIs(n) = 0;
%         MuscleISs(n) = 0;
%     else
%         imgIndex = ones(size(fishSlice,1),size(fishSlice,2),size(fishSlice,3));
%         imgIndex(muscleImgSlicePixelIDs{muscleImgSliceMaxRegion}) = 0;
%         muscleThreshImg(logical(imgIndex)) = 0;
% 
%         B = bwboundaries(muscleThreshImg,'noholes');
%         softImgPerimNumPixels = cellfun(@numel,B);
%         [~,softImgPerimMax] = max(softImgPerimNumPixels);
%         
%         perim = B{softImgPerimMax};
%         imgIndex = ones(size(fishSlice,1),size(fishSlice,2),size(fishSlice,3));
%         imgIndex(sub2ind(size(imgIndex),perim(:,1),perim(:,2))) = 0;
%         scaleSlice = imdilate(~logical(imgIndex),ones(8,8));
%         
%         muscleThreshImg(muscleThreshImg > threshold) = 0;
% 
%         muscleThreshImg(~((muscleThreshImg~=0) - scaleSlice)) = 0;
%         MuscleMIs(n) = mean(nonzeros(muscleThreshImg));
%         MuscleISs(n) = std(nonzeros(muscleThreshImg));
%         
%         muscleThreshImg(muscleThreshImg~=0) = 1;
%         
%         muscleImgSliceRegions = bwconncomp(muscleThreshImg); %regions is a structure containing connected regions
%         muscleImgSlicePixelIDs = muscleImgSliceRegions.PixelIdxList; %compute pixel index list
%         muscleImgSliceNumPixels = cellfun(@numel,muscleImgSlicePixelIDs);
%         [~,muscleImgSliceMaxRegion] = max(muscleImgSliceNumPixels);
%         imgIndex = ones(size(fishSlice,1),size(fishSlice,2),size(fishSlice,3));
%         imgIndex(muscleImgSlicePixelIDs{muscleImgSliceMaxRegion}) = 0;
%         muscleThreshImg(logical(imgIndex)) = 0;
%         
%         fishSlice(~logical(muscleThreshImg)) = 0;
%         MuscleImages{n} = fishSlice;
%         
%         softThreshSliceProps = regionprops(muscleThreshImg,'Eccentricity',...
%             'MinorAxisLength','MajorAxisLength','Area','Perimeter','Orientation');
% 
%         MuscleEccentricity(n) = softThreshSliceProps.Eccentricity;
%         MuscleMinAxLength(n) = softThreshSliceProps.MinorAxisLength;
%         MuscleMajAxLength(n) = softThreshSliceProps.MajorAxisLength;
%         MuscleArea(n) = softThreshSliceProps.Area;
%         MusclePerimiter(n) = softThreshSliceProps.Perimeter;
%         
%     end
    
end

resizefactor = 0.5;
resI = imresize(FishSegFishImage,sqrt(resizefactor),'nearest');
resI = imresize(flip(permute(resI, [1 3 2]),3),sqrt(resizefactor),'nearest');
resI = imresize(flip(permute(resI, [2 1 3]),2),sqrt(resizefactor),'nearest');
resI = flip(permute(resI, [2 1 3]),1);
resI = flip(permute(resI, [1 3 2]),2);

newFishSegFishImage = resI;

resI = imresize(FishSegVertImage,sqrt(resizefactor),'nearest');
resI = imresize(flip(permute(resI, [1 3 2]),3),sqrt(resizefactor),'nearest');
resI = imresize(flip(permute(resI, [2 1 3]),2),sqrt(resizefactor),'nearest');
resI = flip(permute(resI, [2 1 3]),1);
resI = flip(permute(resI, [1 3 2]),2);

newFishSegVertImage = resI;

resI = imresize(grayim2,sqrt(resizefactor),'nearest');
resI = imresize(flip(permute(resI, [1 3 2]),3),sqrt(resizefactor),'nearest');
resI = imresize(flip(permute(resI, [2 1 3]),2),sqrt(resizefactor),'nearest');
resI = flip(permute(resI, [2 1 3]),1);
resI = flip(permute(resI, [1 3 2]),2);

newfullim = resI;

R = newfullim; G = newfullim; B = newfullim;

ColorPattern1 = [1 0 0; 0 1 0; 0 0 1];
ColorPattern2 = [1 0.7812 0.4975 ; 1 1 0; 0 1 1];

%ColorPattern2 = [2/3 1/3 1/2 ; 1/3 1/3 0;0 1 1];

Ri = newFishSegVertImage == 1;
Gi = newFishSegVertImage == 2;
Bi = newFishSegVertImage == 3;

colormipR = squeeze(max(newfullim,[],2));
colormipG = squeeze(max(newfullim,[],2));
colormipB = squeeze(max(newfullim,[],2));

for i = 1:size(VertebralImages,2)
    
    if rem(i,2) == 1
        clmap = ColorPattern1;
    else
        clmap = ColorPattern2;
    end

    R(newFishSegFishImage == i & Ri) = clmap(1,1);
    G(newFishSegFishImage == i & Ri) = clmap(1,2);
    B(newFishSegFishImage == i & Ri) = clmap(1,3);
    
    R(newFishSegFishImage == i & Gi) = clmap(2,1);
    G(newFishSegFishImage == i & Gi) = clmap(2,2);
    B(newFishSegFishImage == i & Gi) = clmap(2,3);

    R(newFishSegFishImage == i & Bi) = clmap(3,1);
    G(newFishSegFishImage == i & Bi) = clmap(3,2);
    B(newFishSegFishImage == i & Bi) = clmap(3,3);
    
    Rind = logical(squeeze(max(newFishSegFishImage == i & Ri,[],2)));
    colormipR(Rind) = clmap(1,1);
    colormipG(Rind) = clmap(1,2);
    colormipB(Rind) = clmap(1,3);

    Gind = logical(squeeze(max(newFishSegFishImage == i & Gi,[],2)));
    colormipR(Gind) = clmap(2,1);
    colormipG(Gind) = clmap(2,2);
    colormipB(Gind) = clmap(2,3);

    Bind = logical(squeeze(max(newFishSegFishImage == i & Bi,[],2)));
    colormipR(Bind) = clmap(3,1);
    colormipG(Bind) = clmap(3,2);
    colormipB(Bind) = clmap(3,3);
end

figure(8)
allmip(:,:,1) = colormipR;
allmip(:,:,2) = colormipG;
allmip(:,:,3) = colormipB;
colorHandle = imshow(allmip);


ColoredSegmentation(:,:,:,1) = R;
ColoredSegmentation(:,:,:,2) = G;
ColoredSegmentation(:,:,:,3) = B;

%Save segmented images
d = date;
listing = dir(inPath);
filenames = {listing.name};
numanalysis = cellfun(@(x) strcmp(x(1:min([length(d) length(x)])), d),filenames);
if sum(numanalysis) ~= 0
    type = 0;
    [~,id] = find(numanalysis);
    for i = 1:length(id)
        fname = filenames{id(i)};
        fname(1:length(d)) = [];
        if isempty(fname)
            type = max(type,1);
        elseif numel(find(fname == '(')) == 1 && numel(find(fname == ')')) == 1
            [~,openid] = find(fname == '(');
            [~,closeid] = find(fname == ')');
            if  closeid - openid > 1
                btwstring = fname(openid+1:closeid-1);
                filterid = isstrprop(btwstring,'digit');
                if sum(filterid) == length(filterid)
                    type = max(type,str2double(btwstring)+1);
                end
            end
        end
    end
    
    if type > 0
        d = [d ' (' int2str(type) ')'];
    end
end
mkdir([inPath d]);

CentrumLength = CentrumLength*Resolution;

VertebralVolumes = VertebralVolumes*(Resolution)^3;
CentrumVolumes = CentrumVolumes*(Resolution)^3;
HaemalVolumes = HaemalVolumes*(Resolution)^3;
NeuralVolumes = NeuralVolumes*(Resolution)^3;
VertebralSAs = VertebralSAs*(Resolution)^2;
CentrumSAs = CentrumSAs*(Resolution)^2;
HaemalSAs = HaemalSAs*(Resolution)^2;
NeuralSAs = NeuralSAs*(Resolution)^2;

ma = 2^15;
VertebralTMDs = (VertebralMIs*2*ma-ma)./4096*281.709-195.402;
CentrumTMDs = (CentrumMIs*2*ma-ma)./4096*281.709-195.402;
HaemalTMDs = (HaemalMIs*ma*2-ma)./4096*281.709-195.402;
NeuralTMDs = (NeuralMIs*ma*2-ma)./4096*281.709-195.402;

VertebralISs = VertebralISs*2*ma./4096*281.709;
CentrumISs = CentrumISs*2*ma./4096*281.709;
HaemalISs = HaemalISs*ma*2./4096*281.709;
NeuralISs = NeuralISs*ma*2./4096*281.709;

% MuscleMinAxLength = MuscleMinAxLength*(21);
% MuscleMajAxLength = MuscleMajAxLength*(21);
% MuscleArea = MuscleMajAxLength*(21)^2;
% MusclePerimiter = MusclePerimiter*(21);

VertebralMTs = VertebralMTs * Resolution;
CentrumMTs = CentrumMTs * Resolution;
HaemalMTs = HaemalMTs * Resolution;
NeuralMTs = NeuralMTs * Resolution;

VertebralTSs = VertebralTSs * Resolution;
CentrumTSs = CentrumTSs * Resolution;
HaemalTSs = HaemalTSs * Resolution;
NeuralTSs = NeuralTSs * Resolution;

PropertyMatrix = [VertebralVolumes CentrumVolumes HaemalVolumes NeuralVolumes...
    VertebralSAs CentrumSAs HaemalSAs NeuralSAs VertebralTMDs CentrumTMDs HaemalTMDs NeuralTMDs...
    VertebralISs CentrumISs HaemalISs NeuralISs VertebralMTs CentrumMTs HaemalMTs...
    NeuralMTs VertebralTSs CentrumTSs HaemalTSs NeuralTSs CentrumLength];

NormalizedMatrix = (PropertyMatrix-repmat(mean(PropertyMatrix),[numBodies,1]))./repmat(std(PropertyMatrix),[numBodies,1]);
VertebralNames = cellfun(@(x) ['Vertebrae' int2str(x)], num2cell(1:n), 'UniformOutput', 0);

PropertyTable = table(VertebralVolumes, CentrumVolumes, HaemalVolumes, NeuralVolumes...
    , VertebralSAs, CentrumSAs, HaemalSAs, NeuralSAs, VertebralTMDs, CentrumTMDs, HaemalTMDs,...
    NeuralTMDs, VertebralISs, CentrumISs, HaemalISs, NeuralISs, VertebralMTs, CentrumMTs,...
    HaemalMTs, NeuralMTs, VertebralTSs, CentrumTSs, HaemalTSs, NeuralTSs, ...
    CentrumLength, 'RowNames', VertebralNames);

fileID = fopen([inPath d f 'VertIDS.txt'],'w');
[nrows,~] = size(vertIDs);
for i = 1:nrows
    nelements = length(vertIDs{i,:});
    format = repmat('%d ',[1 nelements]);
    format(end) = [];
    format = [format '\r\n'];
    fprintf(fileID,format,vertIDs{i,:});
end

fclose(fileID);

fileID = fopen([inPath d f 'Thresholds.txt'],'w');
fprintf(fileID,'Bone threshold: %.5f\r\n',threshold);
fprintf(fileID,'T.Multiplier: %.5f\r\n',thresholdMultiplier);
fprintf(fileID,['Algorithm: ' Algorithm '\r\n']);

fclose(fileID);

fileID = fopen([inPath d f 'FishMetadata.txt'],'w');
fprintf(fileID,['Version: ' Version '\r\n']);
fprintf(fileID,[inFileName '\r\n']);
fprintf(fileID,['Rotation angle: ' num2str(rotAngle)]);

fclose(fileID);

fileID = fopen([inPath d f 'UserInputs.txt'],'w');

fprintf(fileID,['VertebralSegData: ' int2str(vloc) '\r\n']);

for i = 1:size(AllCoords,2)
    data2write = AllCoords{i};
    for j = 1:size(data2write,1)
        fprintf(fileID,'%d %d %d %d\r\n',data2write(j,:));
    end
    fprintf(fileID,'\r\n');
end

fclose(fileID);

h = figure(9);
imshow(NormalizedMatrix, 'InitialMagnification','fit');

colormap(jet)
caxis([min(NormalizedMatrix(:)) max(NormalizedMatrix(:))])
colorbar

for i = 1:size(VertebralImages,2)
    text(1,i,['VB' num2str(i) '      '],'HorizontalAlignment','right');
end

parameters = {'Vert Vol', 'Cent Vol', 'Haem Vol', 'Neur Vol'...
    , 'Vert SA', 'Cent SA', 'Haem SA', 'Neur SA', 'Vert TMD', 'Cent TMD', 'Haem TMD',...
    'Neur TMD', 'Vert IS', 'Cent IS', 'Haem IS', 'Neur IS', 'Vert MT', 'Cent MT', 'Haem MT', 'Neur MT', 'Vert TS', 'Cent TS',...
    'Haem TS', 'Neur TS', 'Cent Len'};

for i = 1:length(parameters)
    text(i,1,['       ' parameters{i}],'Rotation',90,'fontsize',6);
end

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(h, [inPath d f inFileName '.jpg']);

writetable(PropertyTable, [inPath d f inFileName '.txt'], 'Delimiter', '\t', 'WriteRowNames', true);

figure(5)
nsip = imshow(squeeze(max(graynewim,[],2)));
saveas(nsip, [inPath d f inFileName 'MIP_Segmented.jpg']);
figure(6)
mipHandle = imshow(omip,'DisplayRange',[0 numRegions],'Border','tight');
map = colorcube;
map = map(1:end-64,:);
map = flipud(repmat(map,round(numRegions/length(map)),1));
map(1,:) = [0 0 0];
colormap(map);
saveas(mipHandle, [inPath d f inFileName 'CC.jpg']);

figure(7)
sipoHandle = imshow(squeeze(max(grayim2,[],2)));
saveas(sipoHandle, [inPath d f inFileName 'Maximum_Projection.jpg']);

tiffccpath = [inPath d f 'CCStack.tif'];
imwrite(squeeze(uint16(labelImg(:,1,:))),tiffccpath);

for i = 2:size(labelImg,2)
    imwrite(squeeze(uint16(labelImg(:,i,:))),tiffccpath,'WriteMode','append');
end

tiffccpath = [inPath d f 'ColorSegStack.tif'];
imwrite(squeeze(ColoredSegmentation(:,1,:,:)),tiffccpath);

for i = 2:size(ColoredSegmentation,2)
    imwrite(squeeze(ColoredSegmentation(:,i,:,:)),tiffccpath,'WriteMode','append');
end

% tiffccpath = [inPath d f 'MuscleSlices.tif'];
% imwrite(squeeze(uint16(MuscleImages{1}*2^16)),tiffccpath);
% 
% for i = 2:length(MuscleImages)
%     imwrite(squeeze(uint16(MuscleImages{i}*2^16)),tiffccpath,'WriteMode','append');
% end

saveas(colorHandle, [inPath d f inFileName 'Color_Projection.jpg']);