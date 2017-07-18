function FishCuTScriptv1_1(inPath,inFileName,FilterIndex)

%Close the windows from 1 through 9, used by FishCut.
javaaddpath 'C:\Program Files\MATLAB\R2017a\java\jar\mij.jar';

for i = 1:9
    figure(i)
    close(i)
end

[PATHSTR,~,~] = fileparts(mfilename('fullpath'));
fileID = fopen([PATHSTR filesep 'ImagingParameters.txt']);
C = textscan(fileID,'%s %f');
fclose(fileID);

params = C{2};

%%Initialize parameters, 
%Imaging resolution (21 micron or 14 micron used in paper).
Resolution = params(1);
Slope = params(2);
Intercept = params(3);
thresholdMultiplier = params(4);

%Set file separator character based on OS.
f = filesep;
%Version of FishCut saved in metadata file.
Version = 'beta';
%Folder path obtained from homepage.
inPath = [inPath filesep];
%Threshold algorithm used in Miji (Default).
Algorithm = 'Default';
%Manual implementation of rotation.
rotAngle = 0;

%start at first verterbra following webberian apparatus, and count until
%two vertebrae from tail fins 

%%Read in the image
switch FilterIndex
    case 1
        im = squeeze(dicomread([inPath inFileName]));
        dataType = class(im);
    case 2
        dataType = class(imread([inPath inFileName]));
        im_info = imfinfo([inPath inFileName]);
        im(1:im_info(1).Height,1:im_info(1).Width,1:numel(im_info)) = 0;
        for i=1:numel(im_info)
            im(:,:,i) = imread([inPath inFileName],i);
        end
        
end

%Turn the matrix into grayscale values based on data type.
switch dataType
    case 'int16'
        grayim = mat2gray(im,[-2^15 2^15]);
    case 'uint16'
        grayim = mat2gray(im,[0 2^16]);
end

%Store the resulting maximum intensity projection.
sipDummy = max(grayim,[],2);

%Flip the fish to a horizontal axis.
sip = reshape(sipDummy,size(sipDummy,1),size(sipDummy,3));

%Plot the image on a scale from 0 to 1.
figure(1);
imshow(sip,[0 1],'Border','tight');

%Attempt at automatically detecting fish boundary.
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

%Display fish with automatically detected boundaries.
figure(1);
siphandle1 = imshow(outlinedSip);

%Either manually outline fish boundary or proceed with thresholding.
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

%Run Miji if it has not been already run.
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

%Run the default thresholding algorithm (IsoData) in ImageJ.
fishim(fishim==0.5) = 0;
MIJ.createImage('SIP',uint16(fishim*2^16),true);
MIJ.run('Auto Threshold',['method=' Algorithm ' ignore_black ignore_white white']);
threshindxs = MIJ.getCurrentImage;
MIJ.run('Close All')
threshindxs(threshindxs~=0) = 1;
fishim(~logical(threshindxs)) = 0;

switch dataType
    case 'int16'
        threshold = (((2^16*min(fishim(fishim~=0))-2^15)*thresholdMultiplier)+2^15)/2^16;
    case 'uint16'
        threshold = min(fishim(fishim~=0))*thresholdMultiplier;
end


%Turn image into binary image. 
newim = grayim; 
newim(newim < threshold) = 0;
newim(newim ~= 0) = 1; %newim is the original binary image. 

bwImg = newim; %bwImg is the binary image used to compute connectivity.
manualMarks = zeros(size(newim,1),size(newim,2),size(newim,3)); %manualMarks 3D array records the manual marks drawn by the user.

vk = 1; %vk is the container variable for the current vertebrae.
VertebralxCoordsAll = []; %contains the x-coordinates for the vertebral separation lines.
VertebralzCoordsAll = []; %contains the z-coordinates for the vertebral separation lines.
AllCoords = cell(0); %AllCoords stores the record of user input lines.
vloc = 0;
cutspine = 1; %cutspine is a variable for the condition to cut variable.

while cutspine
    
    %compute connectivity
    disp('Computing connectivity...');
    regions = bwconncomp(bwImg,6); %regions is a structure containing connected regions
    voxelIDs = regions.PixelIdxList; %compute voxel index list
    numPixels = cellfun(@numel,voxelIDs); %compute the number of pixels in each component.
    regProps = regionprops(regions,'Centroid'); %regProps is an array of structures with region properties
    centroids = cat(1,regProps.Centroid); %centroids is an array with the centroids of the detected regions, used to sort the locations of vertebrae.
    
    sortedpixelIDs = sortrows([voxelIDs.' num2cell(numPixels.') num2cell(centroids)],2); %sort the voxel index list based on the size of components.
    centroids = cell2mat(sortedpixelIDs(:,3:end)); %
    sortedpixelIDs = sortedpixelIDs(:,1:2);
    
    newregions = regions;
    newregions.PixelIdxList = sortedpixelIDs(:,1).';
    labelImg = labelmatrix(newregions); %compute labelled matrix
    numRegions = regions.NumObjects;

    %compute MIP for label matrix, visualize using custom colormap
    disp('Making labelled projections...'); 
    mipDummy = max(labelImg,[],2); %compute maxima along dimension 2 
    mip = reshape(mipDummy,size(mipDummy,1),size(mipDummy,3)); %reshape mip
    omip = mip; %omip is the image that will be saved near the end of the program.
    %Open a clean figure 3.
    figure(3)
    close(3)
    figure(3)
    map = colorcube(64); %Obtain the default colorcube map values.
    map = map(1:end-36,:); %Remove the similar gradient colors.
    map = [0 0 0; flipud(repmat(map,round(numRegions/length(map)),1))]; %Set the largest component to white, and repeat colors for the number of elements
    
    mipHandle = imshow(mip,'DisplayRange',[0 numRegions],'Border','tight','Colormap',map); %Display the connected components projection.
    
    %compute SIP for segmented image; visualize the image with a bone
    %colormap.
    graynewim = grayim;
    graynewim(~bwImg)=0;
    sipDummy = max(graynewim,[],2);
    sip = reshape(sipDummy,size(sipDummy,1),size(sipDummy,3));
    figure(4)
    close(4)
    figure(4)
    sipHandle = imshow(sip,'DisplayRange',[0 1],'Border','tight');
    colormap(bone);
    
    askinput = 1; %variable for whether input is within constraints.
    while askinput
        remainInLoop1 = input('Enter "s" to segment, or "e" to exit: ','s');
        askinput = ~ismember(remainInLoop1,{'s','e'});
        if askinput
            disp('');
            disp('Make sure to enter either "s" or "e"!')
        end
    end
    
    %Segmentation prompt to user (semi-automatic or manual).
    AutoSegment = 0;
    ManualSegment = 0;
    if remainInLoop1 == 'e'
        break
    else
        askinput = 1;
        while askinput %while loop for correct user input.
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
        bwImg = newim; %work with new image for each entry into semi-automatic mode.
        bwImg(logical(manualMarks)) = 0; %start by applying the user's manual marks.
        
        if vloc ~= 0 %vloc contains the number of user inputs.
            AllCoords(vloc) = [];
        end
        
        vloc = size(AllCoords,2)+1; %vloc contains the number of vertebral locations specified by user.
        
        vk = 0; %vk is iteration variable denoting the vertebrae being used by the program.
        VertebralxCoordsAll = []; %We reset the vertebral coordinate information.
        VertebralzCoordsAll = [];
        %We obtain line coordinates from user
        disp('Enter line segments, press space to finish');        
        remainInLoop2 = 1;
        lineData = []; %lineData is the temporary variable containing a pair of user-defined points.
        numPoints = 0; %numPoints contains the number of formal points entered by user.
        while remainInLoop2 == 1
            figure(4)
            hold on
            [zCoords,xCoords,buttonPress] = ginputc(1,'Color',[1 1 1]); %get x and z line coordinates
            if buttonPress == 1
                numPoints = numPoints+1; %increase the number of points.
                plot(zCoords,xCoords,'*','Color','w','Tag',['Point' num2str(numPoints)]) %plot a white point for user input.
                lineData = [lineData;xCoords zCoords];
                
            elseif buttonPress == 98 %allow user to delete point by pressing "b"
                if size(lineData,1) == 1 %If only one of the user inputs of the pair is entered.
                    lineData = []; %reset the variable containing paired information.
                    tempPoint1 = findobj('Tag',['Point' num2str(vk*2+1)]);
                    delete(tempPoint1) %delete the corresponding white point 
                    numPoints = numPoints - 1; %decreases the formal point count by 1.
                elseif size(lineData,1) == 0 && size(VertebralxCoordsAll,1) > 0 %If lineData is currently empty and there has been a previous entry.
                    VertebralxCoordsAll(vk,:) = []; %reset the stored coordinates of the preceding formal paired points. 
                    VertebralzCoordsAll(vk,:) = [];
                    tempLine = findobj('Tag',['VertebralLine' num2str(vk)]); %delete the white line on the plot.
                    delete(tempLine)
                    tempPoint1 = findobj('Tag',['Point' num2str(vk*2-1)]); %delete the white points on the plot.
                    delete(tempPoint1)
                    tempPoint2 = findobj('Tag',['Point' num2str(vk*2)]);
                    delete(tempPoint2)
                    vk = vk-1; %decreases the vertebral count by 1.
                    numPoints = numPoints - 2; %decreases the formal point count by 2.
                end
                
            elseif buttonPress == 32 %If the user input is space.
                remainInLoop2 = 0; %Escape this loop.
                
            elseif ~ishandle(4) %If the figure has been deleted.
                remainInLoop2 = 0; %Escape this loop.
            
            else %If input is something else, request input again.
                disp('Invalid point section because of extraneous key or mouse press.')
                disp('Please re-enter point or hit "b" to backspace.')
            end

            %If lineData variable has reached this portion of code with
            %size 2,
            if size(lineData,1) == 2
                vk = vk + 1; %increase vertebral count by 1.
                line(lineData(:,2),lineData(:,1),'Tag',['VertebralLine' num2str(vk)],'Color','w'); %draw a white line between the points.
                VertebralxCoordsAll(vk,:) = lineData(:,1); %Formally record the coordinates of the points.
                VertebralzCoordsAll(vk,:) = lineData(:,2);
                lineData = []; %reset the lineData variable.
            end
        end
        
        if size(VertebralxCoordsAll,1) > 1 %If at least one vertebrae has been entered 
            %sort the user-defined lines based on their location.
            ids = VertebralxCoordsAll(:,1) > VertebralxCoordsAll(:,2); 
            VertebralxCoordsAll(ids,:) = fliplr(VertebralxCoordsAll(ids,:));
            VertebralzCoordsAll(ids,:) = fliplr(VertebralzCoordsAll(ids,:));
            %add the user-defined vertebral lines to the global record.
            AllCoords = [AllCoords {[VertebralxCoordsAll VertebralzCoordsAll]}];
            %create the temp variables, fullxs and fullzs.
            fullxs = VertebralxCoordsAll;
            fullzs = VertebralzCoordsAll;
            %take the midpoints of the lines.
            avexCoords = mean(VertebralxCoordsAll,2);
            avezCoords = mean(VertebralzCoordsAll,2);
            
            ylist = zeros(size(VertebralxCoordsAll,1),1);
            %Loop through each user-defined line.
            disp('Running segmentation algorithm. Please wait.');
            for m = 1:size(VertebralxCoordsAll,1)
                [x1,xmin] = min(VertebralxCoordsAll(m,:));%Find the lower point.
                [x2,xmax] = max(VertebralxCoordsAll(m,:));%Find the higher point.

                x1 = floor(x1);%Take the integer indices of input.
                if x1 <= 0
                    x1 = 1;
                end
                x2 = ceil(x2);
                if x2 > size(sipDummy,1)
                    x2 = size(sipDummy,1);
                end
                xBar = floor(avexCoords(m));%Turn the midpoints into integers.
                zBar = floor(avezCoords(m));
                if m == 1%Get the length of the vertebrae and store in dz.
                    dz = avezCoords(2)-avezCoords(1);
                elseif m == size(VertebralxCoordsAll,1)
                    dz = avezCoords(end)-avezCoords(end-1);
                else
                    dz = 0.5*(avezCoords(m+1)-avezCoords(m-1));
                end
                dx = dz;
                dy = dz;
                mdz = VertebralzCoordsAll(m,xmax)-VertebralzCoordsAll(m,xmin); %mdz is the difference in z-coordinates.
                mdx = x2-x1; %mdx is the difference in x-coordinates.
                xcoords = round(x1):round(x2);
                zarray1 = round(mdz/mdx*(xcoords-x1)+VertebralzCoordsAll(m,xmin)); %zarray1 contains the line with the calculated slope and intercept.
                zarray2 = zarray1+1; %zarray2 widens the width of the vertebral boundary line so that adjacent vertebrae will be separate connected components.
                grow = 1; %grow is the boolean as a condition to stay in the while loop.

                while grow 
                    [s1,~,s3] = size(bwImg);
                    %select adjacent vertebrae as separate regions and check boundaries. 
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
                    bwSubImage = bwImg(xel1:xel2,:,zel1:zel2); %bwImg contains the set of two vertebrae.
                    [~,yDummy,~] = ind2sub(size(bwSubImage),find(bwSubImage)); %yDummy contains the y-coordinates of the selected regions.
                    y = floor(mean(yDummy)); %y are the means of the y-coordinates in a single dimension.
                    ylist(m) = mean(yDummy); %ylist will be used later.
                    %make line mask
                    if ~isnan(y)
                        for n = 1:size(bwImg,2)
                            sliceToMask = squeeze(bwImg(:,n,:)); %
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
    %Code to handle the manual segmentation inputs.
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
                numPoints = numPoints+1; %For a click, increase the count of points by 1.
                plot(zCoords,xCoords,'*','Color','w','Tag',['Point' num2str(numPoints)]) 
                lineData = [lineData;xCoords zCoords]; 
                
            elseif buttonPress == 98 %pressing a b, delete one or two points.
                if size(lineData,1) == 1 %if size of the point is 1.
                    lineData = []; %reset the value of lineData.
                    tempPoint1 = findobj('Tag',['Point' num2str(k*2+1)]); %erase the single point for the plot.
                    delete(tempPoint1)
                    numPoints = numPoints - 1; %decrease the count of the formal points by 1.
                elseif size(lineData,1) == 0 && size(xCoordsAll,1) > 0 %If the previous entry contained 2 points and there existed a previous line.
                    xCoordsAll(k,:) = []; %reset the preceding stored data.
                    zCoordsAll(k,:) = [];
                    tempLine = findobj('Tag',['Line' num2str(k)]); %erase the line from the figure.
                    delete(tempLine)
                    tempPoint1 = findobj('Tag',['Point' num2str(k*2-1)]); %erase the two points from the figure.
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
                for n = 1:size(manualMarks,2)
                    sliceToMask = squeeze(manualMarks(:,n,:));
                    sliceToMask([sub2ind(size(sliceToMask),xcoords,zarray1)' sub2ind(size(sliceToMask),xcoords,zarray2)']) = 1;
                    %set the plane of the line equal to 0.
                    manualMarks(:,n,:) = sliceToMask;
                end
                bwImg(logical(manualMarks)) = 0;
            end
        end
    end
end

%Start selecting connected components (CCs) from the 2D CC projection.
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
                disp(['Specifying vertebral body ',num2str(n),'/',num2str(size(vertIDs,1))]);
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
                    for sn = 1:size(manualMarks,2)
                        sliceToMask = squeeze(manualMarks(:,sn,:));
                        sliceToMask([sub2ind(size(sliceToMask),xcoords,zarray1)' sub2ind(size(sliceToMask),xcoords,zarray2)']) = 1;
                        manualMarks(:,sn,:) = sliceToMask;
                    end
                    bwImg(logical(manualMarks)) = 0;
                end

                regions = bwconncomp(bwImg,6); %regions is a structure containing connected regions
                voxelIDs = regions.PixelIdxList; %compute pixel index list
                numPixels = cellfun(@numel,voxelIDs);    
                regProps = regionprops(regions,'Centroid'); %regProps is an array of structures with region properties
                centroids = cat(1,regProps.Centroid);

                sortedpixelIDs = sortrows([voxelIDs.' num2cell(numPixels.') num2cell(centroids)],2);
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
                map = colorcube(64); %Obtain the default colorcube map values.
                map = map(1:end-36,:); %Remove the similar gradient colors.
                map = [0 0 0; flipud(repmat(map,round(numRegions/length(map)),1))]; %Set the largest component to white, and repeat colors for the number of elements
                mipHandle = imshow(mip,'DisplayRange',[0 numRegions],'Border','tight','Colormap',map); %Display the connected components projection.
                set(mipHandle,'cdata',mip)
                break
            end
        end
        set(mipHandle,'cdata',mip)
    end
end

%Tissue Properties
grayim2 = grayim;

FishSegFishImage = zeros(size(im));
FishSegVertImage = zeros(size(im));

%Isolate Vertebral parts
numBodies = size(vertIDs,1);
VertebralImages{numBodies} = [];
CentrumImages{numBodies} = [];
HaemalArchImages{numBodies} = [];
NeuralArchImages{numBodies} = [];
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
CentrumLength = zeros(numBodies,1);

%this for loop goes through each vertebrae and calculates the parameters.
%the loop is structured so that the centrum is delineated by user input, 
%the neural arch is above the centrum, and the haemal arch is the rest of
%the image.
for n = 1:size(VertebralImages,2)

    CentrumLength(n) = sqrt((mean(fullxs(n,:)) - mean(fullxs(n+1,:)))^2+...
        (mean(fullzs(n,:)) - mean(fullzs(n+1,:)))^2+(ylist(n)-ylist(n+1))^2);
    %ylist was saved earlier as the centroids of the centra.
    
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
    if maxbounds1 > size(centrumImage,1)
        maxbounds1 = size(centrumImage,1);
    end
    if maxbounds2 > size(centrumImage,2)
        maxbounds2 = size(centrumImage,2);
    end
    if maxbounds3 > size(centrumImage,3)
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
        if maxbounds1 > size(neuralImage,1)
            maxbounds1 = size(neuralImage,1);
        end
        if maxbounds2 > size(neuralImage,2)
            maxbounds2 = size(neuralImage,2);
        end
        if maxbounds3 > size(neuralImage,3)
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
        if maxbounds1 > size(haemalImage,1)
            maxbounds1 = size(haemalImage,1);
        end
        if maxbounds2 > size(haemalImage,2)
            maxbounds2 = size(haemalImage,2);
        end
        if maxbounds3 > size(haemalImage,3)
            maxbounds3 = size(haemalImage,3);
        end
        HaemalArchImages{n} = haemalImage(bounds(2):maxbounds1,bounds(1):maxbounds2,bounds(3):maxbounds3);
        HaemalVolumes(n) = nnz(haemalImage);
        HaemalSAs(n) = nnz(bwperim(bwHaemalImage, 26));
        HaemalMIs(n) = mean(nonzeros(haemalImage(:)));
        HaemalISs(n) = std(nonzeros(haemalImage(:)));
    end
end

%Since imresize acts in 2-dimensions. Resize each dimension by a factor
%of (rf)^(1/2), so that the color labeling will run faster.
resizefactor = 0.5;
newrf = nthroot(resizefactor,2);

resI = imresize(flip(permute(FishSegFishImage, [2 1 3]),2),newrf,'nearest');
resI = imresize(flip(permute(resI, [1 3 2]),3),newrf,'nearest');
resI = flip(permute(resI, [2 1 3]),2);
resI = imresize(flip(permute(resI, [1 3 2]),3),newrf,'nearest');

resI = flip(permute(resI, [1 3 2]),2);
resI = flip(permute(resI, [2 1 3]),1);
resI = flip(permute(resI, [1 3 2]),2);
resI = flip(permute(resI, [2 1 3]),1);

newFishSegFishImage = resI;

resI = imresize(flip(permute(FishSegVertImage, [2 1 3]),2),newrf,'nearest');
resI = imresize(flip(permute(resI, [1 3 2]),3),newrf,'nearest');
resI = flip(permute(resI, [2 1 3]),2);
resI = imresize(flip(permute(resI, [1 3 2]),3),newrf,'nearest');

resI = flip(permute(resI, [1 3 2]),2);
resI = flip(permute(resI, [2 1 3]),1);
resI = flip(permute(resI, [1 3 2]),2);
resI = flip(permute(resI, [2 1 3]),1);


newFishSegVertImage = resI;

resI = imresize(flip(permute(grayim2, [2 1 3]),2),newrf,'nearest');
resI = imresize(flip(permute(resI, [1 3 2]),3),newrf,'nearest');
resI = flip(permute(resI, [2 1 3]),2);
resI = imresize(flip(permute(resI, [1 3 2]),3),newrf,'nearest');

resI = flip(permute(resI, [1 3 2]),2);
resI = flip(permute(resI, [2 1 3]),1);
resI = flip(permute(resI, [1 3 2]),2);
resI = flip(permute(resI, [2 1 3]),1);

newfullim = resI;

R = newfullim; G = newfullim; B = newfullim;

%These contain the RBG scheme for coloring fish segmentation.
ColorPattern1 = [1 0 0; 0 1 0; 0 0 1];
ColorPattern2 = [1 0.7812 0.4975 ; 1 1 0; 0 1 1];

%ColorPattern2 = [2/3 1/3 1/2 ; 1/3 1/3 0;0 1 1];

%1,2,and 3 label the centra, neural arches, and haemal arches.
Ri = newFishSegVertImage == 1;
Gi = newFishSegVertImage == 2;
Bi = newFishSegVertImage == 3;

%colormip contains the 2D version of the 3D colored fish.
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

%Color is added as the 4th dimension for the color tiff. 
ColoredSegmentation(:,:,:,1) = R;
ColoredSegmentation(:,:,:,2) = G;
ColoredSegmentation(:,:,:,3) = B;

%Create a new folder, under where the current fish is saved, containing the analysis outputs. 
%If the fish has been previously analyzed on the same day, create a folder with (n) attached, where n denotes the (n-1)th same-day analysis. 
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

%Convert the grayscale values to the correct units with the given slope,
%intercept, and resolution.
CentrumLength = CentrumLength*Resolution;

VertebralVolumes = VertebralVolumes*(Resolution)^3;
CentrumVolumes = CentrumVolumes*(Resolution)^3;
HaemalVolumes = HaemalVolumes*(Resolution)^3;
NeuralVolumes = NeuralVolumes*(Resolution)^3;
VertebralSAs = VertebralSAs*(Resolution)^2;
CentrumSAs = CentrumSAs*(Resolution)^2;
HaemalSAs = HaemalSAs*(Resolution)^2;
NeuralSAs = NeuralSAs*(Resolution)^2;

switch dataType
    case 'int16'
        VertebralTMDs = (VertebralMIs*2^16-2^15)./4096*Slope+Intercept;
        CentrumTMDs = (CentrumMIs*2^16-2^15)./4096*Slope+Intercept;
        HaemalTMDs = (HaemalMIs*2^16-2^15)./4096*Slope+Intercept;
        NeuralTMDs = (NeuralMIs*2^16-2^15)./4096*Slope+Intercept;
    case 'uint16'
        VertebralTMDs = (VertebralMIs*2^16)./4096*Slope+Intercept;
        CentrumTMDs = (CentrumMIs*2^16)./4096*Slope+Intercept;
        HaemalTMDs = (HaemalMIs*2^16)./4096*Slope+Intercept;
        NeuralTMDs = (NeuralMIs*2^16)./4096*Slope+Intercept;
end
VertebralISs = VertebralISs*2^16./4096*Slope;
CentrumISs = CentrumISs*2^16./4096*Slope;
HaemalISs = HaemalISs*2^16./4096*Slope;
NeuralISs = NeuralISs*2^16./4096*Slope;

VertebralMTs = VertebralMTs * Resolution;
CentrumMTs = CentrumMTs * Resolution;
HaemalMTs = HaemalMTs * Resolution;
NeuralMTs = NeuralMTs * Resolution;

VertebralTSs = VertebralTSs * Resolution;
CentrumTSs = CentrumTSs * Resolution;
HaemalTSs = HaemalTSs * Resolution;
NeuralTSs = NeuralTSs * Resolution;

%Create a matrix to be plotted 2D containing all of the computed parameters.
PropertyMatrix = [VertebralVolumes CentrumVolumes HaemalVolumes NeuralVolumes...
    VertebralSAs CentrumSAs HaemalSAs NeuralSAs VertebralTMDs CentrumTMDs HaemalTMDs NeuralTMDs...
    VertebralISs CentrumISs HaemalISs NeuralISs VertebralMTs CentrumMTs HaemalMTs...
    NeuralMTs VertebralTSs CentrumTSs HaemalTSs NeuralTSs CentrumLength];

%Create a matrix with the normalized values across each vertebrae (in the
%paper the mean and std are calculated from WT clutchmates).
NormalizedMatrix = (PropertyMatrix-repmat(mean(PropertyMatrix),[numBodies,1]))./repmat(std(PropertyMatrix),[numBodies,1]);
%VertebralNames will be the row-labels on the heatmap.
VertebralNames = cellfun(@(x) ['Vertebrae' int2str(x)], num2cell(1:n), 'UniformOutput', 0);

%PropertyTable will be written to the text file.
PropertyTable = table(VertebralVolumes, CentrumVolumes, HaemalVolumes, NeuralVolumes...
    , VertebralSAs, CentrumSAs, HaemalSAs, NeuralSAs, VertebralTMDs, CentrumTMDs, HaemalTMDs,...
    NeuralTMDs, VertebralISs, CentrumISs, HaemalISs, NeuralISs, VertebralMTs, CentrumMTs,...
    HaemalMTs, NeuralMTs, VertebralTSs, CentrumTSs, HaemalTSs, NeuralTSs, ...
    CentrumLength, 'RowNames', VertebralNames);

%Saved VertIDs is necessary for replicating the segmentation with a future
%version of FishCuT.
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

%We save the thresholding parameters in Thresholds.txt.
fileID = fopen([inPath d f 'Thresholds.txt'],'w');
fprintf(fileID,'Bone threshold: %.5f\r\n',threshold);
fprintf(fileID,'T.Multiplier: %.5f\r\n',thresholdMultiplier);
fprintf(fileID,['Algorithm: ' Algorithm '\r\n']);

fclose(fileID);

fileID = fopen([inPath d f 'FishMetadata.txt'],'w');
fprintf(fileID,['Version: ' Version '\r\n']);
fprintf(fileID,[inFileName '\r\n']);
fprintf(fileID,['Rotation angle: ' num2str(rotAngle)]);
fprintf(fileID,['Resolution: ' num2str(params(1))]);
fprintf(fileID,['Slope: ' num2str(params(2))]);
fprintf(fileID,['Intercept ' num2str(params(3))]);

fclose(fileID);

%Write the coordinates of user-inputted lines to a text file.
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

%Visualize the computed measures on a heatmap.
h = figure(9);
imshow(NormalizedMatrix, 'InitialMagnification','fit','ColorMap',jet);
caxis([min(NormalizedMatrix(:)) max(NormalizedMatrix(:))])
colorbar
%Label the rows.
for i = 1:size(VertebralImages,2)
    text(1,i,['VB' num2str(i) '      '],'HorizontalAlignment','right');
end
%Label the columns.
parameters = {'Vert Vol', 'Cent Vol', 'Haem Vol', 'Neur Vol'...
    , 'Vert SA', 'Cent SA', 'Haem SA', 'Neur SA', 'Vert TMD', 'Cent TMD', 'Haem TMD',...
    'Neur TMD', 'Vert IS', 'Cent IS', 'Haem IS', 'Neur IS', 'Vert MT', 'Cent MT', 'Haem MT', 'Neur MT', 'Vert TS', 'Cent TS',...
    'Haem TS', 'Neur TS', 'Cent Len'};

for i = 1:length(parameters)
    text(i,1,['       ' parameters{i}],'Rotation',90,'fontsize',6);
end

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%Save the heatmap to a jpeg.
saveas(h, [inPath d f inFileName '.jpg']);
%write the measures to a text file.
writetable(PropertyTable, [inPath d f inFileName '.txt'], 'Delimiter', '\t', 'WriteRowNames', true);

%save the thresholded and segmented MIP to a text file.
figure(5)
nsip = imshow(squeeze(max(graynewim,[],2)));
saveas(nsip, [inPath d f inFileName 'MIP_Segmented.jpg']);
figure(6)

map = colorcube(64); %Obtain the default colorcube map values.
map = map(1:end-36,:); %Remove the similar gradient colors.
map = [0 0 0; flipud(repmat(map,round(numRegions/length(map)),1))]; %Set the largest component to white, and repeat colors for the number of elements
mipHandle = imshow(omip,'DisplayRange',[0 numRegions],'Border','tight','Colormap',map); %Display the connected components projection.
%save the MIP containing the colored CC's to a jpeg.
saveas(mipHandle, [inPath d f inFileName 'CC.jpg']);

%save the MIP of the original fish to a jpeg.
figure(7)
sipoHandle = imshow(squeeze(max(grayim2,[],2)));
saveas(sipoHandle, [inPath d f inFileName 'Maximum_Projection.jpg']);

%save the 3D CC image to a tiff stack.
tiffccpath = [inPath d f 'CCStack.tif'];
imwrite(squeeze(uint16(labelImg(:,1,:))),tiffccpath);

for i = 2:size(labelImg,2)
    imwrite(squeeze(uint16(labelImg(:,i,:))),tiffccpath,'WriteMode','append');
end

%save the 3D colored image to a tiff stack.
tiffccpath = [inPath d f 'ColorSegStack.tif'];
imwrite(squeeze(ColoredSegmentation(:,1,:,:)),tiffccpath);

for i = 2:size(ColoredSegmentation,2)
    imwrite(squeeze(ColoredSegmentation(:,i,:,:)),tiffccpath,'WriteMode','append');
end

%save the 2D colored fish to a jpeg.
saveas(colorHandle, [inPath d f inFileName 'Color_Projection.jpg']);