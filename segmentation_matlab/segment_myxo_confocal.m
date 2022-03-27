function segment_myxo_confocal(file_name, total_channel_nb, channel_to_segment)

file_name_pre = strcat(file_name(1:end-4),'_',num2str(channel_to_segment));

% Load all images
fprintf('Reading image file %s\n',file_name_pre);
if contains(file_name,'.tif')
    cube_allchan = [];
    info = imfinfo(file_name);
    num_images = numel(info);
    for k = 1:num_images
        A = imread(file_name, k, 'Info', info);
        cube_allchan = cat(3,cube_allchan,A);
    end
elseif contains(file_name,'.nd2')
    cube_allchan = bfopen(file_name);
    cube_allchan = cat(3,cube_allchan{1}{:,1});
end

% Isolate target channel to segment
cube = double(cube_allchan(:,:,channel_to_segment:total_channel_nb:end));
cubeN = cube;

fprintf('\tProcessing images\n');

% Normalize image intensity
for k = 1:size(cubeN,3)
    tmp = cubeN(:,:,k);
    m = mode(tmp(tmp>1));
    tmp = tmp - m;
    tmp = tmp / prctile(tmp(:),99.99);
    tmp(tmp<0) = 0;
    tmp(tmp>1) = 1;
    cubeN(:,:,k) = tmp;
end
cubeN = cubeN-min(cubeN(:));
cubeN = cubeN/max(cubeN(:));

% Apply 3D hessian filter
[e1,e2,e3] = hessianFilter3D(cubeN,2);

fprintf('\tSegmenting images\n');

% Erode and binarize filtered 3D image
J = -e3-e2;
Je = imerode(J,ones(3,3,3));
cutoff = prctile(Je(:),99);
bw = Je>cutoff;
bw = bwlabeln(bw);

% Open BW 3D image and remove small objects
bw = imopen(bw,ones(2,2,2));
bw = bwareaopen(bw,100);
bw = bwlabeln(bw);
bw = dilateLabels(bw,J,3,cutoff);

fprintf('\tAdjusting segmentation\n');

% Post-processing of segmented 3D objects
cell_stats = regionprops3(bw,'Volume');
volcutoff = median(cell_stats.Volume);

segBW = zeros(size(bw));
objN = 1;

for k = 1:max(bw(:))
    
    obj = bw == k;
    reg = regionprops3(obj);

    if ~isempty(reg)
        % extract object
        bb = round(reg.BoundingBox);
        bb = [max(1,bb(2)-5) min((bb(2)+bb(5)-1)+5,size(bw,1)) ...
            max(1,bb(1)-5) min((bb(1)+bb(4)-1)+5,size(bw,2)) ...
            max(1,bb(3)-5) min((bb(3)+bb(6)-1)+5,size(bw,3))];
        obj = obj(bb(1):bb(2),bb(3):bb(4),bb(5):bb(6));

        if isempty(obj)
            continue;
        end
        
        % find skeleton endpoints
        ep = bwmorph3(bwskel(obj),'endpoints');
        
        % objects with more then 2 endpoints are likely not segmented correctly
        if sum(ep(:)) > 2 && reg.Volume > volcutoff
            epid = find(ep);
            [x,y,z] = ind2sub(size(obj),epid);
            D = pdist([x, y, z]);
            DG = [];
            GDIST = cell(numel(epid)-1,1);
            for i = 1:numel(epid)
                gd = double(bwdistgeodesic(obj,epid(i),'quasi-euclidean'));
                GDIST{i} = round(gd*8)/8;
                DG = [DG gd(epid((i+1:end)))'];
            end
            
            R = D./DG;

            % find straight paths in the skeleton to trim branches and prodduce objects with only 2 endpoints
            Rs = R(R>0.6 & DG>75);
            Rs = sort(Rs,'descend');
            R = squareform(R);
            for j = 1:numel(Rs)
                [a,b] = find(R==Rs(j));
                if ~isempty(a)
                    R(a(1),:) = 0;
                    R(:,b(1)) = 0;
                    R(:,a(1)) = 0;
                    R(b(1),:) = 0;
                    d = GDIST{a(1)} + GDIST{b(1)};
                    objS = abs(d - d(epid(a(1))))<4;
                    mask = objS;
                    [x,y,z] = ind2sub(size(objS),find(mask));
                    x = x + bb(1) - 1;
                    y = y + bb(3) - 1;
                    z = z + bb(5) - 1;
                    ind = sub2ind(size(segBW),x,y,z);
                    segBW(ind) = objN;
                    objN = objN + 1;
                    obj(mask) = 0;
                end
            end
        end
        
        obj = bwlabeln(obj);
        
        % store newly segmented objects
        for j = 1:max(obj(:))
            [x,y,z] = ind2sub(size(obj),find(obj==j));
            x = x + bb(1) - 1;
            y = y + bb(3) - 1;
            z = z + bb(5) - 1;
            ind = sub2ind(size(segBW),x,y,z);
            segBW(ind) = objN;
            objN = objN + 1;
        end
    end
end

% remove small objects
[N,ed] = histcounts(segBW(:),'BinMethod','integers');
ed = ceil(ed);
segBW(ismember(segBW,ed(N<350))) = 0;

% dilate remaining objects
segBW = dilateLabels(segBW,J,3,cutoff);

% calculate object properties
fprintf('\tExtracting cell properties\n');

opts = {
    'Volume'
    'Centroid'
    'EigenValues'
    'EigenVectors'
    'EquivDiameter'
    'Orientation'
    'PrincipalAxisLength'
    'Solidity'
    'SurfaceArea'
    'MaxIntensity'
    'MeanIntensity'
    'MinIntensity'
    'WeightedCentroid'
    };

cell_stats = regionprops3(segBW,cube,opts);
cell_stats.sphericity = pi^(1/3) * (6 * cell_stats.Volume).^(2/3) ./ cell_stats.SurfaceArea;
cell_stats.id = (1:size(cell_stats,1))';
cell_stats = cell_stats(cell_stats.Volume>0,:);
cell_stats.nb_branch_pt = zeros(size(cell_stats,1),1);

% calculate object skeleton properties
segBWskel = zeros(size(segBW));
for i = 1:size(cell_stats,1)
    tmp = segBW == cell_stats.id(i);
    iskel = bwskel(tmp>0);
    skel_branch = bwmorph3(iskel,'branchpoints');
    segBWskel(iskel) = cell_stats.id(i);
    cell_stats.nb_branch_pt(i) = sum(skel_branch(:));
end

segBWskeldist = bwdist(segBW == 0);
segBWskeldist(segBWskel == 0) = 0;

opts = {'Volume', 'MaxIntensity', 'MinIntensity', 'MeanIntensity'};
skel_stats = regionprops3(segBWskel, segBWskeldist, opts);
skel_stats = skel_stats(skel_stats.Volume>0,:);

cell_stats.SkeletonSize = skel_stats.Volume;
cell_stats.MeanWidth    = skel_stats.MeanIntensity;
cell_stats.MaxWidth     = skel_stats.MaxIntensity;
cell_stats.MinWidth     = skel_stats.MinIntensity;


% save processsed and labeled images
fprintf('\tWriting images\n');

outputFileName = strcat(file_name_pre,'_labels.tif');
if exist(outputFileName,'file')
    delete(outputFileName);
end
for K=1:size(segBW,3)
    imwrite(uint16(segBW(:, :, K)), outputFileName, 'WriteMode', 'append');
end

outputFileName = strcat(file_name_pre,'_skeleton.tif');
if exist(outputFileName,'file')
    delete(outputFileName);
end
for K=1:size(segBWskel,3)
    imwrite(uint16(segBWskel(:, :, K)), outputFileName, 'WriteMode', 'append');
end

cubeN = uint16(cubeN*65535);

outputFileName = strcat(file_name_pre,'_norm.tif');
if exist(outputFileName,'file')
    delete(outputFileName);
end
for K=1:size(cubeN,3)
    imwrite(cubeN(:, :, K), outputFileName, 'WriteMode', 'append');
end

J = (-e3-e2)-min(-e3(:)-e2(:));
J = uint16(J/max(J(:))*65535);

outputFileName = strcat(file_name_pre,'_hessian.tif');
if exist(outputFileName,'file')
    delete(outputFileName);
end
for K=1:size(J,3)
    imwrite(J(:, :, K), outputFileName, 'WriteMode', 'append');
end

% calculate signal intensity of each object in the remaining channels
fprintf('\tAnalyzing extra channels\n');

opts = {
    'Volume'
    'MaxIntensity'
    'MeanIntensity'
    'MinIntensity'
    'WeightedCentroid'
    };

extraChan = setdiff(1:total_channel_nb, channel_to_segment);

for k = 1:numel(extraChan)
    
    cubeChan = double(cube_allchan(:,:,extraChan:total_channel_nb:end));
    chan_stats = regionprops3(segBW,cubeChan,opts);
    chan_stats = chan_stats(chan_stats.Volume>0,:);
    cell_stats.(sprintf('channel_%d_MaxIntensity',extraChan(k))) = chan_stats.MaxIntensity;
    cell_stats.(sprintf('channel_%d_MeanIntensity',extraChan(k))) = chan_stats.MeanIntensity;
    cell_stats.(sprintf('channel_%d_MinIntensity',extraChan(k))) = chan_stats.MinIntensity;
    cell_stats.(sprintf('channel_%d_WeightedCentroid',extraChan(k))) = chan_stats.WeightedCentroid;
    
    cubeChan = (1+cubeChan) ./ (1+cube);
    chan_stats = regionprops3(segBW,cubeChan,opts);
    chan_stats = chan_stats(chan_stats.Volume>0,:);
    cell_stats.(sprintf('channel_%d_Ratio_MaxIntensity',extraChan(k))) = chan_stats.MaxIntensity;
    cell_stats.(sprintf('channel_%d_Ratio_MeanIntensity',extraChan(k))) = chan_stats.MeanIntensity;
    cell_stats.(sprintf('channel_%d_Ratio_MinIntensity',extraChan(k))) = chan_stats.MinIntensity;
    cell_stats.(sprintf('channel_%d_Ratio_WeightedCentroid',extraChan(k))) = chan_stats.WeightedCentroid;
    
end

% save summary statistics table
fprintf('\tWriting data table\n');

outputFileName = strcat(file_name_pre,'_cellstats.txt');
if exist(outputFileName,'file')
    delete(outputFileName);
end

writetable(cell_stats,outputFileName);

end