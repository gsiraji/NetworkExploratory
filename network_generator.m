Im = imread('Beads.tif');ptm = (1e-6)/1.61; % pixels to meters
V1 = 160000;
%%
Im1 = Im(:,100:400,:); 
%%
Im1 = Im(:,1903:end,:);
%%
[l,w] = size(Im1);

%     BW2 = imbinarize(Im1,graythresh(Im1));
Sns =0.4;
BW = imbinarize((Im1),'adaptive','ForegroundPolarity','dark','Sensitivity',Sns);
%         imshowpair(Im1,BW,'Montage'); 
    se = strel('disk',11);
    closeBW = imclose(BW,se);
%     openBW = imopen(closeBW,se);
%     imshowpair(Im1,closeBW,'Montage'); 
    se = strel('disk',3);
    openBW = imerode(closeBW,se);
%     figure
%     imshowpair(Im1,openBW,'Montage');
    openInv = ~openBW;
    CC = bwconncomp(openInv);
    S = regionprops(CC, 'Centroid');
    %%
%     mask = zeros(size(Im1));
    centers = ones(CC.NumObjects,2);
    for i = 1:CC.NumObjects
%         mask(CC.PixelIdxList{1,i}) = 1;
        centers(i,:) = S(i,:).Centroid;
    end
    %%
    [vx,vy]  = voronoi(centers(:,1),centers(:,2));
    %%
    vx2 = round(vx, -1);
    vy2 = round(vy, -1);
    %%
    G2 = voronoi2graph(vx2,vy2,w,l);
    %%
    G2 = load('Gvoronoi'); G2 = G2.G2;
%     G = voronoi2graph(vx,vy,w,l);
    %     imshow(mask)
    %% 
    G=G2;
    [bin,binsize] = conncomp(G);%,'Type','weak');
    idx = binsize(bin) == max(binsize);
    G = subgraph(G, idx);
    G.Edges.Lengths = findEdgeLength(G);
    eps1 = 11;
    G = pruneG(G);
    G = redNode(G,eps1);
    %%
    dists = bwdist(imcomplement(BW(:,:,1)));
    %%
    G.Edges.Widths = ptm.*FindEdgeProps(G,dists);
    G.Edges.Widths(G.Edges.Widths < 2e-6) = mean(G.Edges.Widths(G.Edges.Widths > 0));
    bk = boundary(G.Nodes.comx,G.Nodes.comy);
    Xvals = G.Nodes.comx(bk);
    [bk_left,~,~] = find(Xvals < 90);
    [bk_right,~,~] = find(Xvals > (max(Xvals)-90));
    bk_left = bk(bk_left);
    bk_right = bk(bk_right);
    C = 1:numnodes(G);
%     C = setdiff(C,[bk_left',bk_right']);
    G = subgraph(G, unique([bk_left',bk_right',C ],'stable'));
    G.Edges.Lengths = ptm.* G.Edges.Lengths;
    G.Nodes.comx = ptm.* G.Nodes.comx;
    G.Nodes.comy = ptm.* G.Nodes.comy;
    %%
    [bin,binsize] = conncomp(G0);%,'Type','weak');
    idx = binsize(bin) == max(binsize);
    G0 = subgraph(G0, idx);

%%
figure()
h = imshow(Im2); set(h, 'AlphaData', 0.2);hold on
p = plot(G);
p.XData = G.Nodes.comx/ptm;
p.YData =G.Nodes.comy/ptm;
%%
figure()
h = imshow(Im1); set(h, 'AlphaData', 0.2);hold on

p = plot(G2);
p.XData = G2.Nodes.comx;
p.YData =G2.Nodes.comy;
%%
figure()
h = imshow(Im2); set(h, 'AlphaData', 0.5);hold on
%%
p = plot(G);
p.XData = G.Nodes.comx;
p.YData =G.Nodes.comy;
% p.EdgeCData = G.Edges.Flows; 
%%
p.NodeCData = G.Nodes.Potentials; 
p.LineWidth = 5;
% set(gca,  'CLim', [0 5].*1e-6);
%% failed attempts to remove white spots
%%
    n = 1;  
    Idouble = im2double(Im1); 
    avg = mean2(Idouble);
    sigma = std2(Idouble);
    %Adjust the contrast based on the standard deviation.
    J = imadjust(Im1,[avg-n*sigma avg+n*sigma],[]);
    %Display the adjusted image.
    imshowpair(J,Im1,'montage')
    %%

mask = imbinarize(Im3, 0.2);
    % Dilate it a bit to make the spots bigger.
    se = strel('line',1,90);
%     mask = imdilate(mask, se);
    rIm = regionfill(Im3, mask);
    imshow(rIm, [])
    %%
    mask = Im3 > 30;
% Find the areas
    props = regionprops(logical(mask), 'Area');
    allAreas = sort([props.Area]);
    %%
% Extract only blobs larger than 25.
    mask = bwareaopen(mask, 25);
    %%
%     mask = imclearborder(mask);

    maskedGrayImage = Im3; % First initialize.
    maskedGrayImage(~mask) = 0;  % Now mask
    imshow(maskedGrayImage);


    %%
    % create the image skeleton from the binary image
    C = bwskel(BW(:,:,1));
    % calculate the distance to the nearest bead from each pixel
    disto = bwdist(imcomplement(BW(:,:,1)));
    

    subplot(2,2,1);
%     imshow(Im1)
    BW2 = imbinarize(Im1(l-202:end,:,:),'adaptive');
        % use the package 'Skel2Graph' to calculate the graph corresponding to
    % the skeleton
    [G0,Nodes0, Links00] = Skel2Graph3D(C,100);
    G0 = graph(G0);
    G0.Nodes.x = [Nodes0.comy]';
    G0.Nodes.y = [Nodes0.comx]';
    div = 1.2;
ThR=graythresh(Im1)/div;
    BW1 = imbinarize(Im1(1:l-201,:),ThR);
    BW = [BW1;BW2];
    %%
        %%
    subplot(2,2,2);
    imshow(BW)
    subplot(2,2,3);
    imshow(BW1)
    subplot(2,2,4);
    imshow(BW2)
    %%
%     figure();imshowpair(Im1,BW1,'Montage'); 

    labels = 1:93;
    labels = num2str(labels');
    labels = cellstr(labels);
    plot([S.Centroid])
    text([S.Centroid],labels)
    %%

function GV = voronoi2graph(VX,VY,Xfin,Yfin)
GV = graph;
%(zeros(sz1),VNodes);

V1 = [VX(1,:)',VY(1,:)'];
V2 = [VX(2,:)',VY(2,:)'];
% node1 = V1(2,:);
% node2 = V2(2,:);
XY = [[,]; [,]];
% GV = addnode(GV,table(node1(1,1),node1(1,2), 'VariableNames',{'comx','comy'}));
% GV = addnode(GV,table(node2(1,1),node2(1,2), 'VariableNames',{'comx','comy'}));
% if node1_idx ~= node2_idx
%         GV = addedge(GV, 1, 2);
% end
for i = 2:1:length(V1)
    node1 = V1(i,:);
    node2 = V2(i,:);
    if node1(1,1) > Xfin+1 || node1(1,2) > Yfin+1 
        continue
    elseif node2(1,1) > Xfin+1 || node2(1,2) > Yfin+1 
        continue
    elseif sum(node1 < -1) > 0
        continue
    elseif  sum(node2 < -1) > 0
        continue
    end
     if length(XY) > 0  

        [Lia,Locb] = ismember(node1, XY,'rows');
        if Locb == 0
            XY = [XY; node1];
            GV = addnode(GV,table(node1(1,1),node1(1,2), 'VariableNames',{'comx','comy'}));
            node1_idx = numnodes(GV);
        else
            node1_idx = Locb;
        end
        [Lia,Locb] = ismember(node2, XY,'rows');
        if Locb == 0
            XY = [XY; node2];
            GV = addnode(GV,table(node2(1,1),node2(1,2), 'VariableNames',{'comx','comy'}));
            node2_idx = numnodes(GV);
        else
            node2_idx = Locb;
        end
        
     else
        XY = [XY; node1];
        GV = addnode(GV,table(node1(1,1),node1(1,2), 'VariableNames',{'comx','comy'}));
        node1_idx = numnodes(GV);
        XY = [XY; node2];
        GV = addnode(GV,table(node2(1,1),node2(1,2), 'VariableNames',{'comx','comy'}));
        node2_idx = numnodes(GV);
    end
    
    
    if node1_idx ~= node2_idx
        GV = addedge(GV, node1_idx, node2_idx);
    end
end
    end

function edgeLength = findEdgeLength(G) 
        edgeLength = zeros(numedges(G),1);
        XY = [(G.Nodes(:,:).comx),(G.Nodes(:,:).comy)]; %faster w XY 
        for i  = 1:numedges(G)
            nlist = G.Edges(i,:).EndNodes;
            edgeLength(i,:) = pdist( [XY(nlist(1),:); XY(nlist(2),:)]);
%             edgeLength(i,:) = pdist( [table2array(G.Nodes(nlist(1),1:2)); table2array(G.Nodes(nlist(2),1:2))]);
        end
end

function G = pruneG(G)

D = degree(G);
Idx = find(D > 1);
G = subgraph(G,Idx);

end

function G = redNode(G,eps1)
% find edges that are too short
% IdxEdges = find(G.Edges.Lengths < eps1);
% find end nodes of the redundant edges
IdxNodes = G.Edges(find(G.Edges.Lengths < eps1),:).EndNodes;
% replace the redundant connections
newEdgeList = [];
for i = 1:length(IdxNodes)
     % find neighbors of the node being removed
    neiList1 = neighbors(G,IdxNodes(i,1));
    % remove the node that is already connected
    neiList1 = neiList1(neiList1 ~= IdxNodes(i,2));
    neiList2 = neighbors(G,IdxNodes(i,2));
    for k = 1:length(neiList1)
        if ismember(neiList1(k),neiList2) == 0
            newEdge = G.Edges(1,:);
            newEdge.EndNodes(:,:) = [IdxNodes(i,2),neiList1(k)];
            newEdgeList = [newEdgeList;newEdge];
        end
    end
end
G = addedge(G, newEdgeList);
G = rmnode(G,IdxNodes(:,1));

end

