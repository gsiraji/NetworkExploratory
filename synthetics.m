%% synthetics



function [X,GV] = synthetics(disorder, heter, Xfin,Yfin)

X = zeros(Xfin*Yfin,2);

% disorder = 0.1;
% heter = 0.2;
% Xfin = 30;

mm = 0:2:Xfin;
nn = 1:2:Xfin;
counter = 1;
% while n < Xfin
%     if mod(n,2) == 0
for n = nn
        for m = mm
%             if mod(n,2) == 1
            
%             xii = [m,n];
            rndi = rand(1,1);
            rndii = rand(1,1);
            
            if heter > rand(1,1)
                continue
            end
            
            if disorder > rndi
                xi = [n + rndi, m];
            else
                xi = [n,m];
            end
            if m < Yfin
%                 X = [X;  xi];
                X(counter,:) = xi;
                counter = counter + 1;
            end
            
            if disorder > rndii
                xii = [m , n+ rndi];
            else
                xii = [m,n];
            end
            if n < Yfin
%                 X = [X;  xii];
                X(counter,:) = xii;
                counter = counter + 1;
            end
%             X = unique(X);
        end
%         end
%     else 
%         n = n + 1;
    
end
X = X(1:counter-1,:)    ;
% X = X';    
% end
%%
% X = [0 0; 0 1; 0 2; 1 0; 1 2; 2 2; ...
%     1.5 1; 2 0; 2 1; 2.5 1;  ...
%       4 0;4 0.5; 4 1; 4 1.5; 4 2; 5 0; 5 1; ...
%       5 2; 6 0; 6 1; 6 2;...
%       8 0; 9 0; 10 0; 10 0.5; 10 1; 10 2;...
%       8 1; 9 1; 8 1.5; 8 2; 9 2;...
%       12 0; 13 0; 14 0; 14 0.5; 14 1; 14 2;...
%       12 1; 13 1; 12 1.5; 12 2; 13 2];
      
      




[VX,VY] = voronoi(X(:,1),X(:,2));



% h = plot(VX,VY,'-b',X(:,1),X(:,2),'ok'); hold on;
% xlim([-1 Xfin+1])
% ylim([-1 Xfin+1])
%%
% dt = delaunayTriangulation(X);
% [V,R] = voronoiDiagram(dt);
% 
% VV = unique(V(:,1:2), 'rows');
% VV = sort(V);
% 
% VNodes = table;
% VNodes.X = VV(2:end,1);
% VNodes.Y = VV(2:end,2);


GV = graph;
%(zeros(sz1),VNodes);

V1 = [VX(1,:)',VY(1,:)'];
V2 = [VX(2,:)',VY(2,:)'];
% node1 = V1(2,:);
% node2 = V2(2,:);
XY = [[,]; [,]];
% GV = addnode(GV,table(node1(1,1),node1(1,2), 'VariableNames',{'X','Y'}));
% GV = addnode(GV,table(node2(1,1),node2(1,2), 'VariableNames',{'X','Y'}));
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
            GV = addnode(GV,table(node1(1,1),node1(1,2), 'VariableNames',{'X','Y'}));
            node1_idx = numnodes(GV);
        else
            node1_idx = Locb;
        end
        [Lia,Locb] = ismember(node2, XY,'rows');
        if Locb == 0
            XY = [XY; node2];
            GV = addnode(GV,table(node2(1,1),node2(1,2), 'VariableNames',{'X','Y'}));
            node2_idx = numnodes(GV);
        else
            node2_idx = Locb;
        end
        
     else
        XY = [XY; node1];
        GV = addnode(GV,table(node1(1,1),node1(1,2), 'VariableNames',{'X','Y'}));
        node1_idx = numnodes(GV);
        XY = [XY; node2];
        GV = addnode(GV,table(node2(1,1),node2(1,2), 'VariableNames',{'X','Y'}));
        node2_idx = numnodes(GV);
    end
    
    
    if node1_idx ~= node2_idx
        GV = addedge(GV, node1_idx, node2_idx);
    end
end

 


% [V, R] = voronoi(X(1,:), X(2,:)) ;       

% 
% X = [0 1; 0 3;0 5 ; 1 0; ...
%       1 1; 3 1; 5 1;...
%       2 0.5; 4 0.5; 

% dalunaytriangulation