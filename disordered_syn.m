%% community.m


%%

%colormap

% psymap = [255 20 40
% 255 47 146
% 255 138 216
% 255 255 255
% 150 230 80
% 255 251 0
% 55 145 230]./255;
% %55 145 230255 255 255
% psymap1 = [255 20 40
% 55 145 230]./255;
% 
% Rval = interp1(0:10:60,psymap(:,1,1),0:1:60);
% Gval = interp1(0:10:60,psymap(:,2,1),0:1:60);
% Bval = interp1(0:10:60,psymap(:,3,1),0:1:60);
% 
% psymap2 = [Rval(:),Gval(:),Bval(:)];
% psymap3 =[Rval(20:60);Gval(20:60);Bval(20:60)]';



%%
graph_mode = 'syn';
if graph_mode == 'exp'

tic
Im = imread('Beads.tif');Im = Im(:,1903:end,:);
Im2 = imread('beadred.png');Im2 = Im2(:,1903:end,:);
V1 = 155000;
% Im = Im(:,6000:end,:); Im2 = Im2(:,6000:end,:);V1 = 20000;

pixels_to_meters = (1e-6)/1.61;ptm = pixels_to_meters;

w = size(Im,1);l = size(Im,2);boundary_limits = [50, l-54];dyn_vis = 60e-3;depth = 1e-3;p_dmeter = 1e-6;
Nstep = 4000;tot_time = 4;dt = tot_time/Nstep; % dt has to be < 0.002

eros_prob = 0.023;
prob_fun = @(y) y/(exp(2).*exp (- y) + 1);


[G1, Nodes1, Links1, dists] = Im2Graph(Im, 1.5);
G0 = graph(G1);
toc

tic

[bin,binsize] = conncomp(G0);
idx = binsize(bin) == max(binsize);
SG = subgraph(G0, idx);
nodez_SG = Nodes1(idx);
SG.Nodes.comx = [nodez_SG.comy]';
SG.Nodes.comy = [nodez_SG.comx]';

idx_B1 =[];  %%%%%%%%%%% separate the boundary and center nodes %%%%%%%%%%
idx_B2 = [];
idx_C = [];


bk = boundary([nodez_SG.comy]',[nodez_SG.comx]');
% bk_left = bk(154:164);  
bk_left = bk(1:10);
idx_B1 = bk_left';
nodeset = 1:1:numnodes(SG);
% nodeset = nodeset(nodeset ~= bk_left');
for nodenum = nodeset
        if ismember(nodenum,bk_left)
            continue
        elseif nodez_SG(nodenum).comy > boundary_limits(2)
               idx_B2 = [idx_B2 nodenum];
        else       
            if degree(SG,nodenum) > 1 
                idx_C = [idx_C nodenum];
            end
        end
end






G = subgraph(SG, [idx_B1 idx_B2 idx_C]);   
nodez_HSG = nodez_SG([idx_B1 idx_B2 idx_C]);

% [EBC, NBC]=edge_betweenness_bin(adjacency(G));


G.Nodes.comx = [nodez_HSG.comy]';
G.Nodes.comy = [nodez_HSG.comx]';
G.Nodes.ID = (1:1:numnodes(G))';
G.Edges.ID = (1:1:numedges(G))';
limit1 = length(idx_B1);
limit2 = length(idx_B2);
Boundary2.left = G.Nodes.ID(1:limit1);
Boundary2.right = G.Nodes.ID(1+limit1:limit1+limit2);
% G.Nodes.comy = ptm.*(w.*ones(length(G.Nodes.comx),1)-[Nodes1.comx]');

G.Edges.Lengths = ptm.*(G.Edges.Weight);
EWidths = FindEdgeProps(G.Edges,G.Nodes,ptm,dists);
ELengths = G.Edges.Lengths;
G.Nodes.comx = ptm.*[nodez_HSG.comy]';
G.Nodes.comy = ptm.*[nodez_HSG.comx]';

% EWidths = randi([3,15],numedges(G),1).*1e-6;
% ELengths = randi([70,140],numedges(G),1).*1e-6;
% ELengths = 1e-5.*ones(numedges(G),1);
EResistances = (12*dyn_vis.*ELengths./((EWidths.^3).*depth.*(1-0.063.*EWidths./depth))); %limit where w>>h
EConds = 1./EResistances;
G.Edges.Widths = EWidths;
% G.Edges.Lengths = ELengths;
G.Edges.Resistances = EResistances;
I = incidence(G);
% D = diag(degree(G).^(-1/2));
% ws = (1/256).*ones(1,m*n);
% ws = [1, 1/4,1];
ws = EConds;
D = (abs(I)*ws).^(-1/2);
D = diag(D);
W = diag(ws);
L = I*W*I';

lb = length(idx_B1)+length(idx_B2);
lL = length(L);
LBB = L(1:lb,1:lb);
LBC = L(1:lb,lb+1:lL);
LCB = L(lb+1:lL,1:lb);
LCC = L(lb+1:lL,lb+1:lL);

psi_B1 = V1.*ones(length(idx_B1),1); %%% set up the voltage vector for boundaries
psi_B2 = zeros(length(idx_B2),1);

psi_B = [psi_B1; psi_B2];
LCCinv = inv(LCC);
LS = LBB - (LBC*(LCCinv*LCB));
JB = LS*psi_B;

J = zeros(lL,1);
J(1:lb) = JB;
psi_V = L\J;
%psi_C = -LCCinv*(LCB*psi_B); %%the solution for the special case of two b%nodes
CV1 = -psi_V(1) + V1;
testtest=CV1.*ones(lL,1);
psi_V = CV1+psi_V; %% getting the unique solution by a shift

G.Nodes.Potentials = psi_V;
G.Edges.Flows = FindEdgeFlows(G.Edges, G.Nodes);
G.Edges.Open = ones(numedges(G),1);
tau_w = abs(G.Edges.Flows)./((G.Edges.Widths.^3));
G.Edges.Shear = tau_w;
local_pot0 = abs(G.Edges.Flows).*(G.Edges.Resistances);

flow_cr = mean(G.Edges.Flows(:,1))./1e14;
p_cr = mean(local_pot0);





beginning_set = 1:1:limit1;

toc


else



    disorder1 = 0;
    heter1 = 0;
    Xfin1 = 100;
    Yfin1 = 20;
%     Xfin1 = 10;
%     Yfin1 = 5;
    
    [X,Gs] = synthetics(disorder1, heter1, Xfin1, Yfin1);
    
    V1 =8000;
    V1 =200;
    idx_B1 =[];  %%%%%%%%%%% separate the boundary and center nodes %%%%%%%%%%
    idx_B2 = [];
    idx_C = [];
    boundary_limits = [1, Xfin1-1];
    
    
    for nodenum = 1:1:numnodes(Gs)
            if Gs.Nodes(nodenum,:).X < boundary_limits(1)
                idx_B1 = [idx_B1 nodenum];
            elseif Gs.Nodes(nodenum,:).X > boundary_limits(2)
                idx_B2 = [idx_B2 nodenum];
            else       
                if degree(Gs,nodenum) > 1 
                    idx_C = [idx_C nodenum];
                end
            end
    end
    
    G = subgraph(Gs, [idx_B1 idx_B2 idx_C]);   
    
    
    G.Nodes.ID = (1:1:numnodes(G))';
    G.Edges.ID = (1:1:numedges(G))';
    limit1 = length(idx_B1);
    limit2 = length(idx_B2);
    Boundary2.left = G.Nodes.ID(1:limit1);
    Boundary2.right = G.Nodes.ID(1+limit1:limit1+limit2);
    
    

    G.Nodes.X = (1e-6).*G.Nodes.X;
    G.Nodes.Y = (1e-6).*G.Nodes.Y;
    for i= 1:1:numedges(G)
        endnode_v = G.Edges.EndNodes(i,:);
        X_v = G.Nodes.X(endnode_v);
        Y_v =  G.Nodes.Y(endnode_v);
        
        G.Edges.Lengths(i) = sqrt((X_v(1) - X_v(2)).^2+(Y_v(1) - Y_v(2)).^2);
        
    end
    
    
    G.Nodes.comx = G.Nodes.X;
    G.Nodes.comy = G.Nodes.Y;
    mask = ones(numedges(G),1);
    
    



    G.Edges.Widths = 1+ 13.*rand(numedges(G),1);
    G.Edges.Widths = (1e-6).*G.Edges.Widths;
    ELengths = G.Edges.Lengths;
    EWidths = G.Edges.Widths;
    
    % EResistances = (12*dyn_vis.*ELengths./((EWidths.^3).*depth.*(1-0.063.*EWidths./depth))); %limit where w>>h
    
    
    
    p_norm1 = 0.5;
    p_dmeter = 1e-6;
    dep_prob = 1000;
     
    I = incidence(G);
    inedges = find(I(Boundary2.left,:));
    avg_capacity = sum(floor(G.Edges.Widths(inedges,:)./p_dmeter),'all');
%     avg_capacity = mean(floor(G.Edges.Widths(find(I(Boundary2.left,:)),:)./p_dmeter)).*length(inedges);

    dyn_vis = 3e-3; %3mpa-s
    depth = 1e-3; %A = 1mm^2
    
%     EResistances = (8/pi*dyn_vis.*ELengths./(EWidths./2).^4);
    EResistances=(12*60e-3.*ELengths./((EWidths.^3).*1e-3.*(1-0.063.*EWidths./1e-3)));
    EConds = 1./EResistances;
    
    G.Edges.Resistances = EResistances;
    G.Edges.Ci = EConds;
    G01 = G;
    % G = find_flows(G, Boundary2, V0, G);
    G = potSolver(G,V1,Boundary2,G01);
    G.Edges.Flows = FindEdgeFlows(G.Edges, G.Nodes);
    
    G.Edges.Open = ones(numedges(G),1);
    tau_w = abs(G.Edges.Flows)./((G.Edges.Widths.^3));
    G.Edges.Shear = tau_w;


end





 





%%
G01 = G;

% for E = 1:1:numedges(G01)
%     i = G01.Edges.EndNodes(E,1);
%     j = G01.Edges.EndNodes(E,2);
%    G01.Edges.Betweenness(E,1) = EBC(i,j);
% end
flow_cr = mean(G.Edges.Flows)./1e9;
tlcf = flow_cr;
r0 = mean(G.Edges.Widths);
tau_cr =median(G.Edges.Shear) ;%1e+02; %10^4;%3.35;%;%
pot_mean = mean(G01.Edges.Flows.*G01.Edges.Resistances);
local_pot0 = abs(G.Edges.Flows).*(G.Edges.Resistances);

% flow_cr = mean(G.Edges.Flows(:,1))./1e14;
p_cr = 5*mean(local_pot0);
k1 = 1;

%% particle_sim



tot_time = 0.001;
Nstep = 150; %dt has to be < 1e-4
dt = tot_time/Nstep;
test_G2 = G;
test_G1 = G;
inj_step = 100;
Nparticle = 100;
% inj_rate = Nparticle/inj_step;

% tot_Nparticle = ceil(tot_time./inj_step)*Nparticle;
tot_Nparticle = ceil(Nstep./inj_step+1)*Nparticle;
% if inj_step > tot_time
if inj_step > Nstep
    tot_Nparticle = Nparticle;
end
Nparticle_active = Nparticle;

k1 = 1;
particle_set = particle.empty(0,tot_Nparticle);
for p_idx = 1:1:tot_Nparticle
    particle_set(p_idx) = particle();
%     if p_idx < Nparticle + 1
%         particle_set(p_idx) = init_particle(particle_set(p_idx),test_G1.Nodes);
        particle_set(p_idx) = init_particle(particle_set(p_idx),test_G2.Nodes(1:limit1,:));
%     end
end


current_node_list = ones(tot_Nparticle, 1);
current_edge_list = zeros(tot_Nparticle, 1);
% Nstep = 800;
%%
make_vid = 1;
frame_k = 1;
if make_vid == 1
    v = VideoWriter('psim400t1200step_40v10_50injstep1000np8k');
    v.Quality = 100;
    v.FrameRate = 10;
    open(v)
end
passed_time = 0;
deposition_v = zeros(Nstep,1);
avg_r = zeros(Nstep,1);
target1 = Boundary2.right;
%%
for step = 1:1:Nstep
    GID = test_G2.Edges.ID;
    
    if current_node_list ~= 0
        if ismember(current_node_list,target1) 
            'done'
            break
        end
    end
    avg_r(step,1) = sum(test_G2.Edges.Widths)./numedges(test_G1);
    if passed_time > 0
%         if mod(passed_time, inj_step) == 0
        if mod(step, inj_step) == 0    
            Nparticle_active = Nparticle_active + Nparticle;%%
%             for p_idx = (Nparticle_active - Nparticle +1):1:Nparticle_active
%                 particle_set(p_idx) = init_particle(particle_set(p_idx),test_G2.Nodes);
% %                 particle_set(p_idx) = init_particle(particle_set(p_idx),test_G2.Nodes(beginning_set,:));
%             end
            if abs(mean(test_G2.Edges.Flows)) < tlcf
                'broken'
                break
        
            elseif sum(test_G2.Edges.Widths)./numedges(test_G1) < r0/2
                break
            end
            
        end
    end
    

  for p_idx = 1:1:Nparticle_active
    
    [particle_set(p_idx),test_G2,test_G1,mask] = next_action(particle_set(p_idx),incidence(test_G2),test_G2.Edges,target1,depth,dt,p_norm1,dep_prob,test_G2.Nodes,p_dmeter,test_G2,mask,test_G1);
%     pall_xy(step,2*p_idx-1:2*p_idx) =  [particle_set(p_idx).comx; particle_set(p_idx).comy];
    current_node_list(p_idx) = particle_set(1,p_idx).node_num;
    current_edge_list(p_idx) = particle_set(1,p_idx).edge_num;
    deposition_v(step,1) = deposition_v(step,1) + (particle_set(1,p_idx).deposited)./Nparticle_active;%%
%     test_G1.Edges.Open(:,1) = mask;
%     test_G2 = evolve_G(test_G2, test_G1,mask);
  end
    GEdges = test_G2.Edges;

    if size(unique(test_G2.Edges.ID),1) ~= size(unique(GID),1)
        GIDtest = unique(GID);
        GIDtest2 = unique(test_G2.Edges.ID);
        test_G2 = potSolver(test_G2,V1,Boundary2,G01);
        test_G2.Edges.Flows = FindEdgeFlows(test_G2.Edges, test_G2.Nodes);
%         test_G2 = find_flows(test_G2, Boundary1, V1,test_G1);
    elseif unique(test_G2.Edges.ID) ~= unique(GID)
        GIDtest = unique(GID);
        GIDtest2 = unique(test_G2.Edges.ID);
%         test_G2 = find_flows(test_G2, Boundary1, V1,test_G1);
         test_G2 = potSolver(test_G2,V1,Boundary2,G01);
         test_G2.Flows = FindEdgeFlows(test_G2.Edges, test_G2.Nodes);
    end

    blocked_edges = find(~mask);
    %
    
    if ~isempty(blocked_edges) 
        for ei = blocked_edges'
            G2edge = find(test_G2.Edges.ID == ei,1);
            if ~isempty(G2edge) 
%                 Gedge = find(test_G1.Edges.ID == ei,1);
            end_nodes = test_G1.Edges.EndNodes(ei,:);
            ei_id = test_G1.Edges.ID(ei,:);
            nodeid = test_G1.Nodes.ID(end_nodes',:);
            endnode1 = test_G2.Nodes(find(test_G2.Nodes.ID(:,1) == nodeid(1),1),:);
            endnode2 = test_G2.Nodes(find(test_G2.Nodes.ID(:,1) == nodeid(2),1),:);
            
            if abs(endnode1.Potentials-endnode2.Potentials) > p_cr
                p_idxx = find(current_edge_list == ei);
                for idxx = p_idxx'
                    particle_set(1,idxx).deposited = 0;
                end
                test_G1.Edges.Open(ei,1) = 1;
                mask(ei,1) = 1;
                new_edge = test_G1.Edges(ei,:);
%                 end_node1 = new_edge.EndNodes(:,1);
%                         end_node1ID = test_G1.Nodes.ID(end_node1);
%                         end_node2 = new_edge.EndNodes(:,2);
%                         end_node2ID = G1.Nodes.ID(end_node2);
                new_edge.EndNodes(:,1)  =find(test_G2.Nodes.ID == nodeid(1,1),1);
                new_edge.EndNodes(:,2) = find(test_G2.Nodes.ID == nodeid(2,1),1);             
                        test_G2 = addedge(test_G2,new_edge);
            end
            
            end
        end
    end
    
    %


%     res2 = res(test_G2.Edges.Open > 0);
    
 %%%%%%%%%%%% vid
    
    if mod(step,frame_k) == 0
        step
    if make_vid == 1
% %         figure(1);
% %             h = plot(X(:,1),X(:,2),'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 9); hold on;
%             p2 = plot(test_G2,'XData',test_G2.Nodes.comx,'YData',test_G2.Nodes.comy);
%             p2.EdgeCData = test_G2.Edges.Flows;% ones(numedges(test_G2),1) ;
%              p2.Marker = 'none';
%             % p2.NodeCData = test_G2.Nodes.Potential;
%             p2.LineWidth = 3;
%             hold on;
% %             for iii = 1:1:Nparticle_active
% %                 plot((particle_set(1,iii).comx),(particle_set(1,iii).comy),'o',...
% %                     'lineWidth',1.5, 'MarkerFaceColor', 'k', 'markerSize', 2);
% %                 hold on;
% %             end
            ff = figure('units','pixels','position',[0 0 1020 540]);
            h = plot(X(:,1).*1e-6,X(:,2).*1e-6,'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 6); hold on;
            for iii = 1:1:Nparticle_active
                if particle_set(1,iii).out == 0
                    if particle_set(1,iii).deposited == 1
                        p_color = 'r';
                        e_color = 'r';
                    else
                        p_color = 'w';
                        e_color = 'k';
                    end
                        plot((particle_set(1,iii).comx),(particle_set(1,iii).comy),'o',...
                            'lineWidth',1, 'MarkerFaceColor', p_color,...
                            'MarkerEdgeColor',e_color, 'markerSize', 2);hold on;
                end
            end
%             xlim([0 Xfin1])
%             ylim([0 Yfin1])
            axis off
% %             % labeledge(p2,1:numedges(test_G2),test_G2.Edges.Flow)
%               colormap(white)
% % %                colorbar
            set(gca, 'CLim', [1 10]);
% % %             title('{G2}')
% % % % 
            F(k1) = getframe(gcf); 
            % 
            % %      for framez = 1:1:pace
            writeVideo(v,F(k1));
            % %      end
            close all ;
            k1 = k1 +1;
    end
    end
%%%%%%%%%%    
    passed_time = dt*step;
    
    
    
end

if make_vid == 1
    close(v)
end

Ni = numedges(test_G2);
Order_par = (1/(Ni-1)).*(Ni - ((sum((test_G2.Edges.Flows).^2,'all'))^2)/(sum((test_G2.Edges.Flows).^4,'all')));



%% 
% make_vid = 1;
% % [bin,binsize] = conncomp(G01);
% if make_vid == 1
%     v = VideoWriter('community10_medtau_zoom_syn.avi');
%     v.Quality = 100;
%     v.FrameRate = 5; 
%     open(v)
% end
% binsize_vec = ones(200,1);
% maxflow_vec = ones(200,1);
% maxind_vec = ones(200,1);
% [bin,binsize] = conncomp(G);
% tic
% %
% Q = [];
% while k1<190%length(binsize) < 2
% 
%     
%     % total flow
%     %     tic
% %     for x_c = 100.*ptm
% %         Q2 = 0;
% %         for ed=1:1:numedges(G)
% %             twonodes = G.Edges.EndNodes(ed,:);
% %             twonodes_x = [G.Nodes(twonodes(1),:).comx,G.Nodes(twonodes(2),:).comx];
% %             twonodes_y = [G.Nodes(twonodes(1),:).comy,G.Nodes(twonodes(2),:).comy];
% %             Y_x = (G.Nodes(twonodes(1),:).comy-G.Nodes(twonodes(2),:).comy)/(G.Nodes(twonodes(1),:).comx-G.Nodes(twonodes(2),:).comx)*(x_c-G.Nodes(twonodes(1),:).comx)+G.Nodes(twonodes(1),:).comy;
% %             if (min(twonodes_y) < Y_x) && (Y_x < max(twonodes_y))
% %                 Q2 = Q2+G.Edges.Flows(ed);
% %                 
% %             end
% %         end
% %         Q = [Q;Q2];
% %     end
% %     toc
%     
%     
%     
%     
%     
%     
%     binsize_vec(k1) = max(binsize);
%     dep_edges = G.Edges(find(tau_w < tau_cr),:);
% 
%     [max_flow, max_ind] = max(dep_edges.Flows);
%     maxflow_vec(k1) = max(max_flow);
%     maxind_vec(k1) = max(max_ind);
%     % [max_flow, max_ind] = max(dep_edges.Flows);
% 
%     G = rmedge(G, max_ind);
% %     k2 = 0;
% %     while k2<4
% %         [max_flow, max_ind] = max(dep_edges.Flows);
% %         G = rmedge(G, max_ind);
% %         k2 = k2 +1;
% %     end
%     
%     [bin,binsize] = conncomp(G);
%     
%     
%     
%     if mod(k1, 1) == 0
%         
%         if make_vid == 1
%             ff = figure('units','pixels','position',[0 0 1020 1080]);
%             tilez = tiledlayout(ff,3,1);
%             %Tile1
%             tile1 =  nexttile;
%     %       
%             if graph_mode == 'exp'
%                 h= imshow(Im2);set(h, 'AlphaData', 0.2);hold on;
%                 pxdata = (G.Nodes.comx)./ptm;
%                 pydata =(G.Nodes.comy)./ptm;
%                 xlim([0 4200]);
%             else
%                 h = plot(X(:,1),X(:,2),'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 9); hold on;
%                 pxdata = (G.Nodes.comx);
%                 pydata =(G.Nodes.comy);
%             end
%             p = plot(G);p.LineWidth = (G.Edges.Widths).*3e5;p.Marker = 'none';
%             p.EdgeCData = abs(G.Edges.Flows);p.LineWidth = 5;
%            
%             p.XData = pxdata;
%             p.YData =pydata;
%             title(strcat(num2str(V1),'Pa,',num2str(k1), ',' , num2str(tau_cr)))
%             colormap(tile1,psymap2);
%             
%             set(gca, 'CLim', [8e-12 2e-11]);
%             colorbar
%     %         %Tile2
%             tile2 = nexttile;
%             if graph_mode == 'exp'
%                 h= imshow(Im2);set(h, 'AlphaData', 0.2);hold on;
%                 pxdata = (G.Nodes.comx)./ptm;
%                 pydata =(G.Nodes.comy)./ptm;
%                 xlim([0 4200]);
%             else
%                 h = plot(X(:,1),X(:,2),'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 9); hold on;
%                 pxdata = (G.Nodes.comx);
%                 pydata =(G.Nodes.comy);
%             end
%             p = plot(G);p.LineWidth = 5;p.Marker = 'none';
%             p.EdgeCData = abs(G.Edges.Shear);
%             p.XData = pxdata;
%             p.YData =pydata;
% %             xlim([0 4200]);
%             title('shear')
%             colormap(tile2, psymap1);
%     
%             set(gca, 'CLim', [0 tau_cr*2]);
%             colorbar
%             tilez.TileSpacing = 'none';
%             tilez.Padding = 'compact';
%             % Tile 3
%             tile3 = nexttile;
%             if graph_mode == 'exp'
%                 h= imshow(Im2);set(h, 'AlphaData', 0.2);hold on;
%                 pxdata = (G.Nodes.comx)./ptm;
%                 pydata =(G.Nodes.comy)./ptm;
%                 xlim([0 4200]);
%             else
%                 h = plot(X(:,1),X(:,2),'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 9); hold on;
%                 pxdata = (G.Nodes.comx);
%                 pydata =(G.Nodes.comy);
%             end
%             p = plot(G);p.LineWidth = 5;p.Marker = 'none';
%             p.EdgeCData = abs(G.Edges.Flows.*G.Edges.Resistances);
%             p.XData = pxdata;
%             p.YData =pydata;
% %             xlim([0 4200]);
%             title('local potential')
%             colormap(tile3, psymap1);
%     
%             set(gca, 'CLim', [0 pot_mean]);
%             colorbar
%     %         tilze.title(strcat(num2str(V1),'Pa,',num2str(k1), ',' , num2str(tau_cr)))
%     %    saveas(gcf,[pwd strcat('/video_c/tau10_35',num2str(k1),'.tif')]);
%             F(k1) = getframe(gcf); 
% 
%             writeVideo(v,F(k1));
%         end
%     end
%     G = potSolver(G,V1,Boundary2,G01);
%     G.Edges.Flows = FindEdgeFlows(G.Edges, G.Nodes);
%     tau_w = abs(G.Edges.Flows)./((G.Edges.Widths.^3));
%     G.Edges.Shear = tau_w;
%     if mean(G.Edges.Flows(:,1)) < flow_cr %& G.Edges.Flows(randi(numedges(G)),1) < flow_cr)
%         if graph_mode == 'exp'
%             h= imshow(Im2);set(h, 'AlphaData', 0.2);hold on;
%             p = plot(G);p.LineWidth = (G.Edges.Widths).*3e5;p.Marker = 'none';
%             p.EdgeCData = abs(G.Edges.Flows);
%             p.XData = (G.Nodes.comx)./ptm;
%             p.YData =(G.Nodes.comy)./ptm;
%             title(strcat(num2str(V1),'Pa,',num2str(k1), ',' , num2str(tau_cr)))
%             colormap(psymap2);
%             set(gca, 'CLim', [8e-12 2e-11]);
%             colorbar
%             saveas(gcf,[pwd strcat('/figz/tau_median2',num2str(k1),'.tif')]);
%         end
%         k_cr = k1;
%         break
%     end
%    
%    
%    
%    
%    k1 = k1 + 1;
%    
%    close all
% end
% 
% toc
% 
% if make_vid ==1
%     close(v)
% end







%% functions





function edgeQ = FindEdgeFlows(Egraph, nodez_HSG)
El = length(table2array(Egraph));
edgeQ = ones(El,1);

for ed =1:1:El
    twonodes = (Egraph.EndNodes(ed,:));
%     twonodes_x = [nodez_HSG.comx(twonodes(1)),nodez_HSG.comx(twonodes(2))];
    twonodes_x = nodez_HSG(twonodes,:).comx;
%     twonodes_y = [nodez_HSG(twonodes(1)).comy,nodez_HSG(twonodes(2)).comy];
    [~, max_x_idx] = max(twonodes_x);
    [~, min_x_idx] = min(twonodes_x);
    edq  = (nodez_HSG.Potentials(twonodes(min_x_idx)) - nodez_HSG.Potentials(twonodes(max_x_idx)))./Egraph.Resistances(ed);
%         ed_col  = abs((nodez_HSG.Potentials(twonodes(2)) - nodez_HSG.Potentials(twonodes(1))))./Egraph.Resistances(ed);
%         if nodez_HSG.Potentials(twonodes(max_x_idx)) > nodez_HSG.Potentials(twonodes(min_x_idx))
%               ed_col = -1*ed_col; 
%         end
    edgeQ(ed) = edq;
%     lambda = 24./((1-0.351.*edge_thic(ed)./depth).*(1+edge_thic(ed)./depth)).^2;
%     shear_rate(ed) = abs(ed_col*4/(pi*(edge_thic(ed)/2)^3));
%     shear_rate2(ed) = abs(ed_col.^2.*res_vec(ed).*lambda./(8*(edge_thic(ed).*depth)^2));
    

end
end

function [Graph00, Nodes00, Links00, disto] = Im2Graph(Im1, div)


ThR=graythresh(Im1)/div;
BW = imbinarize(Im1,ThR);
C = bwskel(BW(:,:,1));
disto = bwdist(imcomplement(BW(:,:,1)));
% BW1 = BW;
[Graph00,Nodes00, Links00] = Skel2Graph3D(C,100);

end
function [edge_thic] = FindEdgeProps(Egraph,nodez_HSG,ptm,dists)
El = length(table2array(Egraph));
edge_thic = ones(El,1);
% edge_leng = ones(El,1);
% cntr_xs = ones(El,1);
for ed = 1:1:El
    
    twonodes = (Egraph.EndNodes(ed,:));
    cntr_y = round((nodez_HSG.comx(twonodes(1)) + nodez_HSG.comx(twonodes(2)))/(2));%.*ptm .*ptm
    cntr_x = round((nodez_HSG.comy(twonodes(1)) + nodez_HSG.comy(twonodes(2)))/(2));
    edge_thic(ed) = 2*log(2).*dists(cntr_x,cntr_y).*ptm;
    if edge_thic(ed) == 0
        edge_thic(ed) = 11.967*ptm;  %%%% 'eliminating' these edges
    end
    
%     dis = sqrt((nodez_HSG(twonodes(1)).comy - nodez_HSG(twonodes(2)).comy)^2 + (nodez_HSG(twonodes(1)).comx - nodez_HSG(twonodes(2)).comx)^2);
%     edge_leng(ed) = dis.*ptm;
    
end
end

function G00 = potSolver(G00,V00,bdry_lm,G01)

[bin,binsize] = conncomp(G00);%,'Type','weak');
idx = binsize(bin) == max(binsize);
SG00 = subgraph(G00, idx);
nodez_SG00 = G00.Nodes(idx,:);

idx_B100 =[];  %%%%%%%%%%% separate the boundary and center nodes %%%%%%%%%%
idx_B200 = [];
idx_C00 = [];
% boundary_limits = bdry_lm.*ptm; %%%%% sample_bead


% for nodenum = 1:1:numnodes(SG00)
%     if( nodez_SG00.comx(nodenum) < bdry_lm(1)  )
%         idx_B100 = [idx_B100 nodenum];
%     else
%         if nodez_SG00.comx(nodenum) > bdry_lm(2)
%            idx_B200 = [idx_B200 nodenum];
%        else       
%    if degree(SG00,nodenum) > 1 
%         idx_C00 = [idx_C00 nodenum];
%    end
%         end
%     end
% end


    for iii=1:1:numnodes(SG00)
        if ismember(SG00.Nodes.ID(iii,1),bdry_lm.left) 
            idx_B100 = [idx_B100; iii];
        elseif ismember(SG00.Nodes.ID(iii,1),bdry_lm.right) 
            idx_B200 = [idx_B200; iii];
        else
            idx_C00 = [idx_C00; iii];
        end
    end



G00 = subgraph(SG00, [idx_B100; idx_B200; idx_C00]);   
% nodez_HSG00 = nodez_SG00([idx_B100 idx_B200 idx_C00]);

% G.Nodes.comx = ptm.*[nodez_HSG00.comy]';
% G.Nodes.comy = ptm.*[nodez_HSG00.comx]';
% G.Nodes.comy = ptm.*(w.*ones(length(G.Nodes.comx),1)-[Nodes1.comx]');
limit1 = length(idx_B100);
limit2 = length(idx_B200);
% target00 = ((limit1+1):1:limit1+limit2);
% G.Edges.Lengths = ptm.*(G.Edges.Weight);
% EWidths = FindEdgeProps(G.Edges,G.Nodes,ptm,dists);
% ELengths = G00.Edges.Lengths;
% EWidths = G00.Edges.Widths;
% depth = 1e-3;
% EResistances = (12*dyn_vis.*ELengths./((EWidths.^3).*depth.*(1-0.063.*EWidths./depth))); %limit where w>>h

% G.Edges.Widths = EWidths;
% G.Edges.Lengths = ELengths;
EResistances00 = G00.Edges.Resistances ;
EConds00 = 1./EResistances00;
I = incidence(G00);

ws = EConds00;
% D = (abs(I)*ws).^(-1/2);
% D = diag(D);
W = diag(ws);
L = I*W*I';

lb = limit1+limit2;
lL = length(L);
LBB = L(1:lb,1:lb);
LBC = L(1:lb,lb+1:lL);
LCB = L(lb+1:lL,1:lb);
LCC = L(lb+1:lL,lb+1:lL);

psi_B100 = V00.*ones(limit1,1); %%% set up the voltage vector for boundaries
psi_B200 = zeros(limit2,1);

psi_B00 = [psi_B100; psi_B200];
LCCinv = inv(LCC);
LS = LBB - (LBC*(LCCinv*LCB));
JB = LS*psi_B00;

J = zeros(lL,1);
J(1:lb) = JB;
pot_vec = L\J;
%psi_C = -LCCinv*(LCB*psi_B); %%the solution for the special case of two b%nodes
CV1 = -pot_vec(1) + V00;
% testtest=CV1.*ones(lL,1);
pot_vec = CV1+pot_vec; %% getting the unique solution by a shift

G00.Nodes.Potentials = pot_vec;

    for iii = 1:1:numnodes(G01)
        if ismember(G01.Nodes.ID(iii),G00.Nodes.ID)
            continue
        else
            new_node  = G01.Nodes(iii,:);
            new_node.Potentials = 0;
            G00 = addnode(G00,new_node);
        end
    end


end