
chain_mode = 0;

if chain_mode == 0
    disorder1 = 0;
    heter1 = 0;
    Xfin1 = 200;
    Yfin1 = 36;
    
    [X,Gs] = synthetics(disorder1, heter1, Xfin1, Yfin1);
    X(:,2) = 1e-5*sqrt(3).*X(:,2);
    Gs.Nodes.Y =  1e-5*sqrt(3).*Gs.Nodes.Y;
    Gs.Nodes.X =  1e-5.*Gs.Nodes.X;
    V0 = 1e5/2.1378; 
    idx_B1 =[];  %%%%%%%%%%% separate the boundary and center nodes %%%%%%%%%%
    idx_B2 = [];
    idx_C = [];
    boundary_limits = 1e-5.*[1, Xfin1-1];
    
    
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
    
    
    for i= 1:1:numedges(G)
        endnode_v = G.Edges.EndNodes(i,:);
        X_v = G.Nodes.X(endnode_v);
        Y_v =  G.Nodes.Y(endnode_v);
        
        G.Edges.Lengths(i) = sqrt((X_v(1) - X_v(2)).^2+(Y_v(1) - Y_v(2)).^2);
    
    end


    G.Nodes.comx = G.Nodes.X;
    G.Nodes.comy = G.Nodes.Y;
    G.Edges.Widths = 1.2486e-05.*ones(numedges(G),1); % 1+ 13.*rand(numedges(G),1);
    dyn_vis = 3e-3; %3mpa-s
    depth = 1e-3; %A = 1mm^2

else

    
    chain_l = 1000;
    V0 = chain_l;

    s = 1:1:chain_l;
    t = 2:1:chain_l+1;
    weights = 1.*ones(chain_l,1);
    G = graph(s,t, weights);
    G.Nodes.comx = (0:1:chain_l)';
    G.Nodes.comy = zeros(chain_l+1,1);
    G = subgraph(G, [1 chain_l+1 2:1:chain_l]);
    
    G.Nodes.ID = (1:1:numnodes(G))';
    G.Edges.ID = (1:1:numedges(G))';
    
    limit1 = 1;
    limit2 = 1;
    Boundary2.left = G.Nodes.ID(1);
    Boundary2.right = G.Nodes.ID(1+limit1:limit1+limit2);
    G.Edges.Lengths = ones(chain_l, 1);
%     G.Edges.Widths = 1+ 13.*rand(numedges(G),1);
    G.Edges.Widths = 2.*ones(numedges(G),1);
    dyn_vis = pi/8; %3mpa-s
    depth = 1e-3; %A = 1mm^2
end
%%


ELengths = G.Edges.Lengths;
EWidths = G.Edges.Widths;

% EResistances = (12*dyn_vis.*ELengths./((EWidths.^3).*depth.*(1-0.063.*EWidths./depth))); %limit where w>>h


p_norm1 = 0.5;
p_dmeter = 1;
dep_prob = 1000;



EResistances = (8/pi*dyn_vis.*ELengths./(EWidths./2).^4);
EConds = 1./EResistances;

G.Edges.Resistances = EResistances;
G.Edges.Ci = EConds;

G = find_flows(G, Boundary2, V0, G);


G.Edges.Open = ones(numedges(G),1);
tau_w = abs(G.Edges.Flows)./((G.Edges.Widths.^3));
G.Edges.Shear = tau_w;
% tau_cr = 1e-2;
% dep_p = dep_prob.*(tau_cr.*ones(length(tau_w),1) - tau_w).*dt;

G.Edges.pCount = zeros(numedges(G),1);


%%

%%%% vid
k1 = 1;
v = VideoWriter('pchainsim00.avi');
v.Quality = 100;
v.FrameRate = 10;

G2 = G;


Nstep = 30;tot_time = 3;dt = tot_time/Nstep;

beginning_set = 1:1:limit1;
target1 = G.Nodes.ID((limit1+1):limit1+limit2);


% Nii = numedges(G2);
% Order_pari = (1/(Nii-1)).*(Nii - ((sum((G2.Edges.Flows).^2,'all'))^2)/(sum((G2.Edges.Flows).^4,'all')));   



tlcf = mean(abs(G.Edges.Flows))./(1e3);

inj_step = 1;
Nparticle = 100;
% inj_rate = Nparticle/inj_step;

tot_Nparticle = round(tot_time./inj_step)*Nparticle;
if inj_step > tot_time
    tot_Nparticle = Nparticle;
end
Nparticle_active = Nparticle;
particle_set = particle.empty(0,tot_Nparticle);




for p_idx = 1:1:tot_Nparticle
    particle_set(p_idx) = particle();
%     if p_idx < Nparticle + 1
    particle_set(p_idx) = init_particle(particle_set(p_idx),G.Nodes(beginning_set,:));
%     end
end



current_node_list = zeros(tot_Nparticle, 1);
mask = ones(numedges(G2),1);

passed_time = 0;


open(v)
%%%%%%
for step= 1:1:Nstep
    GID = G2.Edges.ID;
    G.Edges.pCount = zeros(numedges(G),1);
    EResistances = G2.Edges.Resistances;
    if current_node_list ~= 0
        if ismember(current_node_list,target1) 
            'done'
            break
        end
    end
    
    
    if passed_time > 0
        if mod(passed_time, inj_step) == 0
            Nparticle_active = Nparticle_active + Nparticle;
%             for p_idx = (Nparticle_active - Nparticle +1):1:Nparticle_active
%                 particle_set(p_idx) = init_particle(particle_set(p_idx),test_G2.Nodes(beginning_set,:));
%             end
        end
    end
    
    
  for p_idx = 1:1:Nparticle_active
    
%     [particle_set(p_idx),G2,mask] = next_action(particle_set(p_idx),incidence(G2),G2.Edges,target1,depth,dt,p_norm1,dep_prob,G2.Nodes,p_dmeter,G2,mask);
     [particle_set(p_idx),G2,G,mask] = next_action(particle_set(p_idx),...
         incidence(G2),G2.Edges,target1,depth,dt,p_norm1,dep_prob,G2.Nodes,p_dmeter,G2,mask,G);
%     pall_xy(step,2*p_idx-1:2*p_idx) =  [particle_set(p_idx).comx; particle_set(p_idx).comy];
    current_node_list(p_idx) = particle_set(1,p_idx).node_num;
    pedge_num = particle_set(1,p_idx).edge_num;
    
    GEdges = G2.Edges;
  end
    

    


 %%%%%%%%%%%% vid
    if mod(step,4) == 0
        step
        
%         tilez = tiledlayout(2,1);
        % Tile 1
%         nexttile;
        f = figure;
        f.Position = [100 100 540 200];
        if chain_mode == 0
            h = plot(X(:,1),X(:,2),'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 9); hold on;
        end
        p = plot(G, 'LineWidth', 5); hold on;
        p.Marker = 'none';
        p.XData = (G.Nodes.comx);
        p.YData =(G.Nodes.comy);

        p.EdgeCData = abs(G.Edges.Flows);
        colormap(autumn)
% %         set(gca, 'CLim', [0 0.5]);
%         set(gca, 'CLim', [8e+02 2e+03]);
%         set(gca, 'CLim', [0.5e+04 3e+04]);
%         for iii = 1:1:Nparticle
%             plot(pall_xy(1:step-1,2*iii-1)./ptm,pall_xy(1:step-1,2*iii)./ptm,'o', 'lineWidth',1.5, 'MarkerFaceColor', 'k', 'markerSize', 5);hold on;
%         end
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
                        'MarkerEdgeColor',e_color, 'markerSize', 3);hold on;
            end
        end
        axis off
        
        if size(unique(G2.Edges.ID),1) ~= size(unique(GID),1) | unique(G2.Edges.ID) ~= unique(GID)
    %         GIDtest = unique(GID);
    %         GIDtest2 = unique(G2.Edges.ID);
    %     if sum(G2.Edges.Resistances ~= EResistances) < numedges(G2)
%             G2 = potSolver(G2,V1,Boundary2,G);
%             G2.Edges.Flows = FindEdgeFlows(G2.Edges, G2.Nodes);
            G2 =  find_flows(G2, Boundary2, V0, G);
            tau_w = abs(G2.Edges.Flows)./((G2.Edges.Widths.^3));
            G2.Edges.Shear = tau_w;
        end
        
        
    %     plot(p100_xy(1:1:step-1,1)./ptm, p100_y(1:1:step-1)./ptm, ':o', 'lineWidth',2, 'color', 'k', 'MarkerFaceColor', 'm', 'markerSize', 5 );
%         hold on;
        
        
%         p.NodeCData = G2.Nodes.Potentials;
% %         labelnode(p,1:numnodes(G2),G2.Nodes.Potentials)
        
        %tile2

%         saveas(gcf,[pwd strcat('/video_2/test80k',num2str(k1),'.tif')]);
%         save([pwd strcat('/video_1/p155k',num2str(k1),'.mat')],'particle_set');
%         saveas(gcf, strcat('p80k',num2str(k1),'.tif'));
        F(k1) = getframe(gcf); 

        writeVideo(v,F(k1));

        close all ;
        k1 = k1 +1;
%         local_pot = abs(G2.Edges.Flows).*(G2.Edges.Resistances);
    end
%%%%%%%%%%

    if abs(mean(G2.Edges.Flows)) < tlcf
        'broken'
        break
    end
    passed_time = dt*step;


end

% if pedge_num > 0
%         edgeInd = find(G.Edges.ID == pedge_num);
%         G.Edges.pCount(edgeInd) = G.Edges.pCount(edgeInd)+1;
%     
% end




close(v)









function G = find_flows(G, Bndry, V, G1)
   
    V1 = V;   
%     G1 = G;
%     G.Edges.Open = ones(numedges(G),1);
    [bin,binsize] = conncomp(G);%,'Type','weak');
    idx = binsize(bin) == max(binsize);
    
%     G.Nodes.X = X_pos';
%     G.Nodes.Y = Y_pos';

   
    G = subgraph(G, idx); 
%     idx_B1 = intersect(Bndry.left,G.Nodes.ID);
%     idx_B2 = intersect(Bndry.right,G.Nodes.ID);
    idx_B1 = [];
    idx_B2 = [];
    idx_C = [];
    for iii=1:1:numnodes(G)
        if ismember(G.Nodes.ID(iii,1),Bndry.left) 
            idx_B1 = [idx_B1; iii];
        elseif ismember(G.Nodes.ID(iii,1),Bndry.right) 
            idx_B2 = [idx_B2; iii];
        else
            idx_C = [idx_C; iii];
        end
    end
            
    
    G = subgraph(G, [idx_B1; idx_B2; idx_C]); 
    D_G = incidence(G);
    
    if isempty(idx_B1)
        error('network broken')
    end
    if isempty(idx_B2)
        error('network broken')
    end
    muu = 1;
    radii = (G.Edges.Widths)./2;
    G.Edges.Ci = pi/(8*muu).*(radii.^4)./(G.Edges.Lengths);
%     cndct_vector = 1./resistances(idx);
    L = D_G*diag(G.Edges.Ci)*D_G'; %ones(length(D_HSG),1)
    
     % L = laplacian(HSG);  %%%%%%%%%%% uniform resistance
    Node_IDs = G.Edges.EndNodes; 
    
     
    

    lb = length(idx_B1)+length(idx_B2);
    lL = length(L);
    LBB = L(1:lb,1:lb);
    LBC = L(1:lb,lb+1:lL);
    LCB = L(lb+1:lL,1:lb);
    LCC = L(lb+1:lL,lb+1:lL);
    % V1 = 32; %
    psi_B1 = V1.*ones(length(idx_B1),1); %%% set up the voltage vector for boundaries
    psi_B2 = zeros(length(idx_B2),1);

    psi_B = [psi_B1; psi_B2];
    LCCinv = inv(LCC);
    LS = LBB - (LBC*(LCCinv*LCB));
    JB = LS*psi_B;

    J = zeros(lL,1);
    J(1:lb) = JB;

    psi_V = L\J;

    CV1 = -psi_V(1) + V1;
    flows = ones(numedges(G),1);
    psi_V = CV1+psi_V; %% getting the unique solution by a shift
%     psi_Vtemp = psi_V;
%     psi_V1 =  psi_V(1:length(idx_B1),1);
%     psi_V2 = psi_V(length(idx_B1)+1:lb,1);
%     psi_Vtemp(idx_B1,1) = psi_V1;
%     psi_Vtemp(idx_B2,1) = psi_V2;
%     psi_Vtemp(length(idx_B1)+1:length(psi_V)-length(idx_B2),1) = psi_V(lb+1:end);
    G.Nodes.Potentials = psi_V;
%     psi_V = psi_Vtemp;
    X_pos = G.Nodes.comx;
    
    for ed =1:1:length(Node_IDs)
        twonodes = Node_IDs(ed,:);
        twonodes_x = [X_pos(twonodes(1)),X_pos(twonodes(2))];
        twonodes_v = psi_V(twonodes);
        [max_x, max_x_idx] = max(twonodes_v);
        [min_x, min_x_idx] = min(twonodes_v);
%         if nodez_HSG(twonodes(1)).comy < nodez_HSG(twonodes(2)).comy
%         ed_col  = (psi_V(twonodes(1)) - psi_V(twonodes(2))).*cndct_vec(ed);
%         else
%             lilq  = -(psi_V(twonodes(max_x_idx)) - psi_V(twonodes(min_x_idx))).*G.Edges.Ci(ed);
            lilq  = (max(twonodes_v) - min(twonodes_v)).*G.Edges.Ci(ed);
%             if psi_V(twonodes(max_x_idx)) > psi_V(twonodes(min_x_idx))
            if twonodes_x(min_x_idx) < twonodes_x(max_x_idx)
                  lilq = -1*lilq; 
            end
%         end
        flows(ed) = lilq;
    end
   G.Edges.Flows = flows;
   
   
   for iii = 1:1:numnodes(G1)
        if ismember(G1.Nodes.ID(iii),G.Nodes.ID)
            continue
        else
            new_node  = G1.Nodes(iii,:);
            new_node.Potentials = 0;
            G = addnode(G,new_node);
        end
    end
   
   
   
   
end

function G = evolve_G(G,G1,mask) %G1-> test_G1
    rmID = G1.Edges.ID(~logical(mask));
    addID = G1.Edges.ID(logical(mask));
    
%     for i = 1:1:numedges(G1)
%         if mask(i,1) == 0
%             rmID = [rmID G1.Edges.ID(i)];
%         else
%             addID = [addID G1.Edges.ID(i)];
%         end
%     end
    
    for iii = addID'
        if nnz(ismember(G.Edges.ID,iii)) > 0 
            continue
        else
            % need to find the nodes that correspond to the node IDs
            new_edge = G1.Edges(find(G1.Edges.ID == iii,1),:);
            end_node1 = new_edge.EndNodes(:,1);
            end_node1ID = G1.Nodes.ID(end_node1);
            end_node2 = new_edge.EndNodes(:,2);
            end_node2ID = G1.Nodes.ID(end_node2);
           
            new_edge.EndNodes(:,1)  = find(G.Nodes.ID == end_node1ID,1);
            new_edge.EndNodes(:,2) = find(G.Nodes.ID == end_node2ID,1);
%             
            G = addedge(G,new_edge);
        end
    end
    
    for ii = rmID'
        if nnz(ismember(G.Edges.ID,ii)) > 0
            testID = find(G.Edges.ID == ii,1);
            G = rmedge(G,testID);
        end
    end

end