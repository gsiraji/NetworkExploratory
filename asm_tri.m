%% directed_asm_tri


N = 71; %only odd #s here
N = 20;
rate_d = 1;

make_vid = 1;
%% diamond grid
% disorder = 0;
% heter = 0;
% [Adj,Xij] = makedigrid(N,N);
% % [X,GV] = synthetics(disorder, heter, N,N);
% G = digraph(Adj);
% G = assignXY(G,N);
% isites = 1:2:N-2;
%% square grid
conn = 2; %1 directed 3 neighbors - 2 undirected 4 neighbors
[G,Adj] = squaregrid(N, conn);
isites = 1:N;

% for i = 1:numnodes(G)
%     
%     if mod(i,2) == 1
%         G.Nodes.comy(i,:) = rem(i,N-1)+1; 
%         G.Nodes.comx(i,:) = floor(i/(N-1))*2+1;
%     else
%         G.Nodes.comy(i,:) = mod(i,N+1)-1;
%         G.Nodes.comx(i,:) = (floor(i/(N+1))+1)*2;
%         if mod(i,N+1) == 0
%             G.Nodes.comy(i,:) = N;
%             G.Nodes.comx(i,:) = (floor(i/(N+1)))*2;
%         end
%         
%     end
% 
% end
%%

% p = plot(G);
% p.XData = (G.Nodes.comx);
% p.YData =(G.Nodes.comy);

%%
v = VideoWriter('pasm08.avi');v.Quality = 100;v.FrameRate = 8;
if make_vid == 1
    %     v = VideoWriter('pasm07.avi');
        open(v)
end

max_n  = 4;
for count3r = 1%1:200
    sites  = table2array(G.Nodes);
    sites(:,3) = zeros(numnodes(G),1);
    % Delta = adjacency(GV);
    Delta = Adj;
    % Delta = diag(max_n.*ones(N,1));
    
    T = 480;
    make_vid = 0;
    
    starting_time = 1;
    pause_t = 0;

    tvec = 1:1:T;
    toppling_vec = zeros(length(tvec),1);
    ava_vec = zeros(length(tvec),1);
    k1 = make_vid;
    
    colormap1 = [0 0 0
        1 0 0];
    sites(20,3) = 10;
    ids_vec = zeros(35,length(tvec));
    for t = tvec
        k2 = 1;
    
    
        indx = isites(randi(length(isites)));
        sites(indx,3) = sites(indx,3) +1;
        if sites(indx,3) >= max_n
            toppling_vec(t,1) = indx; 
            [ids,avalanche_sz,sites,v,k1]  = relax(sites, max_n,indx,N,Delta,v,k1,t,starting_time);
            ava_vec(t) =  avalanche_sz;
            ids_vec(k2:length(ids),t) = ids;
            k2 = k2+length(ids);
        end
        
%         if t > starting_time
%     %         if mod(t,pause_t) == 0
%     %                 t
%                 if make_vid == 1
%                     [v,k1]=Vid(v,k1,sites,max_n);
%                     k1 = k1+1;
%                 end
%     %         end
%         end
    
    end
    
    if make_vid == 1
        close(v)
    end
    
    %% measure s distribution
    u = numnodes(G);
    KU_temp = DistS(ava_vec,u);
    if count3r > 1
        KU = mean([KU, KU_temp],2);
    else
        KU = KU_temp;
    end
end

%% test relax
% 
% [ids,avalanche_sz,sites,v,k1]  = relax(sites, max_n,indx,N,Delta,v,k1,t,starting_time);
% plot_surf(sites,max_n)


%% plot KU
% plot( [KU(1:100);KU(100:25:end)]./(sum(KU,'all')), '*k');set(gca, 'YScale', 'log');set(gca, 'XScale', 'log');xlabel('s')
% ylabel('D(s)'); hold on
% ylim([1e-5 max(KU)])
%%


% for i  = 0:1:N-1
%     for k = 0:1:sites(i+1,1)-1
%         rectangle('Position', [i k 1 1], 'FaceColor', 'c'); hold on; 
% %         ectangle('Position', [i 0 1 sites(i,1)], 'FaceColor', 'c'); hold on; 
%     end
% end
% axis equal
% xlim([0 N])
% ylim([0 max_n])

%     figure(1)
%     p = plot(G);
%     p.XData = (G.Nodes.comx);
%     p.YData =(G.Nodes.comy);
%     above_thr = find(sites(:,3) >= max_n);
%     color_vec = zeros(length(sites),1);
%     color_vec(above_thr) = 1;
%     p.NodeCData =  color_vec;
%     colormap(colormap1)
%     pause



%% functions


function [G,A] = squaregrid(N,connection)
    A = sparse(N^2,N^2); %sparse(i,j,v) v = ones(length(i),1)
    
        if connection == 1
            for i  = 1:N^2
%         if i < N^2
%             A(i, i+1) = 1;
                if i < N*(N-1)+1
                    A(i, i+N) = 1;
                    if mod(i,N) == 1
                        A(i, i+2*N-1) = 1;
                    else
                        A(i, i+N-1) = 1;
                    end
                    if mod(i,N) == 0
                        A(i, i+1) = 1;
                    else
                        A(i, i+N+1) = 1;
                    end
                end
                G = digraph(A);
            end
        elseif connection == 2
            for i  = 1:N^2
%         if i < N^2
%             A(i, i+1) = 1;
                if i < N*(N-1)+1
                    A(i, i+N) = 1;
                end
                if i < N^2
                    A(i, i+1) = 1;
                end
                if i > 1
                    A(i, i-1) = 1;
                    if i > N
                        A(i, i-N) = 1;
                    end
                end
                
            end
            G = graph(A);
        end
%         end
%         if mod(i,N) ~= 1
%             A(i, i-1) = 1;
%         end
    


    G.Nodes.comx = ones(numnodes(G),1);
    G.Nodes.comy = ones(numnodes(G),1);
    
    for i = 1:numnodes(G)
        G.Nodes.comx(i,1) = floor((i-1)/N);
        G.Nodes.comy(i,1) = mod(i-1,N);
    end

end





function G = assignXY(G,N)
    G.Nodes.comx = ones(numnodes(G),1);
    G.Nodes.comy = ones(numnodes(G),1);
    
    for i = 1:numnodes(G)
        
        if mod(i,2) == 1
            G.Nodes.comy(i,:) = rem(i,N-1)+1; 
            G.Nodes.comx(i,:) = floor(i/(N-1))*2+1;
        else
            G.Nodes.comy(i,:) = mod(i,N+1)-1;
            G.Nodes.comx(i,:) = (floor(i/(N+1))+1)*2;
            if mod(i,N+1) == 0
                G.Nodes.comy(i,:) = N;
                G.Nodes.comx(i,:) = (floor(i/(N+1)))*2;
            end
            
        end
    
    end

end




function [A,X] = makedigrid(Xfin,Yfin)
    

    mm = 0:2:Yfin;
    nn = 1:2:Xfin;
    heter = 0;
    disorder = 0;
    max_num = (Xfin-1)*(Xfin+1)/2;
    A = zeros(max_num);
    X = [];
    for i = 1:max_num
        if mod(i,2) == 0 % even
            if mod(i,Xfin+1) == 0 % even, top
                ind_down = (i+2*floor((max_num-i)/(Xfin+1))-1); % adj
                if ind_down<=max_num
                    A(i,ind_down) = 1;
                end
                if i+2 <= max_num % periodic bc
                        A(i,i+2) = 1;
                end
              
            elseif mod(i,Xfin+1) == 2 % even, bottom
                ind_up = (i+2*floor((max_num-i)/(Xfin+1))+1); % adj
                if ind_up<= max_num
                    A(i,ind_up) = 1;
                end
                if ((Xfin)*2+i) <= max_num %periodic bc
                    A(i,Xfin*2+i) = 1;
                end
                
            else % interior
                ind_down = (i+2*floor((max_num-i)/(Xfin+1))-1);
                ind_up = (i+2*floor((max_num-i)/(Xfin+1))+1);
                if ind_down <= max_num
                    A(i,ind_down) = 1;
                    if ind_up <= max_num
                        A(i,ind_up) = 1;
                    end
                end
              
            end
        else % interior, odd
            ind_down = (i+2*floor(i/(Xfin-1))+1); 
            ind_up = (i+2*floor(i/(Xfin-1))+3); 
            if ind_down <=max_num
                A(i,ind_down) = 1;
                if ind_up <=max_num
                    A(i,ind_up) = 1;
                end
            end
            
        end
    end



    for n = nn
            for m = mm
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
                    X = [X;  xi];
                end
                
                if disorder > rndii
                    xii = [m , n+ rndi];
                else
                    xii = [m,n];
                end
                if n < Yfin
                    X = [X;  xii];
                end
    
            end
    
    end
end




function [IDs,a,sites,v,k1] = relax(sites, max_p,inds,N,Delta,v,k1,t,starting_t)
function [v,k1]=Vid(v,k1,sites,max_n,N)
%             plot_recs(N,max_n,sites)
function plot_surf(sites,max_n,N)
% %     h = plot(X(:,1),X(:,2),'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 9); hold on
%     scatter(sites(:,1), sites(:,2),2,"black")
%     hold on;
% %     (sites(:,3)+0.01).*10
%     scatter(sites(:,1), sites(:,2), 50.*sites(:,3)+1,(sites(:,3)+0.01), 'filled')
% %     colormap("winter")
    
    for iik = 0:N-1
        for j = 1:N
    %         indz = find(sites(:,1) == i);
    %         Z(i,j) = sites(indz(find(sites(indz,2) ==  j))) 
            Z(iik+1,j) = sites(iik*N+j,3);
        end
    end
    [X,Y] = meshgrid(1:N,1:N);
    surf(X,Y,Z, 'EdgeColor','none')
    view([-235.139032646288 51.2882514403584])
    set(gca, 'CLim', [0 max_n]);
    colorbar
    axis off  
    
end

    plot_surf(sites,max_n,N);
            F(k1) = getframe(gcf); 
            writeVideo(v,F(k1));
            close all ;
            k1 = k1 +1;
end


    indx = inds;
    IDs = [];
%     a = 0;
    
    while indx > 0
        if (k1>0 & t >= starting_t)
            [v,k1] = Vid(v,k1,sites,max_p,N);
        end  
        IDs = [IDs; indx];
%         a = length(unique(IDs));
        a = length((IDs));
%         a = a + 1;
%         topple_size = randi(max_p);
%         topple_size = sites(indx,3);
        
        if sites(indx,1) < N 
%             & sites(indx,2) < N
            neighbors = find(Delta(indx,:));
            for i = neighbors
                sites(i,3) = sites(i,3) + 1;
            end
            topple_size = sites(indx,3)-1;
            topple_size = sites(indx,3);
            sites(indx,3) = sites(indx,3) - topple_size;
%             neighbors = [];
%             neighbors_temp = find(Delta(indx,:));
%             neighbors_x = sites(neighbors_temp,1);
%             for i = 1:1:length(neighbors_x)
%                 if neighbors_x(i,1) > sites(indx,1)
%                     neighbors = [neighbors; neighbors_temp(i)];
%                 end
%             end
%             indx_nei = randi(length(neighbors));
%             sites(neighbors(indx_nei),3) = sites(neighbors(indx_nei),3) + topple_size;
        end
        indx=find(sites(:,3) >= max_p,1);

    end

end






function plot_recs(N,max_n,sites)

    for i  = 0:1:N-1
        if sites(i+1,1) == max_n-1
            rec_color = 'r';
        else
            rec_color = 'c';
        end
        for k = 0:1:sites(i+1,1)-1
            rectangle('Position', [i k 1 1], 'FaceColor', rec_color); hold on; 
    %         ectangle('Position', [i 0 1 sites(i,1)], 'FaceColor', 'c'); hold on; 
        end
    end
    axis equal
    xlim([0 N])
    ylim([0 max_n])
    
    axis off
end

function plot_surf(sites,max_n)
%     h = plot(X(:,1),X(:,2),'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 9); hold on
    scatter(sites(:,1), sites(:,2),2,"black")
    hold on;
%     (sites(:,3)+0.01).*10
    scatter(sites(:,1), sites(:,2), 100,(sites(:,3)+0.01), 'filled')
    colormap("winter")
    set(gca, 'CLim', [0 max_n]);
    colorbar
    axis off  
    
end



function [v,k1]=Vid(v,k1,sites,max_n)
%             plot_recs(N,max_n,sites)
function plot_surf(sites,max_n)
%     h = plot(X(:,1),X(:,2),'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 9); hold on
    scatter(sites(:,1), sites(:,2),2,"black")
    hold on;
%     (sites(:,3)+0.01).*10
    scatter(sites(:,1), sites(:,2), 100.*sites(:,3),(sites(:,3)+0.01), 'filled')
%     colormap("winter")
    set(gca, 'CLim', [0 max_n]);
    colorbar
    axis off  
    
end

    plot_surf(sites,max_n);
            F(k1) = getframe(gcf); 
            writeVideo(v,F(k1));
            close all ;
            k1 = k1 +1;
end


function D_s =DistS(U,u)
    D_s = zeros(u,1);
%     D_s = ones(max(U),1);
%     Ui = unique(U);
    for i = 1:max(U)
        D_s(i) = length(find(U == i));
    end
    
end

%   %     U = unique(ava_vec);
%     U = (ava_vec);  %     ep = 1e-3;    KU_i = ava_vec(ava_vec < U(i)+ep); 
%         KU_i = KU_i(KU_i > U(i)-ep);   
%         KU(i) = length(KU_i);


function G0 = correlations(ids_vec,ssize)

G0 = zeros(ssize);

for i = 1:ssize
    [~, coli] = find(ids_vec == i);
    for j = 1:ssize
        
        G0(i,j) = length(find(coli == j));


    end
end


end
