%% CMAM chain
make_vid = 1;
if make_vid == 1
    v = VideoWriter('cmam1d02.avi');
    v.Quality = 100;
    v.FrameRate = 10;
    k1 = 1;
    open(v)
end



bc_con = 'periodic';
init_mode = 1; % 1poisson,   2uniform
lamb = 3;
T = 1;Nstep = 128;
dt  = T/Nstep;
N = 50;
N_empty = 0;
s = ones(Nstep,1);

w = 3;
sc = (w + 2 - 2*sqrt(w+1))/w;

sites = initi(N,N_empty,init_mode,lamb);

m_all = sum(sites,"all");
m_vec  = 1:m_all;
% P = zeros(m_all,Nstep);
P = zeros(m_all,1000);
rho = zeros(Nstep,1);
iii = 1;
alpha  = 0;
%%
for step = 1:Nstep
    s(step) = nnz(sites)./N;
    if step > Nstep - 1000
        for m= 1:m_all
            P(m,step) = numel(find(sites == m))/N;
        end
        iii = iii+1;
    end
%     rho(step) = sum((m_vec)'.*P(:,step),"all");
    if make_vid == 1
        [v,k1]=Vid(v,k1,sites); 
    end
    ij_vec = randperm(numel(sites));
    for ij = ij_vec
        m = sites(ij);
        if m ~= 0
            
            ij_nei = [ij+1,ij-1];
            if ij == 1
                ij_nei(2) = N;
            elseif ij == N
                ij_nei(1) = 1;
            end
           
            if w*dt > rand(1,1) & m > 0
                sites(ij) = sites(ij) - 1;
                nei_id =ij_nei(randi(2));
%                 nei_id = ij_nei(1);
                sites(nei_id) = sites(nei_id) + 1;
                m = sites(ij);
            end
            
            if m > 0 & (m^-alpha)*dt > rand(1,1)
                
                nei_id =ij_nei(randi(2));
                sites(nei_id) = sites(nei_id)+sites(ij);
                sites(ij) = 0;
                m = sites(ij);
                
            end
        
            


        end


    
    end


end
% %%
% 
% sites = initi(N,N_empty,init_mode,lamb);
% m_all = sum(sites,"all");
% m_vec  = 1:m_all;
% %%
% tic
% t_now = 0;
% t_fin  = 210;
% iii = 1;
% Pnew = zeros(m_all,200);
% Nstep = 1000000;
% snew = zeros(Nstep./10,1);
% % while t_now < t_fin
% for ai = 1:Nstep
%     [~,all_events]=find(sites);
%     l_e = numel(all_events);
%     
%     if t_now > t_fin - 10 %& iii < t_fin
%         for m= 1:m_all
%             Pnew(m,iii) = numel(find(sites == m))/N;
%         end
%         snew(iii,1) = l_e./N;
%         iii = iii+1;
%     end
% %     all_rates = 2.*l_e.*(w+1);
% %     tau = exprnd(-t_fin*(all_rates));
%     array1 = (1/w).*unifrnd(0,1,l_e,1);
%     array2 = unifrnd(0,1,l_e,1);
%     array = [array1;array2];    array = -log(array);
% %     [tau,I] = min(-log([unifrnd(0,1,l_e,1);(1/w).*unifrnd(0,1,l_e,1)]));
%     [tau,I] = min(array);
%     if I > l_e
%         site_ind = all_events(I-l_e);
%         site_nei = [site_ind+1,site_ind-1];
%         if site_ind == 1
%             site_nei(2) = N;
%         elseif site_ind == N
%             site_nei(1) = 1;
%         end
%         site_nei =site_nei(randi(2));
%         sites(site_nei) = sites(site_nei) + 1;
%         sites(site_ind) = sites(site_ind) - 1;
% 
%     else
%         site_ind = all_events(I);
%         site_nei = [site_ind+1,site_ind-1];
%         if site_ind == 1
%             site_nei(2) = N;
%         elseif site_ind == N
%             site_nei(1) = 1;
%         end
%         site_nei =site_nei(randi(2));
%         sites(site_nei) = sites(site_nei) + sites(site_ind);
%         sites(site_ind) = 0;
%     end
% %229
%     t_now = t_now + tau;
% end
% toc
% %%


%%

if make_vid == 1
    close(v)
end

% P7 = P;
% s7 = mean(s(Nstep-1000:end));
% s7std = std(s(Nstep-1000:end));
% % loglog(m_vec, mean(P2,2), 'kx')%,'MarkerFaceColor', 'k')P(:,Nstep)
% % xlabel('m', 'Interpreter','latex') 
% % ylabel('P(m)', 'Interpreter','latex') 
% % figure(2)
% plot(s)
% 
% %%
% n = 100;
% pl = length(P7);
% clrvec = (1:n)./n;%1-clrvec(l)
% for l = 1:n
%     loglog(mean(P7(:,pl-l:end),2),'o','color', [0 0 clrvec(l)], 'MarkerSize', l);hold on;
% end




%%

function sites = initi(N,N_empty,init_mode,lamb)

if init_mode == 1
    sites = poissrnd(lamb,[1,N]);
elseif init_mode == 2
    sites = lamb.*ones(1,N);
end


empty_idx = ones(1,N);
empty_idx(randperm(numel(empty_idx), N_empty)) = 0;
sites = sites.*empty_idx;
sites = sites';

end




function plot_recs(sites)
N = numel(sites);
for i  = 0:1:N-1
        rec_color = 'c';
    for k = 0:1:sites(i+1,1)-1
        rectangle('Position', [i k 1 1], 'FaceColor', rec_color); hold on; 
%         ectangle('Position', [i 0 1 sites(i,1)], 'FaceColor', 'c'); hold on; 
    end
end
axis equal
xlim([0 N])
% ylim([0 max()])
axis off
end

function [v,k1]=Vid(v,k1,sites)
            plot_recs(sites)
            F(k1) = getframe(gcf); 
            writeVideo(v,F(k1));
            close all ;
            k1 = k1 +1;
end



%     prob_vec = [1 ; w];
%     prob_norm=[0 prob_vec']/sum(prob_vec);
%     %%create cumlative distribution
%     p_dist=cumsum(prob_norm);
%     %%calculate which bin the random number falls into (which edge the particle selects)
%     [~,~,inds] = histcounts(rand,p_dist); 
%     site_i=randi(l_e);
%     site_nei = [site_i+1,site_i-1];
%     if site_i == 1
%         site_nei(2) = N;
%     elseif site_i == N
%         site_nei(1) = 1;
%     end

%     if inds == 1
%         sites(site_nei) = sites(site_nei) + site_i;
%         site_i = 0;
%     else
%         sites(site_nei) = sites(site_nei) + 1;
%         site_i = site_i - 1;
%         sites(all_events(I-l_e)) = site_i;??
%     end