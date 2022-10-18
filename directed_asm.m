%% directed_asm


N = 20;
rate_d = 1;
max_n  = 10;
sites  = zeros(N,1);
Delta = diag(max_n.*ones(N,1));
T = 200;
make_vid = 1;


if make_vid == 1
    v = VideoWriter('pasm00.avi');
    v.Quality = 100;
    v.FrameRate = 10;
    k1 = 1;
    open(v)
end



for i = 1:1:N-1
    Delta(i,i+1) = 1;
end

for t = 1:1:T
    indx = randi(N);
    sites(indx,1) = sites(indx,1) +1;
    if sites(indx,1) == max_n
        sites  = relax(sites, max_n,indx,N);
    end
    
    if make_vid == 1
        [v,k1]=Vid(v,k1,N,max_n,sites); 
    end

end

close(v)

%% plot

% for i  = 0:1:N-1
%     for k = 0:1:sites(i+1,1)-1
%         rectangle('Position', [i k 1 1], 'FaceColor', 'c'); hold on; 
% %         ectangle('Position', [i 0 1 sites(i,1)], 'FaceColor', 'c'); hold on; 
%     end
% end
% axis equal
% xlim([0 N])
% ylim([0 max_n])
%% functions

function sites = relax(sites, max_p,inds,N)
    indx = inds;
    while indx > 0
        topple_size = randi(max_p);
        sites(indx,1) = sites(indx,1) - topple_size;
        if indx < N
            sites(indx+1,1) = sites(indx+1,1) + topple_size;
        end
        indx=find(sites >= max_p,1);
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

function [v,k1]=Vid(v,k1,N,max_n,sites)
            plot_recs(N,max_n,sites)
            F(k1) = getframe(gcf); 
            writeVideo(v,F(k1));
            close all ;
            k1 = k1 +1;
end