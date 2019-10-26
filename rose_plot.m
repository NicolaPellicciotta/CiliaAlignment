%%%% sript to analysie trajectories of beads from multiple video from the
%%%% Post_PIV_long.m

path_fold = '/home/np451/Desktop/ependymal data/November_2017/7.11/analysis'
cd(path_fold)
subdir={'Flow_control','Flow_20X' };

for j=1:1%numel(subdir)

    
start=0;
cd(path_fold);cd(subdir{j}); d=dir('*.movie');
clear xx;clear yy;clear uu;clear vv;

for i=1:size(d,1)  %%%% load all coordinates and velocities in 4 arrays (xx,yy,uu,vv)
    filename= d(i).name;
 
    cd(path_fold);cd(subdir{j});cd(filename);
  
    Apost= load('PIV_post_std.mat');
        um = Apost.uu;
        vm = Apost.vv;

        uu(start+1:start+numel(um)) = um;
        vv(start+1:start+numel(vm)) = vm;

        start= start+ numel(Apost.uu(:));
    end
cd(path_fold);cd(subdir{j});   

%%% correlation 

% fXX= xx;
% fYY= yy;
% fUU= uu;
% fVV= vv;
% 
% 
% 
% dist_x= repmat(fXX(:), [1,numel(fXX(:))]);
% dist_y= repmat(fYY(:), [1,numel(fYY(:))]);
% DR= sqrt((dist_x - dist_x').^2 + (dist_y - dist_y').^2) ;
% 
% %CORR = repmat(a(:),[1,numel(a)]) - (repmat(a(:),[1,numel(a)]))'; 
% CORRu  = repmat(fUU(:), [1,numel(fUU(:))]);
% CORRv  = repmat(fVV(:), [1,numel(fVV(:))]);
% 
% CORR = (CORRu.*CORRu' + CORRv.*CORRv')./( sqrt(CORRu.^2 + CORRv.^2).*sqrt(CORRu'.^2 + CORRv'.^2)); 
% bin_res = 96;
% bins = 0:bin_res:max(DR(:));
% 
% clear cc;
% clear n;
% for i=1:(numel(bins)-1)
% 
%        CORR_bin= (CORR( (DR(:) > bins(i)) & (DR(:) <= bins(i+1))));
%        cc(i)= mean(CORR_bin(~isnan(CORR_bin)));
%        ecc(i) = mean((CORR_bin(~isnan(CORR_bin))- cc(i)).^2);
%        n(i) = numel(~isnan(CORR_bin));
%        ecc(i)=sqrt(ecc(i)/n(i));
% end

M= sqrt(uu.^2 +vv.^2);
nu= uu./M;
nv= vv./M;
pol=sqrt(nanmean(nu)^2+ nanmean(nv)^2 );
UM= nanmean(nu);
VM= nanmean(nv);

figure(); rose(angle(complex(nu,nv)),18); 
title(strcat('Orientation distribution;  mean pol=',num2str(pol))); saveas(gcf,'rose.png');

% 
% figure(1);
% quiver(xx,yy,nu,nv,0.7);
% xlabel('X [um]');
% ylabel(' Y[um]');
% title('Flow pattern without validation');
% fig=figure(1);
% %saveas(fig,'Flow pattern.pdf');
% close all;
% 
% axis image
% 
% save()
% save('corr_Results.mat','cc','ecc','n','xx','yy','uu','vv','bin_res');
save('pol.mat','pol')
cd(path_fold);
end

%%