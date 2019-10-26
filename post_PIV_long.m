%path = '/home/np451/Desktop/ependymal data/10.8/beads/';



path = '/u/homes/np451/Desktop/ependymal data/December_2017/18.12/';
cd(path);
path_dir= strcat(path,'analysis');
store_dir='/u/homes/np451/Desktop/ependymal data/December_2017/18.12/';

cd(path_dir);

%subdir={'v0_1','v0.5_4','v1_4','v1.5_4','v2_4','v0_0','v0.5_5','v1_5','v1.5_5','v2_5'};
%subdir={'40X_Ch2'}%,'40X_Ch1_nr'};
subdir={''};

for jj=1:numel(subdir)
disp(jj);
%exp=subdir{jj};
exp=num2str(subdir(jj));
cd(path_dir);    
exp_path= strcat(path_dir,'/',exp);
exp_store= strcat(store_dir,exp);

cd(exp_store);
d=dir('*.movie');

for nf=1:size(d,1)
    
        filename= d(nf).name;
        cd(exp_store);    %%% moving to store to load frames

        mo=moviereader(filename);
        fs=mo.read();
        cd(exp_path); cd(filename);  %%%% moving to the analysis dit

        load('PIV_result.mat');
        close('all')
        

        %% PostProcessing  removing velocity out of range

        %%%%% remove velocity that do not make sense
        
        ulim= [nanmean(u(:)) - 3*nanstd(u(:)),nanmean(u(:)) + 3*nanstd(u(:))];
        vlim= [nanmean(v(:)) - 3*nanstd(v(:)),nanmean(v(:)) + 3*nanstd(v(:))];
        u1=u;v1=v;
        u1( u<ulim(1) | u>ulim(2))= nan;
        v1( v<vlim(1) | v>vlim(2))= nan;

        fst= std(double(fs),[],3);

        
        bs=30/2;   %%% box size PIV
        clear um;clear vm;
        
        %%%% weighted average over time of the velocity giving much weight
        %%%% (array w) depending on the fluctuation on intensity of the mean inside the box 
        %%%% box has size bs as the PIV
        
        
        %% find value to use for p_val
%         imagesc(fs(:,:,1)); 
%         disp(strcat(num2str(nf),filename));
% %%%% uncomment this part if you want interactive selection of noise and
% %%% pval
%         %%% select noise 
%         disp('Select Noise');
%         R_n =getrect();R_n=floor(R_n); noise= fs(R_n(2):R_n(2)+bs*2,R_n(1):R_n(1)+bs*2,1); v_noise=mean(noise(:));
% 
%         %%% select focused bead
%         disp('Select very focused bead');
%         R_g =getrect();R_g=floor(R_g);  good= fs(R_g(2):R_g(2)+R_g(4),R_g(1):R_g(1)+R_g(3),1); 
%         %[gp,ip]= max(good(:)); [ip1,ip2]=ind2sub(size(good),ip);
%         %good= good(ip1-bs:ip1+bs,ip2-bs:ip2+bs);  
%         v_good=mean(good(:));
% 
% %         %%% select halo
% %         disp('Select halo bead');
% %         R_h =getrect();R_h=floor(R_h);  halo= fs(R_h(2):R_h(2)+R_h(4),R_h(1):R_h(1)+R_h(3),1); 
% %         [gp,ip]= max(halo(:)); [ip1,ip2]=ind2sub(size(halo),ip);
% %         halo= halo(ip1-bs:ip1+bs,ip2-bs:ip2+bs,1);  v_halo=mean(halo(:));
%         
%         close all;
%         subplot(2,2,1);imagesc(noise);title(strcat('noise ',num2str(v_noise)));
%         subplot(2,2,2);imagesc(good);title(strcat('good ',num2str(v_good)));
% %         subplot(2,2,3);imagesc(halo);title(strcat('halo ',num2str(v_halo)));
%         saveas(gcf,'good_halo','pdf'); close all;
%         
        
%         p_val= abs(v_good- v_noise) -abs((v_good - v_noise)/5) % -(abs(v_good-v_halo)/3);
%         if p_val<0.8; p_val=0.8; end; 
        
        %%
        
 
        p_val=250;
         
         for b=1:numel(x);
             y_box=floor(y(b)-bs);
             x_box=floor(x(b)-bs);
             box_fs=fs(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs),:);
             box_m = mean(mean(box_fs,1),2);box_m= box_m(:);
             if mod(numel(box_m),2)==1; box_m=box_m(1:end-1);end;
             box_m=mean(reshape(box_m(:),[2,floor(numel(box_m(:))/2)]),1);    


             %%%% weight array with value proportional with the prominency of the peak 
             [pks,locs,wi,p]=findpeaks(box_m(:));
             w= zeros(size(box_m(:)));
             indx=1:numel(p); indx=indx(p>p_val); 
             for k=1:numel(indx); width= floor(wi(indx(k))/2); w_min =locs(indx(k))-width; 
                 if w_min<1; w_min=1; end;  w(w_min:locs(indx(k))+width) = p(indx(k)); end;            
             w=w(1:numel(box_m));
             [I,J]=ind2sub(size(x),b);
             um(I,J)= nansum( squeeze(u(I,J,:)).*w)/sum(w);
             vm(I,J)= nansum( squeeze(v(I,J,:)).*w)/sum(w);
         end
  %      xx=x(~isnan(um) | ~isnan(vm));yy=y(~isnan(um) | ~isnan(vm));uu= um(~isnan(um) | ~isnan(vm)); vv=vm(~isnan(um) | ~isnan(vm));
        

        
        %%%remove weird vectors
        ulim= [nanmean(um(:)) - 3*nanstd(um(:)),nanmean(um(:)) + 3*nanstd(um(:))];
        vlim= [nanmean(vm(:)) - 3*nanstd(vm(:)),nanmean(vm(:)) + 3*nanstd(vm(:))];
        u1m=um;v1m=vm;
        u1m( um<ulim(1) | um>ulim(2))= nan;
        v1m( vm<vlim(1) | vm>vlim(2))= nan;
        um=u1m;vm=v1m;
        
        xx=x(~isnan(um) | ~isnan(vm));yy=y(~isnan(um) | ~isnan(vm));uu= um(~isnan(um) | ~isnan(vm)); vv=vm(~isnan(um) | ~isnan(vm));
   
        
        figure(1);
        title(strcat(filename(1:end-5),' PostProcessing'));
        imagesc(fst);hold on;
        quiver(xx,yy,uu,vv,2,'r');
        fig=figure(1);
        saveas(fig,strcat(filename(1:end-5),'_PostPro_fst.png'));
      %  close(1);        
        
        
        save('PIV_post_std.mat','ulim','vlim','uu','vv','xx','yy','p_val');

      

        end
end


%% Load data to find mean polarisation over multiple field of view

%%%%%% this script is suitable with PIV_script_maker
%%%% it read in order of date
ppm= 0.14*2

path_fold = '/u/homes/np451/Desktop/ependymal data/December_2017/18.12/analysis/' 
cd(path_fold)
%subdir={'v0_1','v0.5_4','v1_4','v1.5_4','v2_4','v0_0','v0.5_5','v1_5','v1.5_5','v2_5'};
subdir=1:3;

for j=1:numel(subdir)

    
start=0;
cd(path_fold);cd(num2str(subdir(j))); d=dir('*.movie');
clear xx;clear yy;clear uu;clear vv;

for i=1:size(d,1)  %%%% load all coordinates and velocities in 4 arrays (xx,yy,uu,vv)
    filename= d(i).name;

    cd(path_fold);cd(num2str(subdir(j)));cd(filename);
  
    Apost= load('PIV_post_std.mat');
        xn= Apost.xx      %%%%% weird but it work seems rotation of 90 degree anticlockwise? (x,y)-> (y,-x)
        yn= Apost.yy 
        um = Apost.uu;
        vm = Apost.vv;


        xx(start+1:start+numel(xn)) = xn;
        yy(start+1:start+numel(yn)) = yn;

        uu(start+1:start+numel(um)) = um;
        vv(start+1:start+numel(vm)) = vm;

        start= start+ numel(Apost.xx(:));
    end
cd(path_fold);cd(num2str(subdir(j)));   

% %% correlation 
% 
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
close all;

% 
% figure();
% quiver(xx,yy,nu,nv,0.7);
% xlabel('X [um]');
% ylabel(' Y[um]');
% title('Flow pattern without validation');
% axis image
% saveas(gcf,'FlowPattern.pdf');
% close all;



% save('corr_Results.mat','cc','ecc','n','xx','yy','uu','vv','bin_res');
% save('pol.mat','pol');
cd(path_fold);
end



%% plot number box
 bs=30/2;   %%% box size PIV
figure();imagesc(fst);
 for b=1:numel(x);
     y_box=floor(y(b)-bs);
     x_box=floor(x(b)-bs);
     hold on;rectangle('position',[x_box y_box bs*2 bs*2]);
     text(x_box+bs,y_box+bs,num2str(b),'HorizontalAlignment','center');
%                 box_fst=fst(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs));
%                 val= mean(box_fst(:));
%                 if val>0.2
%                     [I,J]=ind2sub(size(x),b);
%                     xx=cat(1,xx,x(b));yy=cat(1,yy,y(b)); uu=cat(1,uu,nanmean(u1(I,J,:)));vv=cat(1,vv,nanmean(v1(I,J,:)));
% 
% 
% 
%                 end
 end
