
%%%%% PIV Mosaic
ppm= 0.0014 ; %%%% mm  (using 4X objective 1 pixel is 1.4 um) 

pos=    {'CC','BC','B2','TC','T2','CB','CT','BT','TB','TT','BB'};
y_shift=floor([  0 ,-1.5,-3  ,1.5 , 3  , 0  , 0  ,-1.5, 1.5, 1.5,-1.5]/ppm); 
x_shift=-floor([  0 , 0  , 0  , 0  , 0  ,-1.5, 1.5, 1.5,-1.5, 1.5,-1.5]/ppm); 


path_fold = '/home/np451/Desktop/Cilia_data'

cd(path_fold) 
suffix='*.movie';
direc = dir([path_fold,'4X_*']);
N_files= size(direc,1);

start=0;


for i=1:N_files    %%%% load all coordinates and velocities in 4 arrays (xx,yy,uu,vv)

    indx=0;    
    exp_name= direc(i).name;
    for j=1:size(pos,2);
        if (direc(i).isdir ==1) & (sum(pos{j} == exp_name(4:5)) == 2); 
        indx=j;
        end
    end
    disp(indx)
    if indx~=0   %%% indx corrispond to the
    cd(strcat(path_fold, exp_name));
    
    load('PIV_result.mat');
    load('PIV_post.mat');
    xn= x-x_shift(indx);
    yn= y-y_shift(indx);
    u1m = mean(u1,3);
    v1m = mean(v1,3);
    um = mean(u,3);
    vm = mean(v,3);
    
    
        
    
    xx(start+1:start+numel(xn)) = xn;
    yy(start+1:start+numel(yn)) = yn;
    uu1(start+1:start+numel(u1m)) = u1m;
    vv1(start+1:start+numel(v1m)) = v1m;
    uu(start+1:start+numel(um)) = um;
    vv(start+1:start+numel(vm)) = vm;
    
    start= start+ numel(x(:));
    
    end
end
    
cd(path_fold);

%% SET new lattice with space dx and do a velocy average over a radius dx/2

dx=96;

px_max=(floor(max(xx(:))/dx)+1)*dx;  %%%% set the new axis limit
py_max=(floor(max(yy(:))/dx)+1)*dx;

XX= repmat((-px_max:dx:px_max),[numel((-py_max:dx:py_max)),1]);
YY= repmat((-py_max:dx:py_max)',[1,numel((-px_max:dx:px_max))]);
UU= zeros(size(XX));
VV= zeros(size(YY));
UU1= zeros(size(XX));
VV1= zeros(size(YY));



min_d=dx;

uu(isnan(uu))=0;
vv(isnan(vv))=0;
uu1(isnan(uu))=0;
vv1(isnan(vv))=0;



for i=1:size(XX,1)
    for j=1:size(XX,2)
        dr= sqrt( (xx(:)-XX(i,j)).^2 + (yy(:)-YY(i,j)).^2 ) ;
        UU1(i,j) = mean(uu1(dr(:)<min_d));
        VV1(i,j) = mean(vv1(dr(:)<min_d));
        UU(i,j) = mean(uu(dr(:)<min_d));
        VV(i,j) = mean(vv(dr(:)<min_d));
    
    
    end
end 


%% figures

close('all')

figure(1);
quiver(XX*ppm,-YY*ppm,UU,-VV,3);
xlabel('[mm]');
ylabel('[mm]');
title('Flow pattern without validation');
axis image


figure(2);
quiver(XX*ppm,-YY*ppm,UU1,-VV1,3);
xlabel('[mm]');
ylabel('[mm]');
title('Flow pattern with validation');
axis image;

figure(3);
quiver(xx*ppm,-yy*ppm,uu,-vv,3);
xlabel('[mm]');
ylabel('[mm]');
title('Flow pattern preAverage without validation');
axis image;


%% orientation correlation
% dist_x= repmat(XX(:), [1,numel(XX(:))]);
% dist_y= repmat(YY(:)', [numel(YY(:))],1);
% DR= sqrt((dist_x - dist_x').^2 + (dist_y - dist_y').^2) ;
%%% take data inside circle
xc= 960;
yc= 600;
circle_r = floor(3.25/ppm);
mask= sqrt((XX-xc).^2 + (YY-yc).^2) < circle_r;
fXX= XX(mask);
fYY= YY(mask);
fUU= UU(mask);
fVV= VV(mask);


dist_x= repmat(fXX(:), [1,numel(fXX(:))]);
dist_y= repmat(fYY(:), [1,numel(fYY(:))]);
DR= sqrt((dist_x - dist_x').^2 + (dist_y - dist_y').^2) ;

%CORR = repmat(a(:),[1,numel(a)]) - (repmat(a(:),[1,numel(a)]))'; 
CORRu  = repmat(fUU(:), [1,numel(fUU(:))]);
CORRv  = repmat(fVV(:), [1,numel(fVV(:))]);

CORR = (CORRu.*CORRu' + CORRv.*CORRv')./( sqrt(CORRu.^2 + CORRv.^2).*sqrt(CORRu'.^2 + CORRv'.^2)); 
bin_res = 96;
bins = 0:bin_res:max(DR(:));

clear cc;
clear n;
for i=1:(numel(bins)-1)

       CORR_bin= (CORR( (DR(:) > bins(i)) & (DR(:) <= bins(i+1))));
       cc(i)= mean(CORR_bin(~isnan(CORR_bin)));
       ecc(i) = mean((CORR_bin(~isnan(CORR_bin))- cc(i)).^2);
       n(i) = numel(~isnan(CORR_bin));
       ecc(i)=sqrt(ecc(i)/n(i));
end
%%
figure(4);
plot(bins(2:end)*ppm,cc,'ro'); hold on;
plot(bins(2:end)*ppm, zeros(numel(bins(2:end))),'k-');
ylabel('(u1*u2)/(|u1||u2|)','FontSize',15);
xlabel('distance [mm]','FontSize',15);
xlim([0,7]);
title('orientation correlation as function of distance ');


%% velocity distribution

subplot(1,2,1)
histogram(fUU(:))
xlabel('velocity [px/frame]','FontSize',15);
title('histogram velocities along X');

subplot(1,2,2)
histogram(fVV(:))
xlabel('velocity [px/frame]','FontSize',15);
title('histogram velocities along Y');


%     XX= repmat(xx,[numel(xx(:)),1]);
%     YY= repmat(xx,[numel(xx(:)),1]);
%     I1 = repmat(1:numel(xx),[numel(xx),1]);
%     I2 = repmat((1:numel(xx))',[1,numel(xx)]);
%         
%     dX= XX- XX';
%     dY= YY -YY';
%     
%     dR= sqrt(dX.^2 +dY.^2);
%     dR(dR==0) =  min_d;
%     
%    
%     
%     for i=1:size(dR,1)
%         [MIN,argmin]= min(dR(i,:));
%         if MIN < min_d
%             xxn(i)=  (xx(i)+ xx(argmin))/2 ;
%             yyn(i)=  (yy(i)+ yy(argmin))/2 ;
%             uu1n(i)= (uu1(i)+ uu1(argmin))/2 ;
%             vv1n(i)= (vv1(i)+ vv1(argmin))/2 ;
%         else
%             xxn(i)=  xx(i);
%             yyn(i)=  yy(i);
%             uu1n(i)=  uu1(i);
%             vv1n(i)= vv1(i) ;
%         end
% end
