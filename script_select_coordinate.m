%%% analisys 

fstd=std(double(frame_stack),[],3);
thresh = multithresh(fstd);
fstd_t= imquantize(fstd,thresh)-1;
std_mask= repmat(fstd_t,[1,1,size(frame_stack,3)]);

% clean trajectory x,y, mean(u), mean(v)

boxs = (s{8,2}/2)-1;
N_coord= numel(x(:));
mask=zeros(size(x(:)));
for j=1:N_coord
    val=fstd_t(floor(y(j)-boxs):floor(y(j)+boxs),floor(x(j)-boxs):floor(x(j)+boxs));
    val= mean(val(:));
    if val>0.3;
        disp('bau')
        mask(j)=1;
    end
end
mask=logical(mask);

um=nanmean(u,3);
vm=nanmean(v,3);
imagesc(fstd_t)
hold on;
quiver(x(mask),y(mask),um(mask),vm(mask),1,'r');