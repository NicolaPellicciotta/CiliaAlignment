function [N_near,df] = find_nearcilia_df(fXX,fYY,F32,d_lim)

dist_x= repmat(fXX(:), [1,numel(fXX(:))]);
dist_y= repmat(fYY(:), [1,numel(fYY(:))]);
FF= repmat(F32(:),[1,numel(F32(:))]);FF_temp=nan(size(FF));
DR= sqrt((dist_x - dist_x').^2 + (dist_y - dist_y').^2) ;
N_near=sum(DR<d_lim,2);
FF_temp(DR<d_lim)= FF(DR<d_lim);
local_f =nanmean(FF_temp,2);
local_df=nanstd(FF_temp,2);
df=local_df./local_f;


end
