function [bins,cc,ecc,n] = corr_orientation_func(x,y,nu,nv,bin_res)
%find the spatial correlation function (x,y) for the normalised orientation (nu,nv)
%   Detailed explanation goes here
    if numel(x)~=numel(nu);
        disp('ERROR: position and velocity should have same size');
    end
    fXX= x;fYY= y;fUU= nu-nanmean(nu(:));fVV= nv-nanmean(nv(:));
    dist_x= repmat(fXX(:), [1,numel(fXX(:))]);
    dist_y= repmat(fYY(:), [1,numel(fYY(:))]);
    DR= sqrt((dist_x - dist_x').^2 + (dist_y - dist_y').^2) ;

    CORRu  = repmat(fUU(:), [1,numel(fUU(:))]);
    CORRv  = repmat(fVV(:), [1,numel(fVV(:))]);

    CORR = (CORRu.*CORRu' + CORRv.*CORRv')./( sqrt(CORRu.^2 + CORRv.^2).*sqrt(CORRu'.^2 + CORRv'.^2)); 
    bins = 0:bin_res:max(DR(:));  %%% in pixels

    clear cc;
    clear n;
    for i=1:(numel(bins)-1)

       CORR_bin= (CORR( (DR(:) > bins(i)) & (DR(:) <= bins(i+1))));
       cc(i)= mean(CORR_bin(~isnan(CORR_bin)));
       ecc(i) = mean((CORR_bin(~isnan(CORR_bin))- cc(i)).^2);
       n(i) = numel(~isnan(CORR_bin));
       ecc(i)=sqrt(ecc(i)/n(i));
    end

end

