function [fp,fq,m_pxx] = average_fft_meanpk(fs,mask,FR,fq_lim,res)
    if nargin < 5 | isempty(res)
        res=2;
    end
    
    if nargin < 4 | isempty(fq_lim)
        fq_lim(1)= 10;
        fq_lim(2)= 30;
    end
    N_frames= size(fs,3);
    fs_roi= fs(repmat(logical(mask),[1,1,N_frames]));
    fs_roi=reshape(fs_roi,[sum(mask(:)),N_frames]);
    %roi=double(fs_roi)- mean(fs_roi,2);
    roi=double(fs_roi)- movmean(fs_roi,60,2);

%%%% periodogram needs time on the row and pixel on the coloun; 
    window = hann(floor(N_frames));
    window= repmat(window,[1,size(roi,1)])';
    n= floor(N_frames/res);
    if mod(n,2)==0; n= n-1;end
    pxx= abs(fft(double(roi).*window,n,2)).^2;
    m_pxx= mean(pxx(:,1:floor(n/2)),1);
    fq= (0:(FR./n):(FR./2-FR./n));
    f_range=fq> (fq_lim(1)) & fq<(fq_lim(2));
    [pks,locs,w,p] = findpeaks(m_pxx(f_range),fq(f_range));%%%% find the peaks frequency in the selected freq range
    [~,ind_sort]= sort(pks);                                           %%%% sort peaks and get an index 
    pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);  %%% order all the variables with the same index
    [~,where] = max(pks);                                       
      
    
    %%%%% asign the frequency as the one with the highest peak, with safety
    %%%%% for harmonics
    if numel(locs)==0; fp=nan;    
    elseif numel(locs)==1; fp= locs(end);
    elseif numel(locs)>1; 
        if abs((locs(end)/2)-locs(end-1))< 0.3 & pks(end-1)> 0.66*pks(end);
            fp= locs(end-1); 
        else fp=locs(end);
        end
    end
end