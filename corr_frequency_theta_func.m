function [idR,idTH,cc,ecc,n] = corr_frequency_theta_func(x,y,F,bin_res,bin_theta)
    %corelletation fucntion of the orientation in 2d
    % bin_res in pixels
    %bin_theta in rad

    if nargin < 5 | isempty(bin_theta)
      %  bin_theta=pi/2;  %%% this is for standard x and y compare
      bin_theta= 2*pi+1; %%%% so you get a single variable
    end

%find the spatial correlation function (x,y) for the normalised orientation (nu,nv)
%   Detailed explanation goes here
    if numel(x)~=numel(F);
        disp('ERROR: position and velocity should have same size');
    end

    fXX= x;fYY= y;ff= F-nanmean(F(:));
    dist_x= repmat(fXX(:), [1,numel(fXX(:))]);
    dist_y= repmat(fYY(:), [1,numel(fYY(:))]);
    DR= ceil(sqrt((dist_x - dist_x').^2 + (dist_y - dist_y').^2)/bin_res) ;

    TH= atan2( dist_y - dist_y',dist_x - dist_x');
    TH=(angle(exp(i*TH)*exp(i*pi/4))); %%% rotation 45degree
    TH(TH<0)=TH(TH<0)+pi;   
    TH=floor(TH/bin_theta);
    [th_max]=max(TH(:));TH(TH==th_max)=th_max-1;

    [G,idR,idTH]=findgroups(DR(:),TH(:));
    
    CORRf  = repmat( ff(:), [1,numel(ff(:))]);
    CORR = (CORRf.*CORRf');%./( sqrt(CORRf.^2).*sqrt(CORRf'.^2) ); 
    cc=accumarray(G,CORR(:),[],@nanmean); 
    cc=cc./cc(idR==0); %%% normalise to 1
    ecc=accumarray(G,CORR(:),[],@std)./sqrt(accumarray(G,CORR(:),[],@numel)) ;
    ecc=ecc./cc(idR==0); %%% normalise error
    n=accumarray(G,CORR(:),[],@numel);
    
    idR=idR*bin_res;
    idTH=idTH*bin_theta;
end

