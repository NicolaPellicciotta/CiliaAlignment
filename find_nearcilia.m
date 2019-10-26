function [N_near] = find_nearcilia(fXX,fYY,d_lim)

dist_x= repmat(fXX(:), [1,numel(fXX(:))]);
dist_y= repmat(fYY(:), [1,numel(fYY(:))]);
DR= sqrt((dist_x - dist_x').^2 + (dist_y - dist_y').^2) ;
N_near=sum(DR<d_lim,2);


end

