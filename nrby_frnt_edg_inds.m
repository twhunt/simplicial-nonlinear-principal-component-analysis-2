function nrby_frnt_edg_inds = nrby_frnt_edg_inds(...
    nrby_vrtx_inds, frnt_edg_inds, edg_vrtx_inds)

if isempty(nrby_vrtx_inds)
    
    nrby_frnt_edg_inds = [];
    
else
    
    is_nrby_frnt_edg = false(size(frnt_edg_inds));
    
    for k=1:numel(is_nrby_frnt_edg)
        
        is_nrby_frnt_edg(k) = ...
            any(nrby_vrtx_inds == edg_vrtx_inds(frnt_edg_inds(k), 1)) ...
            || any(nrby_vrtx_inds == edg_vrtx_inds(frnt_edg_inds(k), 2));
        
    end
    
    nrby_frnt_edg_inds = frnt_edg_inds(is_nrby_frnt_edg);
    
end

