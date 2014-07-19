function nrby_frnt_tri_inds = nrby_frnt_tri_inds(...
    tri_vrtx_inds, edg_vrtx_inds, nrby_frnt_edg_inds)

if isempty(nrby_frnt_edg_inds)
    
    nrby_frnt_tri_inds = [];
    
else
    
    nrby_frnt_tri_inds = zeros(size(nrby_frnt_edg_inds));
    
    tri_contains_frnt_edg = tri_contains_edg(...
        edg_vrtx_inds(nrby_frnt_edg_inds(1), 1), ...
        edg_vrtx_inds(nrby_frnt_edg_inds(1), 2), ...
        tri_vrtx_inds);
    
    frnt_tri_edg_ind = find(tri_contains_frnt_edg);
        
    assert(numel(frnt_tri_edg_ind) == 1);

    nrby_frnt_tri_inds(1) = frnt_tri_edg_ind;
    
    num_unique_tris_found = 1;
    for k=2:numel(nrby_frnt_edg_inds)
        
        tri_contains_frnt_edg = tri_contains_edg(...
            edg_vrtx_inds(nrby_frnt_edg_inds(k), 1), ...
            edg_vrtx_inds(nrby_frnt_edg_inds(k), 2), ...
            tri_vrtx_inds);
        
        frnt_tri_edg_ind = find(tri_contains_frnt_edg);
        
        assert(numel(frnt_tri_edg_ind) == 1);
        
        if ~any(frnt_tri_edg_ind ...
                == nrby_frnt_tri_inds(1:num_unique_tris_found))
            num_unique_tris_found = num_unique_tris_found + 1;

            nrby_frnt_tri_inds(num_unique_tris_found) = frnt_tri_edg_ind;

        end
            
    end
        
        
    end
    nrby_frnt_tri_inds = nrby_frnt_tri_inds(1:num_unique_tris_found);
end

