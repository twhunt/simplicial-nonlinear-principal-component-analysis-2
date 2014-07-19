function [intl_tri_is_vbl] = intl_tris_are_vbl(...
    tri_vrtx_inds, ...
    edg_vrtx_inds, ...
    vrtx_crdnts, ...
    intl_edg_srfc_inds, ...
    intl_vrtx_inds, ...
    srfc_crdnts, ...
    nrbyVrtxDstncTol, ...
    non_adj_tri_dist_tol, ...
    non_adj_tri_dist_tol_ovrlp)

%the vertex coordinates of the intial triangles are indexed by 
%intl_edg_srfc_inds and intl_vrtx_inds
%the number of initial triangles is equal to the number of entries in
%intl_vrtx_inds.
%If there are two triangles, intl_edg_srfc_inds is the common edge (i.e.,
%it holds the indices into srfc_crdnts of the two vertices common to both
%triangles)

%Basic sanity check:
%Vertices of shared edge are distinct
%non-shared vertices when there are 2 initial triangles are distinct



num_intl_tris = numel(intl_vrtx_inds);
assert(num_intl_tris <= 2);

intl_tri_is_vbl = false(1, num_intl_tris);
if isempty(intl_tri_is_vbl)
    return
end


if intl_edg_srfc_inds(1)==intl_edg_srfc_inds(2)
    
    warning('Vertices of initial edge are not distinct')
    intl_tri_is_vbl(1:end) = false;
    return
    
end


if num_intl_tris == 2 && intl_vrtx_inds(1) == intl_vrtx_inds(2)
    
    %the non-shared vertex indices are not distinct, i.e. two initial
    %triangles coincide
    intl_tri_is_vbl(2) = false;
    num_intl_tris = 1;
    
end

%check if candidate vertex coincides with the edge vertex, if so,
%don't bother checking for conflict
for tri_i=1:num_intl_tris
    
    intl_tri_is_vbl(tri_i) = ...
        ~any(intl_vrtx_inds(tri_i) == intl_edg_srfc_inds);
    
end

if all(~intl_tri_is_vbl)
    %all initial triangles are inviable due to non-uniqueness of vertices
    return;
end


if num_intl_tris == 2 && all(intl_tri_is_vbl)
    %test if the two initial triangles conflict with each other
   
    tmp_vrtx_crdnts = ...        
        srfc_crdnts(:, [intl_edg_srfc_inds(:)' , intl_vrtx_inds(2)]);
    
    tmp_cand_tri_vrtx_crdnts = ...
        srfc_crdnts(:, [intl_edg_srfc_inds(:)' , intl_vrtx_inds(1)]);
    
    tmp_cand_tri_vrtx_is_shrd = [true true false];
    
    tmp_tri_vrtx_inds = [1 2 3];

    tmp_exstng_tri_shrd_vrtx_inds = [1 2 0];
    
    tmp_exstng_tri_num_shrd_vrtxs = 2;
    
    [cand_tri_ovrlap, exstng_tri_ovrlap, intl_tri_is_vbl(2)] ...
    = cand_tri_cnflct2(...
    tmp_cand_tri_vrtx_crdnts, ...
    tmp_cand_tri_vrtx_is_shrd, ...
    tmp_tri_vrtx_inds, ...
    tmp_exstng_tri_shrd_vrtx_inds, ...
    tmp_exstng_tri_num_shrd_vrtxs, ...
    tmp_vrtx_crdnts, ...
    non_adj_tri_dist_tol, ...
    non_adj_tri_dist_tol_ovrlp);
    
end


if any(intl_tri_is_vbl);
    
    shrdEdgMdpntCrdnts = mean(srfc_crdnts(:, [intl_edg_srfc_inds(:)']), 2);

    isa_nrby_vrtx = is_nrby_vrtx(...
        vrtx_crdnts, ...
        shrdEdgMdpntCrdnts, ...
        nrbyVrtxDstncTol);
    
    nghbrhd_vrtx_inds = find(isa_nrby_vrtx);
    
    is_nrby_tri = false(size(tri_vrtx_inds, 1), 1);

    for vrtx_i=1:numel(nghbrhd_vrtx_inds)
        
        is_nrby_tri = is_nrby_tri ...
            | nghbrhd_vrtx_inds(vrtx_i) == tri_vrtx_inds(:, 1) ...
            | nghbrhd_vrtx_inds(vrtx_i) == tri_vrtx_inds(:, 2) ...
            | nghbrhd_vrtx_inds(vrtx_i) == tri_vrtx_inds(:, 3);
        
    end

    nghbrhd_frnt_tri_inds = find(is_nrby_tri);

    if ~isempty(nghbrhd_frnt_tri_inds)
        %Check if viable initial triangle conflicts with nearby portion of
        %the triangulation
        for tri_i=1:find(intl_tri_is_vbl)
            
            tmp_cand_tri_vrtx_crdnts = ...
                srfc_crdnts(:, ...
                [intl_edg_srfc_inds(:)' , intl_vrtx_inds(tri_i)]);
            
            tmp_cand_tri_vrtx_is_shrd = [false false false];
            
            tmp_exstng_tri_shrd_vrtx_inds = ...
                zeros(numel(nghbrhd_frnt_tri_inds), 3);
            
            tmp_exstng_tri_num_shrd_vrtxs = ...
                zeros(numel(nghbrhd_frnt_tri_inds), 1);
            
            [cand_tri_ovrlap, exstng_tri_ovrlap, intl_tri_is_vbl(tri_i)]...
                = cand_tri_cnflct2(...
                tmp_cand_tri_vrtx_crdnts, ...
                tmp_cand_tri_vrtx_is_shrd, ...
                tri_vrtx_inds, ...
                tmp_exstng_tri_shrd_vrtx_inds, ...
                tmp_exstng_tri_num_shrd_vrtxs, ...
                vrtx_crdnts, ...
                non_adj_tri_dist_tol, ...
                non_adj_tri_ovrlap_dist_tol);

        end
    end
end

end