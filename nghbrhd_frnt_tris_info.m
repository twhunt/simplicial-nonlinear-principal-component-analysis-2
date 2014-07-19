function [...
    cand_tri_vrtx_is_shrd, ...
    exstng_tri_shrd_vrtx_inds, ...
    exstng_tri_num_shrd_vrtxs]...
    = nghbrhd_frnt_tris_info(...
    cndt_tri_exstng_edg_vrtx_inds, ...
    cndt_vrtx_ind, ...
    nghbrhd_frnt_tri_vrtx_inds)
    
num_nghbrhd_tris = size(nghbrhd_frnt_tri_vrtx_inds, 1);

cand_tri_vrtx_is_shrd     = false(num_nghbrhd_tris, 3);
exstng_tri_shrd_vrtx_inds = zeros(num_nghbrhd_tris, 3);
exstng_tri_num_shrd_vrtxs = zeros(num_nghbrhd_tris, 1);


for cand_vrtx_i=1:2
    
    is_crrnt_vrtx_ind = ...
        cndt_tri_exstng_edg_vrtx_inds(cand_vrtx_i) ...
        == nghbrhd_frnt_tri_vrtx_inds;

    exstng_tri_shrd_vrtx_inds(is_crrnt_vrtx_ind) = cand_vrtx_i;
    
    [rows_with_crrnt_vrtx, cols_with_crrnt_vrtx] = find(is_crrnt_vrtx_ind);
    
    cand_tri_vrtx_is_shrd(rows_with_crrnt_vrtx, cand_vrtx_i) = true;
    
end


if ~isempty(cndt_vrtx_ind)
    
    is_crrnt_vrtx_ind = cndt_vrtx_ind == nghbrhd_frnt_tri_vrtx_inds;

    exstng_tri_shrd_vrtx_inds(is_crrnt_vrtx_ind) = 3;
        
    [rows_with_crrnt_vrtx, cols_with_crrnt_vrtx] = find(is_crrnt_vrtx_ind);
    
    cand_tri_vrtx_is_shrd(rows_with_crrnt_vrtx, 3) = true;
 
end

exstng_tri_num_shrd_vrtxs = sum(cand_tri_vrtx_is_shrd, 2);

%assert(all(sum(cand_tri_vrtx_is_shrd(exstng_tri_num_shrd_vrtxs==2, :),2)));
%assert(all(sum(cand_tri_vrtx_is_shrd(exstng_tri_num_shrd_vrtxs==1, :),2)));
%sum(cand_tri_vrtx_is_shrd(exstng_tri_num_shrd_vrtxs==0, :))