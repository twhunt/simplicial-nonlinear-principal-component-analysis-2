function adjcncy_lst = edge_adjacency_list(edg_vrtx_inds)

num_edgs    = size(edg_vrtx_inds, 1);
adjcncy_lst = cell(num_edgs, 1);

if num_edgs == 0
    return
end

%is_adjcnt = false(num_edgs, 1);

for edg_i=1:num_edgs
    
    is_adjcnt = ...
        edg_vrtx_inds(edg_i, 1) == edg_vrtx_inds(:, 1) ...
        | edg_vrtx_inds(edg_i, 2) == edg_vrtx_inds(:, 1) ...
        | edg_vrtx_inds(edg_i, 1) == edg_vrtx_inds(:, 2) ...
        | edg_vrtx_inds(edg_i, 2) == edg_vrtx_inds(:, 2);
    
    is_adjcnt(edg_i) = false;

    adjcncy_lst{edg_i} = find(is_adjcnt);
    
end

