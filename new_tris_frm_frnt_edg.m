function [triVrtxInds, numTris, edgVrtxInds, numEdgs, ...
    invblEdgInds, numInvblEdgInds] ...
    = new_tris_frm_frnt_edg(...
    srfcCrdnts, ...
    vrtx_ind_to_srfc_pt_ind, ...
    vrtx_crdnts, ...
    edg_fifo, ...
    triVrtxInds, ...
    edgVrtxInds, ...
    numTris, ...
    numEdgs, ...
    P_dmnt_egnvctrs, ...
    P_dmnt_egnvals, ...
    emprcl_drctn_crrltn_eval_bias, ...
    nnzEvals, ...
    emprcl_drctn_egn_sprs_algrthm, ...
    eigsOpts, ...
    invblEdgInds, ...
    numInvblEdgInds, ...
    non_adj_tri_dist_tol, ...
    non_adj_tri_ovrlap_dist_tol, ...
    new_tri_max_edg_lngth, ...
    plot_hndls, ...
    plot_frqncy, ...
    vrtx_crdnts_3D)


ttl_num_new_edgs = 0;
ttl_num_new_tris = 0;

fifo_ind = numel(edg_fifo);

cndt_new_edgs = zeros(2,3);
is_exstng_edg = false(1,2);

triVrtxInds(1:numTris, :) = sort(triVrtxInds(1:numTris, :), 2);

vrtxNghbrhdIsEmpty = false(1,2);

%dstnc_objctvs_decomp = cell(1,2);


%initialize is_frnt_edg to indicate that all edges are front edges
%entries in is_frnt_edg gets get set to correct values at top of while loop
%is_frnt_edg = true(size(edg_vrtx_inds,1), 1);

while fifo_ind ~= 0
    
    
    %dbg_new_tri  = zeros(1,3);
    dbg_new_edgs = zeros(2,3);
    
    crrnt_edg_is_frnt_edg = false;
    crrnt_edg_is_invbl    = true;
    stack_is_empty        = false;
    
    %inefficient to recompute at the top of the loop
    %is_frnt_edg should be updated at the bottom of the loop
    is_frnt_edg ...
        = edge_blngs_to_xctly_one_tri(...
        edgVrtxInds(1:numEdgs, :), triVrtxInds(1:numTris, :));
    
    while (~crrnt_edg_is_frnt_edg || crrnt_edg_is_invbl ...
            || any(vrtxNghbrhdIsEmpty)) && ~stack_is_empty

        %\/ DEBUG \/
        %warning('Setting current front edge to 1798 for debug')
        %prmtn=[fifo_ind 120]
        %edg_fifo(prmtn) = [1798 edg_fifo(fifo_ind)];
        %actv_edg_ind = 1488;
        %/\ DEBUG /\
        
        actv_edg_ind = edg_fifo(fifo_ind);
        
        
        edg_fifo(fifo_ind) = 0; %not strictly neccesary to zero out edge 
        %index in array, but helpful for debugging
        fifo_ind = fifo_ind - 1;
        disp(['fifo_ind: ' num2str(fifo_ind)]);
        
        stack_is_empty        = fifo_ind == 0;
        crrnt_edg_is_invbl    ...
            = any(actv_edg_ind == invblEdgInds(1:numInvblEdgInds));
        crrnt_edg_is_frnt_edg = is_frnt_edg(actv_edg_ind);
        
        %Edges created during seam sewing stage may not have empirical
        %direction covariance matrices computed for their vertices yet
        for vrtx_i=1:2

            vrtxNghbrhdIsEmpty(vrtx_i) = ...
                isempty(P_dmnt_egnvctrs{actv_edg_ind, vrtx_i});
            
            if vrtxNghbrhdIsEmpty(vrtx_i)
               
                actvEdgLngth = ...
                    norm(...
                    vrtx_crdnts(:, edgVrtxInds(actv_edg_ind, 1)) ...
                    - vrtx_crdnts(:, edgVrtxInds(actv_edg_ind, 2)));

                is_in_nghbrhd = in_nghbrhd_crdnts(...
                    srfcCrdnts, ...
                    vrtx_crdnts(:, edgVrtxInds(actv_edg_ind, vrtx_i)), ...
                    actvEdgLngth);
                
                %take out surface data point that coincides with edge
                %vertex, so there isn't a zero vector when summing the
                %direction vector outer products
                is_in_nghbrhd(...
                    vrtx_ind_to_srfc_pt_ind(...
                    edgVrtxInds(actv_edg_ind, vrtx_i))) ...
                    = false;
                
                num_pts_in_ngbrhd          = sum(is_in_nghbrhd);
                vrtxNghbrhdIsEmpty(vrtx_i) = num_pts_in_ngbrhd == 0;
                
                if ~vrtxNghbrhdIsEmpty(vrtx_i)
                    
                    [P_dmnt_egnvctrs{actv_edg_ind, vrtx_i}, ...
                        P_dmnt_egnvals{actv_edg_ind, vrtx_i}] ...
                        = ...
                        pnts_to_egn_dcmp(...
                        vrtx_crdnts(:, edgVrtxInds(actv_edg_ind, vrtx_i)), ...
                        srfcCrdnts(:,is_in_nghbrhd), ...
                        nnzEvals, ...
                        emprcl_drctn_egn_sprs_algrthm, ...
                        eigsOpts);
                    
                end
                
                
            end
            
            
            
        end
                
        if any(vrtxNghbrhdIsEmpty)
            
            if numInvblEdgInds == numel(invblEdgInds)               
                invblEdgInds((end+1):(end+numel(invblEdgInds))) = 0;
            end
            
            numInvblEdgInds = numInvblEdgInds + 1;
            invblEdgInds(numInvblEdgInds) = actv_edg_ind;
            
        end
        
        
    end
    
    if (~crrnt_edg_is_frnt_edg || crrnt_edg_is_invbl ...
            || any(vrtxNghbrhdIsEmpty))
        %break out of outer while loop
        continue;
    end
        
    
    %unique list of all vertices belonging to a front edge
    frnt_vrtx_inds = unique(edgVrtxInds(is_frnt_edg,1:2));
    
    frnt_edg_inds = find(is_frnt_edg);
    
    actv_edg_vrtx_inds = edgVrtxInds(actv_edg_ind, 1:2);
    %actv_edg_vrtx_inds = sort(actv_edg_vrtx_inds);
    
    actv_edg_crdnts = vrtx_crdnts(:, actv_edg_vrtx_inds);
    
    %\/ DEBUG \/
    if exist('dbg_actv_edg_hndl', 'var') && ishandle(dbg_actv_edg_hndl)
        delete(dbg_actv_edg_hndl);
    end
    
    actv_edg_crdnts_3D = vrtx_crdnts_3D(:, actv_edg_vrtx_inds);
    dbg_actv_edg_hndl = plot3(...
        actv_edg_crdnts_3D(1, :), ...
        actv_edg_crdnts_3D(2, :), ...
        actv_edg_crdnts_3D(3, :), 'r-o');    
    %/\ DEBUG /\

    
    for k=1:2
        
        QMtrc{k} = @(x) QMtrcDecomp(x - actv_edg_crdnts(:,k), ...
            P_dmnt_egnvctrs{actv_edg_ind, k}, ...
            1./(P_dmnt_egnvals{actv_edg_ind, k} + emprcl_drctn_crrltn_eval_bias), ...
            1/emprcl_drctn_crrltn_eval_bias);
        
        %dstnc_objctvs_decomp{k} = @(x) dstnc_objctv_decomp(...
        %    x-actv_edg_crdnts(:,k), ...
        %    P_dmnt_egnvctrs{actv_edg_ind, k}, ...
        %    P_dmnt_egnvals{actv_edg_ind, k}, ...
        %    emprcl_drctn_crrltn_eval_bias);
        
    end
        
    
    actv_edg_mdpnt_crdnts = mean(actv_edg_crdnts, 2);
    
    isa_nrby_vrtx = is_nrby_vrtx(...
        vrtx_crdnts, ...
        actv_edg_mdpnt_crdnts, ...
        2.5*new_tri_max_edg_lngth);
    
    nrst_vrtx_inds = find(isa_nrby_vrtx);
        
    nghbrhd_frnt_edg_inds = nrby_frnt_edg_inds(...
        nrst_vrtx_inds, frnt_edg_inds, edgVrtxInds(1:numEdgs, :));
    
    nghbrhd_frnt_tri_inds = nrby_frnt_tri_inds(...
        triVrtxInds(1:numTris, :), ...
        edgVrtxInds(1:numEdgs, :), ...
        nghbrhd_frnt_edg_inds);
    
    %     nghbrhd_frnt_tri_
    %     frnt_
    
    cand_tri_vrtx_crdnts(:, 1:2) = ...
        vrtx_crdnts(:, edgVrtxInds(actv_edg_ind, 1:2));
    
    %Take out vertices that don't belong to a front edge, so we only
    %generate candidate triangles composed of front edge vertices
    %Take out active edge vertices so we don't generate a candidate
    %triangle with two unique vertices
    
    %recompute nearby vertices so a vertex is nearby if and only if it
    %belongs to a candidate triangle whose edges are less than the maximum
    %edge length
    %isa_nrby_vrtx = is_nrby_vrtx(...
    %    vrtx_crdnts, ...
    %    actv_edg_mdpnt_crdnts, ...
    %    new_tri_max_edg_lngth);
    
    %nrst_vrtx_inds = find(isa_nrby_vrtx);

    
    nrst_vrtx_inds = intersect(nrst_vrtx_inds, frnt_vrtx_inds);
    nrst_vrtx_inds = nrst_vrtx_inds(...
        actv_edg_vrtx_inds(1) ~= nrst_vrtx_inds ...
        & actv_edg_vrtx_inds(2) ~= nrst_vrtx_inds);
    
    num_cndt_tris = numel(nrst_vrtx_inds);
    
    if true
        dbg_nrby_vrtx_hndl = plot3(...
            vrtx_crdnts_3D(1, nrst_vrtx_inds), ...
            vrtx_crdnts_3D(2, nrst_vrtx_inds), ...
            vrtx_crdnts_3D(3, nrst_vrtx_inds), ...
            'bo', 'MarkerSize', 10);
        
        delete(dbg_nrby_vrtx_hndl);
    end    
    
    %compute induced and euclidean edge lengths of candidate triangles
    %cndt_tri_edgs_sqrd_indcd_lngth = Inf(num_cndt_tris, 3);
    cndt_tri_edgs_indcd_lngth      = Inf(num_cndt_tris, 3);
    edg_eucldn_lngths              = Inf(num_cndt_tris, 2);
    for cand_tri_i=1:num_cndt_tris        
        
        for edg_i=1:2
                        
            edg_eucldn_lngths(cand_tri_i, edg_i) = norm(...
                actv_edg_crdnts(:, edg_i) ...
                - vrtx_crdnts(:, nrst_vrtx_inds(cand_tri_i)));
            
            %cndt_tri_edgs_sqrd_indcd_lngth(cand_tri_i, edg_i) = ...
            %    dstnc_objctvs_decomp{edg_i}(...
            %    vrtx_crdnts(:, nrst_vrtx_inds(cand_tri_i)));
            
            cndt_tri_edgs_indcd_lngth(cand_tri_i, edg_i) = ...
                QMtrc{edg_i}(vrtx_crdnts(:, nrst_vrtx_inds(cand_tri_i)));
            
        end
        
        %cndt_tri_edgs_sqrd_indcd_lngth(cand_tri_i, 3) = ...
        %    sum(cndt_tri_edgs_sqrd_indcd_lngth(cand_tri_i, 1:2));

        cndt_tri_edgs_indcd_lngth(cand_tri_i, 3) = ...
            sum(cndt_tri_edgs_indcd_lngth(cand_tri_i, 1:2));

    end
    
    %Don't consider any candidate triangle with an edge length longer than
    %the maximum acceptable edge length
    cndt_tri_legs_edg_lgnth_acctble = ...
        edg_eucldn_lngths(:,1) < new_tri_max_edg_lngth ...
        & edg_eucldn_lngths(:,2) < new_tri_max_edg_lngth;
    
    nrst_vrtx_inds = nrst_vrtx_inds(cndt_tri_legs_edg_lgnth_acctble);
    edg_eucldn_lngths = ...
        edg_eucldn_lngths(cndt_tri_legs_edg_lgnth_acctble, :);
    %cndt_tri_edgs_sqrd_indcd_lngth = ...
    %    cndt_tri_edgs_sqrd_indcd_lngth(cndt_tri_legs_edg_lgnth_acctble, :);
    
    cndt_tri_edgs_indcd_lngth = ...
        cndt_tri_edgs_indcd_lngth(cndt_tri_legs_edg_lgnth_acctble, :);
    
    num_cndt_tris = sum(cndt_tri_legs_edg_lgnth_acctble);
    
    if num_cndt_tris == 0
        
        if numInvblEdgInds == numel(invblEdgInds)
            invblEdgInds((end+1):(end+numel(invblEdgInds))) = 0;
        end

        numInvblEdgInds = numInvblEdgInds + 1;
        invblEdgInds(numInvblEdgInds) = actv_edg_ind;
        %there was not a nearby front edge vertex
        continue
        %mark this edge as inviable
    end
    
    cndt_tri_cnflcts = true(1, num_cndt_tris);
    
    %[sortd_dstncs sort_prmtn] = sort(cndt_tri_edgs_sqrd_indcd_lngth(:,3));
    [sortd_dstncs sort_prmtn] = sort(cndt_tri_edgs_indcd_lngth(:,3));
    nrst_vrtx_inds = nrst_vrtx_inds(sort_prmtn);
    edg_eucldn_lngths = edg_eucldn_lngths(sort_prmtn, :);
    %cndt_tri_edgs_sqrd_indcd_lngth = ...
    %    cndt_tri_edgs_sqrd_indcd_lngth(sort_prmtn, :);
    cndt_tri_edgs_indcd_lngth = ...
        cndt_tri_edgs_indcd_lngth(sort_prmtn, :);

    %min_sum_sqrd_dstnc_ind = 0;    
    %min_sum_sqrd_dstnc = Inf;
    min_sum_dstnc     = Inf;
    min_sum_dstnc_ind = 0;    

    for cand_tri_i=1:num_cndt_tris
        
        if exist('nrstv_vrtx_hndl') && ishandle(nrstv_vrtx_hndl)
            delete(nrstv_vrtx_hndl);
        end
    
        nrstv_vrtx_hndl = plot3(...
            vrtx_crdnts_3D(1, nrst_vrtx_inds(cand_tri_i)), ...
            vrtx_crdnts_3D(2, nrst_vrtx_inds(cand_tri_i)), ...
            vrtx_crdnts_3D(3, nrst_vrtx_inds(cand_tri_i)), ...
            'r+', 'MarkerSize', 10);

        
        %         if any(~cndt_tri_cnflcts) ...
        %                 && ...
        %                 cndt_tri_edgs_sqrd_indcd_lngth(cand_tri_i, 3) ...
        %                 > min_sum_sqrd_dstnc
        if any(~cndt_tri_cnflcts) ...
                && ...
                cndt_tri_edgs_indcd_lngth(cand_tri_i, 3) ...
                > min_sum_dstnc

            %Don't check if the current triangle conflicts
            %We've found a triangle that doesn't conflict whose sum of
            %squared leg lengths is less than the current triangle
            continue;
            
        end                    
        
        cand_tri_vrtx_inds = ...
            [actv_edg_vrtx_inds nrst_vrtx_inds(cand_tri_i)];
        
        srtd_cand_tri_vrtx_inds = ...
            sort(cand_tri_vrtx_inds, 'ascend');
        
        %Do not accept the current candidate triangle if it is an existing
        %triangle
        triInd = vrtx_inds_to_tri_ind(...
            srtd_cand_tri_vrtx_inds(1), ...
            srtd_cand_tri_vrtx_inds(2), ...
            srtd_cand_tri_vrtx_inds(3), ...
            edgVrtxInds);
        
        cndt_tri_cnflcts(cand_tri_i) = ~isempty(triInd);                
        
        if ~cndt_tri_cnflcts(cand_tri_i)
            
            cand_tri_vrtx_crdnts(:, 3) = ...
                vrtx_crdnts(:, nrst_vrtx_inds(cand_tri_i));            
            
            [cand_tri_vrtx_is_shrd, ...
                exstng_tri_shrd_vrtx_inds, ...
                exstng_tri_num_shrd_vrtxs]...
                = nghbrhd_frnt_tris_info(...
                actv_edg_vrtx_inds, ...
                nrst_vrtx_inds(cand_tri_i), ...
                triVrtxInds(nghbrhd_frnt_tri_inds, :));
            
            [cand_tri_ovrlap, exstng_tri_ovrlap, tmp_cand_tri_cnflcts] ...
                = cand_tri_cnflct2(...
                cand_tri_vrtx_crdnts, ...
                cand_tri_vrtx_is_shrd, ...
                triVrtxInds(nghbrhd_frnt_tri_inds ,:), ...
                exstng_tri_shrd_vrtx_inds, ...
                exstng_tri_num_shrd_vrtxs, ...
                vrtx_crdnts, ...
                -1*non_adj_tri_dist_tol, ...
                non_adj_tri_ovrlap_dist_tol);
            
            cndt_tri_cnflcts(cand_tri_i) = any(tmp_cand_tri_cnflcts);
            
            %if ~cndt_tri_cnflcts(cand_tri_i) ...
            %        && ...
            %        cndt_tri_edgs_sqrd_indcd_lngth(cand_tri_i, 3) ...
            %        < min_sum_sqrd_dstnc
            if ~cndt_tri_cnflcts(cand_tri_i) ...
                    && ...
                    cndt_tri_edgs_indcd_lngth(cand_tri_i, 3) ...
                    < min_sum_dstnc
                
                %min_sum_sqrd_dstnc = ...
                %    cndt_tri_edgs_sqrd_indcd_lngth(cand_tri_i, 3);
                %min_sum_sqrd_dstnc_ind = cand_tri_i;
                min_sum_dstnc = ...
                    cndt_tri_edgs_indcd_lngth(cand_tri_i, 3);
                min_sum_dstnc_ind = cand_tri_i;
                
            end
            
        end
        
    end
    
    if any(~cndt_tri_cnflcts)
        
        %cand_tri_i = min_sum_sqrd_dstnc_ind;
        cand_tri_i = min_sum_dstnc_ind;
        
        cndt_new_edgs(1, :) = ...
            [sort([actv_edg_vrtx_inds(1) nrst_vrtx_inds(cand_tri_i)])...
            actv_edg_vrtx_inds(2)];
        
        cndt_new_edgs(2, :) = ...
            [sort([actv_edg_vrtx_inds(2) nrst_vrtx_inds(cand_tri_i)])...
            actv_edg_vrtx_inds(1)];
        
        isShrdEdg = false;
        new_tri_exstng_edg_ind = [0 0];
        for k=1:2
            
            edgInds = vrtx_inds_to_edg_inds(...
                cndt_new_edgs(k, 1),...
                cndt_new_edgs(k, 2), ...
                edgVrtxInds);

            is_exstng_edg(k)    = ~isempty(edgInds);
            
            %assert(numel(edgInds) == 0 || numel(edgInds) == 1);
            isShrdEdg = numel(edgInds) > 1;
            if isShrdEdg
                %This is something of a cheat.
                %It should be impossible to add an edge that already
                %belongs to two triangles.
                warning(...
                    ['Tried to add triangle with edge that already ' ....
                    'belongs to two triangles']);
                numInvblEdgInds = numInvblEdgInds + 1;
                invblEdgInds(numInvblEdgInds) = actv_edg_ind;
                
                break; 
            
            end
                
            
            if is_exstng_edg(k)
                new_tri_exstng_edg_ind(k) = edgInds;
            end
                        
        end
        
        if isShrdEdg
            
            %Don't add the triangle
            continue
        
        end
                
        ttl_num_new_tris = ttl_num_new_tris + 1;
                        
        if numTris == size(triVrtxInds, 1)
            
            triVrtxInds = [triVrtxInds; zeros(size(triVrtxInds))];
            
        end
        
        numTris = numTris + 1;
        triVrtxInds(numTris, :) = ...
            sort([actv_edg_vrtx_inds(:)' nrst_vrtx_inds(cand_tri_i)]);
        
        %dbg_new_tri = tri_vrtx_inds(end,:);
                       
        num_new_edgs = 2 - sum(is_exstng_edg);
        
        if numEdgs + num_new_edgs > size(edgVrtxInds, 1)
            
            edgVrtxInds = [edgVrtxInds; zeros(size(edgVrtxInds))];
            
        end
        
        if fifo_ind + num_new_edgs > numel(edg_fifo)
            
            %double size of edge fifo stack, without changing from row to
            %column vector, or vice versa
            edg_fifo((end+1):(end+numel(edg_fifo))) = 0;
            
        end
        
        for k=find(~is_exstng_edg)
            
            numEdgs = numEdgs + 1;
            edgVrtxInds(numEdgs, :) = cndt_new_edgs(k, :);
            
            fifo_ind = fifo_ind + 1;
            %push the newly generated edge onto the stack
            %its index is equal to the numer of edges in
            %edg_vrtx_info
            edg_fifo(fifo_ind) = numEdgs;
            
            dbg_new_edgs(k, :) = cndt_new_edgs(k, :);
            
            %end
            
        end
        
        %end
        ttl_num_new_edgs = ttl_num_new_edgs + num_new_edgs;
        
        %frnt_cycl_edg_inds = frnt_cycl_basis(edg_vrtx_inds, tri_vrtx_inds);
        
        %is_frnt_edg(end:(end+num_new_edgs)) = true;
        
        if mod(numTris, plot_frqncy) == 0
            
            plot_hndls.tris_hndl = plot_tris(...
                triVrtxInds(1:numTris, :), ...
                vrtx_crdnts_3D(1,:), ...
                vrtx_crdnts_3D(2,:), ...
                vrtx_crdnts_3D(3,:), ...
                plot_hndls.tris_hndl);
            
            drawnow;
        end
        
    else
        
        disp('Did not add a triangle')
        if numInvblEdgInds == numel(invblEdgInds)
            invblEdgInds((end+1):(end+numel(invblEdgInds))) = 0;
        end

        numInvblEdgInds = numInvblEdgInds + 1;
        invblEdgInds(numInvblEdgInds) = actv_edg_ind;
        
    end
     
end

disp([...
    'added ' num2str(ttl_num_new_edgs) ...
    ' edges and ' num2str(ttl_num_new_tris) ' triangles'])

if ttl_num_new_tris > 0
    %warning('Marking all edges as valid. Likely inefficient!')
    invblEdgInds    = [];
    numInvblEdgInds = 0;
else
    invblEdgInds = invblEdgInds(1:numInvblEdgInds);
end

end