function [...
    vrtxCrdnts, ...
    triVrtxInds, ...
    edgVrtxInds, ...
    numVrtxs, ...
    numTris, ...
    numEdgs, ...
    srfc_pt_ind_to_vrtx_ind, ...
    vrtx_ind_to_srfc_pt_ind,...
    intl_mnmzr_crdnts, ...
    P_dmnt_egnvctrs, ...
    P_dmnt_egnvals]...
    = advancing_front_main_loop(...
    srfc_crdnts, ...
    edg_fifo, ...
    edg_fifo_ind, ...
    vrtxCrdnts, ...
    triVrtxInds, ...
    edgVrtxInds, ...
    numVrtxs, ...
    numTris, ...
    numEdgs, ...
    srfc_pt_ind_to_vrtx_ind, ...
    vrtx_ind_to_srfc_pt_ind,...
    intl_mnmzr_crdnts, ...
    P_dmnt_egnvctrs, ...
    P_dmnt_egnvals, ...
    eigs_opts, ...
    SNPCA_params, ...
    plot_hndls)



cndt_vrtx_optns             = optimset('fmincon');
%cndt_vrtx_optns.Algorithm   = 'interior-point';
cndt_vrtx_optns.Algorithm   = 'active-set';
%cndt_vrtx_optns.GradObj     = 'on';
cndt_vrtx_optns.GradObj     = 'off';
cndt_vrtx_optns.GradConstr  = 'on';
cndt_vrtx_optns.Display     = 'off';
%cndt_vrtx_optns.Diagnostics = 'off';


%options for x0 unconstrained minimizer
%turn on user supplied gradient
x0_optns             = optimset('fmincon');
%find x0 by solving constrained minimization problem with squared
%differnce of induced distances as the objective, so the function 
%tolerance of the x0 problem should be square of the constraint tolerance
%of the constrained optimization problem that finds the new vertex
%coordinates
x0_optns.TolFun      = cndt_vrtx_optns.TolCon^2;
%x0_optns.Algorithm   = 'interior-point';
x0_optns.Algorithm   = 'active-set';
x0_optns.GradObj     = 'off';
%x0_optns.GradObj     = 'on';
x0_optns.GradConstr  = 'on';
%x0_optns.Diagnostics = 'on';
x0_optns.Display     = 'off';
%x0_optns.Display     = 'iter';


%num_surf_pts = size(srfc_crdnts, 2);
%data_dmnsn = size(srfc_crdnts, 1);

%distance as measured by induced metrics
dstnc_objctvs_decomp = cell(1,2);

is_adj_frnt_edg      = cell(1,3);
is_adj_frnt_edg_inds = cell(1,3);

cand_tri_vrtx_crdnts = zeros(size(vrtxCrdnts, 1), 3);

%plotting variables
ordinate_index = 3; %1 for x, 2 for y, 3 for z. Used for triangle colo
local_color_map = bone();

isNewEdg = false(2,1);

while edg_fifo_ind > 0
    
    %sav_cnt = sav_cnt + 1;
    %save(['saved_data/nt_' num2str(sav_cnt) '.mat'])

    if SNPCA_params.plot_frqncy ~= 0
        delete(plot_hndls.actv_edge_hndl(...
            ishandle(plot_hndls.actv_edge_hndl)));
        delete(plot_hndls.actv_edge_intr_hndl(ishandle(...
            plot_hndls.actv_edge_intr_hndl)));
        %delete(plot_hndls.front_hndls(plot_hndls.ishandle(...
        %    front_hndls)));
        delete(plot_hndls.cand_vert_and_surf_hndl(ishandle(...
            plot_hndls.cand_vert_and_surf_hndl)));
        delete(plot_hndls.x0_hndl(ishandle(...
            plot_hndls.x0_hndl)));                
    end
    
    %tic
    
    is_frnt_edg = ...
        edge_blngs_to_xctly_one_tri(...
        edgVrtxInds(1:numEdgs, :), ...
        triVrtxInds(1:numTris, :));
    
    tmp_frnt_edg_inds = find(is_frnt_edg);
    tmp_num_frnt_edgs = numel(tmp_frnt_edg_inds);
    %     for edg_i=tmp_num_frnt_edgs:-1:1
    %         tmp_frnt_h(edg_i) = plot3(...
    %             vrtxCrdnts(1,edgVrtxInds(tmp_frnt_edg_inds(edg_i), 1:2)),...
    %             vrtxCrdnts(2,edgVrtxInds(tmp_frnt_edg_inds(edg_i), 1:2)),...
    %             vrtxCrdnts(3,edgVrtxInds(tmp_frnt_edg_inds(edg_i), 1:2)),...
    %             'm-');
    %     end
    %     drawnow
    %     for edg_i=tmp_num_frnt_edgs:-1:1
    %         delete(tmp_frnt_h(edg_i));
    %     end
            
    
    %pop edges off the stack until the stack is empty, or popped edge is a
    %front edge
    crrnt_edg_is_frnt_edg = false;
    stack_is_empty = false;
    while ~stack_is_empty > 0 && ~crrnt_edg_is_frnt_edg
       
        actv_edge_ind          = edg_fifo(edg_fifo_ind);
        edg_fifo(edg_fifo_ind) = 0; %zero out edge index to aid debugging
        edg_fifo_ind           = edg_fifo_ind - 1;
        crrnt_edg_is_frnt_edg  = is_frnt_edg(actv_edge_ind);
        stack_is_empty         = edg_fifo_ind == 0;
        
    end
    
    if ~crrnt_edg_is_frnt_edg
        %condition at top of outer while loop will evaluate to false, and
        %advancing front stage will terminate
       continue; 
    end
        
    
    %edg_fifo(1:edg_fifo_ind)
    
    error_struct    = new_error_struct();
    cndt_vrtx_info  = new_cndt_vrtx_info(size(srfc_crdnts,1));
    
    vert_error      = false;
    
    cand_tri_ind    = NaN;
    
    %dbg_prvs_num_tris = dbg_num_tris;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START: get indices and coordinates of active edge's vertices and
    % interior vertex
    actv_vert_inds        = edgVrtxInds(actv_edge_ind, 1:2);
    
    %\/ DEBUG \/
    %     sortd_actv_vert_inds = sort(actv_vert_inds);
    %     if false%sortd_actv_vert_inds(1) == 807 && sortd_actv_vert_inds(2) == 811
    %         filePathStr = ['interleaved_data/problem_edge.mat'];
    %         save(filePathStr)
    %     elseif false
    %         filePathStr = ['interleaved_data/problem_edge.mat'];
    %         load(filePathStr);
    %         SNPCA_params.plot_frqncy = 0
    %         x0_optns.TolFun      = cndt_vrtx_optns.TolCon^2;
    %         x0_optns.GradObj     = 'on';
    %
    %         cndt_vrtx_optns.GradObj     = 'on';
    %         cndt_vrtx_optns.Display = 'iter-detailed';
    %         cndt_vrtx_optns.Diagnostics = 'on';
    %         %cndt_vrtx_optns.ScaleProblem = 'obj-and-constr';
    %         cndt_vrtx_optns.InitTrustRegionRadius = ...
    %             .25*6.5536e-01 %active edge length = 6.5536e-01
    %         cndt_vrtx_optns.InitBarrierParam = 10;
    %         cndt_vrtx_optns.PlotFcns = @optimplotconstrviolation;
    %         cndt_vrtx_optns.DerivativeCheck = 'on';
    %         cndt_vrtx_optns.Algorithm = 'active-set';
    %         cndt_vrtx_optns.MaxFunEvals = 1000;
    %         %P_dmnt_egnvctrs
    %     end
    %/\ DEBUG /\
    actv_edge_vert_coords = vrtxCrdnts(:, actv_vert_inds(1:2));

    actv_tri_ind = edge_ind_to_tri_inds(...
        actv_edge_ind, ...
        edgVrtxInds(1:numEdgs, :), ...
        triVrtxInds(1:numTris, :));
    assert(numel(actv_tri_ind) == 1)
    
    isActvEdgVrtx = ...
        actv_vert_inds(1) == triVrtxInds(actv_tri_ind, :) ...
        | actv_vert_inds(2) == triVrtxInds(actv_tri_ind, :);
    assert(sum(isActvEdgVrtx) == 2);
    
    intr_vert_ind         = triVrtxInds(actv_tri_ind, ~isActvEdgVrtx);
    actv_intr_coords      = vrtxCrdnts(:, intr_vert_ind);
    actv_edg_vctr         = ...
        actv_edge_vert_coords(:,2) - actv_edge_vert_coords(:,1);
    actv_edg_lngth        = norm(actv_edg_vctr);
    actv_edg_mdpt_crdnts  = ...
        actv_edge_vert_coords(:,2) + actv_edge_vert_coords(:, 1);
    actv_edg_mdpt_crdnts  = .5*actv_edg_mdpt_crdnts;
    
    disp(['active edge length: ' num2str(actv_edg_lngth)]);
    % END: get indices and coordinates of active edge's vertices and
    % interior vertex
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START: plot tessellation, front, and active edge
    if mod(numTris, SNPCA_params.plot_frqncy) == 0
        %delete(plot_hndls.actv_edge_hndl(...
        %    ishandle(plot_hndls.actv_edge_hndl)));
        %delete(plot_hndls.actv_edge_intr_hndl(ishandle(...
        %    plot_hndls.actv_edge_intr_hndl)));
        %delete(plot_hndls.front_hndls(plot_hndls.ishandle(...
        %    front_hndls)));
        %delete(plot_hndls.cand_vert_and_surf_hndl(ishandle(...
        %    plot_hndls.cand_vert_and_surf_hndl)));
        
        vrtx_crdnts_3D = ...
            SNPCA_params.rtn_mtrx(:,1:3)'*vrtxCrdnts(:, 1:numVrtxs);
        
        plot_hndls = plot_tris_actv_edg(...
            plot_hndls, ...
            triVrtxInds(1:numTris, :), ...
            edgVrtxInds(1:numEdgs, :), ...
            actv_edge_ind, ...
            vrtx_crdnts_3D);
        
        %color_indices = tri_FaceVertexCData(triVrtxInds(1:numTris, :),...
        %    vrtx_crdnts_3D(ordinate_index, :), ...
        %    size(local_color_map,1), ...
        %    min(vrtx_crdnts_3D(ordinate_index, :)),...
        %    max(vrtx_crdnts_3D(ordinate_index, :)));
        
        %set(plot_hndls.tris_hndl, ...
        %    'CDataMapping', 'direct', ...
        %    'FaceVertexCData', color_indices(:));
        
    end
    % END: plot tessellation, front, and active edge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START: place candidate vertex
    %try
    
    %DEBUG test
    %frnt_edg_adj_tris_driver(...
    %    edgVrtxInds(1:numEdgs, :), triVrtxInds(1:numTris, :), frnt_cycl_edg_inds{frnt_cycl_ind}, actv_edge_ind, ...
    %    vrtx_crds_x, vrtx_crds_y, vrtx_crds_z)
    %end DEBUG test
            
    for vrtx_i=1:2
        
        %calculate new direction empirical correlation matrix
        
        %get surface data points that are in the Euclidean
        %neighborhood of the current edge vertex
        is_in_nghbrhd = in_nghbrhd_crdnts(...
            srfc_crdnts, ...
            vrtxCrdnts(:, edgVrtxInds(actv_edge_ind, vrtx_i)), ...
            SNPCA_params.srch_rad_fac1*actv_edg_lngth);
        
        %dbg_h = plot3(...
        %    srfc_crdnts(1, is_in_nghbrhd), ...
        %    srfc_crdnts(2, is_in_nghbrhd), ....
        %    srfc_crdnts(3, is_in_nghbrhd), ....
        %    'ro');
        %delete(dbg_h);
        
        %take out surface data point that coincides with edge
        %vertex, so there isn't a zero vector when summing the
        %direction vector outer products
        is_in_nghbrhd(...
            vrtx_ind_to_srfc_pt_ind(...
            edgVrtxInds(actv_edge_ind, vrtx_i))) ...
            = false;
                
        vert_error = ~any(is_in_nghbrhd);
        
        if vert_error
            %there are no surface data points in the search sphere            
            warning('No surface data points in search sphere')
            break
            
        else
            
            %compute eigen decomposition of empirical direction
            %covariance matrices
            [P_dmnt_egnvctrs{actv_edge_ind, vrtx_i}, ...
                P_dmnt_egnvals{actv_edge_ind, vrtx_i}] = ...
                pnts_to_egn_dcmp(vrtxCrdnts(:, edgVrtxInds(actv_edge_ind, vrtx_i)), ...
                srfc_crdnts(:,is_in_nghbrhd), ...
                SNPCA_params.nnz_egnvals, ...
                SNPCA_params.emprcl_drctn_egn_sprs_algrthm, ...
                eigs_opts);
            
            if false
                tmp_h = plot3(...
                    srfc_crdnts(1, is_in_nghbrhd), ...
                    srfc_crdnts(2, is_in_nghbrhd), ...
                    srfc_crdnts(3, is_in_nghbrhd), 'ro');
                delete(tmp_h);
            end
            
        end
        
    end
      
    
    if vert_error
        %Problem computing empirical direction covariance matrices
        continue
    end

            
    %Set up metrics induced by perturbed inverse of empirical
    %covariance matrix
    
    %calculate constraint sphere radius as a weighted average of
    %the active edge length and the characteristic length
    cnstrnt_sphr_rds = ...
        (1-SNPCA_params.prfrd_cnstrnt_rds_wght)...
        *SNPCA_params.cnstrnt_rad_fac*actv_edg_lngth...
        + ...
        SNPCA_params.prfrd_cnstrnt_rds_wght...
        *SNPCA_params.prfrd_cnstrnt_rds;
    
    %Generate objective and constraint functions (and their gradients)
    %required by fmincon for generating prenudged new vertex
    %coordinates, and the initial point required by fmincon to generate
    %the new vertex
    [...
        isclsCnstrnt, isclsObjctv, ...
        sphrCnstrnt, ...
        Q1objctv, Q2objctv, ...
        Q1MtrxVctrPrdct, Q2MtrxVctrPrdct, ...
        ineqCnstrntLHS, ineqCnstrntRHS] ...
        = optmztnCnstrntsObjctv(...
        P_dmnt_egnvctrs{actv_edge_ind, 1}, ...
        P_dmnt_egnvctrs{actv_edge_ind, 2}, ...
        P_dmnt_egnvals{actv_edge_ind, 1}, ...
        P_dmnt_egnvals{actv_edge_ind, 2}, ...
        SNPCA_params.emprcl_drctn_crrltn_eval_bias, ...
        SNPCA_params.emprcl_drctn_crrltn_eval_bias, ...
        actv_edge_vert_coords(:, 1), ...
        actv_edge_vert_coords(:, 2), ...
        actv_intr_coords, ...
        cnstrnt_sphr_rds);
        
    
    %Find initial point satisfying constraints
    [x0, eqlty_cnstrnt_val, x0_exit_flag] = cnstrnt_x0( ...
        Q1MtrxVctrPrdct, ...
        isclsObjctv, ...
        sphrCnstrnt, ...
        actv_edg_mdpt_crdnts, ...
        cnstrnt_sphr_rds, ...
        ineqCnstrntLHS, ...
        ineqCnstrntRHS, ...
        actv_edge_vert_coords(:, 1), ...
        actv_edge_vert_coords(:, 2), ...
        actv_intr_coords, ...
        x0_optns);
        
    
    if x0_exit_flag <= 0;
        warning(...
            ['Trouble finding initial point satisfying ' ...
            'equality constraints. fmincon exit flag: %d']...
            , x0_exit_flag);
        
        continue;
    end
    
    
    %find minimizer of constrained minimization problem
    
    [cand_tri_vert_coords, Q1ObjctvVal, mnmzr_exit_flag] ...
        = gen_tri_vert_coords3_decomp(...
        Q1objctv, ...
        isclsCnstrnt, ...
        sphrCnstrnt, ...
        ineqCnstrntLHS, ...
        ineqCnstrntRHS, ...
        x0, ...
        cndt_vrtx_optns);
        
    if mnmzr_exit_flag <= 0;
        warning(...
            'Problem finding minimizer. fmincon exit flag: %d', ...
            mnmzr_exit_flag);
        
        continue;
        
    end
    
    if mod(numTris, SNPCA_params.plot_frqncy) == 0
        plot_hndls.cand_vert_and_surf_hndl = plot3(...
            cand_tri_vert_coords(1), ...
            cand_tri_vert_coords(2), ...
            cand_tri_vert_coords(3), ...
            'rx');
    end
    
    x0ObjctvVal    = Q1objctv(x0);
    mnmzrObjctvVal = Q1objctv(cand_tri_vert_coords); 
    disp(...
        ['X0 objective value = ' num2str(x0ObjctvVal) ....
        ' minimizer objective value = ' num2str(mnmzrObjctvVal)]);
    
    %save(...
    %    'pancake_subspace_info.mat', ...
    %    'srfc_crdnts', 'is_in_nghbrhd', ...
    %    'tmp_P_dmnt_egnvctrs', 'tmp_P_dmnt_egnvals', ...
    %    'actv_edge_vert_coords', 'actv_intr_coords', ...
    %    'NLPCA_params', ...
    %    'cand_tri_vert_coords');
        
            
        
    %candidate vertex was generated succesfuly
    %nudge the vertex to the nearest surface data point, and reject if
    %the distance between the original and nudged candidate vertices is
    %too great
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %START: nudge candidate vertex to nearest surface data point
    isInMaxNudgeNghbrhd = in_nghbrhd_crdnts(...
        srfc_crdnts, ...
        cand_tri_vert_coords, ...
        SNPCA_params.cand_vert_max_nudge_dist);
    
    if ~any(isInMaxNudgeNghbrhd)
       
        warning('No points within max nudge distance of minimizer')
        continue;   
        
    end
    
    
    mnmzrActvEdgMdpntDsplcmnt = ...
        cand_tri_vert_coords - actv_edg_mdpt_crdnts;
    mnmzrActvEdgMdpntDsplcmntDstnc = norm(mnmzrActvEdgMdpntDsplcmnt);
    
    is_in_nghbrhd = in_nghbrhd_crdnts(...
        srfc_crdnts, ...
        cand_tri_vert_coords, ...
        mnmzrActvEdgMdpntDsplcmntDstnc);
    
    if ~any(is_in_nghbrhd)

        warning('No surface points near enough minimizer to build Q metric')
        continue
        
    end
    
    
    %\/ build metric based on surface data points near initial
    %candidate vertex \/
    [cand_vert_P_dmnt_egnvctrs, cand_vert_P_dmnt_egnvals] = ...
        pnts_to_egn_dcmp(...
        cand_tri_vert_coords, ...
        srfc_crdnts(:, is_in_nghbrhd), ...
        SNPCA_params.nnz_egnvals, ...
        SNPCA_params.emprcl_drctn_egn_sprs_algrthm, ...
        eigs_opts);
    
    dstnc_objctv_decomp_cand_vert = @(x) QMtrcDecomp(...
        x-cand_tri_vert_coords, ...
        cand_vert_P_dmnt_egnvctrs, ...
        1./(cand_vert_P_dmnt_egnvals ...
        + SNPCA_params.emprcl_drctn_crrltn_eval_bias),...
        1/SNPCA_params.emprcl_drctn_crrltn_eval_bias);
    %/\ build metric based on surface data points near initial
    %candidate vertex /\
    
    %find nearest vertex in Q metric associated with initial candidate
    %vertex
    [minQDstnc, indrctNrstInd] = nearest_srfc_pt_Q2(...
        srfc_crdnts(:,isInMaxNudgeNghbrhd), dstnc_objctv_decomp_cand_vert);
    
    maxNudgeNghbrhdVrtxInds = find(isInMaxNudgeNghbrhd);
    nearest_surf_pt_ind     = maxNudgeNghbrhdVrtxInds(indrctNrstInd);
        
    cndt_vrtx_info.srfc_pt_ind = nearest_surf_pt_ind;
    
    cand_vert_nudge_dist = norm(...
        srfc_crdnts(:, nearest_surf_pt_ind) - cand_tri_vert_coords);
    
    %     %Was the minimizer nudged too far?
    %     vert_error = ...
    %         cand_vert_nudge_dist > SNPCA_params.cand_vert_max_nudge_dist;
    %
    %     if vert_error
    %
    %         warning(...
    %             ['Candidate vertex too far from surface data, ' ...
    %             'max distance "%e", actual distance "%e"'], ...
    %             SNPCA_params.cand_vert_max_nudge_dist, ...
    %             cand_vert_nudge_dist);
    %
    %         %indicate candidate vertex placement error
    %         error_struct.err_id = 'NLPCA:too_much_nudge';
    %
    %         continue;
    %
    %     end
    
    
    vert_error = ...
        srfc_pt_ind_to_vrtx_ind.isKey(nearest_surf_pt_ind) ...
        && (srfc_pt_ind_to_vrtx_ind(nearest_surf_pt_ind) ...
        == actv_vert_inds(1) ...
        || srfc_pt_ind_to_vrtx_ind(nearest_surf_pt_ind) ...
        == actv_vert_inds(2));
    
    if vert_error
        warning(['Nearest surface data point is a vertx of the ' ...
            'active edge.'])
        
        continue;
        
    end        
    %END: nudge candidate vertex to nearest surface data point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START: get tessellation info about candidate vertex
    %if ~vert_error
    frnt_edg_inds = find(is_frnt_edg);
    
    [cndt_vrtx_info.is_exstng_vrtx, ...
        cndt_vrtx_info.is_frnt_edg_vrtx,  ...
        cndt_vrtx_info.is_actv_edg_vrtx, ...
        cndt_vrtx_info.vrtx_ind] ...
        = srfc_pt_edg_info(...
        nearest_surf_pt_ind, actv_edge_ind, ...
        {frnt_edg_inds}, ...
        edgVrtxInds(1:numEdgs, :), srfc_pt_ind_to_vrtx_ind);
    
    cndt_vrtx_info.srfc_pt_ind      = nearest_surf_pt_ind;
    cndt_vrtx_info.crds             = ...
        srfc_crdnts(:, nearest_surf_pt_ind);
    % END: get tessellation info about candidate vertex
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %START: plot candidate vertex and nearest surface data point
    if mod(numTris, SNPCA_params.plot_frqncy) == 0
        
        cand_tri_vert_coords_3D = ...
            SNPCA_params.rtn_mtrx(:,1:3).'*cand_tri_vert_coords;
        cndt_vrtx_info_crdnts_3D = ...
            SNPCA_params.rtn_mtrx(:,1:3).'*cndt_vrtx_info.crds(:);
        
        %disp('Supressing plot of minimizer')
        delete(plot_hndls.cand_vert_and_surf_hndl(ishandle(plot_hndls.cand_vert_and_surf_hndl)));
        plot_hndls.cand_vert_and_surf_hndl = plot3(...
            [cand_tri_vert_coords_3D(1) cndt_vrtx_info_crdnts_3D(1)], ...
            [cand_tri_vert_coords_3D(2) cndt_vrtx_info_crdnts_3D(2)], ...
            [cand_tri_vert_coords_3D(3) cndt_vrtx_info_crdnts_3D(3)], ...
            'x-r');
        plot_hndls.x0_hndl = plot3(x0(1), x0(2), x0(3), 'ro');
        
        drawnow
        
        %\/ animation \/
        %set(gcf, 'renderer', 'painter'); drawnow;
        %print('-depsc', ['advancing_front_animation/cand_vert_' num2str(sav_cnt)])
        if false
            disp('Generating advancing front animation')
            if ~exist('advancing_front_plot_count', 'var')
                advancing_front_plot_count = 0;
            else
                advancing_front_plot_count = ...
                    advancing_front_plot_count + 1;
            end
            print(...
                '-dpng', '-r150', ...
                ['animation/advancing_front/' ...
                num2str(advancing_front_plot_count)]);
        end
        %/\ animation /\
    end
    %END: plot candidate vertex and nearest surface data point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % \/ Get list of non-shared vertices belonging to a front edge
    %adjacent to the active edge, and within the distance tolerance
    %specified by NLPCA_params.adj_vert_max_nudge_dist \/
    adjcnt_frnt_edg_cndt_vrtx_inds = cndt_adjcnt_vrtx_indxs(...
        actv_edge_ind, frnt_edg_inds, edgVrtxInds(1:numEdgs, :), ...
        cndt_vrtx_info.crds ,vrtxCrdnts(:, 1:numVrtxs), ...
        SNPCA_params.adj_vert_max_nudge_dist);
    % /\ Get list of non-shared vertices belonging to a front edge
    %adjacent to the active edge, and within the distance tolerance
    %specified by NLPCA_params.adj_vert_max_nudge_dist /\
    
    %candidate triangles are defined by the surface point nearest the
    %minimizer of the constrained optimization problem, and the
    %non-shared adjacent front edge vertices in
    %adjcnt_frnt_edg_cndt_vrtx_inds
    num_adj_cndt_vrtxs = numel(adjcnt_frnt_edg_cndt_vrtx_inds);
    num_cand_tris      = 1 + num_adj_cndt_vrtxs;
    
    
    %build list of info structs for candidate vertices
    %call cand_vert_error sequentially on entries of cand_tri_vrtx_info
    %so the nonadjacent candidate vertex is always checked last
    %add the adjacent existing candidate vertices in arbitrary order
    %(prefer smaller angles between adjacent active edges?)
    cand_tri_vrtx_info(num_adj_cndt_vrtxs+1) = cndt_vrtx_info;
    for k=1:num_adj_cndt_vrtxs
        
        cand_tri_vrtx_info(k) = cndt_vrtx_info;
        cand_tri_vrtx_info(k).is_exstng_vrtx   = true;
        cand_tri_vrtx_info(k).is_frnt_edg_vrtx = true;
        cand_tri_vrtx_info(k).is_actv_edg_vrtx = false;
        cand_tri_vrtx_info(k).vrtx_ind         = ...
            adjcnt_frnt_edg_cndt_vrtx_inds(k);
        cand_tri_vrtx_info(k).crds             = ...
            vrtxCrdnts(:, adjcnt_frnt_edg_cndt_vrtx_inds(k));
        
    end
    
    
    % END: place candidate vertex
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %if error_struct.actv_edg_is_vbl
    
    %Build list of front edges that share a vertex with the active edge
    isa_nrby_vrtx = is_nrby_vrtx(...
        vrtxCrdnts(:, 1:numVrtxs), ...
        cndt_vrtx_info.crds, ...
        SNPCA_params.nrby_vrtx_dstnc_tol);
    
    nghbrhd_vrtx_inds = find(isa_nrby_vrtx);
    
    is_nrby_tri = false(numTris, 1);
    
    for vrtx_i=1:numel(nghbrhd_vrtx_inds)
        
        is_nrby_tri = is_nrby_tri ...
            | nghbrhd_vrtx_inds(vrtx_i) == triVrtxInds(1:numTris, 1)...
            | nghbrhd_vrtx_inds(vrtx_i) == triVrtxInds(1:numTris, 2)...
            | nghbrhd_vrtx_inds(vrtx_i) == triVrtxInds(1:numTris, 3);
        
    end
    
    nghbrhd_frnt_tri_inds = find(is_nrby_tri);
    
    num_nghbrhd_frnt_tris = sum(is_nrby_tri);
    disp(['num nghbrhd tris = ' num2str(num_nghbrhd_frnt_tris) ' num tris = ' num2str(numTris)])
    
    cand_tri_vrtx_crdnts(:, 1:2) = ...
        vrtxCrdnts(:, edgVrtxInds(actv_edge_ind, 1:2));
    
    for cand_tri_i=1:num_cand_tris
        
        %check if the candidate triangle is an existing triangle
        if cand_tri_vrtx_info(cand_tri_i).is_exstng_vrtx
            
            triInd = vrtx_inds_to_tri_ind(...
                actv_edge_vert_coords(1), ...
                actv_edge_vert_coords(2), ...
                cand_tri_vrtx_info(cand_tri_i).vrtx_ind, ...
                triVrtxInds(1:numTris, :));
            
            vert_error = ~isempty(triInd);
            if vert_error
                %the candidate triangle is an existing triangle.
                continue;
            end
            
        end        
        
        cand_tri_vrtx_crdnts(:, 3) = ...
            cand_tri_vrtx_info(cand_tri_i).crds;
        
        [cand_tri_vrtx_is_shrd, ...
            exstng_tri_shrd_vrtx_inds, ...
            exstng_tri_num_shrd_vrtxs]...
            = nghbrhd_frnt_tris_info(...
            edgVrtxInds(actv_edge_ind, 1:2), ...
            cand_tri_vrtx_info(cand_tri_i).vrtx_ind, ...
            triVrtxInds(nghbrhd_frnt_tri_inds, :));
        
        [cand_tri_ovrlap, exstng_tri_ovrlap, cand_tri_cnflcts] ...
            = cand_tri_cnflct2(...
            cand_tri_vrtx_crdnts, ...
            cand_tri_vrtx_is_shrd, ...
            triVrtxInds(nghbrhd_frnt_tri_inds ,:), ...
            exstng_tri_shrd_vrtx_inds, ...
            exstng_tri_num_shrd_vrtxs, ...
            vrtxCrdnts(:, 1:numVrtxs), ...
            SNPCA_params.non_adj_tri_dist_tol, ...
            SNPCA_params.non_adj_tri_ovrlap_dist_tol);
        
        
        vert_error = any(cand_tri_cnflcts);
        if ~vert_error
            %Candidate triangle does NOT conflict with an existing
            %triangle.
            %Don't check any more triangles.
            break;
            
        end
        
    end
    
    if vert_error
        %all candidate triangles conflicted with the existing triangulation
        continue;
        
    end
        
        
    %Double available storage for triangles and vertex coordinates if
    %neccessary    
    if  ~cand_tri_vrtx_info(cand_tri_i).is_exstng_vrtx ...
            &&  numVrtxs == size(vrtxCrdnts, 2)
        
        %double size of matrix of vertex coordinates
        vrtxCrdnts = [vrtxCrdnts, zeros(size(vrtxCrdnts))];
        
    end
    
    if numTris == size(triVrtxInds, 1)
        
        %double the number of triangles that can be held
        triVrtxInds = [triVrtxInds; zeros(size(triVrtxInds))];
        
    end

    %\/ REFACTORED UPDATE EDGE AND TRIANGLE DATA STRUCTURES \/
    %Add the vertex of the new triangle that doesn't belong to the active
    %edge if it's new
    if cand_tri_vrtx_info(cand_tri_i).is_exstng_vrtx

        triNonActvEdgVrtxInd = cand_tri_vrtx_info(cand_tri_i).vrtx_ind;
        
    else

        numVrtxs             = numVrtxs + 1;
        triNonActvEdgVrtxInd = numVrtxs;
        
        vrtxCrdnts(:, triNonActvEdgVrtxInd) = ....
            cand_tri_vrtx_info(cand_tri_i).crds(:);
    
        srfc_pt_ind_to_vrtx_ind(nearest_surf_pt_ind) = ...
            triNonActvEdgVrtxInd;
        
        vrtx_ind_to_srfc_pt_ind(triNonActvEdgVrtxInd) = ...
            nearest_surf_pt_ind;
        
    end    

    %Add the new triangle
    numTris                = numTris + 1;
    triVrtxInds(numTris,:) = ...
        sort([edgVrtxInds(actv_edge_ind, 1:2), triNonActvEdgVrtxInd]);
    
    num_new_edges = 0;
    for edg_i=1:2
        
        %don't search for an existing edge if
        %cand_tri_vrtx_info(cand_tri_i).is_exstng_vrtx is false
        isExstngEdge = ...
            cand_tri_vrtx_info(cand_tri_i).is_exstng_vrtx ...
            &&...
            (any(...
            actv_vert_inds(edg_i) == edgVrtxInds(1:numVrtxs, 1) ...
            & triNonActvEdgVrtxInd == edgVrtxInds(1:numVrtxs, 2)) ...
            || ...
            any(...
            actv_vert_inds(edg_i) == edgVrtxInds(1:numVrtxs, 2) ...
            & triNonActvEdgVrtxInd == edgVrtxInds(1:numVrtxs, 1)));
        
        if ~isExstngEdge
           
            num_new_edges = num_new_edges + 1;
            
            if edg_i == 1
                tmp_i = 2;
            else
                tmp_i = 1;
            end
            
            numEdgs = numEdgs + 1;
            
            if numEdgs + num_new_edges > size(edgVrtxInds, 1)
                
                %double the number of edges that can be held
                edgVrtxInds = [edgVrtxInds; zeros(size(edgVrtxInds))];
                
                %Matlab pads the cell with []
                P_dmnt_egnvctrs{(numEdgs+size(P_dmnt_egnvctrs, 1)), 2} = [];
                P_dmnt_egnvals{(numEdgs+size(P_dmnt_egnvctrs, 1)), 2} = [];
            end

            
            edgVrtxInds(numEdgs,:) = [...
                sort([actv_vert_inds(edg_i) triNonActvEdgVrtxInd]) ...
                actv_vert_inds(tmp_i)];
            
                edg_fifo(edg_fifo_ind+1) = numEdgs;
                    edg_fifo_ind = edg_fifo_ind + 1;


        end
        
    end
    %/\ REFACTORED UPDATE EDGE AND TRIANGLE DATA STRUCTURES /\
    
%     if ~cand_tri_vrtx_info(cand_tri_i).is_frnt_edg_vrtx
%         %the candidate vertex is not an existing vertex and is
%         %not a vertex of the front
%         %one new vertex
%         %two new edges
%         
%         num_new_edges   = 2;
%         isNewEdg(1:end) = true;
%         
%         numVrtxs     = numVrtxs + 1;
%         new_vrtx_ind = numVrtxs;
%         
%         vrtxCrdnts(:, new_vrtx_ind) = ....
%             cand_tri_vrtx_info(cand_tri_i).crds(:);
%         
%         %associate the surface data point with the new vertex
%         srfc_pt_ind_to_vrtx_ind(nearest_surf_pt_ind) = ...
%             new_vrtx_ind;
%         
%         vrtx_ind_to_srfc_pt_ind(new_vrtx_ind) = nearest_surf_pt_ind;
%         
%         triNonActvEdgVrtxInd = new_vrtx_ind;
%         
%         %the 2 new edges
%         new_edg_vrtx_inds(1,:) = ...
%             [actv_vert_inds(1) new_vrtx_ind  actv_vert_inds(2)];
%         
%         new_edg_vrtx_inds(2,:) = ...
%             [actv_vert_inds(2) new_vrtx_ind  actv_vert_inds(1)];
%         
%     else
%         %the candidate vertex is an existing vertex in the
%         %front
%         %0, 1, or 2 new edges
% 
%         new_edg_vrtx_inds(1, :) = ...
%             [actv_vert_inds(1), ...
%             cand_tri_vrtx_info(cand_tri_i).vrtx_ind, ...
%             actv_vert_inds(2)];
% 
%         new_edg_vrtx_inds(2, :) = ...
%             [actv_vert_inds(2), ...
%             cand_tri_vrtx_info(cand_tri_i).vrtx_ind, ...
%             actv_vert_inds(1)];
% 
%         for edg_i=1:2
% 
%             edgInds = vrtx_inds_to_edg_inds(...
%                 new_edg_vrtx_inds(edg_i,1), ...
%                 new_edg_vrtx_inds(edg_i,2), ...
%                 edgVrtxInds(1:numEdgs, 1:2));
%             
%             assert(numel(edgInds) <= 1);
%             
%             isNewEdg(edg_i) = ~isempty(edgInds);
%         end
%                         
%         num_new_edges = sum(isNewEdg);
%                         
%         triNonActvEdgVrtxInd = cand_tri_vrtx_info(cand_tri_i).vrtx_ind;
%     end
%         
%     if numEdgs + num_new_edges > size(edgVrtxInds, 1)
%         
%         %double the number of edges that can be held
%         edgVrtxInds = [edgVrtxInds; zeros(size(edgVrtxInds))];
%         
%         %Matlab pads the cell with []
%         P_dmnt_egnvctrs{(numEdgs+size(P_dmnt_egnvctrs, 1)), 2} = [];        
%         P_dmnt_egnvals{(numEdgs+size(P_dmnt_egnvctrs, 1)), 2} = [];
%     end
%     
%     numTris                = numTris + 1;
%     triVrtxInds(numTris,:) = ...
%         [edgVrtxInds(actv_edge_ind, 1:2), triNonActvEdgVrtxInd];    
%     
%     edgVrtxInds(numEdgs+(1:num_new_edges), :) ...
%         = new_edg_vrtx_inds(isNewEdg, :);
%     
%     edg_fifo(edg_fifo_ind+(1:num_new_edges)) ...
%         = numEdgs+(1:num_new_edges);
%     
%     edg_fifo_ind = edg_fifo_ind + num_new_edges;
%     
%     numEdgs = numEdgs + num_new_edges;
    
       
    if SNPCA_params.plot_frqncy > 0
        
        vrtx_crdnts_3D = ...
            SNPCA_params.rtn_mtrx(:,1:3)'*vrtxCrdnts(:, 1:numVrtxs);
        plot_hndls.tris_hndl = plot_tris(triVrtxInds(1:numTris, :), ...
            vrtx_crdnts_3D(1,:), ...
            vrtx_crdnts_3D(2,:), ...
            vrtx_crdnts_3D(3,:), plot_hndls.tris_hndl);
        
    end
    
    %toc
    
end
