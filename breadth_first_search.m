function vrtx_q = breadth_first_search(vrtx_adjcncy_lst, intl_vrtx)

num_vrtxs = numel(vrtx_adjcncy_lst);

vrtx_q = zeros(num_vrtxs, 1);

if num_vrtxs == 0
    return
end

vrtx_unvstd = true(num_vrtxs,1);

vrtx_q_pop_pstn = 1;
unvstd_vrtx_ind = intl_vrtx;

unvstd_vrtx_exsts = true;

while unvstd_vrtx_exsts

    vrtx_q(vrtx_q_pop_pstn) = unvstd_vrtx_ind;
    vrtx_q_enq_pstn = vrtx_q_pop_pstn + 1;
    vrtx_unvstd(unvstd_vrtx_ind) = false;
    
    while vrtx_q_pop_pstn < vrtx_q_enq_pstn
        
        crrnt_vrtx =  vrtx_q(vrtx_q_pop_pstn);
        vrtx_q_pop_pstn = vrtx_q_pop_pstn + 1;
        
        adjcnt_vrtx_unvstd = vrtx_unvstd(vrtx_adjcncy_lst{crrnt_vrtx});
        num_unvstd_adjcnt_vrtx = sum(adjcnt_vrtx_unvstd);
        if num_unvstd_adjcnt_vrtx ~= 0
            
            vrtx_unvstd(...
                vrtx_adjcncy_lst{crrnt_vrtx}(adjcnt_vrtx_unvstd)) = false;
            
            vrtx_q(...
                vrtx_q_enq_pstn:(vrtx_q_enq_pstn+num_unvstd_adjcnt_vrtx-1))...
                = vrtx_adjcncy_lst{crrnt_vrtx}(adjcnt_vrtx_unvstd);
                                    
            vrtx_q_enq_pstn = vrtx_q_enq_pstn + num_unvstd_adjcnt_vrtx;
        end        
        
    end

    unvstd_vrtx_ind = find(vrtx_unvstd(:),1, 'first');
    unvstd_vrtx_exsts = ~isempty(unvstd_vrtx_ind);    
    
end