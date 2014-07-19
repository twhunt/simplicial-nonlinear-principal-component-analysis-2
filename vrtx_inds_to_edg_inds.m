function edgInds = vrtx_inds_to_edg_inds(vrtxInd1, vrtxInd2, edgVrtxInds)

vrtx1inEdg = any(vrtxInd1 == edgVrtxInds(:,1:2), 2);
vrtx1EdgInds = find(vrtx1inEdg);

vrtx12inEdg = any(vrtxInd2 == edgVrtxInds(vrtx1inEdg,1:2), 2);

edgInds = vrtx1EdgInds(vrtx12inEdg);