function triInd = vrtx_inds_to_tri_ind(...
    vrtxInd1, vrtxInd2, vrtxInd3, triVrtxInds)


vrtx1InTri     = any(vrtxInd1 == triVrtxInds, 2);
vrtx1InTriInds = find(vrtx1InTri);

vrtx12InTri     = any(vrtxInd2 == triVrtxInds(vrtx1InTri, :), 2);
vrtx12InTriInds = vrtx1InTriInds(vrtx12InTri);

vrtx123InTri = any(vrtxInd3 == triVrtxInds(vrtx12InTriInds, :), 2);
triInd       = vrtx12InTriInds(vrtx123InTri);