function [...
    isclsCnstrnt, isclsObjctv, ...
    sphrCnstrnt, ...
    Q1objctv, Q2objctv, ...
    Q1MtrxVctrPrdct, Q2MtrxVctrPrdct, ...
    ineqCnstrntLHS, ineqCnstrntRHS] ...
    = optmztnCnstrntsObjctv(...
    P1DmntEgnVctrs, P2DmntEgnVctrs, ...
    P1DmntEgnVals, P2DmntEgnVals, ...
    evalOffset1, evalOffset2, ...
    actvEdgVrtx1Crdnts, actvEdgVrtx2Crdnts, actvEdgIntrVrtxCrdnts, ...
    sphrRads)

Q1WeakEvals = 1./(P1DmntEgnVals + evalOffset1);
Q2WeakEvals = 1./(P2DmntEgnVals + evalOffset2);

Q1DmntEval = 1/evalOffset1;
Q2DmntEval = 1/evalOffset2;

Q1objctv = @(x) QMtrcDecomp(x - actvEdgVrtx1Crdnts(:), ...
    P1DmntEgnVctrs, Q1WeakEvals, Q1DmntEval);

Q2objctv = @(x) QMtrcDecomp(x - actvEdgVrtx2Crdnts(:), ...
    P2DmntEgnVctrs, Q2WeakEvals, Q2DmntEval);

Q1MtrxVctrPrdct = @(x) QMtrxVctrPrdctDecomp(x, ...
    P1DmntEgnVctrs, Q1WeakEvals, Q1DmntEval);

Q2MtrxVctrPrdct = @(x) QMtrxVctrPrdctDecomp(x, ...
    P2DmntEgnVctrs, Q2WeakEvals, Q2DmntEval);

isclsCnstrnt = @(x) isoscelesConstraint(x, Q1objctv, Q2objctv);

isclsObjctv  = @(x) isoscelesObjctv(x, Q1objctv, Q2objctv);

actvEdgMdpntCrdnts = .5*(actvEdgVrtx1Crdnts(:) + actvEdgVrtx2Crdnts(:));

sphrCnstrnt = @(x) sphereConstraint(x - actvEdgMdpntCrdnts, sphrRads);

actvEdgDsplcmnt = actvEdgVrtx2Crdnts(:) - actvEdgVrtx1Crdnts(:);
actvEdgMdpntIntrDsplcmnt = actvEdgMdpntCrdnts - actvEdgIntrVrtxCrdnts(:);

actvEdgDsplcmntNormSqrd = actvEdgDsplcmnt'*actvEdgDsplcmnt;

ineqCnstrntLHS = ...
    actvEdgMdpntIntrDsplcmnt ...
    - ((actvEdgDsplcmnt'*actvEdgMdpntIntrDsplcmnt)/actvEdgDsplcmntNormSqrd)...
    *actvEdgDsplcmnt;

ineqCnstrntLHS = (-1/norm(ineqCnstrntLHS))*ineqCnstrntLHS(:)';
ineqCnstrntRHS = ineqCnstrntLHS*actvEdgMdpntCrdnts;
    
% orgn = zeros(3,1);
% Q1objctv(orgn)
% Q2objctv(orgn)
% isclsCnstrnt(orgn)
% isclsObjctv(orgn)
