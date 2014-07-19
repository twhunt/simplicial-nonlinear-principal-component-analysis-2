function [x0, IsclsObjctvVal, exit_flag] = cnstrnt_x0( ...
    Q1MtrxVctrPrdct, ...
    isclsObjctv, ...
    sphrCnstrnt, ...
    CnstrntSphrCntrCrdnts, ...
    CnstrntSphrRds, ...
    ineqCnstrntLHS, ...
    ineqCnstrntRHS, ...
    actvEdgVrtxCrdnts1, ...
    actvEdgVrtxCrdnts2, ...
    actvEdgIntrCrdnts, ...
    x0_optns)

%generate x0 that:
%lies on the constraint sphere,
%is Q-orthogonal to the active edge
%lies on the correct side of the active edge 

% triBasis = [...
%     actvEdgVrtxCrdnts1(:) - actvEdgIntrCrdnts(:), ...
%     actvEdgVrtxCrdnts2(:) - actvEdgIntrCrdnts(:)];
% 
% 
% affineIneqCnstrntLHS = ineqCnstrntLHS*triBasis;
% affineIneqCnstrntRHS = ...
%     ineqCnstrntRHS - ineqCnstrntLHS*actvEdgIntrCrdnts(:);
% 
% affineSphrCnstrnt = @(y) sphrCnstrnt(...
%     affineCombo(y, triBasis, CnstrntSphrCntrCrdnts));
% 
% affineIsclsObjctv = @(y) isclsObjctv(...
%     affineCombo(y, triBasis, actvEdgIntrCrdnts(:)));
% 
% y0 = (CnstrntSphrRds/norm(triBasis(:, 1) + triBasis(:, 2)))*[1; 1];
% 
% % problem_decomp = struct(...
% %     'x0', y0, ...
% %     'Aineq', affineIneqCnstrntLHS, ...
% %     'bineq', affineIneqCnstrntRHS, ...
% %     'lb', [], ...
% %     'ub', [], ...
% %     'nonlcon', affineSphrCnstrnt , ...
% %     'objective', affineIsclsObjctv, ...
% %     'solver', 'fmincon', ...
% %     'options', x0_optns);
% % 
% % % problem_decomp.objective = dstnc_objctvs_decomp{1};
% % % problem_decomp.nonlcon   = sphr_isos_cnstrnt;
% % % problem_decomp.x0        = x0;
% % 
% % %get the vertex coordinates that solve the constrained minimization problem
% % [y0, affineIsclsObjctvVal, exit_flag] = fmincon(problem_decomp);
% 
% 
% x0 = actvEdgIntrCrdnts(:) + triBasis*y0;
% IsclsObjctvVal = isclsObjctv(x0);
% exit_flag = 1;
% 
% dbg_h = plot3(x0(1), x0(2), x0(3), 'mo');
% delete(dbg_h);
% 
% end

triBasis = [...
    actvEdgVrtxCrdnts1(:) - actvEdgIntrCrdnts(:), ...
    actvEdgVrtxCrdnts2(:) - actvEdgIntrCrdnts(:)];


[QtriBasis RtriBasis] = qr(triBasis, 0);
actv_edg_vctr = actvEdgVrtxCrdnts2(:) - actvEdgVrtxCrdnts1(:);

Q_actv_edg_vctr = Q1MtrxVctrPrdct(actv_edg_vctr);

Q_actv_edg_vctr = ...
    (1/norm(Q_actv_edg_vctr))*Q_actv_edg_vctr;

%x0 is Q orhtogonal to the active edge
intrVrtxSprhrCntrDsplcmnt = CnstrntSphrCntrCrdnts(:) - actvEdgIntrCrdnts(:);
x0 = intrVrtxSprhrCntrDsplcmnt...
    - (Q_actv_edg_vctr'*intrVrtxSprhrCntrDsplcmnt)*Q_actv_edg_vctr;

x0 = QtriBasis*(QtriBasis'*x0);

%scale tmp_x0 so it has length equal to the radius of the
%constraint sphere
x0 = (CnstrntSphrRds/norm(x0))*x0;

if x0'*intrVrtxSprhrCntrDsplcmnt < 0
%x0 should point away from interior of the active triangle when
%anchored at the active edge midpoint
    x0 = -x0;
end

x0 = CnstrntSphrCntrCrdnts + x0;

dbg_h = plot3(x0(1), x0(2), x0(3), 'mo');
delete(dbg_h);

fminconSphrCnstrnt = @(x) x0NonlconCnstrnt(x, sphrCnstrnt);
problem_decomp = struct(...
    'x0', x0, ...
    'Aineq', ineqCnstrntLHS, ...
    'bineq', ineqCnstrntRHS, ...
    'lb', [], ...
    'ub', [], ...
    'nonlcon', fminconSphrCnstrnt , ...
    'objective', isclsObjctv, ...
    'solver', 'fmincon', ...
    'options', x0_optns);

% problem_decomp.objective = dstnc_objctvs_decomp{1};
% problem_decomp.nonlcon   = sphr_isos_cnstrnt;
% problem_decomp.x0        = x0;

%get the vertex coordinates that solve the constrained minimization problem
[x0, IsclsObjctvVal, exit_flag, output] = fmincon(problem_decomp);

dbg_h = plot3(x0(1), x0(2), x0(3), 'mo');
delete(dbg_h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return sphere constraint value and gradient in form that fmincon requires
function [...
    x0NonlconIneqCnstrnt, x0NonlconEqCnstrnt, ...
    x0NonlconIneqGrad, x0NonlconEqGrad ...
    ] ...
    = ...
    x0NonlconCnstrnt(x, sphrCnstrnt)

x0NonlconIneqCnstrnt                  = [];
x0NonlconIneqGrad                     = [];
[x0NonlconEqCnstrnt, x0NonlconEqGrad] = sphrCnstrnt(x);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
function xAffine = affineCombo(basisCrdnts, basis, offset)

xAffine = offset(:) + basis*basisCrdnts(:);

end