function [isoscelesObjctvVal, isoscelesObjctvGrad] ...
    = isoscelesObjctv(x, Q1Mtrc, Q2Mtrc)

[isoscelesConstraintVal, isoscelesConstraintGrad] ...
    = isoscelesConstraint(x, Q1Mtrc, Q2Mtrc);

isoscelesObjctvVal  = isoscelesConstraintVal^2;
isoscelesObjctvGrad = 2*isoscelesConstraintVal*isoscelesConstraintGrad;