function [isoscelesConstraintVal, isoscelesConstraintGrad] ...
    = isoscelesConstraint(x, Q1Mtrc, Q2Mtrc)

[Q1Val, Q1Grad] = Q1Mtrc(x);
[Q2Val, Q2Grad] = Q2Mtrc(x);

Q1ValQ2ValSum = Q1Val + Q2Val;

isoscelesConstraintVal  = (Q1Val - Q2Val)/Q1ValQ2ValSum;
isoscelesConstraintGrad = ...
    (2/Q1ValQ2ValSum^2)*(Q2Val*Q1Grad - Q1Val*Q2Grad);
