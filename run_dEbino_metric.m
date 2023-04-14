% Copyright (c) Meta Platforms, Inc. and affiliates.
% All rights reserved.
%
% This source code is licensed under the license found in the
% LICENSE file in the root directory of this source tree.

Lab_1 = [[40;0;0],[40;0;0],[40;10;-50],[40;10;-50]]; 
Lab_2 = [[50;10;0],[50;10;0],[50;20;-30],[50;20;-30]]; 
L_bg = [100, 45, 100, 45]; 

% calculate dEbino for each case (a color pair and background L*)
dEbino = deltaEbino(Lab_1, Lab_2, L_bg)

% calculate dEbino for each case (ignore background L* and B factor)
dEbino_noLbg = deltaEbino(Lab_1, Lab_2)

