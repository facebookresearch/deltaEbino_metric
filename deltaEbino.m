% Copyright (c) Meta Platforms, Inc. and affiliates.
% All rights reserved.
%
% This source code is licensed under the license found in the
% LICENSE file in the root directory of this source tree.

function   [ dEbino ] = deltaEbino(Lab1, Lab2, Lbg, params)
% this function calculates color difference metric for binocular rivalry, delta E bino.

if exist('Lbg','var')==0
    Lbg = []; 
end

if exist('params','var')==0 || isempty(params)
    % maximize AUC, using derivation data
    params.kL = 0.92;
    params.kC = 0.70;
    params.kH = 1;
    params.B_max = 28.6; % max B factor
    params.B_transition = 6.1; % transition amount for B factor to take full effect
    params.B_transition_offset = 5.3; % transition amount for B factor to take full effect
    params.dLmin = 14.5; % minimum deltaL* (L2-L1) to take full effect of B max; this is to make B factor continous from when L1=L2
end

%Parametric factors
kL = params.kL;
kC = params.kC;
kH = params.kH;


%CIELAB Chroma
C1 = sqrt(Lab1(2,:).^2+Lab1(3,:).^2);
C2 = sqrt(Lab2(2,:).^2+Lab2(3,:).^2);

%Lab Prime
mC = (C1+C2)./2;
G=0.5*(1-sqrt((mC.^7)./((mC.^7)+(25.^7))));
LabP1 = [Lab1(1,:) ; Lab1(2,:).*(1+G) ; Lab1(3,:)];
LabP2 = [Lab2(1,:) ; Lab2(2,:).*(1+G) ; Lab2(3,:)];

%Chroma
CP1 = sqrt(LabP1(2,:).^2+LabP1(3,:).^2);
CP2 = sqrt(LabP2(2,:).^2+LabP2(3,:).^2);

%Hue Angle
hP1t = atan2Deg(LabP1(3,:),LabP1(2,:));
hP2t = atan2Deg(LabP2(3,:),LabP2(2,:));

%Add in 360 to the smaller hue angle if absolute value of difference is > 180
hP1 = hP1t + ((hP1t<hP2t)&(abs(hP1t-hP2t)>180)).*360;
hP2 = hP2t + ((hP1t>hP2t)&(abs(hP1t-hP2t)>180)).*360;

%Delta Values
DLP = LabP1(1,:) - LabP2(1,:);
DCP = CP1 - CP2;
DhP = hP1 - hP2;
DHP = 2*(CP1.*CP2).^(1/2).*sinDeg(DhP./2);

%Arithmetic mean of LCh' values
mLP = (LabP1(1,:)+LabP2(1,:))./2;
mCP = (CP1+CP2)./2;
mhP = (hP1+hP2)./2;

%Weighting Functions
SL = 1+(0.015.*(mLP-50).^2)./sqrt(20+(mLP-50).^2);
% SL = 1;
SC = 1+0.045.*mCP;
T = 1-0.17.*cosDeg(mhP-30)+0.24.*cosDeg(2.*mhP)+0.32.*cosDeg(3.*mhP+6)-0.2.*cosDeg(4.*mhP-63);
SH = 1+0.015.*mCP.*T;

%Rotation function
RC = 2.*sqrt((mCP.^7)./((mCP.^7)+25.^7));
DTheta = 30.*exp(-((mhP-275)./25).^2);
RT = -sinDeg(2.*DTheta).*RC;

% calculate B factor
if isempty(Lbg) % ignore B factor if Lbg is not supplied
    B = 0;
else
    [ B ] = calcBfactor(Lab1(1,:), Lab2(1,:), Lbg, ...
        params.B_max, params.B_transition, params.B_transition_offset, params.dLmin);
end

dEbino = ((DLP./kL./SL).^2+(DCP./kC./SC).^2+(DHP./kH./SH).^2+(RT.*(DCP./kC./SC).*(DHP./kH./SH))).^(1/2)+B;

end


function out = atan2Deg(inY,inX)
out = atan2(inY,inX).*180./pi;
out = out+(out<0).*360;
end


function out = sinDeg(in)
out = sin(in.*pi./180);
end


function out = cosDeg(in)
out = cos(in.*pi./180);
end


function [ B_all ] = calcBfactor(L1_all, L2_all, Lbg_all, ...
    B_max, B_transition, B_transition_offset, dLmin)
% calculate B factor (background effect)
% if Lbg is between L1 and L2, L1 and L2 would most likely appear rivarlous

% B_max = 20; % max B factor
% B_transition = 5; % transition amount for B factor to take full effect
% B_transition_offset = 5; % transition offset for B factor to start taking effect
% dLmin = 10; % minimum deltaL* (L2-L1) to take full effect of B max; this is to make B factor continous from when L1=L2

L_range = [0,100];
n_sample = length(Lbg_all);
B_all = nan(size(Lbg_all));

for k = 1:n_sample
    L1 = L1_all(k);
    L2 = L2_all(k);
    Lbg = Lbg_all(k);

    %     if (L1-Lbg)*(L2-Lbg)<0
    %         min_diff = min([abs(L1-Lbg),abs(L2-Lbg)]);
    %
    %         if min_diff>=B_transition
    %             B = B_max;
    %         else
    %             B = min_diff/B_transition*B_max;
    %         end
    %     else
    %         B = 0;
    %     end

    % check which is larger; L1 or L2
    B = Inf;
    if L2<L1
        Lmin = L2; Lmax = L1;
    elseif L1<L2
        Lmin = L1; Lmax = L2;
    else
        B = 0;
    end

    if B>0
        % adjust effective B max value based on dL
        dL = Lmax-Lmin;
        if dLmin<=dL
            B_effective = B_max;
        else
            B_effective = B_max*dL/dLmin;
        end

        % offset Lmin and Lmax
        Lmin = max(Lmin-B_transition_offset, L_range(1));
        Lmax = min(Lmax+B_transition_offset, L_range(2));

        if (Lmin-Lbg)*(Lmax-Lbg)<0
            min_diff = min([abs(Lmin-Lbg),abs(Lmax-Lbg)]);

            if min_diff>=B_transition
                B = B_effective;
            else
                B = min_diff/B_transition*B_effective;
            end
        else
            B = 0;
        end
    end
    B_all(k) = B;
end

end

