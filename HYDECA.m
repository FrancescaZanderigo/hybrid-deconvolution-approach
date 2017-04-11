
%HYDECA.m runs the Hybrid Deconvolution Approach (Zanderigo F, Ogden T, Mann J. Hybrid deconvolution approach for model-free estimation of non-specific binding in PET studies without requiring a reference region. XXVIIth International Symposium on Cerebral Blood Flow, Metabolism and Function & XIIth International Conference. June 2015, Vancouver, Canada)

%HYDECA.m runs the Hybrid Deconvolution Approach for one PET subject - to run the algorithm on multiple subjects, considering calling HYDECA.m within a for-end-loop

%HYDECA.m requires definition of two tuning parameters, gamma and beta, that should be optimized for the PET radiotracer at hand - if they are not supplied by the user, this version of HYDECA runs with parameters optimized for the PET radiotracer [11C]DASB

%HYDECA.m accepts one object as input (Obj_HYDECA_in); Obj_HYDECA_in needs to contain the following fields:

    %thresholds_grid_SVD (optional): a row vector of scalars between 0 (excluded) and 1 for the data-driven automatic selection of the Singular Value Decomposition (SVD) threshold (for details, please see Zanderigo F, Parsey RV, Ogden RT. Model-free quantification of dynamic PET data using nonparametric deconvolution. J Cereb Blood Flow Metab. 2015 Aug;35(8):1368-79); a default grid is used if not supplied by the user
    
    %grid_Vnd_values (optional): a row vector of scalars between 0 (excluded) and infinite that are in the range of potential values of VND for the PET radiotracer at hand; ; a default grid is used if not supplied by the user

    %timeAIF: a column vector with the sampling time points of collected blood samples for the arterial input function
    
    %aif: a column vector with the plasma radioactivity values, corrected for metabolites, in correspondance of the sampling time points

    %tissue_original: a N x R matrix of values, with N the number of time frames on which the brain Time Activity Curves have been rebinned, and R the total number of considered brain regions - Time Activity Curves need to be already corrected for vascular contribution
    
    %timeTISSUE: a column vector with the mid-time points of the frames on which the brain Time Activity Curves have been rebinned

    %timeGap (optional): a scalar indicating the minimum time interval step (expressed in minutes) used to interpolate the blood and brain tissue data to an uniformly sampled grid (necessary for SVD); a default value of 0.5 minutes will be used if not specified by the user

    %GAMMAminutes (optional): a scalar indicating how many minutes after tracer injection need to be considered in the first component of the HYDECA cost function; a default value optimized for [11C]DASB will be used if not specified by the user

    %BETA (optional): a scalar indicating how to weight the second component of the HYDECA cost function relative to the first component; a default value optimized for [11C]DASB will be used if not specified by the user

%%HYDECA.m returns:

    %HYDECA_Vnd: value of Vnd estimated for the subject by HYDECA

    %HYDECA_cost_function: a vector with the values of the HYDECA cost function in correspondance of the Vnd grid values 

    %Vnd_grid: a vector with the used Vnd values grid - this is returned only if not supplied by the user
    
    %Residue_functions: residue functions R(t) estimated via SVD in each region of interest

    %Kvalues: Ki value estimated via SVD in each region of interest
    
    %convolvedcurves: convolved Time Activity Curves
    
    %time_uniform: a grid of time point, uniformly spaced, in correspondance of which Residue_functions and convolvedcurves are reported

    
    
function [Obj_HYDECA_OUT] = HYDECA (Obj_HYDECA_in)

if isfield(Obj_HYDECA_in,'thresholds_grid_SVD')
    thresh_grid = Obj_HYDECA_in.thresholds_grid_SVD;
else
    thresh_grid = [0.02: 0.01 :0.4];
end
thresh_g_long = length(thresh_grid);

if isfield(Obj_HYDECA_in,'grid_Vnd_values')
    Vnd = Obj_HYDECA_in.grid_Vnd_values;
else
    Vnd =[0.1:0.1:30]; 
end
length_Vnd_grid = length(Vnd);

if isfield(Obj_HYDECA_in,'timeGap')
    time_uniform = [Obj_HYDECA_in.timeAIF(1):Obj_HYDECA_in.timeGap:Obj_HYDECA_in.timeAIF(end)]';
else
    time_uniform = [Obj_HYDECA_in.timeAIF(1):0.5:Obj_HYDECA_in.timeAIF(end)]';
end
length_uniform = length(time_uniform);

tissue_original = Obj_HYDECA_in.tissue_original;
CUtimePs = size(tissue_original,1);
NUMroi = size(tissue_original,2);

aif = interp1(Obj_HYDECA_in.timeAIF,Obj_HYDECA_in.aif,time_uniform,'pchip');

candidatesPointsCorr = [];
countCandidates = 0;
for findTIME=1:CUtimePs
    candidatesTIME = find(Obj_HYDECA_in.timeTISSUE(findTIME)== time_uniform);
    if isempty(candidatesTIME)
    else
        countCandidates = countCandidates + 1;
        candidatesPoints(countCandidates,1) = candidatesTIME;
        candidatesPointsCorr(countCandidates,1) = findTIME;
    end
end

if isempty(candidatesPointsCorr)
    display('An error occurred in matching the original time points to the uniformly sampled new grid of time samples - please check the time points you entered for the brain tissue data and consider reducing the value of the entered timeGap variable so that it is smaller that the shortest time gap between original time samples')
else
end

for l1out = 1:length(candidatesPointsCorr)
    if candidatesPointsCorr(l1out,1) == 1
        tissue_original_l1out = tissue_original(2:end, :);
        timeTISSUE_l1out = Obj_HYDECA_in.timeTISSUE(2:end, 1);
        display(['running SVD with excluded time-point: ',num2str(l1out)])
        eval([' WBSobjtemp_' num2str(l1out) '= SVD_WBS_autothresh(Obj_HYDECA_in.timeAIF, aif, timeTISSUE_l1out, tissue_original_l1out, time_uniform, thresh_grid);']);
    else
        if candidatesPointsCorr(l1out,1) == CUtimePs
            tissue_original_l1out = tissue_original(1:end-1, :);
            timeTISSUE_l1out = Obj_HYDECA_in.timeTISSUE(1:end-1, 1);
            display(['running SVD with excluded time-point: ',num2str(l1out)])
            eval([' WBSobjtemp_' num2str(l1out) '= SVD_WBS_autothresh(Obj_HYDECA_in.timeAIF, aif, timeTISSUE_l1out, tissue_original_l1out, time_uniform, thresh_grid);']); 
        else
            tissue_original_l1out = tissue_original(1:candidatesPointsCorr(l1out,1)-1, :);
            tissue_original_l1out(candidatesPointsCorr(l1out,1):CUtimePs-1, :) = tissue_original(candidatesPointsCorr(l1out,1)+1:end, :);
            timeTISSUE_l1out = Obj_HYDECA_in.timeTISSUE(1:candidatesPointsCorr(l1out,1)-1, 1);
            timeTISSUE_l1out(candidatesPointsCorr(l1out,1):CUtimePs-1, 1) = Obj_HYDECA_in.timeTISSUE(candidatesPointsCorr(l1out,1)+1:end, 1);
            display(['running SVD with excluded time-point: ',num2str(l1out)])
            eval([' WBSobjtemp_' num2str(l1out) '= SVD_WBS_autothresh(Obj_HYDECA_in.timeAIF, aif, timeTISSUE_l1out, tissue_original_l1out, time_uniform, thresh_grid);']); 
        end
    end
end

for l1out = 1:length(candidatesPointsCorr)
    for numROI=1:NUMroi
        for euthresh = 1:thresh_g_long
            eval([ 'residuals_2_l1out(numROI, euthresh, l1out) = WBSobjtemp_' num2str(l1out) '.residuals_2(candidatesPoints(l1out,1),numROI,euthresh);' ]);
        end
    end
end
residuals_2_l1out_abs = abs(residuals_2_l1out);
[residuals_2_l1out_min, vvv] = min(residuals_2_l1out_abs,[],2);
vvv_acrossTIME = mean(vvv,3);
vvv_acrossTIME_round = round(vvv_acrossTIME);

for numROI=1:NUMroi
    eval([' WBSobjFinal_ROI' num2str(numROI) '= SVD_WBS_autothresh_wThresh(Obj_HYDECA_in.timeAIF, aif, Obj_HYDECA_in.timeTISSUE, tissue_original(:,numROI),thresh_grid(vvv_acrossTIME_round(numROI)), time_uniform);']);
end

for bqis=1:NUMroi
    eval([' All_riconv(:,bqis) = WBSobjFinal_ROI' num2str(bqis) '.riconv; '])
    eval([' AllR(:,bqis) = WBSobjFinal_ROI' num2str(bqis) '.R; '])
    eval([' All_K1(1,bqis) = WBSobjFinal_ROI' num2str(bqis) '.CBF_rel; '])
end

%shift residue curves R(t) so that max = 1 coincides with t=0 for following computation
R_shift = zeros(length_uniform, NUMroi);
where_R_shift = zeros(NUMroi,1);
for bqis=1:NUMroi
    [maxR, wheremaxR] = max(AllR(:,bqis));
    R_shift(1:length_uniform-wheremaxR+1,bqis) = AllR(wheremaxR:end,bqis);
    where_R_shift(bqis,1) = length_uniform-wheremaxR+1;
    clear maxR wheremaxR
end

%Residual Sum of Squares of the total residue curve R(t) versus interp1(timeTISSUE,Rnd(:,numROI,rstu), time_uniform) in the first GAMMA minutes
for rstu = 1:length_Vnd_grid
    for numROI=1:NUMroi
        alphaFZ(numROI,rstu) = All_K1(1,numROI)/Vnd(rstu);
        Rnd_time_uniform(:,numROI,rstu) = exp(-time_uniform.*alphaFZ(numROI,rstu));
    end
end
for rstu = 1:length_Vnd_grid
    for numROI=1:NUMroi
        for timefz = 1:length_uniform
            residuals_Rnd(numROI, rstu, timefz) = Rnd_time_uniform(timefz,numROI,rstu) - R_shift(timefz,numROI);
        end
    end
end
square_residuals_Rnd = residuals_Rnd.^2;
if isfield(Obj_HYDECA_in,'GAMMAminutes')
    if isfield(Obj_HYDECA_in,'timeGap')
        GAMMAindex = 1+Obj_HYDECA_in.GAMMAminutes/Obj_HYDECA_in.timeGap;
    else
        GAMMAindex = 1+Obj_HYDECA_in.GAMMAminutes/0.5;
    end
else
    GAMMAminutes = 10; %optimized for [11C]DASB
    %GAMMAminutes = 11; %optimized for [11C]CUMI-101
    GAMMAindex = 1+GAMMAminutes/0.5;
end
square_residuals_Rnd_OI = square_residuals_Rnd(:,:,1:GAMMAindex); 
RSS_Rnd = sum(square_residuals_Rnd_OI,3);
RSS_Rnd_sum_GAMMAmin = sum(RSS_Rnd);
    
%Area Under the Curve under the x-axis (negativity) of the corresponding specific component S(t) of R(t)
Rs_corresponding_flipped = zeros(NUMroi, length_Vnd_grid, length_uniform);
for rstu = 1:length_Vnd_grid
    for numROI=1:NUMroi
        for timefz = 1:length_uniform
            Rs_corresponding(numROI, rstu, timefz) = - Rnd_time_uniform(timefz,numROI,rstu) + R_shift(timefz,numROI);
            if Rs_corresponding(numROI, rstu, timefz)<0
                Rs_corresponding_flipped(numROI, rstu, timefz) = -Rs_corresponding(numROI, rstu, timefz);
            else
            end
        end
    end
end
for rstu = 1:length_Vnd_grid
    for numROI=1:NUMroi
        AUC_negative(numROI, rstu) = trapz(time_uniform(1:where_R_shift(numROI,1)),Rs_corresponding_flipped(numROI, rstu,1 :where_R_shift(numROI,1)));
    end
end
AUC_total_negative = sum(AUC_negative);
if isfield(Obj_HYDECA_in,'BETA')
    cost_RSSinGAMMAmin_betaAUC = RSS_Rnd_sum_GAMMAmin + Obj_HYDECA_in.BETA*AUC_total_negative;
else
    cost_RSSinGAMMAmin_betaAUC = RSS_Rnd_sum_GAMMAmin + 3.5*AUC_total_negative;%optimized for [11C]DASB
    %cost_RSSinGAMMAmin_betaAUC = RSS_Rnd_sum_GAMMAmin + 5*AUC_total_negative;%optimized for [11C]CUMI-101
end
[min_cost_RSSinGAMMAmin_betaAUC, where_cost_RSSinGAMMAmin_betaAUC] = min(cost_RSSinGAMMAmin_betaAUC);
Vnd_cost_RSSinGAMMAmin_betaAUC = Vnd(where_cost_RSSinGAMMAmin_betaAUC);

Obj_HYDECA_OUT.HYDECA_Vnd = Vnd_cost_RSSinGAMMAmin_betaAUC;
Obj_HYDECA_OUT.HYDECA_cost_function =cost_RSSinGAMMAmin_betaAUC;
if isfield(Obj_HYDECA_in,'grid_Vnd_values')
else
    Obj_HYDECA_OUT.Vnd_grid =Vnd;
end
Obj_HYDECA_OUT.Residue_functions = AllR;
Obj_HYDECA_OUT.Kvalues = All_K1;
Obj_HYDECA_OUT.convolvedcurves  = All_riconv;
Obj_HYDECA_OUT.time_uniform  = time_uniform; 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                    FUNCTIONS CALLED BY HYDECA                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% SVD_WBS_autothresh.m

function WBSobj = SVD_WBS_autothresh(timeAIF, aifIN, timeTISSUEIN, tissue_originalIN, time_uniformIN, thresh_gridIN)

ns = length(time_uniformIN);
nnnumROI = size(tissue_originalIN,2);
tissue_interp = zeros(ns,nnnumROI);
for numROI=1:nnnumROI
    tissue_interp(:,numROI) = interp1(timeTISSUEIN,tissue_originalIN(:,numROI),time_uniformIN,'pchip');
end
aif_OK = aifIN;

colonna(1)=aif_OK(1);
colonna(ns)=aif_OK(ns);
for i=2:ns-1
    colonna(i)=(aif_OK(i-1)+4*aif_OK(i)+aif_OK(i+1))/6; %for details, please see Ostergaard L, Sorensen AG, Kwong KK, Weisskoff RM, Gyldensted C, Rosen BR. High resolution measurement of cerebral blood flow using intravascular tracer bolus passages. Part II: Experimental comparison and preliminary results. Magn Reson Med. 1996 Nov;36(5):726-36.
end
riga(1)=aif_OK(1);
for i=2:ns
   riga(i)=0;
end

C_in=toeplitz(colonna,riga);
[U,S,V]=svd(C_in); 
U_trasp=U';
xam=max(max(S));
thresh_g_long = length(thresh_gridIN);

for euthresh = 1:thresh_g_long
    for i=1:ns
        for j=1:ns
            if i==j
                if S(i,j)>=xam*thresh_gridIN(euthresh); 
                    singular_values(i,euthresh)=S(i,j);
                else
                    singular_values(i,euthresh)=inf;
                end
            end
        end
    end
end
for euthresh = 1:thresh_g_long
    for k=1:ns
        if singular_values(k,euthresh)==inf
            inv_sing_val(k,euthresh)=0;
        else
            inv_sing_val(k,euthresh)=1/singular_values(k,euthresh);
        end
    end
end
for euthresh = 1:thresh_g_long
    core(:,:,euthresh)=diag(inv_sing_val(:,euthresh));
end
for euthresh = 1:thresh_g_long
    C_in_inversa(:,:,euthresh)=V*core(:,:,euthresh)*U_trasp;
end

CBFdeltatR = zeros(length(aif_OK),nnnumROI,thresh_g_long);
for euthresh = 1:thresh_g_long
    for rOLD=1:nnnumROI
        CBFdeltatR(:,rOLD,euthresh)=C_in_inversa(:,:,euthresh)*tissue_interp(:,rOLD);
    end
end
for euthresh = 1:thresh_g_long
    for rOLD=1:nnnumROI
        CBFcontempo(rOLD,1,euthresh)=max(CBFdeltatR(:,rOLD, euthresh));
        CBF_rel(rOLD,1, euthresh)=CBFcontempo(rOLD,1,euthresh)/(time_uniformIN(2)-time_uniformIN(1));
        riconv(:,rOLD,euthresh)=C_in*CBFdeltatR(:,rOLD,euthresh);
    end
end
WBSobj.riconv = riconv;

residuals_2 = zeros(size(tissue_interp,1),nnnumROI,thresh_g_long);
for euthresh = 1:thresh_g_long
    for numROI=1:nnnumROI
        residuals_2(:,numROI,euthresh) = tissue_interp(:,numROI) - riconv(:,numROI,euthresh);
    end
end
WBSobj.residuals_2 = residuals_2;

for rOLD=1:nnnumROI
    for euthresh = 1:thresh_g_long
        R_pointANY(:,rOLD,euthresh) = (CBFdeltatR(:,rOLD,euthresh)./(time_uniformIN(2)-time_uniformIN(1)))./CBF_rel(rOLD,1,euthresh);
    end
end
WBSobj.R = R_pointANY;
WBSobj.CBF_rel = CBF_rel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%SVD_WBS_autothresh_wThresh

function WBSobj = SVD_WBS_autothresh_wThresh(timeAIF, aifIN, timeTISSUEIN, tissue_originalIN, optThresh, time_uniformIN)

ns = length(time_uniformIN);
nnnumROI = size(tissue_originalIN,2);
tissue_interp = zeros(ns,nnnumROI);

for numROI=1:nnnumROI
    tissue_interp(:,numROI) = interp1(timeTISSUEIN,tissue_originalIN(:,numROI),time_uniformIN,'pchip');
end
aif_OK = aifIN;

WBSobj.tissue_interp = tissue_interp;
WBSobj.aif_OK = aif_OK;

colonna(1)=aif_OK(1);
colonna(ns)=aif_OK(ns);
for i=2:ns-1
    colonna(i)=(aif_OK(i-1)+4*aif_OK(i)+aif_OK(i+1))/6; %for details, please see Ostergaard L, Sorensen AG, Kwong KK, Weisskoff RM, Gyldensted C, Rosen BR. High resolution measurement of cerebral blood flow using intravascular tracer bolus passages. Part II: Experimental comparison and preliminary results. Magn Reson Med. 1996 Nov;36(5):726-36.
end
riga(1)=aif_OK(1);
for i=2:ns
   riga(i)=0;
end

C_in=toeplitz(colonna,riga);
[U,S,V]=svd(C_in); 
U_trasp=U';
xam=max(max(S));

for i=1:ns
    for j=1:ns
        if i==j
            if S(i,j)>=xam*optThresh;
                singular_values(i)=S(i,j);
            else
                singular_values(i)=inf;
            end
        end
    end
end
for k=1:ns
    if singular_values(k)==inf
        inv_sing_val(k)=0;
    else
        inv_sing_val(k)=1/singular_values(k);
    end
end
core=diag(inv_sing_val);
C_in_inversa=V*core*U_trasp;

CBFdeltatR = zeros(length(aif_OK),nnnumROI);
for rOLD=1:nnnumROI
    CBFdeltatR(:,rOLD)=C_in_inversa*tissue_interp(:,rOLD);
end

for rOLD=1:nnnumROI
    CBFcontempo(rOLD,1)=max(CBFdeltatR(:,rOLD));
    CBF_rel(rOLD,1)=CBFcontempo(rOLD,1)/(time_uniformIN(2)-time_uniformIN(1));
    riconv(:,rOLD)=C_in*CBFdeltatR(:,rOLD);
end
WBSobj.riconv = riconv;

residuals_2 = zeros(size(tissue_interp,1),nnnumROI);
for numROI=1:nnnumROI
    residuals_2(:,numROI) = tissue_interp(:,numROI) - riconv(:,numROI);
end
WBSobj.residuals_2 = residuals_2;

for rOLD=1:nnnumROI
    R_pointANY(:,rOLD) = (CBFdeltatR(:,rOLD)./(time_uniformIN(2)-time_uniformIN(1)))./CBF_rel(rOLD,1);
end
WBSobj.R = R_pointANY;
WBSobj.CBF_rel = CBF_rel;
