% ----- INPUT: 
%
% cellbonddata is the filename of Matlab workspace generated 
% from Tissue Analyser (imported bonds.csv and cells.csv files).

load('cellbonddata');
evolver=0; % Set 1 for pseudodata analysis of surface evolver

%% OUTPUT: WHAT DOES THIS PROGRAM COMPUTE/GENERATE
% PRIMARY::
% 1. COMPUTES TENSION AND PRESSURE BY BAYESIAN INFERENCE
%    THERE ARE TWO CHOICES - USUAL BOX COAR-GRAINING OR CELL LEVEL
% 2. ORGANELLE POLARITY CALCULATION
%    POSITION OF AN ORGANELLE W.R.T. NUCLEUS. VECTOR: NUCLEYS -> ORGANELLE
%    RELATIVE ANGLE BETWEEN POLARITY VECTOR AND OVERALL MOVEMENT VECTOR
% 3. COMPUTE FORCE POLARIZATION AT CELL LEVEL
%    NEEDS TO BE FIXED; REQUIRES CELL LEVEL BATCHERLOR STRESS
% 4. COMPURE GEOMETRIC POLARIZATION
%    POLYGONAL SHAPE ANISOTROPY => MAJOR AXIS, MINOR AXIS, ORIENTATION
% 5. 

%% SYSTEM PARAMETERS
SquareSize = 80; % Box size for coarse-grained batchelor stress
nRadius = 5.5;   % Circle radius for cellular level batchelor stress
pixScale = 1;    % pix to actual distance, if requires. Default set to 1

%% CHOICE PARAMETERS
isDelFreeEdge = 0; % Do I delete the free edge i.e. set zero tension?
showLLL= 0;        % Display Log likelihood?
showDist = 1;      % Show distance of each cell from the leading edge?
showLeaderDist = 0;  % Show distance of each cell from the leader cell?
leaderPos = [982 510]; %[1006 482]; % if isLeaderDist = 1, set the leader cell position
showCellForce = 1; % Comppute cell level batchleor stress coarse-grained over a circle with nRadius
showGraphOutput = 1; % Show generated Graph Network Output?
isRandClass = 0 % Randomize class labels
excludePolAng = 1
%% INPUT NUCLEUS AND ORGANELLE IMAGES

% Nuclear channel => Read the file with DAPI staining
[nFlname, nPathname] = uigetfile('*.tif','Select the file with DAPI staining');
nFlpathname = strcat(nPathname,nFlname);
nI = imread(nFlpathname);
if size(nI,3)==3
    nI = rgb2gray(nI);
end
nI = double(nI);
[n_xLen n_yLen] = size(nI);

% Organelle channel => Read the file with Organelle staining
[orgFlname, orgPathname] = uigetfile('*.tif','Select the file with Organelle staining');
orgFlpathname = strcat(orgPathname,orgFlname);
orgI = imread(orgFlpathname);
if size(orgI,3)==3
    orgI = rgb2gray(orgI);
end
orgI = double(orgI);

%% TRANSLATE TISSUE ANALYZER GEOMETRY TO MATLAB
tic;
disp('Generating geometry from Tisse Analyzer data...')
TissueAnalyzerToMatlab;

%% REMOVE FREE EDGES FROM T-P CALCULATION
% Find edges of each cell (CellEdges)
Ne=length(E); Nc=length(C);
for i=1:Nc % for each cell
    clear e
    nv=length(C{i}); % how many vertices (and edges) in the cell
    for j=1:nv-1
        e(j,:)=[C{i}(j) C{i}(j+1)]; % each couple of adjacent vertices is an edge
        for k=1:Ne % find which edge
            if (E(k,1)== e(j,1) && E(k,2)== e(j,2)) || (E(k,1)== e(j,2) && E(k,2)== e(j,1))
                CellEdges{i}(j)=k;
            end
        end
    end
    e(nv,:)=[C{i}(nv) C{i}(1)];
    for k=1:Ne
        if (E(k,1)== e(nv,1) && E(k,2)== e(nv,2)) || (E(k,1)== e(nv,2) && E(k,2)== e(nv,1))
            CellEdges{i}(nv)=k; % List of edges (index in E) in Cell i
        end
    end
    
    CellEdges{i}(CellEdges{i}==0)=[]; 
    % some zeros can appear because of image boundaries, remove them !
    % (edges connecting two vertices on the image boundary 
    % are not in E, because tension is not calculated for these edges)
end

% Remove the edges associated with only one cell => Cells at the edge
cellpEdge = zeros(Ne,1);

for c=1:Nc
    numE = length(CellEdges{c});
    if numE<12                  % cell free contributes to cell initially
        for k=1:numE
            cE = CellEdges{c}(k);
            cellpEdge(cE) = cellpEdge(cE) + 1;
        end
    end
end

edgeToDel = find(cellpEdge == 1);
verEdge = EDGES(edgeToDel,:); verEdgeN = length(verEdge)*2;
verEdgeList = reshape(verEdge,[verEdgeN,1]);
verEdge = unique(verEdgeList); verEdgeN = length(verEdge);
vEx = real(VERTICES(verEdge)); vEy = imag(VERTICES(verEdge));
vE = [vEx vEy];

if isDelFreeEdge == 1
    E(edgeToDel,:) = []; E_tot = length(E);
end

%% COMPUTATION OF THE COEFFICIENTS OF MATRIX A 

[A,B,g]=subprog_Calcul_Matrix_A(V, E, C, V_ent, E_tot, C_tot,evolver);
disp(['Done (' num2str(round(toc,2)) ' seconds).'])

%% MAXIMUM LIKEHOOD ESTIMATION (L)
% This is by far the longest step for large tissues, so the range of mu
% tested should be as small as possible...

tic;
disp('Estimating maximum likelihood...')
murange=0.4:0.1:3; 
% range of mu tested. mu is typically slightly smaller than 1
% look at logL(mu), it should have a maximum !
cpt=1;
for mu=murange  
    logL(cpt)=subprog_logLikelihood_Estimation(mu,A,B,g);
    if cpt>4 && logL(cpt)<logL(cpt-1) && logL(cpt-1)<logL(cpt-2) && logL(cpt-2)<logL(cpt-3)
        break %if a clear maximum has been found, break (to save time)
    end
    cpt=cpt+1;
end
disp(['Done (' num2str(round(toc,2)) ' seconds).'])
[~,I]=max(logL);
mu_bayes=murange(I);
mu_bayes

if showLLL == 1
    figure; set(gcf,'color','w');
    plot(murange(1:min(cpt,length(murange))),logL)
    ylabel('Log Likelihood')
    xlabel('Mu')
end
%saveas(gcf,'likelihood.fig');

%% INVERSE PROBLEM

tic;
disp('Solving the inverse problem...')
[INFERENCE] = subprog_inverse_matrix_A(mu_bayes, A, B, g, E_tot);
disp(['Done (' num2str(round(toc,2)) ' seconds).'])

%% FIGURE / TENSIONS
% this should be customized depending on user's needs

figure; set(gcf,'color','w'); couleur = jet(101); hold on

tmin=1-3.5*std(INFERENCE.TENSIONS); % 3*std : good compromise for color limits
tmax=1+2.5*std(INFERENCE.TENSIONS);
Ti = (INFERENCE.TENSIONS-tmin)./(tmax-tmin);
Ti(Ti<0)=0;
Ti(Ti>1)=1;
Ti = round(Ti*100)+1; %scale to jet101
%color scaled from tmin to tmax
colormap jet(101)
cbr = colorbar;
cbr.Ticks=[0 1];
cbr.TickLabels={num2str(round(tmin,1)) num2str(round(tmax,1))};

for cpt = 1:E_tot,
    plot(V(E(cpt,:),1),V(E(cpt,:),2),...
        'Color',couleur(Ti(cpt),:),'LineWidth',2)
end
axis equal; axis tight; xlabel('Bayesian Inference Tension Map')
axis off
set(gca,'Ydir','reverse');
%saveas(gcf,'tensions.fig');

%% FIGURE / PRESSURES
% this should be customized depending on user's needs

figure; set(gcf,'color','w'); couleur = jet(101); hold on;
% Change colormap for representation (jet) and computing (gray)

c_len = min(length(INFERENCE.PRESSURES),C_ent);
P_plot=INFERENCE.PRESSURES(1:c_len);
Pmax=3*std(P_plot); % 3*std : good compromise for color limits
Pmin=-Pmax;
Pi = (P_plot-Pmin)./(Pmax-Pmin);
Pi(Pi<0)=0;
Pi(Pi>1)=1;
Pi = round(Pi*100)+1; %scale to jet101
colormap(couleur)
cbr = colorbar;
cbr.Ticks=[0 1 2];
cbr.TickLabels={num2str(round(Pmin,2)) '0' num2str(round(Pmax,2))};
cellOutline = DrawCellOutline(E_tot,E,V,1);
for c = 1:c_len
    v_in_f = C{c};
    patch(V(v_in_f, 1), V(v_in_f, 2), 1,'FaceColor',couleur(Pi(c),:))
end
axis equal; axis tight; xlabel('Bayesian Inference Pressure Map')
axis off
set(gca,'Ydir','reverse');
%saveas(gcf,'pressures.fig');

%% COMPUTE MEAN CELL DIAMETERES
cellArea = zeros(c_len,1);
for i=1:c_len
    CellArea(i)=polyarea(V(C{i},1),V(C{i},2));
end
meanCArea = mean(CellArea);
meanDia = 2*((meanCArea/pi)^0.5);

%% OPTION: COMPUTE STRESS TENSOR (BATCHELOR)

if (isDelFreeEdge == 0)&&(showCellForce == 0)
    tic;
    disp('Computing the stress tensor...')
    [N,S,P,BoxCenter]=batchelor(INFERENCE.TENSIONS,INFERENCE.PRESSURES,C,E,V, SquareSize);
    %[N,S,P,BoxCenter]=batchelor(ones(1,length(INFERENCE.TENSIONS)),zeros(1,length(INFERENCE.PRESSURES)),C,E,V, SquareSize);
    disp(['Done (' num2str(round(toc,2)) ' seconds).'])
%     saveas(gcf,'stress.fig');
    Bayesian_Inference.Batchelor.N=N;
    Bayesian_Inference.Batchelor.S=S;
    Bayesian_Inference.Batchelor.P=P;
    Bayesian_Inference.Batchelor.BoxCenter=BoxCenter;
    set(gca,'Ydir','reverse');
    
    lenN = length(N);
    sAvg = zeros(1,lenN); Xgrid = zeros(1,lenN); Ygrid = zeros(1,lenN);
    for i=1:lenN
        sAvg(i)=(N{i}(1,1) + N{i}(2,2))/2;
        Xgrid(i) = BoxCenter{i}(1); Ygrid(i) = BoxCenter{i}(2);
    end
    Xg = linspace(min(Xgrid), max(Xgrid), 100);
    Yg = linspace(min(Ygrid), max(Ygrid), 100);
    [xi, yi] = meshgrid(Xg, Yg);
    F = scatteredInterpolant(Xgrid',Ygrid',sAvg');
    zi = F(xi,yi);
    figure;
    surf(xi,yi,zi, 'EdgeAlpha', 0)
    colormap jet
elseif (isDelFreeEdge == 0)&&(showCellForce == 1)
    tic;
    disp('Computing the stress tensor at cell centroids...')
    [N,S,PI,theta1N,theta2N,lambdasN,CellCenter]=...
        cellbatchelor(INFERENCE.TENSIONS,INFERENCE.PRESSURES,C,E,V,C_ent,CellEdges,nRadius,n_xLen,n_yLen);
    disp(['Done (' num2str(round(toc,2)) ' seconds).'])
else
    disp('Not computing Batchelor stresses')
end

%% OUTPUT of Force Inference Calculations

Bayesian_Inference.Infered_Tensions=INFERENCE.TENSIONS;
Bayesian_Inference.Infered_Pressures=INFERENCE.PRESSURES;
Bayesian_Inference.Geometry.Cells=C;
Bayesian_Inference.Geometry.Vertices=V;
Bayesian_Inference.Geometry.Edges=E;
Bayesian_Inference.mu_bayes=mu_bayes;
Cells = C;

%% Organelle Polarity Calculations
figure; 
tic
disp('Calculating the organelle polarity...')
hold on
pcolor = cool(181); colormap(pcolor);  
cbr = colorbar; cbr.Ticks=[0 1 2]; 
cbr.TickLabels={'Parallel' 'Perpendicular' 'Anti-Parallel'};
amp_fac = 3;
cellOutline = DrawCellOutline(E_tot,E,V,1);
% for e = 1:E_tot,
%     plot(V(E(e,:),1),V(E(e,:),2),'k','LineWidth',1)
% end
ThetaInDegrees = zeros(c_len,1); distEdge = zeros(c_len,1);
cellCentroid_X = zeros(c_len,1); cellCentroid_Y = zeros(c_len,1);
normPolrity = zeros(c_len,1); 
realOrgPol = zeros(c_len,1);
orgVect = zeros(c_len,2);
nucPos = zeros(c_len,1);
for c = 1:c_len
    v_in_f = C{c};
    xV = V(v_in_f,1)'; yV = V(v_in_f,2)';
    BW = poly2mask(xV,yV,n_xLen,n_yLen);
    xg = 1:n_xLen; yg = 1:n_yLen; [X,Y] = meshgrid(xg,yg);
     
    nCentroid_X = mean(mean(X.*BW.*nI))/mean(mean(BW.*nI)); 
    nCentroid_Y = mean(mean(Y.*BW.*nI))/mean(mean(BW.*nI));
    
    orgCentroid_X = mean(mean(X.*BW.*orgI))/mean(mean(BW.*orgI));
    orgCentroid_Y = mean(mean(Y.*BW.*orgI))/mean(mean(BW.*orgI));
    
    cellCentroid_X(c) = mean(mean(X.*BW))/mean(mean(BW)); 
    cellCentroid_Y(c) = mean(mean(Y.*BW))/mean(mean(BW));
    orgVectX = (orgCentroid_X-nCentroid_X);
    orgVectY = (orgCentroid_Y-nCentroid_Y);
    
    realOrgPol(c) = atan2d(orgVectY,orgVectX);
    orgVect(c,:) = [orgVectX orgVectY];
    
    migDir = [1 0]; vecDir = [orgVectX orgVectY];
    
    normPolrity(c) = norm(vecDir)/meanDia;
    CosTheta = dot(migDir,vecDir)/(norm(migDir)*norm(vecDir));
    ThetaInDegrees(c) = acosd(CosTheta); 
    ThetaInDegrees_r = max(round(ThetaInDegrees(c)),1);
    % ThetaInDegrees_r = max(round(realOrgPol(c)+180),1);
 
    patch(V(v_in_f, 1), V(v_in_f, 2), 1,...
        'FaceColor',pcolor(ThetaInDegrees_r,:))
    hold on

    q = quiver(cellCentroid_X(c),cellCentroid_Y(c),...
        amp_fac*orgVectX,amp_fac*orgVectY,'color',[0 0 0]);
    q.MaxHeadSize=2;
    
end
hold off
disp(['Done (' num2str(round(toc,2)) ' seconds).'])
axis equal; axis tight; xlabel('Organelle Polarity Map')
axis off
set(gca,'Ydir','reverse');

%% CELL DISTANCE FROM EDGE
whiteIm = zeros(n_xLen,n_yLen);
distEdge = zeros(c_len,1);

% Creating the mask image
for c = 1:C_tot
    v_in_f = C{c};
    cArea=polyarea(V(C{c},1),V(C{c},2));
    
    if cArea < 3*meanCArea
        xV = V(v_in_f,1)'; yV = V(v_in_f,2)';
        BW = poly2mask(xV,yV,n_xLen,n_yLen);
        whiteIm = whiteIm + BW;
    end
end 

for c = 1:c_len
    currRow = whiteIm(round(cellCentroid_Y(c)),:);
    nonMask = find(currRow == 0);
    Pt = find(nonMask>=round(cellCentroid_X(c)),1); edgePt = nonMask(Pt);
    distEdge(c) = (edgePt-cellCentroid_X(c))/meanDia;
end

%% POLARITY VS DISTANCE PLOT
% DOES THE AVERAGE POLARITY VARY SYSTEMATICALLY WITH THE INCREASING
% DISTANCE FROM THE LEADING EDGE?
figure;
distEdge_r = round(distEdge);
[distEdge_bin,theta_mean] = consolidator(distEdge_r, ThetaInDegrees,@mean);

plot(round(distEdge), ThetaInDegrees,'s','color',[0 0.5 1]);
ylabel('Polrity in degrees'); xlabel('Distance from edge (# of cells)')
hold on
plot(distEdge_bin,theta_mean,'o-', 'color',[0.6 0 0]);
hold off

%% SHOW DISTANCE OF EACH CELL FROM FREE EDGE: FOR TESTING
if showDist == 1
    figure; 
    maxDist = max(distEdge_r);
    pcolor = jet(maxDist); colormap(pcolor);
    for c=1:c_len
        v_in_f = C{c};
        cellColor = distEdge_r(c);
        if cellColor > 0
            patch(V(v_in_f, 1), V(v_in_f, 2), 1,...
                'FaceColor',pcolor(cellColor,:))
            hold on
        else
            patch(V(v_in_f, 1), V(v_in_f, 2), 1,...
                'FaceColor',[0 0 1])
            hold on
        end
    end
    hold off
    axis equal; axis tight; xlabel('Distance From Wound Edge Map')
    axis off
    set(gca,'Ydir','reverse');
end

%% CALCULATE DISTANCE FROM LEADER CELLS
distLead_norm = zeros(c_len,1);
leadX = leaderPos(1); leadY = leaderPos(2);
for c=1:c_len
    distLead = ((cellCentroid_X(c)-leadX)^2 + ...
        (cellCentroid_Y(c)-leadY)^2)^0.5;
    distLead_norm(c)=floor(distLead/meanDia);
end
    
if showLeaderDist == 1
    figure;
    maxDistLead = max(distLead_norm);
    pcolor = jet(maxDistLead); colormap(pcolor);
    for c=1:c_len
        v_in_f = C{c};
        cellColor = distLead_norm(c);
        if cellColor > 0
            patch(V(v_in_f, 1), V(v_in_f, 2), 1,...
                'FaceColor',pcolor(cellColor,:))
            hold on
        end
    end
    hold off
    axis equal; axis tight; xlabel('Distance From Leader Cell Map')
    axis off
    set(gca,'Ydir','reverse');
    
    figure;
    [distLead_bin,theta_mean] = consolidator(distLead_norm,...
        ThetaInDegrees,@mean);
    plot(distLead_norm, ThetaInDegrees,'s','color',[0 0.5 1]);
    ylabel('Polrity in degrees'); xlabel('Distance from Leader Cell')
    hold on
    plot(distLead_bin,theta_mean,'o-', 'color',[0.6 0 0]);
    hold off
end

%% CALCULATING FORCE POLARIZATION
figure; tic
hold on
pcolor = autumn(91); colormap(pcolor);  
cbr = colorbar; cbr.Ticks=[0 1 2]; cbr.TickLabels={'0' '45' '90'};
cellOutline = DrawCellOutline(E_tot,E,V,1);

if showCellForce == 1
    pseudoforceAng = zeros(c_len,1); orgPol360 = realOrgPol; 
    dotfV12 = zeros(c_len,1); thetaF = zeros(c_len,2); 
    relAngFP = zeros(c_len,1);
    
    for c=1:c_len
        if realOrgPol(c)<0
            orgPol360(c) = 360 + realOrgPol(c);
        end
        oV = [cosd(orgPol360(c)) sind(orgPol360(c))];
 
        thetad1 = 180*theta2N{c}/pi;
        if thetad1<0
            thetad1 = 360 + thetad1;
        end
        thetad2 = rem((thetad1+180),360);
        thetaF(c,:) = [thetad1  thetad2];
        
        fV1 = [cosd(thetad1) sind(thetad1)];
        fV2 = [cosd(thetad2) sind(thetad2)];
        dotfV12(c) = dot(fV1,fV2);
        
        dotP1 = dot(oV,fV1); dotP2 = dot(oV,fV2);

        if dotP1>=dotP2
            relAngFP(c) = acosd(dotP1); 
        else
            relAngFP(c) = acosd(dotP2);
        end
        pseudoforceAng(c) = ThetaInDegrees(c) + relAngFP(c);
        
        v_in_f = C{c};
        cellColor = max(round(relAngFP(c)),1);
        patch(V(v_in_f, 1), V(v_in_f, 2), 1,...
            'FaceColor',pcolor(cellColor,:))
        hold on
    end
    axis equal; axis tight; xlabel('Relative Force Polarization Map')
    axis off
    set(gca,'Ydir','reverse');
end

%% Plot stress tensor
% principal directions and amplitudes of the stress tensor N
% 
figure; set(gcf,'color','w'); hold on
for e = 1:Ne
    plot([V(E(e,1),1) V(E(e,2),1)],[V(E(e,1),2) V(E(e,2),2)],'color',[0.7 0.7 0.7]);
end
forceAniso = zeros(c_len,1);
for b=1:c_len
    lambdasN{b}=1000*lambdasN{b}; % for plot
    xC = CellCenter{b}(1); yC = CellCenter{b}(2);
    % plot( [xC-0.5*lambdasN{b}(1)*cos(theta1N{b}) xC+0.5*lambdasN{b}(1)*cos(theta1N{b})] , [yC-0.5*lambdasN{b}(1)*sin(theta1N{b}) yC+0.5*lambdasN{b}(1)*sin(theta1N{b})],'linewidth',2,'color',[0.12 0.56 0.24])
    plot( [xC-0.5*lambdasN{b}(2)*cos(theta2N{b}) xC+0.5*lambdasN{b}(2)*cos(theta2N{b})] , [yC-0.5*lambdasN{b}(2)*sin(theta2N{b}) yC+0.5*lambdasN{b}(2)*sin(theta2N{b})],'linewidth',2,'color',[0.86 0.08 0.24])
    forceAniso(b) = (lambdasN{b}(2)-lambdasN{b}(1))/(lambdasN{b}(2)+lambdasN{b}(1));
end
q = quiver(cellCentroid_X,cellCentroid_Y,amp_fac*orgVect(:,1),amp_fac*orgVect(:,2),'linewidth',2,'color',[0 0.2 0.6]);
q.MaxHeadSize=2;
axis equal; axis tight; xlabel('Stress Tensor Map'); axis off
set(gca,'Ydir','reverse');

%% CALCULATING GEOMETRIC POLARIZATION
relAng = zeros(c_len,1);
figure; 
tic
disp('Calculating the geometric polarization ...')
hold on
scale_fac = 0.5; aniso = zeros(c_len,1); 
cellCircularity = zeros(c_len,1); shapeIndex = zeros(c_len,1);
geoVect = zeros(c_len,2); startPt = zeros(c_len,2);
numVertex = zeros(c_len,1);
for c=1:c_len
    v_in_f = C{c}; v_len = length(v_in_f);
    xV = V(v_in_f,1)'; yV = V(v_in_f,2)';
    BW = poly2mask(xV,yV,n_xLen,n_yLen);
    s = regionprops(BW,'Centroid','MajorAxisLength',...
        'MinorAxisLength','Orientation','Perimeter');
   
    majAx = s.MajorAxisLength; minAx = s.MinorAxisLength; 
    phi(c) = s.Orientation; xcent = s.Centroid(1); ycent = s.Centroid(2);
    cellCircularity(c) = 4*pi*CellArea(c)/(s.Perimeter^2);
    shapeIndex(c) = s.Perimeter/(CellArea(c))^0.5; 
    
    x1 = xcent - 0.5*scale_fac*majAx*cosd(phi(c)); 
    x2 = xcent + 0.5*scale_fac*majAx*cosd(phi(c));
    y1 = ycent - 0.5*scale_fac*majAx*sind(phi(c)); 
    y2 = ycent + 0.5*scale_fac*majAx*sind(phi(c));
    startPt(c,:)=[x1 y1];
    geoVectX = (x2-x1); geoVectY = (y2-y1); 
    geoVect(c,:) = [geoVectX geoVectY];
    CosTheta = dot(geoVect(c,:),orgVect(c,:))/...
        (norm(geoVect(c,:))*norm(orgVect(c,:)));
    relAng(c) = round(acosd(CosTheta));
    
    aniso(c)=(majAx-minAx)/(majAx+minAx);
    q = quiver(x1,y1,geoVectX,geoVectY,'color',[0 0 1]);
    
    % number of verties and coordination number
    numVertex(c) = v_len;
    
end

sF = cellCircularity; % Selecting the appropriate shape factor

sFmax = max(sF); sFmin = min(sF); sFmean = mean(sF);
gP = (sF - sFmin)/(sFmax - sFmin);
gP = round(gP*180)+1; %scale to autumn 181

pcolor = autumn(181); colormap(pcolor);
cbr = colorbar; cbr.Ticks=[0 1 2];
cbr.TickLabels={num2str(round(sFmin,2)) num2str(round(sFmean,2)) ...
    num2str(round(sFmax,2))};
cellOutline = DrawCellOutline(E_tot,E,V,1);
for c = 1:c_len
    v_in_f = C{c};
    patch(V(v_in_f, 1), V(v_in_f, 2), 1,'FaceColor',pcolor(gP(c),:))
end
hold on

q = quiver(startPt(:,1),startPt(:,2),geoVect(:,1),geoVect(:,2),'color',[0 0 1])
q.MaxHeadSize=0;

hold off
disp(['Done (' num2str(round(toc,2)) ' seconds).'])
axis equal; axis tight; xlabel('Geometric Polarity Map')
axis off
set(gca,'Ydir','reverse');

%% OUTPUT A GRAPH NETWORK WITH ADJACENCY, FEATURES, AND CLASSIFICATION
% FEATURES = ORGANELLE POLARITY, DISTANCE FROM LEADING EDGE, AREA, ...
% BATCHELOR STRESS, GEOMETRIC POLARITY
% PURPOSE :: USE THIS REPRESENTATION FOR GRAPH NEURAL NETWORK

% Generate Adjacency matrix -----------------------------------------------
numNode = c_len; gn.numNode = numNode;
AdjMat = zeros(numNode,numNode);
edgeCnt = zeros(e,2);

for c=1:numNode
    numE = length(CellEdges{c});
    for k=1:numE
        cE = CellEdges{c}(k);
        if edgeCnt(cE,1)==0
            edgeCnt(cE,1)=c;
        else
            edgeCnt(cE,2)=c;
        end
    end
end

edgeCntClone = edgeCnt;
for k=1:e
    c1 = edgeCnt(k,1); c2 = edgeCnt(k,2); p = c1*c2;
    if p>0
        AdjMat(c1,c2)=1;
    end
end
% imagesc(AdjMat); colormap cool
rowsToDelete = any((edgeCnt(:,1).*edgeCnt(:,2))==0,2);
edgeCntClone(rowsToDelete,:) = [];

gn.source = edgeCntClone(:,1); gn.target = edgeCntClone(:,2);

% [X Y] node coordinates and node labels
gn.nodeXY = [cellCentroid_X cellCentroid_Y]; gn.nodeLabel = (1:numNode)';

% Features
% Degree and Magnitude of polarity ----------------------------------------
num_features = 1; feature_matrix = zeros(numNode,num_features);
nodePolAngle = orgPol360; feature_matrix = nodePolAngle;

nodePolMag = normPolrity; num_features = num_features + 1; 
feature_matrix = [feature_matrix,nodePolMag];

% Distance from leading edge and Distance from leader ---------------------
nodeEdgeDist = distEdge; num_features = num_features + 1; 
feature_matrix = [feature_matrix, nodeEdgeDist];

nodeDistLeader = distLead_norm; num_features = num_features + 1; 
feature_matrix = [feature_matrix, nodeDistLeader];

% Area, Circularity (=4*pi*A/P^2), and Shape Index (P/A^0.5) --------------
nodeArea = CellArea'; num_features = num_features + 1; 
feature_matrix = [feature_matrix, nodeArea];

nodeCircularity = cellCircularity; num_features = num_features + 1; 
feature_matrix = [feature_matrix, nodeCircularity];

nodeShapeIndex = shapeIndex; num_features = num_features + 1; 
feature_matrix = [feature_matrix, nodeShapeIndex];

nodeVertexNum = numVertex; num_features = num_features + 1;
feature_matrix = [feature_matrix, nodeVertexNum];

% Geometeric Orientation and Anisotropy -----------------------------------
nodeGeoPolarity = phi'; num_features = num_features + 1; 
feature_matrix = [feature_matrix, nodeGeoPolarity];

nodeGeoAniso = aniso; num_features = num_features + 1; 
feature_matrix = [feature_matrix, nodeGeoAniso];

% Force Orientation ------------------------------
nodeForceAniso = forceAniso; num_features = num_features + 1;
feature_matrix = [feature_matrix, nodeForceAniso];

nodeForceAng = pseudoforceAng; num_features = num_features + 1;
feature_matrix = [feature_matrix, nodeForceAng];


feature_matrix_n = colnormalize(feature_matrix);

% ----------------------------------------------------------
if excludePolAng == 1
    feature_matrix_n(:,1) = []; 
    num_features = num_features - 1;
end

% Def Classification :: ---------------------------------------------------
% (Parallel<=45)=2; (45<Perpendicular<135)=1;(Anti-parallel>=135)=0
gn.numClasses = 3;
gn.nodeClass = zeros(numNode,1);
for k=1:numNode
    if ThetaInDegrees(k)<=45
        gn.nodeClass(k)=2;
    elseif ThetaInDegrees(k)>=135
        gn.nodeClass(k)=0;
    else
        gn.nodeClass(k)=1;
    end
end

randLabel = floor(gn.numClasses*rand(numNode,1));
if isRandClass == 1
    gn.nodeClass = randLabel;
end

% Show Graph Construction -------------------------------------------------
if showGraphOutput == 1
    figure; set(gcf,'color','w');
    hold on
    % Outline all cells
    for e = 1:E_tot,
        plot(V(E(e,:),1),V(E(e,:),2),'k','LineWidth',1.5)
    end
    % Mark the nodes
    plot(gn.nodeXY(:,1),gn.nodeXY(:,2),'b*');
    % textLabel = num2str(gn.nodeLabel);
    % text(gn.nodeXY(:,1),gn.nodeXY(:,2),textLabel,'VerticalAlignment','top');
    % Draw Graph
    for k=1:e
        c1 = edgeCnt(k,1); c2 = edgeCnt(k,2); p = c1*c2;
        Xs = zeros(2,1); Ys = zeros(2,1);
        if p>0
            Xs(1) = gn.nodeXY(c1,1);  Xs(2) = gn.nodeXY(c2,1);
            Ys(1) = gn.nodeXY(c1,2);  Ys(2) = gn.nodeXY(c2,2);
            plot(Xs,Ys,'r','LineWidth',1.5)
        end
    end
    hold off
end
axis equal; axis tight; xlabel('Graph Network Representation')
axis off
set(gca,'Ydir','reverse');

% show node classification
figure; set(gcf,'color','w');
hold on
cellOutline = DrawCellOutline(E_tot,E,V,1);

for c = 1:c_len
    v_in_f = C{c};
    currClass = gn.nodeClass(c);
    if currClass==2
        patch(V(v_in_f, 1), V(v_in_f, 2), 1,'FaceColor',[0, 0.4470, 0.7410]);
    elseif currClass == 1
        patch(V(v_in_f, 1), V(v_in_f, 2), 1,'FaceColor',[0.9290, 0.6940, 0.1250]);
    else
        patch(V(v_in_f, 1), V(v_in_f, 2), 1,'FaceColor',[0.6350, 0.0780, 0.1840]);
    end
end
hold on
axis equal; axis tight; xlabel('Train Graph Node Classification')
axis off
set(gca,'Ydir','reverse');

save graphCollective AdjMat gn num_features feature_matrix_n
%% END