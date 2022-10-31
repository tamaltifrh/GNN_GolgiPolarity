function [N,S,PI,theta1N,theta2N,lambdasN,CellCenter]=cellbatchelor(Tensions,Pressures,C,E,V,C_ent,CellEdges,nRadius,xLen,yLen)
% Computes the stress tensor N and its deviatoric part S at locations stored in BoxCenter
% using the Batchelor formula (Ishihara 2012)
%
% BoxCenter{i}=(xi,yi);
% N{i}=stress tensor at position xi,yi
% S{i}=deviator of N at xi,yi (junctional stress : cf Guirao eLife 2016)

% a=nRadius; % size of supracellular regions used to compute averaged stress
c_len = C_ent;
defAngle = 0:15:345;
%% Boxes for averaging the stress


%% Cell centers and areas computation
xcent = zeros(c_len,1); ycent = zeros(c_len,1);

for c = 1:c_len
    v_in_f = C{c};
    xV = V(v_in_f,1)'; yV = V(v_in_f,2)';
    BW = poly2mask(xV,yV,xLen,yLen);
    s = regionprops(BW,'Centroid','Area','Perimeter');
    xcent(c) = s.Centroid(1); ycent(c) = s.Centroid(2);
    CellArea(c) = s.Area; 
    CellCenter{c}=[xcent(c) ycent(c)];
end
avgCellRad = (mean(CellArea)/pi)^0.5;
%% Find edges of each cell (CellEdges)

%% Compute edge lengths
Ne=length(E);
for i=1:Ne
    EdgeLength(i)=sqrt((V(E(i,1),1)-V(E(i,2),1))^2 +(V(E(i,1),2)-V(E(i,2),2))^2);
    EdgeLx(i)=V(E(i,1),1)-V(E(i,2),1);
    EdgeLy(i)=V(E(i,1),2)-V(E(i,2),2);
end

%% Compute N and S in each cell

for c = 1:c_len
    CellRad = (CellArea(c)/pi)^0.5;
    nRadPix = nRadius*CellRad;
    X_cir = xcent(c) + round(nRadPix*cosd(defAngle));
    X_cir = max(1,X_cir);X_cir = min(xLen,X_cir);
    
    Y_cir = ycent(c) + round(nRadPix*sind(defAngle));
    Y_cir = max(1,Y_cir);Y_cir = min(yLen,Y_cir);
    
    Cell_in = inpolygon(xcent,ycent,X_cir,Y_cir);
    BoxCells=find(Cell_in); % cells in the box 
    BoxEdges=[]; % edges in the box
    
    for k=1:length(BoxCells)
        BoxEdges=[BoxEdges CellEdges{BoxCells(k)}];
    end
    BoxEdges=unique(BoxEdges);
    
    Nxx=( -sum(Pressures(BoxCells).*CellArea(BoxCells)) + ...
        sum(Tensions(BoxEdges).*(EdgeLx(BoxEdges).^2)./EdgeLength(BoxEdges)) )...
        / sum(CellArea(BoxCells));
    
    Nxy=( 0 + sum(Tensions(BoxEdges).*EdgeLx(BoxEdges).*EdgeLy(BoxEdges)./EdgeLength(BoxEdges)) )...
        / sum(CellArea(BoxCells));
    
    Nyy=( -sum(Pressures(BoxCells).*CellArea(BoxCells)) + ...
        sum(Tensions(BoxEdges).*(EdgeLy(BoxEdges).^2)./EdgeLength(BoxEdges)) )...
        / sum(CellArea(BoxCells));
    
    N{c}=([Nxx Nxy;Nxy Nyy]);
    PI(c)= 0.5*(Nxx+Nyy);
    
    % Deviator (S)
    Sxx=Nxx-PI(c); Sxy=Nxy; Syy=Nyy-PI(c);
    
    S{c}=([Sxx Sxy;Sxy Syy]);
        
    % eigenvalues and eigenvectors for N
    if sum(sum(isnan(N{c})))==0
        [EigVec,EigVal]=eig(N{c}); % highest eigenvalue is the second one
        theta2N{c}=atan(EigVec(2,2)/EigVec(1,2)); %direction of second eigenvector (principal direction)
        theta1N{c}=atan(EigVec(2,1)/EigVec(1,1)); %direction of first eigenvector
        lambdasN{c}=diag(EigVal);
    else
        theta1N{c}=NaN;
        theta2N{c}=NaN;
        lambdasN{c}=[NaN;NaN];
    end 
    
end

end