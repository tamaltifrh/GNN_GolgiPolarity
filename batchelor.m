function [N,S,PI,BoxCenter]=batchelor(Tensions,Pressures,C,E,V, SquareSize)
% Computes the stress tensor N and its deviatoric part S at locations stored in BoxCenter
% using the Batchelor formula (Ishihara 2012)
%
% BoxCenter{i}=(xi,yi);
% N{i}=stress tensor at position xi,yi
% S{i}=deviator of N at xi,yi (junctional stress : cf Guirao eLife 2016)

a=SquareSize; % size of supracellular regions used to compute averaged stress

%% Boxes for averaging the stress

Xmin=min(V(:,1)); %min and max X and Ys of the image
Xmax=max(V(:,1));
Ymin=min(V(:,2));
Ymax=max(V(:,2));

X=Xmin-0.1:a:Xmax+0.1; % vectors of box boundaries
Y=Ymin-0.1:a:Ymax+0.1;

if X(end)>Xmax-a/2
    X(end)=Xmax;
else
    X(length(X)+1)=X(end)+a;
end
if Y(end)>Ymax-a/2
    Y(end)=Ymax;
else
    Y(length(Y)+1)=Y(end)+a;
end

%% Cell centers and areas computation

Nc=length(C);
for i=1:Nc
    CellCenter(i,1)=mean(V(C{i},1));
    CellCenter(i,2)=mean(V(C{i},2));
    CellArea(i)=polyarea(V(C{i},1),V(C{i},2));
end

%% Find edges of each cell (CellEdges)

Ne=length(E);
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

        
%% Compute edge lengths

for i=1:Ne
    EdgeLength(i)=sqrt((V(E(i,1),1)-V(E(i,2),1))^2 +(V(E(i,1),2)-V(E(i,2),2))^2);
    EdgeLx(i)=V(E(i,1),1)-V(E(i,2),1);
    EdgeLy(i)=V(E(i,1),2)-V(E(i,2),2);
end

%% Compute N and S in each box

cpt=1;
for i=1:length(X)-1
    for j=1:length(Y)-1
        
        BoxCenter{cpt}=[(X(i)+X(i+1))/2 (Y(j)+Y(j+1))/2];
        
        Cell_in = inpolygon(CellCenter(:,1),CellCenter(:,2),[X(i) X(i+1) X(i+1) X(i)],[Y(j) Y(j) Y(j+1) Y(j+1)]);
        BoxCells=find(Cell_in); % cells in the box n�cpt
        
        BoxEdges=[]; % edges in the box n�cpt
        for k=1:length(BoxCells)
            BoxEdges=[BoxEdges CellEdges{BoxCells(k)}];
        end
        BoxEdges=unique(BoxEdges);
        
        %stress tensor (N=S+PI*Id)
        Nxx=( -sum(Pressures(BoxCells).*CellArea(BoxCells)) + sum(Tensions(BoxEdges).*(EdgeLx(BoxEdges).^2)./EdgeLength(BoxEdges)) ) / sum(CellArea(BoxCells));
        Nxy=( 0 + sum(Tensions(BoxEdges).*EdgeLx(BoxEdges).*EdgeLy(BoxEdges)./EdgeLength(BoxEdges)) ) / sum(CellArea(BoxCells));
        Nyy=( -sum(Pressures(BoxCells).*CellArea(BoxCells)) + sum(Tensions(BoxEdges).*(EdgeLy(BoxEdges).^2)./EdgeLength(BoxEdges)) ) / sum(CellArea(BoxCells));
        
        N{cpt}=([Nxx Nxy;Nxy Nyy]);
        
        %isotropic part (PI*Id)
        PI(cpt)= 0.5*(Nxx+Nyy);
        
        %deviator (S)
        Sxx=Nxx-PI(cpt);
        Sxy=Nxy;
        Syy=Nyy-PI(cpt);
%         check : explicit expression of S directly calculated from Batchelor formula yields the same results
%         Sxx=( 0.5*sum(Tensions(BoxEdges).*(EdgeLx(BoxEdges).^2-EdgeLy(BoxEdges).^2)./EdgeLength(BoxEdges)) ) / sum(CellArea(BoxCells));
%         Sxy=( sum(Tensions(BoxEdges).*EdgeLx(BoxEdges).*EdgeLy(BoxEdges)./EdgeLength(BoxEdges)) ) / sum(CellArea(BoxCells));
%         Syy=( 0.5*sum(Tensions(BoxEdges).*(EdgeLy(BoxEdges).^2-EdgeLx(BoxEdges).^2)./EdgeLength(BoxEdges)) ) / sum(CellArea(BoxCells));
%         
        S{cpt}=([Sxx Sxy;Sxy Syy]);
        
        % eigenvalues and eigenvectors for N
        if sum(sum(isnan(N{cpt})))==0
            [EigVec,EigVal]=eig(N{cpt}); % highest eigenvalue is the second one
            theta2N{cpt}=atan(EigVec(2,2)/EigVec(1,2)); %direction of second eigenvector (principal direction)
            theta1N{cpt}=atan(EigVec(2,1)/EigVec(1,1)); %direction of first eigenvector
            lambdasN{cpt}=diag(EigVal);
        else
            theta1N{cpt}=NaN;
            theta2N{cpt}=NaN;
            lambdasN{cpt}=[NaN;NaN];
        end 
        
        cpt=cpt+1;
    end
end

%% Plot stress tensor
% principal directions and amplitudes of the stress tensor N

figure; set(gcf,'color','w'); hold on
for e = 1:Ne
    plot([V(E(e,1),1) V(E(e,2),1)],[V(E(e,1),2) V(E(e,2),2)],'color',[0.7 0.7 0.7]);
end
for b=1:length(BoxCenter)
    lambdasN{b}=1000*lambdasN{b}; % for plot
    plot( [BoxCenter{b}(1)-0.5*lambdasN{b}(1)*cos(theta1N{b}) BoxCenter{b}(1)+0.5*lambdasN{b}(1)*cos(theta1N{b})] , [BoxCenter{b}(2)-0.5*lambdasN{b}(1)*sin(theta1N{b}) BoxCenter{b}(2)+0.5*lambdasN{b}(1)*sin(theta1N{b})],'linewidth',3,'color','r')
    plot( [BoxCenter{b}(1)-0.5*lambdasN{b}(2)*cos(theta2N{b}) BoxCenter{b}(1)+0.5*lambdasN{b}(2)*cos(theta2N{b})] , [BoxCenter{b}(2)-0.5*lambdasN{b}(2)*sin(theta2N{b}) BoxCenter{b}(2)+0.5*lambdasN{b}(2)*sin(theta2N{b})],'linewidth',3,'color','r')
end
    
axis equal; axis tight; xlabel('Stress Tensor Map'); axis off

end