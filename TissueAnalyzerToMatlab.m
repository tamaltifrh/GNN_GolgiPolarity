EDGES0=[vx_1_x vx_1_y vx_2_x vx_2_y];

EDGES1(:,1)=EDGES0(:,1)+1i*EDGES0(:,2);
EDGES1(:,2)=EDGES0(:,3)+1i*EDGES0(:,4);
%create complex numbers x1+i*y1 and x2+i*y2 to ease some manipulations

VERTICES=unique(EDGES1); %find vertices that are not on the borders of image
sizeX=max(real(VERTICES));
sizeY=max(imag(VERTICES));
BorderV=(mod(real(VERTICES),sizeX-1)==1)+(mod(imag(VERTICES),sizeY-1)==1); 
BorderV=BorderV>0;
[~,ord]=sort(BorderV);
VERTICES=VERTICES(ord);
V_ent=sum(BorderV==0);

for c=1:length(EDGES1)%find edges that are not on the border
    EDGES(c,1)=find(EDGES1(c,1)==VERTICES);
    EDGES(c,2)=find(EDGES1(c,2)==VERTICES);
end
BorderE = (EDGES(:,1)>V_ent) + (EDGES(:,2)>V_ent);
BorderE = BorderE == 2;
EDGES(BorderE,:)=[];
EDGES1(BorderE,:)=[];

aaa=setdiff(1:length(VERTICES),unique(EDGES)); %remove vertices unconnected to edges
VERTICES(aaa)=[];

for c=1:length(EDGES1)
    EDGES(c,1)=find(EDGES1(c,1)==VERTICES);
    EDGES(c,2)=find(EDGES1(c,2)==VERTICES);
end % reset edges accordingly

for c=1:length(vx_coords_cells) %get coordinates of vertices constituting cells
    FACES000=strsplit(vx_coords_cells{c},'#');
    FACES0{c}=[];
    for j=1:length(FACES000)
        FACES00=strsplit(FACES000{j},':');
        for k=1:length(FACES00)
            FACES0{c}=[FACES0{c} str2num(FACES00{k})];
        end
    end
end

for cpt = 1 :length(is_border_cell)
    isbordercellvect(cpt)  = strcmp(is_border_cell{cpt},'true');
end
C_ent = sum(isbordercellvect==0);
[~,ord]=sort(isbordercellvect);

for c1=1:length(FACES0) %get vertices constituting cells instead of their coordinates
    C{c1} = [];
    for c2=2:2:length(FACES0{ord(c1)})
        value=find((FACES0{ord(c1)}(c2-1)+1i*FACES0{ord(c1)}(c2))==VERTICES);
        if length(value)==1
            C{c1} = [C{c1} value] ;
        end
    end
end
C = C(~cellfun('isempty', C)); %remove empty cells

V = [real(VERTICES) imag(VERTICES)];
E = EDGES;

V_tot = length(V); 
E_tot = length(E);
C_tot = length(C);

% figure; set(gcf,'color','w'); hold on;
% coul=jet(C_tot);
% for c=1:C_tot
%     h=patch(V(C{c},1),V(C{c},2),'r');
%     set(h,'facecolor',coul(c,:));
% end
% plot(V(1:V_ent,1),V(1:V_ent,2),'ok')
% plot(V(V_ent+1:end,1),V(V_ent+1:end,2),'ow')
% for c=1:E_tot
%     plot(V(E(c,:),1),V(E(c,:),2),'g')
% end
