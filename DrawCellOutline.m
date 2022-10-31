function isDone = DrawCellOutline(E_tot,E,V,lw)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
isDone = 0;
hold on
for e = 1:E_tot,
    plot(V(E(e,:),1),V(E(e,:),2),'k','LineWidth',lw)
end
hold off
isDone = 1;
end


