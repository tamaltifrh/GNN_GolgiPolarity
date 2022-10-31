% show node classification or labels
figure; set(gcf,'color','w');
hold on
cellOutline = DrawCellOutline(E_tot,E,V,1);

for c = 1:c_len
    v_in_f = C{c};
    currClass = nodeClass(c);
    if currClass==2
        patch(V(v_in_f, 1), V(v_in_f, 2), 1,'FaceColor',[0, 0.4470, 0.7410]);
    elseif currClass == 1
        patch(V(v_in_f, 1), V(v_in_f, 2), 1,'FaceColor',[0.9290, 0.6940, 0.1250]);
    else
        patch(V(v_in_f, 1), V(v_in_f, 2), 1,'FaceColor',[0.6350, 0.0780, 0.1840]);
    end
end
hold on
axis equal; axis tight; xlabel('Predicted Graph Node Labels')
axis off
set(gca,'Ydir','reverse');