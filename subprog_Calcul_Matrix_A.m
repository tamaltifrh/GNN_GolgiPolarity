function [A,B,g]=subprog_Calcul_Matrix_A(V,E,C,V_ent,E_tot,C_tot,evolver)
    %   - V: 1st row row: x position | 2nd row: y pos 
    %   - E: 1st  and 2nd row: vertices (rank) involved
    %   - C: array of cells: 1st to nth row: vertices (rank) involved.
    %            cells are oriented (vertices rank are ccw)
    %   - V_tot: nbr tot of vertices, 
    %   - V_ent: nbr of vertices belonging to entire cells, 
    %   - E_tot: nbr tot of edges,
    %   - C_tot: nbr tot of cells


    %% FORCE TENSOR A
    n = 2*V_ent; % nbr of equations  
    m = E_tot+C_tot; % nbr of unknown parameters (T+P)
    A = zeros(n,m);

    %% ADD TENSION COEFFICIENTS in A
    for cptE=1:E_tot,
        v1 = E(cptE,1);  v2 = E(cptE,2);
        x1 = V(v1,1); x2 = V(v2,1);
        y1 = V(v1,2); y2 = V(v2,2);
        if v1<=V_ent,
            A(2*v1-1,cptE) = (x2-x1)/sqrt((x2-x1)^2+(y2-y1)^2);
            A(2*v1  ,cptE) = (y2-y1)/sqrt((x2-x1)^2+(y2-y1)^2);
        end
        if v2<=V_ent,
            A(2*v2-1,cptE) = -(x2-x1)/sqrt((x2-x1)^2+(y2-y1)^2);
            A(2*v2  ,cptE) = -(y2-y1)/sqrt((x2-x1)^2+(y2-y1)^2);
        end
    end

    %% ADD PRESSURE COEFFICIENTS in A
    for cptC=1:C_tot,
        if evolver==0 %ordering of vertices is different in evolver and TA
        v_in_C = C{cptC}([end 1:end 1]);
        else
        v_in_C = C{cptC}([1 end:-1:1 end]);
        end
        for cptV = 2:length(v_in_C)-1,
            if v_in_C(cptV)<=V_ent,
              A(2*v_in_C(cptV)  ,cptC+E_tot) = ... %% y's
                 +(V(v_in_C(cptV-1),1)-V(v_in_C(cptV+1),1))/2;
              A(2*v_in_C(cptV)-1,cptC+E_tot) = ... %% x's
                 +(V(v_in_C(cptV+1),2)-V(v_in_C(cptV-1),2))/2;
            end
        end
    end
    
    %% PRIOR MATRIX : NORMAL DITRIBUTION of TENSION  ||Bp-g||² = SUM((Ti-1)²)
    B = zeros(m,m);
    for cpt=1:E_tot,
        B(cpt,cpt) = 1;
    end
    
    g(1:E_tot,1) = 1;
    g(1+E_tot:m,1) = 0;