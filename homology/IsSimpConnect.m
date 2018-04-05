function [flag_Simp] = IsSimpConnect(adj_matrix)
    no_node = size(adj_matrix, 1);
    flag = zeros(no_node, 1);
    node_father = zeros(no_node, 1);
    node_hop = zeros(no_node, 1);

    queue = [1]; 
    flag(1) = 1;

    while(~isempty(queue))
        vertex = queue(1);
        queue(1) = [];

        hop_temp = node_hop(vertex);

        for k = 1 : no_node
            if adj_matrix(vertex, k) == 1 && flag(k) == 0
                node_father(k) = vertex;
                node_hop(k) = hop_temp + 1;
                flag(k) = 1;
                queue  = [queue, k];
            end
        end
    end

    flag_Simp = 1; 
    flag_found = 0;
    for m = 2 : no_node - 1
        for n = m+1 : no_node
            if adj_matrix(m, n) == 1
                if m == node_father(n) || n == node_father(m) || node_father(m) == node_father(n) || adj_matrix(m, node_father(n)) == 1 || adj_matrix(node_father(m), n) == 1
                    continue;
                else
                    m_father = node_father(m);
                    n_father = node_father(n);

                    inter1 = adj_matrix(m, :) & adj_matrix(n, :);
                    inter2 = adj_matrix(m_father, :) & adj_matrix(n_father, :);
                    inter = inter1 &  inter2;

                    if ~isempty(find(inter == 1, 1))
                        continue;
                    else
                        flag_Simp = 0;
                        flag_found = 1;
                        break;
                    end
                end
            end
        end
        if flag_found == 1
            break;
        end
    end
    
    return;
end