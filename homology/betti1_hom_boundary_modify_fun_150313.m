function [boundary_set_hom, boundary_set_hom2] = betti1_hom_boundary_modify_fun_150313(node_coor)

% node_x = node_coor(1, :);
% node_y = node_coor(2, :);

% to find neighbours in order to build 1-simplex
for i=1: length(node_coor)
    node(i).simp1 = [];
    for j=1: length(node_coor)
        if (j==i) continue;
        else
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node(i).simp1 = [node(i).simp1,j];
                % line([node_x(i),node_x(j)],[node_y(i),node_y(j)]);
            end
        end
    end
end

% to find 2-hop nodes
for i = 1: length(node_coor)
    node(i).hop2 = [];
    if ~isempty(node(i).simp1)
        for j = 1: length(node(i).simp1)
            node_hop1 = node(i).simp1(j);
            node(i).hop2 = union(node(i).hop2, node(node_hop1).simp1);
        end
        node(i).hop2 = setdiff(node(i).hop2, node(i).simp1);
        node(i).hop2 = setdiff(node(i).hop2, i);
    end
end

%to find neighbours in order to build 2-simplex
for i = 1: length(node_coor)
    node(i).simp2 = [];
    if length(node(i).simp1) <= 1 continue;
    end
    for j = 1: length(node(i).simp1)-1
        for k = (j+1): length(node(i).simp1)
            if find(node(node(i).simp1(j)).simp1 == node(i).simp1(k))
                node(i).simp2 = [node(i).simp2; node(i).simp1(j), node(i).simp1(k)];
            end
        end
    end
end

% set the fence_flag for all nodes
for i = 1 : length(node_coor)
    if i > length(node_coor) - 20 
        node(i).fence_flag = 1;
    else
        node(i).fence_flag = 0;
    end
end

% for each node, construct the adjacent matrix
for i = 1: length(node_coor)
    if isempty(node(i).simp1) 
        node(i).adj_matrix = [];
    elseif node(i).fence_flag == 1
        simp1_temp = [node(i).simp1, i];
        neighbor_set = sort(simp1_temp);
        size_matrix = length(neighbor_set);
        node(i).adj_matrix = zeros(size_matrix, size_matrix);
        node_fence_index = find(neighbor_set == i);

        node(i).adj_matrix(node_fence_index, :) = ones(1, size_matrix);
        node(i).adj_matrix(:, node_fence_index) = ones(size_matrix, 1);
        
        if size((node(i).simp2),1)
            size_simp2 = size((node(i).simp2),1);
            for index_simp2 = 1: size_simp2
                row_index = find(neighbor_set == node(i).simp2(index_simp2, 1));
                col_index = find(neighbor_set == node(i).simp2(index_simp2, 2));
                node(i).adj_matrix(row_index, col_index) = 1;
                node(i).adj_matrix(col_index, row_index) = 1;
            end
        end
    else
        neighbor_set = sort(node(i).simp1);
        size_matrix = length(neighbor_set);
        node(i).adj_matrix = zeros(size_matrix, size_matrix);
        
        if size((node(i).simp2),1)
            size_simp2 = size((node(i).simp2),1);
            for index_simp2 = 1: size_simp2
                row_index = find(neighbor_set == node(i).simp2(index_simp2, 1));
                col_index = find(neighbor_set == node(i).simp2(index_simp2, 2));
                node(i).adj_matrix(row_index, col_index) = 1;
                node(i).adj_matrix(col_index, row_index) = 1;
            end
        end
    end
end

% determining hole boundary nodes
for i = 1 : length(node_coor)
    if isempty(node(i).adj_matrix)
        node(i).boundary_flag = 1;
        % plot(node_coor(1,i), node_coor(2,i), 'sr'); hold on;
    else
        degree_flag = 0;
        
        sum_degree = sum(node(i).adj_matrix, 2);
        for no_degree = 1: length(sum_degree)
            if sum_degree(no_degree) < 2 
                degree_flag = 1;
                node(i).boundary_flag = 1;
                % plot(node_coor(1,i), node_coor(2,i), 'sr'); hold on;
                break;
            end
        end
        
        if degree_flag == 1 
            continue;
        else
            size_matrix = size(node(i).adj_matrix, 1);
            visit = zeros(size_matrix, 1);
            cycle = zeros(size_matrix, 1);
            visit(1) = 1;
            cycle(1) = 1;

            search_flag = search_matrix (node(i).adj_matrix, visit, cycle, size_matrix, 1);
        end
    
        if search_flag == 1
            node(i).boundary_flag = 0;
        else
            node(i).boundary_flag = 1;
            % plot(node_coor(1,i), node_coor(2,i), 'sr'); hold on;
        end
    end
end

boundary_set_hom = [];
for i = 1 : length(node_coor)
    if node(i).boundary_flag == 1
        boundary_set_hom = [boundary_set_hom, i];
    end
end

while(1)    
    % for each node, find its hole boundary neighbors
    for i = 1: length(node_coor)
        node(i).hole_neighbor = [];
        if isempty(node(i).simp1) 
            continue;
        else
            for j = 1: length(node(i).simp1)
                if node(node(i).simp1(j)).boundary_flag == 1
                    node(i).hole_neighbor = [node(i).hole_neighbor, node(i).simp1(j)];
                end
            end
        end
    end

    % construct the adjacence matrix of hole boundary neighbors for each node
    for i = 1: length(node_coor)
        no_hole_neighbor = length(node(i).hole_neighbor); % hole_neighbor has been sorted
        if no_hole_neighbor >= 2 
            node(i).hole_adj_matrix = zeros(no_hole_neighbor, no_hole_neighbor);

    %         to search the neighborhood information in the adj_matrix
            if node(i).fence_flag == 1
                neighbor_set = sort([node(i).simp1, i]);
            else
                neighbor_set = sort(node(i).simp1);
            end

            for j = 1: no_hole_neighbor - 1
                for k = j+1: no_hole_neighbor
                    hole_node_ID1 = node(i).hole_neighbor(j);
                    hole_node_ID2 = node(i).hole_neighbor(k);

                    row_index = find(neighbor_set == hole_node_ID1);
                    col_index = find(neighbor_set == hole_node_ID2);

                    node(i).hole_adj_matrix(j, k) = node(i).adj_matrix(row_index, col_index);
                    node(i).hole_adj_matrix(k, j) = node(i).hole_adj_matrix(j, k);
                end
            end        
        end
    end

    flag_new_boundary = 0;
    for i = 1 : length(node_coor)
        if node(i).boundary_flag == 1
            flag_new_neigh = 0; % to indicate whether to choose a neighbor as a boundary node
            no_hole_neigh = length(node(i).hole_neighbor);

            if no_hole_neigh == 1
                flag_new_neigh = 1;
            else
                no_cc = no_concom(node(i).hole_adj_matrix);

                if no_cc == 1
                    no_node = size(node(i).hole_adj_matrix, 1);
                    flag_temp = 0;
                    for j = 1 : no_node
                        if sum(node(i).hole_adj_matrix(j, :)) ~= no_node - 1
                            flag_temp = 1;
                            break;
                        end
                    end

                    if flag_temp == 0
                        flag_new_neigh = 1;
                    end
                end
            end

            if flag_new_neigh == 0
                continue;
            else
                neigh_set_temp = [];
                common_set = node(node(i).hole_neighbor(1)).simp1;
                for k = 1 : no_hole_neigh
                    hole_ID = node(i).hole_neighbor(k);
                    common_set = intersect(common_set, node(hole_ID).simp1);
                end
                neigh_set_temp = setdiff(node(i).simp1, common_set);
                neigh_set_temp = setdiff(neigh_set_temp, node(i).hole_neighbor);

                neigh_set_nofence = [];
                for no = 1 : length(neigh_set_temp) % delete fence neighbors
                    neigh_ID = neigh_set_temp(no);
                    if node(neigh_ID).fence_flag == 1
                        continue;
                    else
                        neigh_set_nofence = [neigh_set_nofence, neigh_ID];
                    end
                end

                neigh_set_temp = neigh_set_nofence;

                if ~isempty(neigh_set_temp)
                    degree_set = [];
                    for no = 1 : length(neigh_set_temp)
                        neigh_ID = neigh_set_temp(no);
                        degree_set(no) = length(node(neigh_ID).simp1);
                    end

                    [degree ind] = sort(degree_set);
                    new_neigh_set = neigh_set_temp(ind);

                    flag_found = 0;
                    for no = 1 : length(new_neigh_set)
                        neigh_ID = new_neigh_set(no);
                        no_hole_node = length(node(neigh_ID).hole_neighbor);
                        if no_hole_node >= 2 && ...
                            sum(ismember(node(neigh_ID).hole_neighbor, node(i).hole_neighbor)) == 0
                            flag_found = 1;
                            flag_new_boundary = 1;
                            node(neigh_ID).boundary_flag = 1;
                            % plot(node_coor(1,neigh_ID), node_coor(2,neigh_ID), 'sg'); hold on;
                            break;
                        end
                    end
                    if flag_found == 0
                        neigh_ID = new_neigh_set(1);
                        flag_new_boundary = 1;
                        node(neigh_ID).boundary_flag = 1;
                        % plot(node_coor(1,neigh_ID), node_coor(2,neigh_ID), 'sg'); hold on;
                    end
                end
            end
        end
    end

    if flag_new_boundary == 0
        break;
    end
end

boundary_set_hom2 = [];
for i = 1 : length(node_coor)
    if node(i).boundary_flag == 1
        boundary_set_hom2 = [boundary_set_hom2, i];
    end
end