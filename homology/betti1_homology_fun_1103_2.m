% This version is the newest and best one up to 13/07/2012
% In this version, we only consider the edges and triangles for each node
% The steps of the algorithm are as follows:
% 1. Each node constructs its 1- and 2-simplex (edges and triangles) and
%    the neighbours of simplices
% 2. Compute the weight of each node, it is 0, 1 or 2
% 3. Delete nodes: if the node has the highest weight among his neighbours 
%    and it can be deleted according to the BNP rules (Betti Number
%    Preserving). The BNP rules are that for each node, if its neighour
%    graph is connected and has at least two nodes and any cycle in the
%    graph can be triangulated, then this node can be deleted.
% 4. Delete edges: delete some special edges in the remained graph, such 
%    edge belongs to only one triangular and after deleting this edge, one 
%    of the end nodes of the edge has a hamilton cycle in its neighbor
%    graph. Repeat steps 3 and 4 until no node or edge can be deleted
% 5. Delete some edges connecting non-boundary nodes and boundary nodes in
%    order to find more boundary edges
% 6. Delete some edges connecting boundary nodes in order to find all
%    boundary edges.
% 7. After step 5 and 6, some special cases may happen. These cases usually
%    occur when there are crossing boundary edges, so it is necessary to
%    delete crossing edges
% 8. Primary cycles detection. Common edges between neighbouring holes are
%    chosen as initiators first. And normal edges are chosen later. In this
%    way, all primary cycles can be found.
% 9. Minimization. Try to minimize primary cycles to discover all
%    fine-grained cycles.


% *************************************************************************
% *************       MODIFY THE WAY TO FIND PRIMARY CYCLES     ***********
% *************************************************************************

% set the degree of boundary edges to be 1 or 2, if the edge has no
% neighbour, its degree is 2 which means the edge belongs to two cycles; if
% the edge has one neighbour, then its degree is 1

% *************************************************************************
% ****** Find all cycles sequentially                           ***********
% *************************************************************************

function [cycle_min, betti1_jplex] = betti1_homology_fun_1103_2(node_coor)

% to find neighbours in order to build 1-simplex
for i=1: length(node_coor)
    node(i).neighbors = [];
    for j=1: length(node_coor)
        if (j==i) continue;
        else
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node(i).neighbors = [node(i).neighbors, j];
                % line([node_x(i),node_x(j)],[node_y(i),node_y(j)]);
            end
        end
    end
end

% construct simplices
for i = 1 : length(node_coor)
%     construct 1-simplex (edge)
    no_neighb = length(node(i).neighbors);
    if no_neighb > 0
        for j = 1 : no_neighb
            new_node = node(i).neighbors(j);
            vert_set_temp = sort([i, new_node]);
            neighb_set_temp = intersect(node(i).neighbors, node(new_node).neighbors);
            simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);

            node(i).simp{1}(j) = simp_temp;
        end
    else
        node(i).simp{1} = [];
    end
    
%     construct 2-simplex (triangle)
    if ~isempty(node(i).simp{1})
        no_edge = size(node(i).simp{1}, 2);
        for j = 1 : no_edge
            vert_set = node(i).simp{1}(j).vert;
            neighb_set = node(i).simp{1}(j).neighb;
            if ~isempty(neighb_set)
                no_neighb = length(neighb_set);
                for k = 1 : no_neighb
                    new_node = neighb_set(k);
                    if new_node < max(setdiff(vert_set, i))
                        continue;
                    else
                        vert_set_temp = union(vert_set, new_node);
                        neighb_set_temp = intersect(neighb_set, node(new_node).neighbors);
                        simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);

                        if size(node(i).simp, 2) == 1 %this is the first new higher simplex
                            node(i).simp{2}(1) = simp_temp;
                        else
                            no_new_simp = size(node(i).simp{2}, 2);
                            node(i).simp{2}(no_new_simp+1) = simp_temp;
                        end
                    end
                end
            end
        end
    end
end

% *************************************************************************
% use Jplex to compute betti1
% *************************************************************************
startJPlex;
s1 = ExplicitStream;

% add vertices to simplical complex
for vertice = 1: length(node_coor)
    s1.add([vertice],[0]);
end

% add edges to simplical complex
for vertice = 1: length(node_coor)
    for i = 1 : length(node(vertice).neighbors)
        neigh_ID = node(vertice).neighbors(i);
        if neigh_ID > vertice
            s1.add([vertice, neigh_ID],[0]);
        end
    end    
end

% add triangulars to simplicial complex
for vertice = 1: length(node_coor)
    if size(node(vertice).simp, 2) == 1
        continue;
    else
        no_tri = size(node(vertice).simp{2}, 2);
        for i = 1 : no_tri
            node_set = node(vertice).simp{2}(i).vert;
            if vertice == min(node_set)
                s1.add(node_set, [0]);
            end
        end
    end
end 
    
s1.dump(0).C;
s1.dump(1).C;

s1.close;
intervals = Plex.Persistence.computeIntervals(s1);
result_temp = Plex.FilterInfinite(intervals);

betti0 = betti2vect(result_temp);

betti1_jplex = betti0(2);

if betti0(1) >= 2 
    cycle_min = [];
    return;
end

% for each node, construct the adjacent matrix
for i = 1: length(node_coor)
    if isempty(node(i).neighbors) 
        node(i).adj_matrix = [];
    else
        neighbor_set = sort(node(i).neighbors);
        size_matrix = length(neighbor_set);
        node(i).adj_matrix = zeros(size_matrix, size_matrix);
        
        if size(node(i).simp, 2) >= 2
            no_tri = size(node(i).simp{2}, 2);
            for index = 1: no_tri
                node_tri = setdiff(node(i).simp{2}(index).vert, i);
                row_index = find(neighbor_set == node_tri(1));
                col_index = find(neighbor_set == node_tri(2));
                node(i).adj_matrix(row_index, col_index) = 1;
                node(i).adj_matrix(col_index, row_index) = 1;
            end
        end
    end
end

% set some flags
for i = 1 : length(node_coor)
    node(i).deleted = 0;
    node(i).neighbors_temp = node(i).neighbors;
    node(i).simp_temp = node(i).simp;
    if i >= length(node_coor) - 19 % fence nodes
        node(i).degree = 0;
        node(i).fence_flag = 1;
        node(i).deletable = 0;
    else
        node(i).fence_flag = 0;
    end
end

node_left_seq = 1: length(node_coor);

% to recursively delete nodes or edges, after this, most boundary edges can
% be found, but there also exists some special edges not found
while(1)
    node_deleted_seq = [];
    % delete nodes according to the rule
    while(1)
        %decide the degree of left nodes
        for i = 1 : length(node_left_seq)
            node_ID = node_left_seq(i);
            if node(node_ID).fence_flag == 1
                continue;
            end

            no_edge = size(node(node_ID).simp_temp{1}, 2);
            flag_edge = 0;
            flag_triangle = 0;
            if no_edge > 0
                for j = 1 : no_edge
                    neighbor_set = node(node_ID).simp_temp{1}(j).neighb;
                    if isempty(neighbor_set)
                        flag_edge = 1;
                        break;
                    elseif length(neighbor_set) == 1
                        flag_triangle = 1;
                    end
                end

                if flag_edge == 1
                    node(node_ID).degree = 0;
                elseif flag_triangle == 1
                    node(node_ID).degree = 2;
                end
            end

            if flag_edge == 1 || flag_triangle == 1
                node(node_ID).deletable = 0;
                continue;
            else
                no_cc = no_concom(node(node_ID).adj_matrix);
                if no_cc > 1 % its neighbors are not connected, it can not be deleted
                    node(node_ID).deletable = 0;
                    node(node_ID).degree = 0;
                else
                    flag_Simp = IsSimpConnect(node(node_ID).adj_matrix);
                    if flag_Simp == 0
                        node(node_ID).deletable = 0;
                        node(node_ID).degree = 3; % in this case, the degree of this node is no smaller than 2, set it to be 3 for convience
                    else
                        node(node_ID).deletable = 1;
                        flag_triangle = 0;
                        no_triangle = size(node(node_ID).simp_temp{2}, 2);
                        for j = 1 : no_triangle
                            neighb_set = node(node_ID).simp_temp{2}(j).neighb;
                            if isempty(neighb_set)
                                flag_triangle = 1;
                                break;
                            end
                        end

                        if flag_triangle == 1
                            node(node_ID).degree = 2;
                        else
                            node(node_ID).degree = 3;
                        end
                    end
                end
            end
        end

        % decide the nodes that can be deleted
        for i = 1 : length(node_left_seq)
            node_ID = node_left_seq(i);
            if (node(node_ID).deletable == 1) && (node(node_ID).degree > 2)
                deg_neigh = zeros(1, length(node(node_ID).neighbors_temp));
                for j = 1 : length(node(node_ID).neighbors_temp)
                    neigh_ID = node(node_ID).neighbors_temp(j);
                    if node(neigh_ID).deletable == 1
                        deg_neigh(j) = node(neigh_ID).degree;
                    else
                        deg_neigh(j) = 0;
                    end
                end

                deg = node(node_ID).degree;
                if deg > max(deg_neigh) % choose the node with highest degree
                    node(node_ID).deleted = 1;
                    % plot(node_coor(1,node_ID), node_coor(2,node_ID), 'sr'); hold on;
                elseif deg == max(deg_neigh) % if two nodes have the same degrees, choose the node with lower ID
                    node_index = find(deg_neigh == deg, 1); 
                    if node_ID < node(node_ID).neighbors_temp(node_index)
                       node(node_ID).deleted = 1;
                       % plot(node_coor(1,node_ID), node_coor(2,node_ID), 'sr'); hold on;
                    end
                end

                if node(node_ID).deleted == 1
                    node_deleted_seq = [node_deleted_seq, node_ID];
                end
            end
        end

        node_left_seq = setdiff(node_left_seq, node_deleted_seq); % the left nodes

        if isempty(node_deleted_seq)
            break;
        else
            while ~isempty(node_deleted_seq)
                node_del_ID = node_deleted_seq(1);
                node_deleted_seq(1) = [];

                if ~isempty(node(node_del_ID).neighbors_temp)
                    for i = 1 : length(node(node_del_ID).neighbors_temp)
                        node_cur_ID = node(node_del_ID).neighbors_temp(i);

                        % delete the node_del_ID from the neighbors of
                        % node_cur_ID
                        position = find(node(node_cur_ID).neighbors_temp == node_del_ID);
                        node(node_cur_ID).adj_matrix(position, :) = [];
                        node(node_cur_ID).adj_matrix(:, position) = [];
                        node(node_cur_ID).neighbors_temp = setdiff(node(node_cur_ID).neighbors_temp, node_del_ID);
                        
                        dim_max = size(node(node_cur_ID).simp_temp, 2);

%                         delete the node_del_ID from the simplices of node_cur_ID
                        for dim_simp = 1 : dim_max
                            no_simp = size(node(node_cur_ID).simp_temp{dim_simp}, 2);
                            simp_del_seq = [];
                            for m = 1 : no_simp
                                if ismember(node_del_ID, node(node_cur_ID).simp_temp{dim_simp}(m).vert)
                                     simp_del_seq = [simp_del_seq, m];
                                elseif ismember(node_del_ID, node(node_cur_ID).simp_temp{dim_simp}(m).neighb)
                                    node(node_cur_ID).simp_temp{dim_simp}(m).neighb = setdiff(node(node_cur_ID).simp_temp{dim_simp}(m).neighb, node_del_ID);
                                end
                            end

                            if ~isempty(simp_del_seq)
                                simp_del_seq = sort(simp_del_seq);
                                for ind = 1 : length(simp_del_seq)
                                    post = simp_del_seq(ind);
                                    node(node_cur_ID).simp_temp{dim_simp}(post-ind+1) = [];
                                end
                            end
                        end

            %             if some certain dimensional simplex is empty, delete that cell structure
                        while dim_simp >= 1
                            if isempty(node(node_cur_ID).simp_temp{dim_simp})
                                node(node_cur_ID).simp_temp(dim_simp) = [];
                                dim_simp = dim_simp - 1;
                            else
                                break;
                            end
                        end
                    end
                end
            end
        end
    end

    % find the edges associated with no or only one triangle
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).boundary_edge = []; 
        no_edge = size(node(node_ID).simp_temp{1}, 2);
        for j = 1 : no_edge
            vert_set = node(node_ID).simp_temp{1}(j).vert;
            neighb_set = node(node_ID).simp_temp{1}(j).neighb;
            neighbor_ID = setdiff(vert_set, node_ID);
            
            if node(node_ID).fence_flag == 1 && node(neighbor_ID).fence_flag == 1
                if isempty(neighb_set)
                    node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                    % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
                end
            else
                if length(neighb_set) <= 1
                    node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                    % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
                end
            end
        end       
    end

    % find the node which has only one boundary edge
    node_special_set = [];
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).boundary_flag = 0;
        if ~isempty(node(node_ID).boundary_edge)
            node(node_ID).boundary_flag = 1;
        end

        if length(node(node_ID).boundary_edge) == 1
            neigh_ID = node(node_ID).boundary_edge(1);
            post = find(node(node_ID).neighbors_temp == neigh_ID);
            common_neigh = node(node_ID).simp_temp{1}(post).neighb;
            if ~isempty(common_neigh)
                node_special_set = [node_special_set, node_ID];
            end
        end
    end

    if isempty(node_special_set)
        break;
    end

    % try to delete the boundary edge one of whose end nodes has only one
    % boundary edge
    flag_edge_deleted = 0; % to indicate whether there is any such edge deleted
    for i = 1 : length(node_special_set)
        node_ID = node_special_set(i);
        if isempty(node(node_ID).boundary_edge)
            node(node_ID).boundary_flag = 0;
            continue;
        end

        adj_matrix_temp = node(node_ID).adj_matrix;
        size_matrix = size(node(node_ID).adj_matrix, 1);
        if size_matrix == 1
            flag_edge_deleted = 1;
            neighbor_edge_ID = node(node_ID).neighbors_temp(1);
            % plot(node_coor(1,node_ID), node_coor(2,node_ID), 'sr'); hold on;
            % line([node_x(node_ID),node_x(neighbor_edge_ID)],[node_y(node_ID),node_y(neighbor_edge_ID)], 'Color','k','LineWidth',2);
            
            node(node_ID).adj_matrix = [];
            node(node_ID).neighbors_temp = [];
            node(node_ID).boundary_edge = [];
            node(node_ID).boundary_flag = 0;

            post2 = find(node(neighbor_edge_ID).neighbors_temp == node_ID);
            node(neighbor_edge_ID).adj_matrix(post2, :) = [];
            node(neighbor_edge_ID).adj_matrix(:, post2) = [];
            node(neighbor_edge_ID).neighbors_temp = setdiff(node(neighbor_edge_ID).neighbors_temp, node_ID);
            node(neighbor_edge_ID).boundary_edge = setdiff(node(neighbor_edge_ID).boundary_edge, node_ID);
%             node(neighbor_edge_ID).boundary_neighbor = setdiff(node(neighbor_edge_ID).boundary_neighbor, node_ID);
        elseif size_matrix < 4
            continue;
        elseif no_concom(adj_matrix_temp) > 1
            continue;
        else
            neighbor_edge_ID = node(node_ID).boundary_edge(1);
            post = find(node(node_ID).neighbors_temp == neighbor_edge_ID);
            common_neigh_ID = node(node_ID).simp_temp{1}(post).neighb;
            if isempty(common_neigh_ID)
                continue;
            end

            post_common_neigh = find(node(node_ID).neighbors_temp == common_neigh_ID);
            common_set = node(node_ID).simp_temp{1}(post_common_neigh).neighb;
            
            % if deleting the edge can make the end node has no new
            % boundary edge, the edge can be deleted
            if length(common_set) > 2
                adj_matrix_temp(post, :) = [];
                adj_matrix_temp(:, post) = [];
                
                flag_edge_deleted = 1;
                % line([node_x(node_ID),node_x(neighbor_edge_ID)],[node_y(node_ID),node_y(neighbor_edge_ID)], 'Color','k','LineWidth',2);

                node(node_ID).adj_matrix = adj_matrix_temp;
                node(node_ID).neighbors_temp = setdiff(node(node_ID).neighbors_temp, neighbor_edge_ID);
                node(node_ID).boundary_edge = [];
                node(node_ID).boundary_flag = 0;
                node(node_ID).simp_temp{1}(post) = [];
                post_common = find(node(node_ID).neighbors_temp == common_neigh_ID);
                node(node_ID).simp_temp{1}(post_common).neighb = setdiff(node(node_ID).simp_temp{1}(post_common).neighb, neighbor_edge_ID);
                no_triangle = size(node(node_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(neighbor_edge_ID, node(node_ID).simp_temp{2}(m).vert)
                        node(node_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(node_ID).simp_temp{2})
                    node(node_ID).simp_temp(2) = [];
                end

                post2 = find(node(neighbor_edge_ID).neighbors_temp == node_ID);
                node(neighbor_edge_ID).adj_matrix(post2, :) = [];
                node(neighbor_edge_ID).adj_matrix(:, post2) = [];
                node(neighbor_edge_ID).neighbors_temp = setdiff(node(neighbor_edge_ID).neighbors_temp, node_ID);
                node(neighbor_edge_ID).boundary_edge = setdiff(node(neighbor_edge_ID).boundary_edge, node_ID);
                node(neighbor_edge_ID).simp_temp{1}(post2) = [];
                post_common2 = find(node(neighbor_edge_ID).neighbors_temp == common_neigh_ID);
                node(neighbor_edge_ID).simp_temp{1}(post_common2).neighb = setdiff(node(neighbor_edge_ID).simp_temp{1}(post_common2).neighb, node_ID);
                no_triangle = size(node(neighbor_edge_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(neighbor_edge_ID).simp_temp{2}(m).vert)
                        node(neighbor_edge_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(neighbor_edge_ID).simp_temp{2})
                    node(neighbor_edge_ID).simp_temp(2) = [];
                end

                % common_neigh_ID = intersect(node(node_ID).neighbors_temp, node(neighbor_edge_ID).neighbors_temp);
                post_temp1 = find(node(common_neigh_ID).neighbors_temp == node_ID);
                post_temp2 = find(node(common_neigh_ID).neighbors_temp == neighbor_edge_ID);
                node(common_neigh_ID).adj_matrix(post_temp1, post_temp2) = 0;
                node(common_neigh_ID).adj_matrix(post_temp2, post_temp1) = 0;
                node(common_neigh_ID).simp_temp{1}(post_temp1).neighb = setdiff(node(common_neigh_ID).simp_temp{1}(post_temp1).neighb, neighbor_edge_ID);
                node(common_neigh_ID).simp_temp{1}(post_temp2).neighb = setdiff(node(common_neigh_ID).simp_temp{1}(post_temp2).neighb, node_ID);
                no_triangle = size(node(common_neigh_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(common_neigh_ID).simp_temp{2}(m).vert) && ismember(neighbor_edge_ID, node(common_neigh_ID).simp_temp{2}(m).vert)
                        node(common_neigh_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(common_neigh_ID).simp_temp{2})
                    node(common_neigh_ID).simp_temp(2) = [];
                end
            end
        end            
    end

    if flag_edge_deleted == 0
        break;
    end
end

% find the edges associated with only one triangle
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_edge = []; 
    no_edge = size(node(node_ID).simp_temp{1}, 2);
    for j = 1 : no_edge
        vert_set = node(node_ID).simp_temp{1}(j).vert;
        neighb_set = node(node_ID).simp_temp{1}(j).neighb;
        neighbor_ID = setdiff(vert_set, node_ID);

        if node(node_ID).fence_flag == 1 && node(neighbor_ID).fence_flag == 1
            if isempty(neighb_set)
                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
            end
        else
            if length(neighb_set) <= 1
                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
            end
        end
    end       
end

% set the boundary_flag for each node 
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_flag = 0;
    if ~isempty(node(node_ID).boundary_edge)
        node(node_ID).boundary_flag = 1;
    end
end

% find the boundary neighbors of each left node
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_neighbor = [];
    for j = 1 : length(node(node_ID).neighbors_temp)
        neighbor_ID = node(node_ID).neighbors_temp(j);
        if node(neighbor_ID).boundary_flag == 1
            node(node_ID).boundary_neighbor = [node(node_ID).boundary_neighbor, neighbor_ID];
        end
    end
end

% try to delete some edges connecting non-boundary node and boundary node,
% after this, some new boundary edges may be found
while(1)
    flag_edge_deleted = 0;
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        if ~isempty(node(node_ID).boundary_neighbor) && node(node_ID).boundary_flag == 0
            boundary_deleted_set = [];
            for j = 1 : length(node(node_ID).neighbors_temp) 
                boundary_ID = node(node_ID).neighbors_temp(j);
                if node(node_ID).fence_flag == 1 && node(boundary_ID).fence_flag == 1
                    continue;
                else
                    common_neigh = node(node_ID).simp_temp{1}(j).neighb;
                    if length(common_neigh) < 2 
                        continue;
                    end
                    post_set = [];
                    adj_matrix_temp1 = [];
                    for index = 1 : length(common_neigh)
                        post_set(index) = find(node(node_ID).neighbors_temp == common_neigh(index));
                        adj_matrix_temp1(index, :) = node(node_ID).adj_matrix(post_set(index), :);
                    end
                    adj_matrix_temp = [];
                    for index_col = 1 : length(post_set)
                        post_temp = post_set(index_col);
                        adj_matrix_temp(:, index_col) = adj_matrix_temp1(:, post_temp);
                    end
                    
                    no_cc = no_concom(adj_matrix_temp);
                    if no_cc > 1
                        continue;
                    else
                        flag_Simp = IsSimpConnect(adj_matrix_temp);
                        if flag_Simp == 0
                            continue;
                        else
                            flag_edge_deleted = 1;
                            boundary_deleted_set = [boundary_deleted_set, j];
                            % line([node_x(node_ID),node_x(boundary_ID)],[node_y(node_ID),node_y(boundary_ID)], 'Color','g','LineWidth',2);
                            
%                             node(node_ID).adj_matrix(post_boundary, :) = [];
%                             node(node_ID).adj_matrix(:, post_boundary) = [];
                            % node(node_ID).neighbors_temp = setdiff(node(node_ID).neighbors_temp, boundary_ID);
                            
                            % modify 1-simplex and their neighbors
                            no_edge = length(node(node_ID).neighbors_temp);
                            for index_edge = 1 : no_edge
                                neighb_temp = node(node_ID).simp_temp{1}(index_edge).neighb;
                                if ismember(boundary_ID, neighb_temp)
                                    node(node_ID).simp_temp{1}(index_edge).neighb = setdiff(node(node_ID).simp_temp{1}(index_edge).neighb, boundary_ID);
                                end
                            end
                            % node(node_ID).simp_temp{1}(j) = [];
                            % modify 2-simplex and their neighors
                            no_tri = size(node(node_ID).simp_temp{2}, 2);
                            delete_tri_set = [];
                            for index_tri = 1 : no_tri
                                if ismember(boundary_ID, node(node_ID).simp_temp{2}(index_tri).vert)
                                    delete_tri_set = [delete_tri_set, index_tri];
                                end
                            end
                            if ~isempty(delete_tri_set)
                                delete_tri_set = sort(delete_tri_set);
                                for ind = 1 : length(delete_tri_set);
                                    post = delete_tri_set(ind);
                                    node(node_ID).simp_temp{2}(post-ind+1) = [];
                                end
                            end
                            no_tri = size(node(node_ID).simp_temp{2}, 2);
                            for index_tri = 1 : no_tri
                                if ismember(boundary_ID, node(node_ID).simp_temp{2}(index_tri).neighb)
                                    node(node_ID).simp_temp{2}(index_tri).neighb = setdiff(node(node_ID).simp_temp{2}(index_tri).neighb, boundary_ID);
                                end
                            end
                            
                            post_node = find(node(boundary_ID).neighbors_temp == node_ID);
                            node(boundary_ID).adj_matrix(post_node, :) = [];
                            node(boundary_ID).adj_matrix(:, post_node) = [];
                            node(boundary_ID).neighbors_temp = setdiff(node(boundary_ID).neighbors_temp, node_ID);
                            node(boundary_ID).simp_temp{1}(post_node) = [];
                            % modify 1-simplex and their neighbors
                            no_edge = length(node(boundary_ID).neighbors_temp);
                            for index_edge = 1 : no_edge
                                neighb_temp = node(boundary_ID).simp_temp{1}(index_edge).neighb;
                                if ismember(node_ID, neighb_temp)
                                    node(boundary_ID).simp_temp{1}(index_edge).neighb = setdiff(node(boundary_ID).simp_temp{1}(index_edge).neighb, node_ID);
                                end
                            end
                            % modify 2-simplex and their neighors
                            no_tri = size(node(boundary_ID).simp_temp{2}, 2);
                            delete_tri_set = [];
                            for index_tri = 1 : no_tri
                                if ismember(node_ID, node(boundary_ID).simp_temp{2}(index_tri).vert)
                                    delete_tri_set = [delete_tri_set, index_tri];
                                end
                            end
                            if ~isempty(delete_tri_set)
                                delete_tri_set = sort(delete_tri_set);
                                for ind = 1 : length(delete_tri_set);
                                    post = delete_tri_set(ind);
                                    node(boundary_ID).simp_temp{2}(post-ind+1) = [];
                                end
                            end
                            no_tri = size(node(boundary_ID).simp_temp{2}, 2);
                            for index_tri = 1 : no_tri
                                if ismember(node_ID, node(boundary_ID).simp_temp{2}(index_tri).neighb)
                                    node(boundary_ID).simp_temp{2}(index_tri).neighb = setdiff(node(boundary_ID).simp_temp{2}(index_tri).neighb, node_ID);
                                end
                            end
                            
                            for no_common_neigh = 1 : length(common_neigh)
                                neigh_ID = common_neigh(no_common_neigh);
                                post_temp1 = find(node(neigh_ID).neighbors_temp == node_ID);
                                post_temp2 = find(node(neigh_ID).neighbors_temp == boundary_ID);
                                node(neigh_ID).adj_matrix(post_temp1, post_temp2) = 0;
                                node(neigh_ID).adj_matrix(post_temp2, post_temp1) = 0;
                                node(neigh_ID).simp_temp{1}(post_temp1).neighb = setdiff(node(neigh_ID).simp_temp{1}(post_temp1).neighb, boundary_ID);
                                node(neigh_ID).simp_temp{1}(post_temp2).neighb = setdiff(node(neigh_ID).simp_temp{1}(post_temp2).neighb, node_ID);
                                % modify 2-simplex and their neighors
                                no_tri = size(node(neigh_ID).simp_temp{2}, 2);
                                for index_tri = 1 : no_tri
                                    if ismember(node_ID, node(neigh_ID).simp_temp{2}(index_tri).vert) && ismember(boundary_ID, node(neigh_ID).simp_temp{2}(index_tri).vert)
                                        node(neigh_ID).simp_temp{2}(index_tri) = [];
                                        break;
                                    end
                                end
                                no_tri = size(node(neigh_ID).simp_temp{2}, 2);
                                for index_tri = 1 : no_tri
                                    if ismember(node_ID, node(neigh_ID).simp_temp{2}(index_tri).vert) && ismember(boundary_ID, node(neigh_ID).simp_temp{2}(index_tri).neighb)
                                        node(neigh_ID).simp_temp{2}(index_tri).neighb = setdiff(node(neigh_ID).simp_temp{2}(index_tri).neighb, boundary_ID);
                                    elseif ismember(boundary_ID, node(neigh_ID).simp_temp{2}(index_tri).vert) && ismember(node_ID, node(neigh_ID).simp_temp{2}(index_tri).neighb)
                                        node(neigh_ID).simp_temp{2}(index_tri).neighb = setdiff(node(neigh_ID).simp_temp{2}(index_tri).neighb, node_ID);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if ~isempty(boundary_deleted_set)
                boundary_deleted_set = sort(boundary_deleted_set);
                for ind = 1 : length(boundary_deleted_set);
                    post = boundary_deleted_set(ind);
                    node(node_ID).neighbors_temp(post-ind+1) = [];
                    node(node_ID).adj_matrix(post-ind+1, :) = [];
                    node(node_ID).adj_matrix(:, post-ind+1) = [];
                    node(node_ID).simp_temp{1}(post-ind+1) = [];
                end
            end
        end
    end
    if flag_edge_deleted == 0
        break;
    end
end

% after last step, it is possible some node can further be deleted
while(1)
    % find the edges associated with only one triangle
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).boundary_edge = []; 
        no_edge = size(node(node_ID).simp_temp{1}, 2);
        for j = 1 : no_edge
            vert_set = node(node_ID).simp_temp{1}(j).vert;
            neighb_set = node(node_ID).simp_temp{1}(j).neighb;
            neighbor_ID = setdiff(vert_set, node_ID);

            if node(node_ID).fence_flag == 1 && node(neighbor_ID).fence_flag == 1
                if isempty(neighb_set)
                    node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                    % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
                end
            else
                if length(neighb_set) <= 1
                    node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                    % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
                end
            end
        end       
    end

    % set the boundary_flag for each node
     node_special_set = [];
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).boundary_flag = 0;
        if ~isempty(node(node_ID).boundary_edge)
            node(node_ID).boundary_flag = 1;
        end

        if length(node(node_ID).boundary_edge) == 1
            neigh_ID = node(node_ID).boundary_edge(1);
            post = find(node(node_ID).neighbors_temp == neigh_ID);
            common_neigh = node(node_ID).simp_temp{1}(post).neighb;
            if ~isempty(common_neigh)
                node_special_set = [node_special_set, node_ID];
            end
        end
    end

    % try to delete the boundary edge one of whose end nodes has only one
    % boundary edge
    flag_edge_deleted = 0; % to indicate whether there is any such edge deleted
    for i = 1 : length(node_special_set)
        node_ID = node_special_set(i);
        if isempty(node(node_ID).boundary_edge)
            node(node_ID).boundary_flag = 0;
            continue;
        end

        adj_matrix_temp = node(node_ID).adj_matrix;
        size_matrix = size(node(node_ID).adj_matrix, 1);
        if size_matrix == 1
            flag_edge_deleted = 1;
            neighbor_edge_ID = node(node_ID).neighbors_temp(1);
            % plot(node_coor(1,node_ID), node_coor(2,node_ID), 'sr'); hold on;
            % line([node_x(node_ID),node_x(neighbor_edge_ID)],[node_y(node_ID),node_y(neighbor_edge_ID)], 'Color','k','LineWidth',2);
            
            node(node_ID).adj_matrix = [];
            node(node_ID).neighbors_temp = [];
            node(node_ID).boundary_edge = [];
            node(node_ID).boundary_flag = 0;

            post2 = find(node(neighbor_edge_ID).neighbors_temp == node_ID);
            node(neighbor_edge_ID).adj_matrix(post2, :) = [];
            node(neighbor_edge_ID).adj_matrix(:, post2) = [];
            node(neighbor_edge_ID).neighbors_temp = setdiff(node(neighbor_edge_ID).neighbors_temp, node_ID);
            node(neighbor_edge_ID).boundary_edge = setdiff(node(neighbor_edge_ID).boundary_edge, node_ID);
%             node(neighbor_edge_ID).boundary_neighbor = setdiff(node(neighbor_edge_ID).boundary_neighbor, node_ID);
        elseif size_matrix < 4
            continue;
        elseif no_concom(adj_matrix_temp) > 1
            continue;
        else
            neighbor_edge_ID = node(node_ID).boundary_edge(1);
            post = find(node(node_ID).neighbors_temp == neighbor_edge_ID);
            common_neigh_ID = node(node_ID).simp_temp{1}(post).neighb;
            if isempty(common_neigh_ID)
                continue;
            end

            post_common_neigh = find(node(node_ID).neighbors_temp == common_neigh_ID);
            common_set = node(node_ID).simp_temp{1}(post_common_neigh).neighb;
            
            % if deleting the edge can make the end node has no new
            % boundary edge, the edge can be deleted
            if length(common_set) > 2
                adj_matrix_temp(post, :) = [];
                adj_matrix_temp(:, post) = [];
                
                flag_edge_deleted = 1;
                % line([node_x(node_ID),node_x(neighbor_edge_ID)],[node_y(node_ID),node_y(neighbor_edge_ID)], 'Color','k','LineWidth',2);

                node(node_ID).adj_matrix = adj_matrix_temp;
                node(node_ID).neighbors_temp = setdiff(node(node_ID).neighbors_temp, neighbor_edge_ID);
                node(node_ID).boundary_edge = [];
                node(node_ID).boundary_flag = 0;
                node(node_ID).simp_temp{1}(post) = [];
                post_common = find(node(node_ID).neighbors_temp == common_neigh_ID);
                node(node_ID).simp_temp{1}(post_common).neighb = setdiff(node(node_ID).simp_temp{1}(post_common).neighb, neighbor_edge_ID);
                no_triangle = size(node(node_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(neighbor_edge_ID, node(node_ID).simp_temp{2}(m).vert)
                        node(node_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(node_ID).simp_temp{2})
                    node(node_ID).simp_temp(2) = [];
                end

                post2 = find(node(neighbor_edge_ID).neighbors_temp == node_ID);
                node(neighbor_edge_ID).adj_matrix(post2, :) = [];
                node(neighbor_edge_ID).adj_matrix(:, post2) = [];
                node(neighbor_edge_ID).neighbors_temp = setdiff(node(neighbor_edge_ID).neighbors_temp, node_ID);
                node(neighbor_edge_ID).boundary_edge = setdiff(node(neighbor_edge_ID).boundary_edge, node_ID);
                node(neighbor_edge_ID).simp_temp{1}(post2) = [];
                post_common2 = find(node(neighbor_edge_ID).neighbors_temp == common_neigh_ID);
                node(neighbor_edge_ID).simp_temp{1}(post_common2).neighb = setdiff(node(neighbor_edge_ID).simp_temp{1}(post_common2).neighb, node_ID);
                no_triangle = size(node(neighbor_edge_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(neighbor_edge_ID).simp_temp{2}(m).vert)
                        node(neighbor_edge_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(neighbor_edge_ID).simp_temp{2})
                    node(neighbor_edge_ID).simp_temp(2) = [];
                end

                % common_neigh_ID = intersect(node(node_ID).neighbors_temp, node(neighbor_edge_ID).neighbors_temp);
                post_temp1 = find(node(common_neigh_ID).neighbors_temp == node_ID);
                post_temp2 = find(node(common_neigh_ID).neighbors_temp == neighbor_edge_ID);
                node(common_neigh_ID).adj_matrix(post_temp1, post_temp2) = 0;
                node(common_neigh_ID).adj_matrix(post_temp2, post_temp1) = 0;
                node(common_neigh_ID).simp_temp{1}(post_temp1).neighb = setdiff(node(common_neigh_ID).simp_temp{1}(post_temp1).neighb, neighbor_edge_ID);
                node(common_neigh_ID).simp_temp{1}(post_temp2).neighb = setdiff(node(common_neigh_ID).simp_temp{1}(post_temp2).neighb, node_ID);
                no_triangle = size(node(common_neigh_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(common_neigh_ID).simp_temp{2}(m).vert) && ismember(neighbor_edge_ID, node(common_neigh_ID).simp_temp{2}(m).vert)
                        node(common_neigh_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(common_neigh_ID).simp_temp{2})
                    node(common_neigh_ID).simp_temp(2) = [];
                end
            end
        end            
    end

    if flag_edge_deleted == 0
        break;
    end
    
    % find the edges associated with only one triangle
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).boundary_edge = []; 
        no_edge = size(node(node_ID).simp_temp{1}, 2);
        for j = 1 : no_edge
            vert_set = node(node_ID).simp_temp{1}(j).vert;
            neighb_set = node(node_ID).simp_temp{1}(j).neighb;
            neighbor_ID = setdiff(vert_set, node_ID);

            if node(node_ID).fence_flag == 1 && node(neighbor_ID).fence_flag == 1
                if isempty(neighb_set)
                    node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                    % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
                end
            else
                if length(neighb_set) <= 1
                    node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                    % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
                end
            end
        end       
    end
    
    % set boundary_flag
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).boundary_flag = 0;
        if ~isempty(node(node_ID).boundary_edge)
            node(node_ID).boundary_flag = 1;
        end
    end
    
    % find the boundary neighbors of each left node
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).boundary_neighbor = [];
        for j = 1 : length(node(node_ID).neighbors_temp)
            neighbor_ID = node(node_ID).neighbors_temp(j);
            if node(neighbor_ID).boundary_flag == 1
                node(node_ID).boundary_neighbor = [node(node_ID).boundary_neighbor, neighbor_ID];
            end
        end
    end
    
    % try to delete some edges connecting non-boundary node and boundary node
    while(1)
        flag_edge_deleted = 0;
        for i = 1 : length(node_left_seq)
            node_ID = node_left_seq(i);
            if ~isempty(node(node_ID).boundary_neighbor) && node(node_ID).boundary_flag == 0
                boundary_deleted_set = [];
                for j = 1 : length(node(node_ID).boundary_neighbor) 
                    boundary_ID = node(node_ID).boundary_neighbor(j);
                    if node(node_ID).fence_flag == 1 && node(boundary_ID).fence_flag == 1
                        continue;
                    else
                        neigh_temp = node(node_ID).neighbors_temp;
                        post = find(neigh_temp == boundary_ID);
                        common_neigh = node(node_ID).simp_temp{1}(post).neighb;
                        if length(common_neigh) < 2
                            continue;
                        end
                        post_set = [];
                        adj_matrix_temp1 = [];
                        for index = 1 : length(common_neigh)
                            post_set(index) = find(neigh_temp == common_neigh(index));
                            adj_matrix_temp1(index, :) = node(node_ID).adj_matrix(post_set(index), :);
                        end
                        adj_matrix_temp = [];
                        for index_col = 1 : length(post_set)
                            post_temp = post_set(index_col);
                            adj_matrix_temp(:, index_col) = adj_matrix_temp1(:, post_temp);
                        end

                        no_cc = no_concom(adj_matrix_temp);
                        if no_cc > 1
                            continue;
                        else
                            flag_Simp = IsSimpConnect(adj_matrix_temp);
                            if flag_Simp == 0
                                continue;
                            else
                                flag_edge_deleted = 1;
                                boundary_deleted_set = [boundary_deleted_set, j];
                                % line([node_x(node_ID),node_x(boundary_ID)],[node_y(node_ID),node_y(boundary_ID)], 'Color','g','LineWidth',2);

                                post_boundary = find(neigh_temp == boundary_ID);
                                node(node_ID).adj_matrix(post_boundary, :) = [];
                                node(node_ID).adj_matrix(:, post_boundary) = [];
                                node(node_ID).neighbors_temp = setdiff(node(node_ID).neighbors_temp, boundary_ID);
                                node(node_ID).simp_temp{1}(post_boundary) = [];
                                % modify 1-simplex and their neighbors
                                no_edge = length(node(node_ID).neighbors_temp);
                                for index_edge = 1 : no_edge
                                    neighb_temp = node(node_ID).simp_temp{1}(index_edge).neighb;
                                    if ismember(boundary_ID, neighb_temp)
                                        node(node_ID).simp_temp{1}(index_edge).neighb = setdiff(node(node_ID).simp_temp{1}(index_edge).neighb, boundary_ID);
                                    end
                                end
                                % modify 2-simplex and their neighors
                                no_tri = size(node(node_ID).simp_temp{2}, 2);
                                delete_tri_set = [];
                                for index_tri = 1 : no_tri
                                    if ismember(boundary_ID, node(node_ID).simp_temp{2}(index_tri).vert)
                                        delete_tri_set = [delete_tri_set, index_tri];
                                    end
                                end
                                if ~isempty(delete_tri_set)
                                    delete_tri_set = sort(delete_tri_set);
                                    for ind = 1 : length(delete_tri_set);
                                        post = delete_tri_set(ind);
                                        node(node_ID).simp_temp{2}(post-ind+1) = [];
                                    end
                                end
                                no_tri = size(node(node_ID).simp_temp{2}, 2);
                                for index_tri = 1 : no_tri
                                    if ismember(boundary_ID, node(node_ID).simp_temp{2}(index_tri).neighb)
                                        node(node_ID).simp_temp{2}(index_tri).neighb = setdiff(node(node_ID).simp_temp{2}(index_tri).neighb, boundary_ID);
                                    end
                                end

                                post_node = find(node(boundary_ID).neighbors_temp == node_ID);
                                node(boundary_ID).adj_matrix(post_node, :) = [];
                                node(boundary_ID).adj_matrix(:, post_node) = [];
                                node(boundary_ID).neighbors_temp = setdiff(node(boundary_ID).neighbors_temp, node_ID);
                                node(boundary_ID).simp_temp{1}(post_node) = [];
                                % modify 1-simplex and their neighbors
                                no_edge = length(node(boundary_ID).neighbors_temp);
                                for index_edge = 1 : no_edge
                                    neighb_temp = node(boundary_ID).simp_temp{1}(index_edge).neighb;
                                    if ismember(node_ID, neighb_temp)
                                        node(boundary_ID).simp_temp{1}(index_edge).neighb = setdiff(node(boundary_ID).simp_temp{1}(index_edge).neighb, node_ID);
                                    end
                                end
                                % modify 2-simplex and their neighors
                                no_tri = size(node(boundary_ID).simp_temp{2}, 2);
                                delete_tri_set = [];
                                for index_tri = 1 : no_tri
                                    if ismember(node_ID, node(boundary_ID).simp_temp{2}(index_tri).vert)
                                        delete_tri_set = [delete_tri_set, index_tri];
                                    end
                                end
                                if ~isempty(delete_tri_set)
                                    delete_tri_set = sort(delete_tri_set);
                                    for ind = 1 : length(delete_tri_set);
                                        post = delete_tri_set(ind);
                                        node(boundary_ID).simp_temp{2}(post-ind+1) = [];
                                    end
                                end
                                no_tri = size(node(boundary_ID).simp_temp{2}, 2);
                                for index_tri = 1 : no_tri
                                    if ismember(node_ID, node(boundary_ID).simp_temp{2}(index_tri).neighb)
                                        node(boundary_ID).simp_temp{2}(index_tri).neighb = setdiff(node(boundary_ID).simp_temp{2}(index_tri).neighb, node_ID);
                                    end
                                end

                                for no_common_neigh = 1 : length(common_neigh)
                                    neigh_ID = common_neigh(no_common_neigh);
                                    post_temp1 = find(node(neigh_ID).neighbors_temp == node_ID);
                                    post_temp2 = find(node(neigh_ID).neighbors_temp == boundary_ID);
                                    node(neigh_ID).adj_matrix(post_temp1, post_temp2) = 0;
                                    node(neigh_ID).adj_matrix(post_temp2, post_temp1) = 0;
                                    node(neigh_ID).simp_temp{1}(post_temp1).neighb = setdiff(node(neigh_ID).simp_temp{1}(post_temp1).neighb, boundary_ID);
                                    node(neigh_ID).simp_temp{1}(post_temp2).neighb = setdiff(node(neigh_ID).simp_temp{1}(post_temp2).neighb, node_ID);
                                    % modify 2-simplex and their neighors
                                    no_tri = size(node(neigh_ID).simp_temp{2}, 2);
                                    for index_tri = 1 : no_tri
                                        if ismember(node_ID, node(neigh_ID).simp_temp{2}(index_tri).vert) && ismember(boundary_ID, node(neigh_ID).simp_temp{2}(index_tri).vert)
                                            node(neigh_ID).simp_temp{2}(index_tri) = [];
                                            break;
                                        end
                                    end
                                    no_tri = size(node(neigh_ID).simp_temp{2}, 2);
                                    for index_tri = 1 : no_tri
                                        if ismember(node_ID, node(neigh_ID).simp_temp{2}(index_tri).vert) && ismember(boundary_ID, node(neigh_ID).simp_temp{2}(index_tri).neighb)
                                            node(neigh_ID).simp_temp{2}(index_tri).neighb = setdiff(node(neigh_ID).simp_temp{2}(index_tri).neighb, boundary_ID);
                                        elseif ismember(boundary_ID, node(neigh_ID).simp_temp{2}(index_tri).vert) && ismember(node_ID, node(neigh_ID).simp_temp{2}(index_tri).neighb)
                                            node(neigh_ID).simp_temp{2}(index_tri).neighb = setdiff(node(neigh_ID).simp_temp{2}(index_tri).neighb, node_ID);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                if ~isempty(boundary_deleted_set)
                    boundary_deleted_set = sort(boundary_deleted_set);
                    for ind = 1 : length(boundary_deleted_set);
                        post = boundary_deleted_set(ind);
                        node(node_ID).boundary_neighbor(post-ind+1) = [];
                    end
                end
            end
        end
        if flag_edge_deleted == 0
            break;
        end
    end

end

% in the following, we try to delete edges connecting boundary edges

% find the edges associated with only one triangle
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_edge = []; 
    no_edge = size(node(node_ID).simp_temp{1}, 2);
    for j = 1 : no_edge
        vert_set = node(node_ID).simp_temp{1}(j).vert;
        neighb_set = node(node_ID).simp_temp{1}(j).neighb;
        neighbor_ID = setdiff(vert_set, node_ID);

        if node(node_ID).fence_flag == 1 && node(neighbor_ID).fence_flag == 1
            if isempty(neighb_set)
                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
            end
        else
            if length(neighb_set) <= 1
                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
            end
        end
    end       
end

% set the boundary_flag for each node 
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_flag = 0;
    if ~isempty(node(node_ID).boundary_edge)
        node(node_ID).boundary_flag = 1;
    end
end

% find the boundary neighbors of each left node
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_neighbor = [];
    for j = 1 : length(node(node_ID).neighbors_temp)
        neighbor_ID = node(node_ID).neighbors_temp(j);
        if node(neighbor_ID).boundary_flag == 1
            node(node_ID).boundary_neighbor = [node(node_ID).boundary_neighbor, neighbor_ID];
        end
    end
end

% try to delete some edges connecting boundary nodes
flag_plot = 0;
while(1)
    flag_edge_deleted = 0;
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        if ~isempty(node(node_ID).boundary_neighbor) && node(node_ID).boundary_flag == 1
            boundary_deleted_set = [];
            for j = 1 : length(node(node_ID).boundary_neighbor) 
                boundary_ID = node(node_ID).boundary_neighbor(j);
                if node(node_ID).fence_flag == 1 && node(boundary_ID).fence_flag == 1
                    continue;
                else
                    neigh_temp = node(node_ID).neighbors_temp;
                    post = find(neigh_temp == boundary_ID);
                    common_neigh = node(node_ID).simp_temp{1}(post).neighb;
                    if length(common_neigh) < 2
                        continue;
                    end
                    post_set = [];
                    adj_matrix_temp1 = [];
                    for index = 1 : length(common_neigh)
                        post_set(index) = find(neigh_temp == common_neigh(index));
                        adj_matrix_temp1(index, :) = node(node_ID).adj_matrix(post_set(index), :);
                    end
                    adj_matrix_temp = [];
                    for index_col = 1 : length(post_set)
                        post_temp = post_set(index_col);
                        adj_matrix_temp(:, index_col) = adj_matrix_temp1(:, post_temp);
                    end
                    
                    no_cc = no_concom(adj_matrix_temp);
                    if no_cc > 1
                        continue;
                    else
                        flag_Simp = IsSimpConnect(adj_matrix_temp);
                        if flag_Simp == 0
                            continue;
                        else
                            flag_plot = 1;
                            flag_edge_deleted = 1;
                            boundary_deleted_set = [boundary_deleted_set, j];
                            % line([node_x(node_ID),node_x(boundary_ID)],[node_y(node_ID),node_y(boundary_ID)], 'Color','g','LineWidth',2);
                            
                            post_boundary = find(neigh_temp == boundary_ID);
                            node(node_ID).adj_matrix(post_boundary, :) = [];
                            node(node_ID).adj_matrix(:, post_boundary) = [];
                            node(node_ID).neighbors_temp = setdiff(node(node_ID).neighbors_temp, boundary_ID);
                            node(node_ID).simp_temp{1}(post_boundary) = [];
                            % modify 1-simplex and their neighbors
                            no_edge = length(node(node_ID).neighbors_temp);
                            for index_edge = 1 : no_edge
                                neighb_temp = node(node_ID).simp_temp{1}(index_edge).neighb;
                                if ismember(boundary_ID, neighb_temp)
                                    node(node_ID).simp_temp{1}(index_edge).neighb = setdiff(node(node_ID).simp_temp{1}(index_edge).neighb, boundary_ID);
                                end
                            end
                            % modify 2-simplex and their neighors
                            no_tri = size(node(node_ID).simp_temp{2}, 2);
                            delete_tri_set = [];
                            for index_tri = 1 : no_tri
                                if ismember(boundary_ID, node(node_ID).simp_temp{2}(index_tri).vert)
                                    delete_tri_set = [delete_tri_set, index_tri];
                                end
                            end
                            if ~isempty(delete_tri_set)
                                delete_tri_set = sort(delete_tri_set);
                                for ind = 1 : length(delete_tri_set);
                                    post = delete_tri_set(ind);
                                    node(node_ID).simp_temp{2}(post-ind+1) = [];
                                end
                            end
                            no_tri = size(node(node_ID).simp_temp{2}, 2);
                            for index_tri = 1 : no_tri
                                if ismember(boundary_ID, node(node_ID).simp_temp{2}(index_tri).neighb)
                                    node(node_ID).simp_temp{2}(index_tri).neighb = setdiff(node(node_ID).simp_temp{2}(index_tri).neighb, boundary_ID);
                                end
                            end
                            
                            post_node = find(node(boundary_ID).neighbors_temp == node_ID);
                            node(boundary_ID).adj_matrix(post_node, :) = [];
                            node(boundary_ID).adj_matrix(:, post_node) = [];
                            node(boundary_ID).neighbors_temp = setdiff(node(boundary_ID).neighbors_temp, node_ID);
                            node(boundary_ID).boundary_neighbor = setdiff(node(boundary_ID).boundary_neighbor, node_ID);
                            node(boundary_ID).simp_temp{1}(post_node) = [];
                            % modify 1-simplex and their neighbors
                            no_edge = length(node(boundary_ID).neighbors_temp);
                            for index_edge = 1 : no_edge
                                neighb_temp = node(boundary_ID).simp_temp{1}(index_edge).neighb;
                                if ismember(node_ID, neighb_temp)
                                    node(boundary_ID).simp_temp{1}(index_edge).neighb = setdiff(node(boundary_ID).simp_temp{1}(index_edge).neighb, node_ID);
                                end
                            end
                            % modify 2-simplex and their neighors
                            no_tri = size(node(boundary_ID).simp_temp{2}, 2);
                            delete_tri_set = [];
                            for index_tri = 1 : no_tri
                                if ismember(node_ID, node(boundary_ID).simp_temp{2}(index_tri).vert)
                                    delete_tri_set = [delete_tri_set, index_tri];
                                end
                            end
                            if ~isempty(delete_tri_set)
                                delete_tri_set = sort(delete_tri_set);
                                for ind = 1 : length(delete_tri_set);
                                    post = delete_tri_set(ind);
                                    node(boundary_ID).simp_temp{2}(post-ind+1) = [];
                                end
                            end
                            no_tri = size(node(boundary_ID).simp_temp{2}, 2);
                            for index_tri = 1 : no_tri
                                if ismember(node_ID, node(boundary_ID).simp_temp{2}(index_tri).neighb)
                                    node(boundary_ID).simp_temp{2}(index_tri).neighb = setdiff(node(boundary_ID).simp_temp{2}(index_tri).neighb, node_ID);
                                end
                            end
                            
                            for no_common_neigh = 1 : length(common_neigh)
                                neigh_ID = common_neigh(no_common_neigh);
                                post_temp1 = find(node(neigh_ID).neighbors_temp == node_ID);
                                post_temp2 = find(node(neigh_ID).neighbors_temp == boundary_ID);
                                node(neigh_ID).adj_matrix(post_temp1, post_temp2) = 0;
                                node(neigh_ID).adj_matrix(post_temp2, post_temp1) = 0;
                                node(neigh_ID).simp_temp{1}(post_temp1).neighb = setdiff(node(neigh_ID).simp_temp{1}(post_temp1).neighb, boundary_ID);
                                node(neigh_ID).simp_temp{1}(post_temp2).neighb = setdiff(node(neigh_ID).simp_temp{1}(post_temp2).neighb, node_ID);
                                % modify 2-simplex and their neighors
                                no_tri = size(node(neigh_ID).simp_temp{2}, 2);
                                for index_tri = 1 : no_tri
                                    if ismember(node_ID, node(neigh_ID).simp_temp{2}(index_tri).vert) && ismember(boundary_ID, node(neigh_ID).simp_temp{2}(index_tri).vert)
                                        node(neigh_ID).simp_temp{2}(index_tri) = [];
                                        break;
                                    end
                                end
                                no_tri = size(node(neigh_ID).simp_temp{2}, 2);
                                for index_tri = 1 : no_tri
                                    if ismember(node_ID, node(neigh_ID).simp_temp{2}(index_tri).vert) && ismember(boundary_ID, node(neigh_ID).simp_temp{2}(index_tri).neighb)
                                        node(neigh_ID).simp_temp{2}(index_tri).neighb = setdiff(node(neigh_ID).simp_temp{2}(index_tri).neighb, boundary_ID);
                                    elseif ismember(boundary_ID, node(neigh_ID).simp_temp{2}(index_tri).vert) && ismember(node_ID, node(neigh_ID).simp_temp{2}(index_tri).neighb)
                                        node(neigh_ID).simp_temp{2}(index_tri).neighb = setdiff(node(neigh_ID).simp_temp{2}(index_tri).neighb, node_ID);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if ~isempty(boundary_deleted_set)
                boundary_deleted_set = sort(boundary_deleted_set);
                for ind = 1 : length(boundary_deleted_set);
                    post = boundary_deleted_set(ind);
                    node(node_ID).boundary_neighbor(post-ind+1) = [];
                end
            end
        end
    end
    if flag_edge_deleted == 0
        break;
    end
end

if flag_plot == 1
    % find the edges associated with only one triangle
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).boundary_edge = []; 
        no_edge = size(node(node_ID).simp_temp{1}, 2);
        for j = 1 : no_edge
            vert_set = node(node_ID).simp_temp{1}(j).vert;
            neighb_set = node(node_ID).simp_temp{1}(j).neighb;
            neighbor_ID = setdiff(vert_set, node_ID);

            if node(node_ID).fence_flag == 1 && node(neighbor_ID).fence_flag == 1
                if isempty(neighb_set)
                    node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                    % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
                end
            else
                if length(neighb_set) <= 1
                    node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                    % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
                end
            end
        end       
    end
    
    % set the boundary_flag for each node 
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).boundary_flag = 0;
        if ~isempty(node(node_ID).boundary_edge)
            node(node_ID).boundary_flag = 1;
        end
    end

    % find the boundary neighbors of each left node
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).boundary_neighbor = [];
        for j = 1 : length(node(node_ID).neighbors_temp)
            neighbor_ID = node(node_ID).neighbors_temp(j);
            if node(neighbor_ID).boundary_flag == 1
                node(node_ID).boundary_neighbor = [node(node_ID).boundary_neighbor, neighbor_ID];
            end
        end
    end
end

% try to delete some special edges, if its two neighboring edges are not
% boundary edges after the deletion of this edge
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    if ~isempty(node(node_ID).boundary_edge)
        no_edge = length(node(node_ID).boundary_edge);
        boundary_deleted = [];
        for j = 1 : no_edge
            neighbor_ID = node(node_ID).boundary_edge(j);
            
            if neighbor_ID < node_ID
                continue;
            end
            
            post = find(node(node_ID).neighbors_temp == neighbor_ID);
            common_neigh = node(node_ID).simp_temp{1}(post).neighb;
            
            if isempty(common_neigh)
                continue;
            else
                post1 = find(node(common_neigh).neighbors_temp == node_ID);
                post2 = find(node(common_neigh).neighbors_temp == neighbor_ID);
                
                common_neigh_set1 = node(common_neigh).simp_temp{1}(post1).neighb;
                common_neigh_set2 = node(common_neigh).simp_temp{1}(post2).neighb;
                if length(common_neigh_set1) >= 3 && length(common_neigh_set2) >= 3 % this edge can be deleted
                    % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','k','LineWidth',2);
                    boundary_deleted = [boundary_deleted, j];
                    node(node_ID).adj_matrix(post, :) = [];
                    node(node_ID).adj_matrix(:, post) = [];
                    node(node_ID).neighbors_temp = setdiff(node(node_ID).neighbors_temp, neighbor_ID);
                    node(node_ID).boundary_neighbor = setdiff(node(node_ID).boundary_neighbor, neighbor_ID);
                    node(node_ID).simp_temp{1}(post) = [];
                    post_common = find(node(node_ID).neighbors_temp == common_neigh);
                    node(node_ID).simp_temp{1}(post_common).neighb = setdiff(node(node_ID).simp_temp{1}(post_common).neighb, neighbor_ID);
                    no_triangle = size(node(node_ID).simp_temp{2}, 2);
                    for m = 1 : no_triangle % delete the triangle from its simplices set
                        if ismember(neighbor_ID, node(node_ID).simp_temp{2}(m).vert)
                            node(node_ID).simp_temp{2}(m) = [];
                            break;
                        end
                    end
                    if isempty(node(node_ID).simp_temp{2})
                        node(node_ID).simp_temp(2) = [];
                    end

                    post2 = find(node(neighbor_ID).neighbors_temp == node_ID);
                    node(neighbor_ID).adj_matrix(post2, :) = [];
                    node(neighbor_ID).adj_matrix(:, post2) = [];
                    node(neighbor_ID).neighbors_temp = setdiff(node(neighbor_ID).neighbors_temp, node_ID);
                    node(neighbor_ID).boundary_edge = setdiff(node(neighbor_ID).boundary_edge, node_ID);
                    node(neighbor_ID).boundary_neighbor = setdiff(node(neighbor_ID).boundary_neighbor, node_ID);
                    if isempty(node(neighbor_ID).boundary_edge)
                        node(neighbor_ID).boundary_flag = 0;
                    end
                    node(neighbor_ID).simp_temp{1}(post2) = [];
                    post_common2 = find(node(neighbor_ID).neighbors_temp == common_neigh);
                    node(neighbor_ID).simp_temp{1}(post_common2).neighb = setdiff(node(neighbor_ID).simp_temp{1}(post_common2).neighb, node_ID);
                    no_triangle = size(node(neighbor_ID).simp_temp{2}, 2);
                    for m = 1 : no_triangle % delete the triangle from its simplices set
                        if ismember(node_ID, node(neighbor_ID).simp_temp{2}(m).vert)
                            node(neighbor_ID).simp_temp{2}(m) = [];
                            break;
                        end
                    end
                    if isempty(node(neighbor_ID).simp_temp{2})
                        node(neighbor_ID).simp_temp(2) = [];
                    end

                    post_temp1 = find(node(common_neigh).neighbors_temp == node_ID);
                    post_temp2 = find(node(common_neigh).neighbors_temp == neighbor_ID);
                    node(common_neigh).adj_matrix(post_temp1, post_temp2) = 0;
                    node(common_neigh).adj_matrix(post_temp2, post_temp1) = 0;
                    node(common_neigh).simp_temp{1}(post_temp1).neighb = setdiff(node(common_neigh).simp_temp{1}(post_temp1).neighb, neighbor_ID);
                    node(common_neigh).simp_temp{1}(post_temp2).neighb = setdiff(node(common_neigh).simp_temp{1}(post_temp2).neighb, node_ID);
                    no_triangle = size(node(common_neigh).simp_temp{2}, 2);
                    for m = 1 : no_triangle % delete the triangle from its simplices set
                        if ismember(node_ID, node(common_neigh).simp_temp{2}(m).vert) && ismember(neighbor_ID, node(common_neigh).simp_temp{2}(m).vert)
                            node(common_neigh).simp_temp{2}(m) = [];
                            break;
                        end
                    end
                    if isempty(node(common_neigh).simp_temp{2})
                        node(common_neigh).simp_temp(2) = [];
                    end
                end
            end
        end
        if ~isempty(boundary_deleted)
            boundary_deleted = sort(boundary_deleted);
            for ind = 1 : length(boundary_deleted);
                post = boundary_deleted(ind);
                node(node_ID).boundary_edge(post-ind+1) = [];
            end
        end
        if isempty(node(node_ID).boundary_edge)
            node(node_ID).boundary_flag = 0;
        end
    else
       node(node_ID).boundary_flag = 0; 
    end
end

% for any node which has only one boundary edge, consider the boundary edge
% connecting its boundary neighbors
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    if length(node(node_ID).boundary_edge) >= 1 && length(node(node_ID).boundary_neighbor) > 2
        neighb_temp_set = setdiff(node(node_ID).boundary_neighbor, node(node_ID).boundary_edge);
        no_boundary_neighbor = length(neighb_temp_set);
        for j = 1 : no_boundary_neighbor - 1
            for k = j+1 : no_boundary_neighbor
                neigh_ID1 = neighb_temp_set(j);
                neigh_ID2 = neighb_temp_set(k);

                post1 = find(node(node_ID).neighbors_temp == neigh_ID1);
                post2 = find(node(node_ID).neighbors_temp == neigh_ID2);
                if node(node_ID).adj_matrix(post1, post2) == 1 && ismember(neigh_ID1, node(neigh_ID2).boundary_edge)
                    common_neigh_set1 = node(node_ID).simp_temp{1}(post1).neighb;
                    common_neigh_set2 = node(node_ID).simp_temp{1}(post2).neighb;

                    if length(common_neigh_set1) == 2 || length(common_neigh_set2) == 2
                        % line([node_x(neigh_ID1),node_x(neigh_ID2)],[node_y(neigh_ID1),node_y(neigh_ID2)], 'Color','k','LineWidth',2);

                        node(node_ID).adj_matrix(post1, post2) = 0;
                        node(node_ID).adj_matrix(post2, post1) = 0;
                        node(node_ID).simp_temp{1}(post1).neighb = setdiff(node(node_ID).simp_temp{1}(post1).neighb, neigh_ID2);
                        node(node_ID).simp_temp{1}(post2).neighb = setdiff(node(node_ID).simp_temp{1}(post2).neighb, neigh_ID1);
                        no_triangle = size(node(node_ID).simp_temp{2}, 2);
                        for m = 1 : no_triangle % delete the triangle from its simplices set
                            if ismember(neigh_ID1, node(node_ID).simp_temp{2}(m).vert) && ismember(neigh_ID2, node(node_ID).simp_temp{2}(m).vert)
                                node(node_ID).simp_temp{2}(m) = [];
                                break;
                            end
                        end
                        if isempty(node(node_ID).simp_temp{2})
                            node(node_ID).simp_temp(2) = [];
                        end

                        post_temp1 = find(node(neigh_ID1).neighbors_temp == neigh_ID2);
                        node(neigh_ID1).adj_matrix(post_temp1, :) = [];
                        node(neigh_ID1).adj_matrix(:, post_temp1) = [];
                        node(neigh_ID1).neighbors_temp = setdiff(node(neigh_ID1).neighbors_temp, neigh_ID2);
                        node(neigh_ID1).boundary_edge = setdiff(node(neigh_ID1).boundary_edge, neigh_ID2);
                        node(neigh_ID1).boundary_neighbor = setdiff(node(neigh_ID1).boundary_neighbor, neigh_ID2);
                        node(neigh_ID1).simp_temp{1}(post_temp1) = [];
                        post_common = find(node(neigh_ID1).neighbors_temp == node_ID);
                        node(neigh_ID1).simp_temp{1}(post_common).neighb = setdiff(node(neigh_ID1).simp_temp{1}(post_common).neighb, neigh_ID2);
                        no_triangle = size(node(neigh_ID1).simp_temp{2}, 2);
                        for m = 1 : no_triangle % delete the triangle from its simplices set
                            if ismember(neigh_ID2, node(neigh_ID1).simp_temp{2}(m).vert)
                                node(neigh_ID1).simp_temp{2}(m) = [];
                                break;
                            end
                        end
                        if isempty(node(neigh_ID1).simp_temp{2})
                            node(neigh_ID1).simp_temp(2) = [];
                        end

                        post_temp2 = find(node(neigh_ID2).neighbors_temp == neigh_ID1);
                        node(neigh_ID2).adj_matrix(post_temp2, :) = [];
                        node(neigh_ID2).adj_matrix(:, post_temp2) = [];
                        node(neigh_ID2).neighbors_temp = setdiff(node(neigh_ID2).neighbors_temp, neigh_ID1);
                        node(neigh_ID2).boundary_edge = setdiff(node(neigh_ID2).boundary_edge, neigh_ID1);
                        node(neigh_ID2).boundary_neighbor = setdiff(node(neigh_ID2).boundary_neighbor, neigh_ID1);
                        node(neigh_ID2).simp_temp{1}(post_temp2) = [];
                        post_common = find(node(neigh_ID2).neighbors_temp == node_ID);
                        node(neigh_ID2).simp_temp{1}(post_common).neighb = setdiff(node(neigh_ID2).simp_temp{1}(post_common).neighb, neigh_ID1);
                        no_triangle = size(node(neigh_ID2).simp_temp{2}, 2);
                        for m = 1 : no_triangle % delete the triangle from its simplices set
                            if ismember(neigh_ID1, node(neigh_ID2).simp_temp{2}(m).vert)
                                node(neigh_ID2).simp_temp{2}(m) = [];
                                break;
                            end
                        end
                        if isempty(node(neigh_ID2).simp_temp{2})
                            node(neigh_ID2).simp_temp(2) = [];
                        end

                        if length(common_neigh_set1) == 2 && length(common_neigh_set2) == 2
                            node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neigh_ID1, neigh_ID2];
                            node(node_ID).boundary_edge = sort(node(node_ID).boundary_edge);
                            % line([node_x(node_ID),node_x(neigh_ID1)],[node_y(node_ID),node_y(neigh_ID1)], 'Color','r','LineWidth',2);
                            % line([node_x(node_ID),node_x(neigh_ID2)],[node_y(node_ID),node_y(neigh_ID2)], 'Color','r','LineWidth',2);

                            node(neigh_ID1).boundary_edge = [node(neigh_ID1).boundary_edge, node_ID];
                            node(neigh_ID1).boundary_edge = sort(node(neigh_ID1).boundary_edge);

                            node(neigh_ID2).boundary_edge = [node(neigh_ID2).boundary_edge, node_ID];
                            node(neigh_ID2).boundary_edge = sort(node(neigh_ID2).boundary_edge);
                        elseif length(common_neigh_set1) == 2
                            node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neigh_ID1];
                            node(node_ID).boundary_edge = sort(node(node_ID).boundary_edge);
                            % line([node_x(node_ID),node_x(neigh_ID1)],[node_y(node_ID),node_y(neigh_ID1)], 'Color','r','LineWidth',2);

                            node(neigh_ID1).boundary_edge = [node(neigh_ID1).boundary_edge, node_ID];
                            node(neigh_ID1).boundary_edge = sort(node(neigh_ID1).boundary_edge);
                        elseif length(common_neigh_set2) == 2
                            node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neigh_ID2];
                            node(node_ID).boundary_edge = sort(node(node_ID).boundary_edge);
                            % line([node_x(node_ID),node_x(neigh_ID2)],[node_y(node_ID),node_y(neigh_ID2)], 'Color','r','LineWidth',2);

                            node(neigh_ID2).boundary_edge = [node(neigh_ID2).boundary_edge, node_ID];
                            node(neigh_ID2).boundary_edge = sort(node(neigh_ID2).boundary_edge);
                        end
                    end 
                end     
            end     
        end
    end
end

% try to delete crossing edges
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    if isempty(node(node_ID).boundary_edge)
        continue;
    end
    boundary_neigh_set = setdiff(node(node_ID).boundary_neighbor, node(node_ID).boundary_edge);
    for j = 1 : length(boundary_neigh_set)
        neigh_ID = boundary_neigh_set(j);
        common_boundary = intersect(node(node_ID).boundary_edge, node(neigh_ID).boundary_edge);
        if length(common_boundary) >= 2
            for k = 1 : length(common_boundary)
                neigh_delete_ID = common_boundary(k);
                post = find(node(node_ID).neighbors_temp == neigh_delete_ID);
                % line([node_x(node_ID),node_x(neigh_delete_ID)],[node_y(node_ID),node_y(neigh_delete_ID)], 'Color','k','LineWidth',2);
                node(node_ID).adj_matrix(post, :) = [];
                node(node_ID).adj_matrix(:, post) = [];
                node(node_ID).neighbors_temp = setdiff(node(node_ID).neighbors_temp, neigh_delete_ID);
                node(node_ID).boundary_edge = setdiff(node(node_ID).boundary_edge, neigh_delete_ID);
                node(node_ID).boundary_neighbor = setdiff(node(node_ID).boundary_neighbor, neigh_delete_ID);
                node(node_ID).simp_temp{1}(post) = [];
                post_common = find(node(node_ID).neighbors_temp == neigh_ID);
                node(node_ID).simp_temp{1}(post_common).neighb = setdiff(node(node_ID).simp_temp{1}(post_common).neighb, neigh_delete_ID);
                no_triangle = size(node(node_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(neigh_delete_ID, node(node_ID).simp_temp{2}(m).vert)
                        node(node_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(node_ID).simp_temp{2})
                    node(node_ID).simp_temp(2) = [];
                end
                
                post2 = find(node(neigh_delete_ID).neighbors_temp == node_ID);
                node(neigh_delete_ID).adj_matrix(post2, :) = [];
                node(neigh_delete_ID).adj_matrix(:, post2) = [];
                node(neigh_delete_ID).neighbors_temp = setdiff(node(neigh_delete_ID).neighbors_temp, node_ID);
                node(neigh_delete_ID).boundary_edge = setdiff(node(neigh_delete_ID).boundary_edge, node_ID);
                node(neigh_delete_ID).boundary_neighbor = setdiff(node(neigh_delete_ID).boundary_neighbor, node_ID);
                node(neigh_delete_ID).simp_temp{1}(post2) = [];
                post_common2 = find(node(neigh_delete_ID).neighbors_temp == neigh_ID);
                node(neigh_delete_ID).simp_temp{1}(post_common2).neighb = setdiff(node(neigh_delete_ID).simp_temp{1}(post_common2).neighb, node_ID);
                no_triangle = size(node(neigh_delete_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(neigh_delete_ID).simp_temp{2}(m).vert)
                        node(neigh_delete_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(neigh_delete_ID).simp_temp{2})
                    node(neigh_delete_ID).simp_temp(2) = [];
                end
                
                post_temp1 = find(node(neigh_ID).neighbors_temp == node_ID);
                post_temp2 = find(node(neigh_ID).neighbors_temp == neigh_delete_ID);
                node(neigh_ID).adj_matrix(post_temp1, post_temp2) = 0;
                node(neigh_ID).adj_matrix(post_temp2, post_temp1) = 0;
                node(neigh_ID).simp_temp{1}(post_temp1).neighb = setdiff(node(neigh_ID).simp_temp{1}(post_temp1).neighb, neigh_delete_ID);
                node(neigh_ID).simp_temp{1}(post_temp2).neighb = setdiff(node(neigh_ID).simp_temp{1}(post_temp2).neighb, node_ID);
                no_triangle = size(node(neigh_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(neigh_ID).simp_temp{2}(m).vert) && ismember(neigh_delete_ID, node(neigh_ID).simp_temp{2}(m).vert)
                        node(neigh_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(neigh_ID).simp_temp{2})
                    node(neigh_ID).simp_temp(2) = [];
                end
            end
            
            post_neigh = find(node(node_ID).neighbors_temp == neigh_ID);
            if length(node(node_ID).simp_temp{1}(post_neigh).neighb) <= 1
                % line([node_x(node_ID),node_x(neigh_ID)],[node_y(node_ID),node_y(neigh_ID)], 'Color','m','LineWidth',2);
                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neigh_ID];
                node(neigh_ID).boundary_edge = [node(neigh_ID).boundary_edge, node_ID];
            end
        end
    end
end 

% try to delete crossing edges
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    boundary_edge_deleted = [];
    for j = 1 : length(node(node_ID).boundary_edge)
        boundary_ID = node(node_ID).boundary_edge(j);
        neigh_temp = node(node_ID).neighbors_temp;
        post = find(neigh_temp == boundary_ID);
        common_neigh_ID = node(node_ID).simp_temp{1}(post).neighb;
        
        if isempty(common_neigh_ID) || ~ismember(common_neigh_ID, node(node_ID).boundary_neighbor)
            continue;
        elseif ismember(common_neigh_ID, node(node_ID).boundary_edge) || ismember(common_neigh_ID, node(boundary_ID).boundary_edge)
            other_neighbor = setdiff(node(common_neigh_ID).boundary_edge, boundary_ID);
            common_node = intersect(node(node_ID).boundary_neighbor, other_neighbor);
            if ~isempty(common_node) % this edge can be deleted
                boundary_edge_deleted = [boundary_edge_deleted, j];
                % line([node_x(node_ID),node_x(boundary_ID)],[node_y(node_ID),node_y(boundary_ID)], 'Color','k','LineWidth',2);
                node(node_ID).adj_matrix(post, :) = [];
                node(node_ID).adj_matrix(:, post) = [];
                node(node_ID).neighbors_temp = setdiff(node(node_ID).neighbors_temp, boundary_ID);
                % node(node_ID).boundary_edge = setdiff(node(node_ID).boundary_edge, boundary_ID);
                node(node_ID).boundary_neighbor = setdiff(node(node_ID).boundary_neighbor, boundary_ID);
                node(node_ID).simp_temp{1}(post) = [];
                post_common = find(node(node_ID).neighbors_temp == common_neigh_ID);
                node(node_ID).simp_temp{1}(post_common).neighb = setdiff(node(node_ID).simp_temp{1}(post_common).neighb, boundary_ID);
                if length(node(node_ID).simp_temp{1}(post_common).neighb) <= 1
                    % line([node_x(node_ID),node_x(common_neigh_ID)],[node_y(node_ID),node_y(common_neigh_ID)], 'Color','m','LineWidth',2);
                    node(node_ID).boundary_edge = [node(node_ID).boundary_edge, common_neigh_ID];
                end
                no_triangle = size(node(node_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(boundary_ID, node(node_ID).simp_temp{2}(m).vert)
                        node(node_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(node_ID).simp_temp{2})
                    node(node_ID).simp_temp(2) = [];
                end

                post2 = find(node(boundary_ID).neighbors_temp == node_ID);
                node(boundary_ID).adj_matrix(post2, :) = [];
                node(boundary_ID).adj_matrix(:, post2) = [];
                node(boundary_ID).neighbors_temp = setdiff(node(boundary_ID).neighbors_temp, node_ID);
                node(boundary_ID).boundary_edge = setdiff(node(boundary_ID).boundary_edge, node_ID);
                node(boundary_ID).boundary_neighbor = setdiff(node(boundary_ID).boundary_neighbor, node_ID);
                node(boundary_ID).simp_temp{1}(post2) = [];
                post_common2 = find(node(boundary_ID).neighbors_temp == common_neigh_ID);
                node(boundary_ID).simp_temp{1}(post_common2).neighb = setdiff(node(boundary_ID).simp_temp{1}(post_common2).neighb, node_ID);
                if length(node(boundary_ID).simp_temp{1}(post_common2).neighb) <= 1
                    % line([node_x(boundary_ID),node_x(common_neigh_ID)],[node_y(boundary_ID),node_y(common_neigh_ID)], 'Color','m','LineWidth',2);
                    node(boundary_ID).boundary_edge = [node(boundary_ID).boundary_edge, common_neigh_ID];
                end
                node(boundary_ID).boundary_edge = sort(node(boundary_ID).boundary_edge);
                no_triangle = size(node(boundary_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(boundary_ID).simp_temp{2}(m).vert)
                        node(boundary_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(boundary_ID).simp_temp{2})
                    node(boundary_ID).simp_temp(2) = [];
                end

                % common_neigh_ID = intersect(node(node_ID).neighbors_temp, node(neighbor_edge_ID).neighbors_temp);
                post_temp1 = find(node(common_neigh_ID).neighbors_temp == node_ID);
                post_temp2 = find(node(common_neigh_ID).neighbors_temp == boundary_ID);
                node(common_neigh_ID).adj_matrix(post_temp1, post_temp2) = 0;
                node(common_neigh_ID).adj_matrix(post_temp2, post_temp1) = 0;
                node(common_neigh_ID).simp_temp{1}(post_temp1).neighb = setdiff(node(common_neigh_ID).simp_temp{1}(post_temp1).neighb, boundary_ID);
                if ~ismember(node_ID, node(common_neigh_ID).boundary_edge) && length(node(common_neigh_ID).simp_temp{1}(post_temp1).neighb) <= 1
                    node(common_neigh_ID).boundary_edge = [node(common_neigh_ID).boundary_edge, node_ID];
                end
                node(common_neigh_ID).simp_temp{1}(post_temp2).neighb = setdiff(node(common_neigh_ID).simp_temp{1}(post_temp2).neighb, node_ID);
                if ~ismember(boundary_ID, node(common_neigh_ID).boundary_edge) && length(node(common_neigh_ID).simp_temp{1}(post_temp2).neighb) <= 1
                    node(common_neigh_ID).boundary_edge = [node(common_neigh_ID).boundary_edge, boundary_ID];
                end
                node(common_neigh_ID).boundary_edge = sort(node(common_neigh_ID).boundary_edge);
                no_triangle = size(node(common_neigh_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(common_neigh_ID).simp_temp{2}(m).vert) && ismember(boundary_ID, node(common_neigh_ID).simp_temp{2}(m).vert)
                        node(common_neigh_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(common_neigh_ID).simp_temp{2})
                    node(common_neigh_ID).simp_temp(2) = [];
                end
            end
        end
    end
    if ~isempty(boundary_edge_deleted)
        boundary_edge_deleted = sort(boundary_edge_deleted);
        for ind = 1 : length(boundary_edge_deleted);
            post = boundary_edge_deleted(ind);
            node(node_ID).boundary_edge(post-ind+1) = [];
        end
        node(node_ID).boundary_edge = sort(node(node_ID).boundary_edge);
    end
end

% try to delete crossing edges
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    boundary_edge_deleted = [];
    for j = 1 : length(node(node_ID).boundary_edge)
        boundary_ID = node(node_ID).boundary_edge(j);
        neigh_temp = node(node_ID).neighbors_temp;
        post = find(neigh_temp == boundary_ID);
        common_neigh_ID = node(node_ID).simp_temp{1}(post).neighb;
        
        if isempty(common_neigh_ID) || ~ismember(common_neigh_ID, node(node_ID).boundary_neighbor)
            continue;
        elseif ~ismember(common_neigh_ID, node(node_ID).boundary_edge) && ~ismember(common_neigh_ID, node(boundary_ID).boundary_edge)
%             other_neighbor = setdiff(node(common_neigh_ID).boundary_neighbor, boundary_ID);
%             common_node = intersect(node(node_ID).boundary_neighbor, other_neighbor);
            common_node = intersect(node(node_ID).boundary_neighbor, node(common_neigh_ID).boundary_edge);
            if ~isempty(common_node) % this edge can be deleted
                boundary_edge_deleted = [boundary_edge_deleted, j];
                % line([node_x(node_ID),node_x(boundary_ID)],[node_y(node_ID),node_y(boundary_ID)], 'Color','k','LineWidth',2);
                node(node_ID).adj_matrix(post, :) = [];
                node(node_ID).adj_matrix(:, post) = [];
                node(node_ID).neighbors_temp = setdiff(node(node_ID).neighbors_temp, boundary_ID);
                % node(node_ID).boundary_edge = setdiff(node(node_ID).boundary_edge, boundary_ID);
                node(node_ID).boundary_neighbor = setdiff(node(node_ID).boundary_neighbor, boundary_ID);
                node(node_ID).simp_temp{1}(post) = [];
                post_common = find(node(node_ID).neighbors_temp == common_neigh_ID);
                node(node_ID).simp_temp{1}(post_common).neighb = setdiff(node(node_ID).simp_temp{1}(post_common).neighb, boundary_ID);
                if length(node(node_ID).simp_temp{1}(post_common).neighb) <= 1
                    % line([node_x(node_ID),node_x(common_neigh_ID)],[node_y(node_ID),node_y(common_neigh_ID)], 'Color','m','LineWidth',2);
                    node(node_ID).boundary_edge = [node(node_ID).boundary_edge, common_neigh_ID];
                end
                no_triangle = size(node(node_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(boundary_ID, node(node_ID).simp_temp{2}(m).vert)
                        node(node_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(node_ID).simp_temp{2})
                    node(node_ID).simp_temp(2) = [];
                end

                post2 = find(node(boundary_ID).neighbors_temp == node_ID);
                node(boundary_ID).adj_matrix(post2, :) = [];
                node(boundary_ID).adj_matrix(:, post2) = [];
                node(boundary_ID).neighbors_temp = setdiff(node(boundary_ID).neighbors_temp, node_ID);
                node(boundary_ID).boundary_edge = setdiff(node(boundary_ID).boundary_edge, node_ID);
                node(boundary_ID).boundary_neighbor = setdiff(node(boundary_ID).boundary_neighbor, node_ID);
                node(boundary_ID).simp_temp{1}(post2) = [];
                post_common2 = find(node(boundary_ID).neighbors_temp == common_neigh_ID);
                node(boundary_ID).simp_temp{1}(post_common2).neighb = setdiff(node(boundary_ID).simp_temp{1}(post_common2).neighb, node_ID);
                if length(node(boundary_ID).simp_temp{1}(post_common2).neighb) <= 1
                    % line([node_x(boundary_ID),node_x(common_neigh_ID)],[node_y(boundary_ID),node_y(common_neigh_ID)], 'Color','m','LineWidth',2);
                    node(boundary_ID).boundary_edge = [node(boundary_ID).boundary_edge, common_neigh_ID];
                end
                node(boundary_ID).boundary_edge = sort(node(boundary_ID).boundary_edge);
                no_triangle = size(node(boundary_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(boundary_ID).simp_temp{2}(m).vert)
                        node(boundary_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(boundary_ID).simp_temp{2})
                    node(boundary_ID).simp_temp(2) = [];
                end

                % common_neigh_ID = intersect(node(node_ID).neighbors_temp, node(neighbor_edge_ID).neighbors_temp);
                post_temp1 = find(node(common_neigh_ID).neighbors_temp == node_ID);
                post_temp2 = find(node(common_neigh_ID).neighbors_temp == boundary_ID);
                node(common_neigh_ID).adj_matrix(post_temp1, post_temp2) = 0;
                node(common_neigh_ID).adj_matrix(post_temp2, post_temp1) = 0;
                node(common_neigh_ID).simp_temp{1}(post_temp1).neighb = setdiff(node(common_neigh_ID).simp_temp{1}(post_temp1).neighb, boundary_ID);
                if length(node(common_neigh_ID).simp_temp{1}(post_temp1).neighb) <= 1
                    node(common_neigh_ID).boundary_edge = [node(common_neigh_ID).boundary_edge, node_ID];
                end
                node(common_neigh_ID).simp_temp{1}(post_temp2).neighb = setdiff(node(common_neigh_ID).simp_temp{1}(post_temp2).neighb, node_ID);
                if length(node(common_neigh_ID).simp_temp{1}(post_temp2).neighb) <= 1
                    node(common_neigh_ID).boundary_edge = [node(common_neigh_ID).boundary_edge, boundary_ID];
                end
                node(common_neigh_ID).boundary_edge = sort(node(common_neigh_ID).boundary_edge);
                no_triangle = size(node(common_neigh_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(common_neigh_ID).simp_temp{2}(m).vert) && ismember(boundary_ID, node(common_neigh_ID).simp_temp{2}(m).vert)
                        node(common_neigh_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(common_neigh_ID).simp_temp{2})
                    node(common_neigh_ID).simp_temp(2) = [];
                end
            end
        end
    end
    if ~isempty(boundary_edge_deleted)
        boundary_edge_deleted = sort(boundary_edge_deleted);
        for ind = 1 : length(boundary_edge_deleted);
            post = boundary_edge_deleted(ind);
            node(node_ID).boundary_edge(post-ind+1) = [];
        end
        node(node_ID).boundary_edge = sort(node(node_ID).boundary_edge);
    end
end    

% find the edges associated with only one triangle
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_edge = []; 
    no_edge = size(node(node_ID).simp_temp{1}, 2);
    for j = 1 : no_edge
        vert_set = node(node_ID).simp_temp{1}(j).vert;
        neighb_set = node(node_ID).simp_temp{1}(j).neighb;
        neighbor_ID = setdiff(vert_set, node_ID);

        if node(node_ID).fence_flag == 1 && node(neighbor_ID).fence_flag == 1
            if isempty(neighb_set)
                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
            end
        else
            if length(neighb_set) <= 1
                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
            end
        end
    end       
end

% set the boundary_flag for each node 
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_flag = 0;
    if ~isempty(node(node_ID).boundary_edge)
        node(node_ID).boundary_flag = 1;
    end
end

% find the boundary neighbors of each left node
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_neighbor = [];
    for j = 1 : length(node(node_ID).neighbors_temp)
        neighbor_ID = node(node_ID).neighbors_temp(j);
        if node(neighbor_ID).boundary_flag == 1
            node(node_ID).boundary_neighbor = [node(node_ID).boundary_neighbor, neighbor_ID];
        end
    end
end

% delete some nodes who have only one boundary neighbor or two connected
% boundary neighbors
while(1)
    flag_node_deleted = 0;
    for i = 1 : length(node_left_seq)
        node_ID = node_left_seq(i);
        while(length(node(node_ID).boundary_edge) == 1 && length(node(node_ID).neighbors_temp) == 1)
            % plot(node_coor(1,node_ID), node_coor(2,node_ID), 'dk'); hold on;
            neigh_ID = node(node_ID).boundary_neighbor(1);
            % line([node_x(node_ID),node_x(neigh_ID)],[node_y(node_ID),node_y(neigh_ID)], 'Color','k','LineWidth',2);
            node(node_ID).boundary_edge = [];
            node(node_ID).boundary_neighbor = [];
            node(node_ID).boundary_flag = 0;
            node(node_ID).simp_temp(1) = [];
            node(node_ID).neighbors_temp = [];
            node(node_ID).adj_matrix = [];
            post = find(node(neigh_ID).neighbors_temp == node_ID);
            node(neigh_ID).adj_matrix(post, :) = [];
            node(neigh_ID).adj_matrix(:, post) = [];
            node(neigh_ID).neighbors_temp = setdiff(node(neigh_ID).neighbors_temp, node_ID);
            node(neigh_ID).boundary_edge = setdiff(node(neigh_ID).boundary_edge, node_ID);
            node(neigh_ID).boundary_neighbor = setdiff(node(neigh_ID).boundary_neighbor, node_ID);
            node(neigh_ID).simp_temp{1}(post) = [];
            node_ID = neigh_ID;
        end

        if length(node(node_ID).boundary_edge) == 2 && length(node(node_ID).neighbors_temp) == 2
            neigh_ID1 = node(node_ID).boundary_edge(1);
            neigh_ID2 = node(node_ID).boundary_edge(2);
            post1 = find(node(node_ID).neighbors_temp == neigh_ID1);
            post2 = find(node(node_ID).neighbors_temp == neigh_ID2);
            if node(node_ID).adj_matrix(post1, post2) == 1
                flag_node_deleted = 1;
                node(node_ID).boundary_edge = [];
                node(node_ID).boundary_neighbor = [];
                node(node_ID).boundary_flag = 0;
                node(node_ID).neighbors_temp = [];
                node(node_ID).simp_temp(2) = [];
                node(node_ID).simp_temp(1) = [];
                % plot(node_coor(1,node_ID), node_coor(2,node_ID), 'dk'); hold on;
                % line([node_x(node_ID),node_x(neigh_ID1)],[node_y(node_ID),node_y(neigh_ID1)], 'Color','k','LineWidth',2);
                % line([node_x(node_ID),node_x(neigh_ID2)],[node_y(node_ID),node_y(neigh_ID2)], 'Color','k','LineWidth',2);

                post_node1 = find(node(neigh_ID1).neighbors_temp == node_ID);
                node(neigh_ID1).simp_temp{1}(post_node1) = [];
                node(neigh_ID1).adj_matrix(post_node1, :) = [];
                node(neigh_ID1).adj_matrix(:, post_node1) = [];
                node(neigh_ID1).neighbors_temp = setdiff(node(neigh_ID1).neighbors_temp, node_ID);
                node(neigh_ID1).boundary_edge = setdiff(node(neigh_ID1).boundary_edge, node_ID);
                node(neigh_ID1).boundary_neighbor = setdiff(node(neigh_ID1).boundary_neighbor, node_ID);
                post_neigh2 = find(node(neigh_ID1).neighbors_temp == neigh_ID2);
                node(neigh_ID1).simp_temp{1}(post_neigh2).neighb = setdiff(node(neigh_ID1).simp_temp{1}(post_neigh2).neighb, node_ID);
                no_triangle = size(node(neigh_ID1).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(neigh_ID1).simp_temp{2}(m).vert)
                        node(neigh_ID1).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(neigh_ID1).simp_temp{2})
                    node(neigh_ID1).simp_temp(2) = [];
                end

                post_node2 = find(node(neigh_ID2).neighbors_temp == node_ID);
                node(neigh_ID2).simp_temp{1}(post_node2) = [];
                node(neigh_ID2).adj_matrix(post_node2, :) = [];
                node(neigh_ID2).adj_matrix(:, post_node2) = [];
                node(neigh_ID2).neighbors_temp = setdiff(node(neigh_ID2).neighbors_temp, node_ID);
                node(neigh_ID2).boundary_edge = setdiff(node(neigh_ID2).boundary_edge, node_ID);
                node(neigh_ID2).boundary_neighbor = setdiff(node(neigh_ID2).boundary_neighbor, node_ID);
                post_neigh1 = find(node(neigh_ID2).neighbors_temp == neigh_ID1);
                node(neigh_ID2).simp_temp{1}(post_neigh1).neighb = setdiff(node(neigh_ID2).simp_temp{1}(post_neigh1).neighb, node_ID);
                no_triangle = size(node(neigh_ID2).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(neigh_ID2).simp_temp{2}(m).vert)
                        node(neigh_ID2).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(neigh_ID2).simp_temp{2})
                    node(neigh_ID2).simp_temp(2) = [];
                end

                if ~ismember(neigh_ID1, node(neigh_ID2).boundary_edge) && length(node(neigh_ID1).simp_temp{1}(post_neigh2).neighb) <= 1
                    node(neigh_ID1).boundary_edge = [node(neigh_ID1).boundary_edge, neigh_ID2];
                    node(neigh_ID2).boundary_edge = [node(neigh_ID2).boundary_edge, neigh_ID1];
                    node(neigh_ID1).boundary_edge = sort(node(neigh_ID1).boundary_edge);
                    node(neigh_ID2).boundary_edge = sort(node(neigh_ID2).boundary_edge);
                    % line([node_x(neigh_ID1),node_x(neigh_ID2)],[node_y(neigh_ID1),node_y(neigh_ID2)], 'Color','r','LineWidth',2);
                end
            end
        end
    end
    if flag_node_deleted == 0
        break;
    end
end

for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    if length(node(node_ID).boundary_edge) == 1
        neighbor_edge_ID = node(node_ID).boundary_edge(1);
        post = find(node(node_ID).neighbors_temp == neighbor_edge_ID);
        common_neigh_ID = node(node_ID).simp_temp{1}(post).neighb;
        if ~isempty(common_neigh_ID)
            adj_matrix_temp = node(node_ID).adj_matrix;
            if size(adj_matrix_temp, 2) < 4 || no_concom(adj_matrix_temp) > 1
                continue;
            end
            
            post_common_neigh = find(node(node_ID).neighbors_temp == common_neigh_ID);
            common_set = node(node_ID).simp_temp{1}(post_common_neigh).neighb;
            
            % if deleting the edge can make the end node has no new
            % boundary edge, the edge can be deleted
            if length(common_set) > 2
                adj_matrix_temp(post, :) = [];
                adj_matrix_temp(:, post) = [];
                
                % line([node_x(node_ID),node_x(neighbor_edge_ID)],[node_y(node_ID),node_y(neighbor_edge_ID)], 'Color','k','LineWidth',2);

                node(node_ID).adj_matrix = adj_matrix_temp;
                node(node_ID).neighbors_temp = setdiff(node(node_ID).neighbors_temp, neighbor_edge_ID);
                node(node_ID).boundary_edge = [];
                node(node_ID).boundary_flag = 0;
                node(node_ID).simp_temp{1}(post) = [];
                post_common = find(node(node_ID).neighbors_temp == common_neigh_ID);
                node(node_ID).simp_temp{1}(post_common).neighb = setdiff(node(node_ID).simp_temp{1}(post_common).neighb, neighbor_edge_ID);
                no_triangle = size(node(node_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(neighbor_edge_ID, node(node_ID).simp_temp{2}(m).vert)
                        node(node_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(node_ID).simp_temp{2})
                    node(node_ID).simp_temp(2) = [];
                end

                post2 = find(node(neighbor_edge_ID).neighbors_temp == node_ID);
                node(neighbor_edge_ID).adj_matrix(post2, :) = [];
                node(neighbor_edge_ID).adj_matrix(:, post2) = [];
                node(neighbor_edge_ID).neighbors_temp = setdiff(node(neighbor_edge_ID).neighbors_temp, node_ID);
                node(neighbor_edge_ID).boundary_edge = setdiff(node(neighbor_edge_ID).boundary_edge, node_ID);
                node(neighbor_edge_ID).simp_temp{1}(post2) = [];
                post_common2 = find(node(neighbor_edge_ID).neighbors_temp == common_neigh_ID);
                node(neighbor_edge_ID).simp_temp{1}(post_common2).neighb = setdiff(node(neighbor_edge_ID).simp_temp{1}(post_common2).neighb, node_ID);
                no_triangle = size(node(neighbor_edge_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(neighbor_edge_ID).simp_temp{2}(m).vert)
                        node(neighbor_edge_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(neighbor_edge_ID).simp_temp{2})
                    node(neighbor_edge_ID).simp_temp(2) = [];
                end

                % common_neigh_ID = intersect(node(node_ID).neighbors_temp, node(neighbor_edge_ID).neighbors_temp);
                post_temp1 = find(node(common_neigh_ID).neighbors_temp == node_ID);
                post_temp2 = find(node(common_neigh_ID).neighbors_temp == neighbor_edge_ID);
                node(common_neigh_ID).adj_matrix(post_temp1, post_temp2) = 0;
                node(common_neigh_ID).adj_matrix(post_temp2, post_temp1) = 0;
                node(common_neigh_ID).simp_temp{1}(post_temp1).neighb = setdiff(node(common_neigh_ID).simp_temp{1}(post_temp1).neighb, neighbor_edge_ID);
                node(common_neigh_ID).simp_temp{1}(post_temp2).neighb = setdiff(node(common_neigh_ID).simp_temp{1}(post_temp2).neighb, node_ID);
                no_triangle = size(node(common_neigh_ID).simp_temp{2}, 2);
                for m = 1 : no_triangle % delete the triangle from its simplices set
                    if ismember(node_ID, node(common_neigh_ID).simp_temp{2}(m).vert) && ismember(neighbor_edge_ID, node(common_neigh_ID).simp_temp{2}(m).vert)
                        node(common_neigh_ID).simp_temp{2}(m) = [];
                        break;
                    end
                end
                if isempty(node(common_neigh_ID).simp_temp{2})
                    node(common_neigh_ID).simp_temp(2) = [];
                end
                continue;
            end
        end
        if length(node(node_ID).boundary_neighbor) >= 2
            neighb_temp_set = setdiff(node(node_ID).boundary_neighbor, node(node_ID).boundary_edge);
            no_boundary_neighbor = length(neighb_temp_set);
            for j = 1 : no_boundary_neighbor - 1
                for k = j+1 : no_boundary_neighbor
                    neigh_ID1 = neighb_temp_set(j);
                    neigh_ID2 = neighb_temp_set(k);
                    
                    post1 = find(node(node_ID).neighbors_temp == neigh_ID1);
                    post2 = find(node(node_ID).neighbors_temp == neigh_ID2);
                    if node(node_ID).adj_matrix(post1, post2) == 1 && ismember(neigh_ID1, node(neigh_ID2).boundary_edge)
                        common_neigh_set1 = node(node_ID).simp_temp{1}(post1).neighb;
                        common_neigh_set2 = node(node_ID).simp_temp{1}(post2).neighb;
                        
                        if length(common_neigh_set1) == 2 || length(common_neigh_set2) == 2
                            % line([node_x(neigh_ID1),node_x(neigh_ID2)],[node_y(neigh_ID1),node_y(neigh_ID2)], 'Color','k','LineWidth',2);
                            
                            node(node_ID).adj_matrix(post1, post2) = 0;
                            node(node_ID).adj_matrix(post2, post1) = 0;
                            node(node_ID).simp_temp{1}(post1).neighb = setdiff(node(node_ID).simp_temp{1}(post1).neighb, neigh_ID2);
                            node(node_ID).simp_temp{1}(post2).neighb = setdiff(node(node_ID).simp_temp{1}(post2).neighb, neigh_ID1);
                            no_triangle = size(node(node_ID).simp_temp{2}, 2);
                            for m = 1 : no_triangle % delete the triangle from its simplices set
                                if ismember(neigh_ID1, node(node_ID).simp_temp{2}(m).vert) && ismember(neigh_ID2, node(node_ID).simp_temp{2}(m).vert)
                                    node(node_ID).simp_temp{2}(m) = [];
                                    break;
                                end
                            end
                            if isempty(node(node_ID).simp_temp{2})
                                node(node_ID).simp_temp(2) = [];
                            end
                            
                            post_temp1 = find(node(neigh_ID1).neighbors_temp == neigh_ID2);
                            node(neigh_ID1).adj_matrix(post_temp1, :) = [];
                            node(neigh_ID1).adj_matrix(:, post_temp1) = [];
                            node(neigh_ID1).neighbors_temp = setdiff(node(neigh_ID1).neighbors_temp, neigh_ID2);
                            node(neigh_ID1).boundary_edge = setdiff(node(neigh_ID1).boundary_edge, neigh_ID2);
                            node(neigh_ID1).boundary_neighbor = setdiff(node(neigh_ID1).boundary_neighbor, neigh_ID2);
                            node(neigh_ID1).simp_temp{1}(post_temp1) = [];
                            post_common = find(node(neigh_ID1).neighbors_temp == node_ID);
                            node(neigh_ID1).simp_temp{1}(post_common).neighb = setdiff(node(neigh_ID1).simp_temp{1}(post_common).neighb, neigh_ID2);
                            no_triangle = size(node(neigh_ID1).simp_temp{2}, 2);
                            for m = 1 : no_triangle % delete the triangle from its simplices set
                                if ismember(neigh_ID2, node(neigh_ID1).simp_temp{2}(m).vert)
                                    node(neigh_ID1).simp_temp{2}(m) = [];
                                    break;
                                end
                            end
                            if isempty(node(neigh_ID1).simp_temp{2})
                                node(neigh_ID1).simp_temp(2) = [];
                            end
                            
                            post_temp2 = find(node(neigh_ID2).neighbors_temp == neigh_ID1);
                            node(neigh_ID2).adj_matrix(post_temp2, :) = [];
                            node(neigh_ID2).adj_matrix(:, post_temp2) = [];
                            node(neigh_ID2).neighbors_temp = setdiff(node(neigh_ID2).neighbors_temp, neigh_ID1);
                            node(neigh_ID2).boundary_edge = setdiff(node(neigh_ID2).boundary_edge, neigh_ID1);
                            node(neigh_ID2).boundary_neighbor = setdiff(node(neigh_ID2).boundary_neighbor, neigh_ID1);
                            node(neigh_ID2).simp_temp{1}(post_temp2) = [];
                            post_common = find(node(neigh_ID2).neighbors_temp == node_ID);
                            node(neigh_ID2).simp_temp{1}(post_common).neighb = setdiff(node(neigh_ID2).simp_temp{1}(post_common).neighb, neigh_ID1);
                            no_triangle = size(node(neigh_ID2).simp_temp{2}, 2);
                            for m = 1 : no_triangle % delete the triangle from its simplices set
                                if ismember(neigh_ID1, node(neigh_ID2).simp_temp{2}(m).vert)
                                    node(neigh_ID2).simp_temp{2}(m) = [];
                                    break;
                                end
                            end
                            if isempty(node(neigh_ID2).simp_temp{2})
                                node(neigh_ID2).simp_temp(2) = [];
                            end
                            
                            if length(common_neigh_set1) == 2 && length(common_neigh_set2) == 2
                                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neigh_ID1, neigh_ID2];
                                node(node_ID).boundary_edge = sort(node(node_ID).boundary_edge);
                                % line([node_x(node_ID),node_x(neigh_ID1)],[node_y(node_ID),node_y(neigh_ID1)], 'Color','r','LineWidth',2);
                                % line([node_x(node_ID),node_x(neigh_ID2)],[node_y(node_ID),node_y(neigh_ID2)], 'Color','r','LineWidth',2);
                                
                                node(neigh_ID1).boundary_edge = [node(neigh_ID1).boundary_edge, node_ID];
                                node(neigh_ID1).boundary_edge = sort(node(neigh_ID1).boundary_edge);
                                
                                node(neigh_ID2).boundary_edge = [node(neigh_ID2).boundary_edge, node_ID];
                                node(neigh_ID2).boundary_edge = sort(node(neigh_ID2).boundary_edge);
                            elseif length(common_neigh_set1) == 2
                                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neigh_ID1];
                                node(node_ID).boundary_edge = sort(node(node_ID).boundary_edge);
                                % line([node_x(node_ID),node_x(neigh_ID1)],[node_y(node_ID),node_y(neigh_ID1)], 'Color','r','LineWidth',2);
                                
                                node(neigh_ID1).boundary_edge = [node(neigh_ID1).boundary_edge, node_ID];
                                node(neigh_ID1).boundary_edge = sort(node(neigh_ID1).boundary_edge);
                            elseif length(common_neigh_set2) == 2
                                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neigh_ID2];
                                node(node_ID).boundary_edge = sort(node(node_ID).boundary_edge);
                                % line([node_x(node_ID),node_x(neigh_ID2)],[node_y(node_ID),node_y(neigh_ID2)], 'Color','r','LineWidth',2);
                                
                                node(neigh_ID2).boundary_edge = [node(neigh_ID2).boundary_edge, node_ID];
                                node(neigh_ID2).boundary_edge = sort(node(neigh_ID2).boundary_edge);
                            end
                        end 
                    end     
                end     
            end
        end
    end
end    

% find the edges associated with only one triangle
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_edge = [];
    if isempty(node(node_ID).neighbors_temp)
        continue;
    else
        no_edge = size(node(node_ID).simp_temp{1}, 2);
    end
    for j = 1 : no_edge
        vert_set = node(node_ID).simp_temp{1}(j).vert;
        neighb_set = node(node_ID).simp_temp{1}(j).neighb;
        neighbor_ID = setdiff(vert_set, node_ID);

        if node(node_ID).fence_flag == 1 && node(neighbor_ID).fence_flag == 1
            if isempty(neighb_set)
                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
            end
        else
            if length(neighb_set) <= 1
                node(node_ID).boundary_edge = [node(node_ID).boundary_edge, neighbor_ID];
                % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
            end
        end
    end       
end

boundary_node_set = [];
% plot boundary edges with red color
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    if ~isempty(node(node_ID).boundary_edge)
        boundary_node_set = [boundary_node_set, node_ID];
        no_boundary_edge = length(node(node_ID).boundary_edge);
        node(node_ID).no_cycle = zeros(1, no_boundary_edge);
        node(node_ID).cycle_seq = [];
        for j = 1 : no_boundary_edge
            neighbor_ID = node(node_ID).boundary_edge(j);
            % set the degree for each boundary edge
            if node(node_ID).fence_flag == 1 && node(neighbor_ID).fence_flag == 1
                node(node_ID).boundary_degree(j) = 1;
            else
                post = find(node(node_ID).neighbors_temp == neighbor_ID);
                if isempty(node(node_ID).simp_temp{1}(post).neighb)
                    node(node_ID).boundary_degree(j) = 2;
                else
                    node(node_ID).boundary_degree(j) = 1;
                end
            end
            
            if neighbor_ID < node_ID
                continue;
            else
                % line([node_x(node_ID),node_x(neighbor_ID)],[node_y(node_ID),node_y(neighbor_ID)], 'Color','r','LineWidth',2);
            end
        end
        
        node_hop2 = [];
        for k = 1 : length(node(node_ID).neighbors_temp)
            neigh_ID = node(node_ID).neighbors_temp(k);
            node_hop2 = union(node_hop2, node(neigh_ID).neighbors_temp);
        end
        node_hop2 =  setdiff(node_hop2, node(node_ID).neighbors_temp);
        node_hop2 = setdiff(node_hop2, node_ID);
        node(node_ID).neigh_hop2 = node_hop2;
    end
end

% set the boundary_flag for each node 
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_flag = 0;
    if ~isempty(node(node_ID).boundary_edge)
        node(node_ID).boundary_flag = 1;
    end
end

% find the boundary neighbors of each left node
for i = 1 : length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).boundary_neighbor = [];
    for j = 1 : length(node(node_ID).neighbors_temp)
        neighbor_ID = node(node_ID).neighbors_temp(j);
        if ~isempty(node(neighbor_ID).boundary_edge)
            node(node_ID).boundary_neighbor = [node(node_ID).boundary_neighbor, neighbor_ID];
        end
    end
end

% *************************************************************************
% find all cycles
% *************************************************************************

cycle = [];
list_initiator = [];
list_init_neighbor = [];

% initialization
for i = 1: length(node_left_seq)
    node_ID = node_left_seq(i);
    node(node_ID).init_flag = 0;
    node(node_ID).init = [];
    node(node_ID).father = [];
    node(node_ID).init_temp = [];
    node(node_ID).father_temp = [];
end

% choose the common edge of two neighbouring holes as initiators
for i = 1 : length(boundary_node_set)
    node_ID = boundary_node_set(i);
    no_boundary_edge = length(node(node_ID).boundary_edge);
    % do not choose fence nodes
    if node(node_ID).fence_flag == 1
        continue;
    end
    
    if no_boundary_edge == 3
        for j = 1 : 3
            neigh_ID = node(node_ID).boundary_edge(j);
            common_set = intersect(node(node_ID).neighbors_temp, node(neigh_ID).neighbors_temp);
            if isempty(common_set)
                if length(node(neigh_ID).boundary_edge) ~= 3
                    node(node_ID).init_flag = 1;
                    list_initiator = [list_initiator, node_ID];
                    list_init_neighbor = [list_init_neighbor, neigh_ID];
                    break;
                else
                    if node_ID < neigh_ID
                        node(node_ID).init_flag = 1;
                        list_initiator = [list_initiator, node_ID];
                        list_init_neighbor = [list_init_neighbor, neigh_ID];
                        break;
                    end
                end
            end
        end
    elseif no_boundary_edge == 2
        neigh_ID1 = node(node_ID).boundary_edge(1);
        neigh_ID2 = node(node_ID).boundary_edge(2);

        no_boundary_edge1 = length(node(neigh_ID1).boundary_edge);
        no_boundary_edge2 = length(node(neigh_ID2).boundary_edge);
        
        common_set1 = intersect(node(node_ID).neighbors_temp, node(neigh_ID1).neighbors_temp);
        common_set2 = intersect(node(node_ID).neighbors_temp, node(neigh_ID2).neighbors_temp);
        
        if isempty(common_set1) && isempty(common_set2)
            if no_boundary_edge1 ~= 2 && no_boundary_edge2 ~= 2
                node(node_ID).init_flag = 1;
                list_initiator = [list_initiator, node_ID];
                list_init_neighbor = [list_init_neighbor, neigh_ID1];
            elseif no_boundary_edge1 == 2 && no_boundary_edge2 == 2
                if node_ID < min(neigh_ID1, neigh_ID2)
                    node(node_ID).init_flag = 1;
                    list_initiator = [list_initiator, node_ID];
                    list_init_neighbor = [list_init_neighbor, neigh_ID1];
                end
            elseif no_boundary_edge1 == 2 && node_ID < neigh_ID1
                node(node_ID).init_flag = 1;
                list_initiator = [list_initiator, node_ID];
                list_init_neighbor = [list_init_neighbor, neigh_ID1];
            elseif no_boundary_edge2 == 2 && node_ID < neigh_ID2
                node(node_ID).init_flag = 1;
                list_initiator = [list_initiator, node_ID];
                list_init_neighbor = [list_init_neighbor, neigh_ID1];
            end
        end
    end
end
    
% If no common edges, choose normal edges
if isempty(list_initiator)
    for i = 1 : length(boundary_node_set)
        node_ID = boundary_node_set(i);
        no_boundary_edge = length(node(node_ID).boundary_edge);

        if no_boundary_edge == 2
            neigh_ID1 = node(node_ID).boundary_edge(1);
            neigh_ID2 = node(node_ID).boundary_edge(2);

            no_boundary_edge1 = length(node(neigh_ID1).boundary_edge);
            no_boundary_edge2 = length(node(neigh_ID2).boundary_edge);

            if no_boundary_edge1 ~= 2 && no_boundary_edge2 ~= 2
                node(node_ID).init_flag = 1;
                list_initiator = [list_initiator, node_ID];
                list_init_neighbor = [list_init_neighbor, neigh_ID1];
            elseif no_boundary_edge1 == 2 && no_boundary_edge2 == 2
                if node_ID < min(neigh_ID1, neigh_ID2)
                    node(node_ID).init_flag = 1;
                    list_initiator = [list_initiator, node_ID];
                    list_init_neighbor = [list_init_neighbor, neigh_ID1];
                end
            elseif no_boundary_edge1 == 2
                if no_boundary_edge2 == 1
                    node(node_ID).init_flag = 1;
                    list_initiator = [list_initiator, node_ID];
                    list_init_neighbor = [list_init_neighbor, neigh_ID2];
                elseif node_ID < neigh_ID1
                    node(node_ID).init_flag = 1;
                    list_initiator = [list_initiator, node_ID];
                    list_init_neighbor = [list_init_neighbor, neigh_ID1];
                end
            elseif no_boundary_edge2 == 2
                if no_boundary_edge1 == 1 || node_ID < neigh_ID2
                    node(node_ID).init_flag = 1;
                    list_initiator = [list_initiator, node_ID];
                    list_init_neighbor = [list_init_neighbor, neigh_ID1];
                end
            end
        end
    end
end

while(1)
    list_broadcastor = [];

    stop_node = [];
    stop_init = [];
    for no = 1 : length(list_initiator)
        node_ID = list_initiator(no);

        no_init = size(node(node_ID).init, 2);
        node(node_ID).init(no_init+1) = node_ID;
        node(node_ID).father(no_init+1) = 0;

        neigh_ID1 = list_init_neighbor(no);
        no_init_temp = size(node(neigh_ID1).init_temp, 2);
        node(neigh_ID1).init_temp(no_init_temp+1) = node_ID;
        node(neigh_ID1).father_temp(no_init_temp+1) = node_ID;

        list_broadcastor = [list_broadcastor, neigh_ID1];

        % plot(node_coor(1,neigh_ID1), node_coor(2,neigh_ID1), 'dg'); hold on;

        while(~isempty(list_broadcastor))
            no_broadcastor = length(list_broadcastor);

            for i = 1 : no_broadcastor
                broadcastor_ID = list_broadcastor(i);

                no_init_temp = size(node(broadcastor_ID).init_temp, 2);

                for j = 1: no_init_temp
                    no_init = size(node(broadcastor_ID).init, 2);

                    init_ID_temp = node(broadcastor_ID).init_temp(j);
                    father_ID_temp = node(broadcastor_ID).father_temp(j);

                    if isAinB(init_ID_temp, node(broadcastor_ID).init)
                        continue;
                    elseif init_ID_temp ~= father_ID_temp && isAinB(init_ID_temp, node(broadcastor_ID).boundary_edge) %find a cycle
                        node(broadcastor_ID).init(no_init+1) = init_ID_temp;
                        node(broadcastor_ID).father(no_init+1) = father_ID_temp;

                        stop_node = [stop_node, broadcastor_ID];
                        stop_init = [stop_init, init_ID_temp];
                    else
                        node(broadcastor_ID).init(no_init+1) = init_ID_temp;
                        node(broadcastor_ID).father(no_init+1) = father_ID_temp;

                        if node(father_ID_temp).boundary_flag == 0 && node(broadcastor_ID).boundary_flag == 0
                            neighbor_ID_set = [];
                        else
                            if ~isempty(node(broadcastor_ID).boundary_edge)
                                neighbor_ID_set = setdiff(node(broadcastor_ID).boundary_edge, father_ID_temp);
                                if isempty(neighbor_ID_set) 
                                    neighbor_ID_set = setdiff(node(broadcastor_ID).neighbors_temp, father_ID_temp);
                                end
                            else
                                neighbor_ID_set = setdiff(node(broadcastor_ID).boundary_neighbor, father_ID_temp);
                                if isempty(neighbor_ID_set)
                                    neighbor_ID_set = setdiff(node(broadcastor_ID).neighbors_temp, father_ID_temp);
                                end
                            end

                            no_left_neighbor = length(neighbor_ID_set);

                            for no_neighbor = 1: no_left_neighbor
                                neighbor_ID = neighbor_ID_set(no_neighbor);

                                node(neighbor_ID).init_temp = [node(neighbor_ID).init_temp, init_ID_temp];
                                node(neighbor_ID).father_temp = [node(neighbor_ID).father_temp, broadcastor_ID];
                            end
                        end
                    end
                end

                node(broadcastor_ID).init_temp = [];
                node(broadcastor_ID).father_temp = [];
            end

            list_broadcastor = [];
            for i = 1: length(node_left_seq)
                if ~isempty(node(node_left_seq(i)).init_temp)
                    list_broadcastor = [list_broadcastor, node_left_seq(i)];
                end
            end
        end
    end

    % construct the cycles
    cycle_seq = [];
    for i = 1 : length(stop_node)
        node_temp = stop_node(i);
        init_ID = stop_init(i);

        no_cycle = size(cycle_seq, 2);

        cycle_seq{no_cycle+1} = [node_temp];

        while(1)
            init_post = find(node(node_temp).init == init_ID);
            father_temp = node(node_temp).father(init_post);

            if ~father_temp
                break;
            else
                cycle_temp = [father_temp, cycle_seq{no_cycle+1}];
                cycle_seq{no_cycle+1} = cycle_temp;
                node_temp = father_temp;
            end
        end
    end

    % for false cycles, delete them
    delete_cycle = [];
    for i = 1 : size(cycle_seq, 2)
        cycle_temp = cycle_seq{i};
        
        if length(cycle_temp) == 3
            delete_cycle = [delete_cycle, i];
            continue;
        end

        if length(cycle_temp) <= 5
            common_set = node(cycle_temp(1)).neighbors_temp;
            for j = 2 : length(cycle_temp)
                common_set = intersect(common_set, node(cycle_temp(j)).neighbors_temp);
            end
            if ~isempty(common_set)
                delete_cycle = [delete_cycle, i];
                continue;
            else
                for j = 1 : length(cycle_temp)
                    flag_found = 0;
                    common_neigh_set = intersect(node(cycle_temp(j)).neighbors_temp, node(cycle_temp(mod(j, length(cycle_temp))+1)).neighbors_temp);
                    if j == 1
                        common_neigh_set = intersect(common_neigh_set, node(cycle_temp(length(cycle_temp))).neighbors_temp);
                    else
                        common_neigh_set = intersect(common_neigh_set, node(cycle_temp(j-1)).neighbors_temp);
                    end
                    if ~isempty(common_neigh_set)
                        for k = 1 : length(common_neigh_set)
                            node_temp_ID = common_neigh_set(k);
                            common_set = node(node_temp_ID).neighbors_temp;
                            cycle_temp_temp = cycle_temp;
                            cycle_temp_temp(j) = [];
                            for n = 1 : length(cycle_temp_temp)
                                common_set = intersect(common_set, node(cycle_temp_temp(n)).neighbors_temp);
                            end
                            if ~isempty(common_set)
                                delete_cycle = [delete_cycle, i];
                                flag_found = 1;
                                break;
                            end
                        end
                    end
                    if flag_found == 1
                        break;
                    end
                end
                
                if flag_found == 0        
                    %   check for the first node
                    j = 1;
                    neighbor_set_hop1 = [cycle_temp(length(cycle_temp)), cycle_temp(1), cycle_temp(2)];
                    other_node_hop1 = setdiff(cycle_temp, neighbor_set_hop1);
                    intersect_hop1 = intersect(other_node_hop1, node(cycle_temp(j)).boundary_neighbor);
                    if ~isempty(intersect_hop1)
                        delete_cycle = [delete_cycle, i];
                        continue;
                    end

                    %   check for other nodes    
                    for j = 2 : length(cycle_temp)
                        neighbor_set_hop1 = [cycle_temp(j-1), cycle_temp(j), cycle_temp(mod(j, length(cycle_temp))+1)];
                        other_node_hop1 = setdiff(cycle_temp, neighbor_set_hop1);
                        intersect_hop1 = intersect(other_node_hop1, node(cycle_temp(j)).boundary_neighbor);
                        if ~isempty(intersect_hop1)
                            delete_cycle = [delete_cycle, i];
                            break;
                        end
                    end
                end
            end
        end
    end
    
    delete_cycle = sort(delete_cycle);
    if ~isempty(delete_cycle)
        for no_seq = 1: length(delete_cycle)
             no_seq_delete = delete_cycle(no_seq);
             cycle_seq(no_seq_delete - no_seq + 1) = [];
        end
    end
    
    no_cycle_temp = size(cycle_seq, 2);

    % for same cycles, only keep  one 
    delete_cycle = [];
    for i = 1 : no_cycle_temp-1
%         if isAinB(i, delete_cycle)
%             continue;
%         end
        
        for j = i+1 : no_cycle_temp
            cycle1 = cycle_seq{i};
            cycle2 = cycle_seq{j};

            cycle_inter = intersect(cycle1, cycle2);

            if length(cycle_inter) == length(cycle1) && length(cycle_inter) == length(cycle2)
                if ~isAinB(j, delete_cycle)
                    delete_cycle = [delete_cycle, j];
                end
            elseif length(cycle_inter) >= 4 && cycle1(1) == cycle2(1)
                if length(cycle1) > length(cycle2)
                    if ~isAinB(i, delete_cycle)
                        delete_cycle = [delete_cycle, i];
                        break;
                    end
                else
                    if ~isAinB(j, delete_cycle)
                        delete_cycle = [delete_cycle, j];
                    end
                end
            end
        end
    end

    cycle_seq_delete = cycle_seq;
    delete_cycle = sort(delete_cycle);
    if ~isempty(delete_cycle)
        for no_seq = 1: length(delete_cycle)
             no_seq_delete = delete_cycle(no_seq);
             cycle_seq_delete(no_seq_delete - no_seq + 1) = [];
        end
    end

    no_cycle_left = size(cycle_seq_delete, 2);
    
    for j = 1 : no_cycle_left
        cycle_temp = cycle_seq_delete{j};
        for k = 1 : length(cycle_temp)
            node_ID = cycle_temp(k);
            if isempty(node(node_ID).boundary_edge)
                continue;
            end
            no_node_cycle = size(node(node_ID).cycle_seq, 2);
            node(node_ID).cycle_seq{no_node_cycle+1} = cycle_temp;
            if k == 1
                node_before_ID = cycle_temp(length(cycle_temp));
            else
                node_before_ID = cycle_temp(k-1);
            end
            if k == length(cycle_temp)
                node_after_ID = cycle_temp(1);
            else
                node_after_ID = cycle_temp(k+1);
            end
            
            post_before = find(node(node_ID).boundary_edge == node_before_ID);
            post_after = find(node(node_ID).boundary_edge == node_after_ID);
            if ~isempty(post_before)
                node(node_ID).no_cycle(post_before) = node(node_ID).no_cycle(post_before) + 1;
            end
            if ~isempty(post_after)
                node(node_ID).no_cycle(post_after) = node(node_ID).no_cycle(post_after) + 1;
            end
        end
    end
    
    delete_cycle = [];
    for i = 1 : no_cycle_left - 1
        for j = i+1 : no_cycle_left
            cycle1 = cycle_seq_delete{i};
            cycle2 = cycle_seq_delete{j};

            cycle_inter = intersect(cycle1, cycle2);
            if length(cycle_inter)/length(cycle1) >= 1/2 && length(cycle_inter)/length(cycle2) >= 1/2
                no_false1 = 0;
                no_false2 = 0;

                for k = 1 : length(cycle1)
                    node_ID = cycle1(k);
                    if isempty(node(node_ID).boundary_edge)
                        continue;
                    end
                    node_after_ID = cycle1(mod(k, length(cycle1)) + 1);                    

                    post_after = find(node(node_ID).boundary_edge == node_after_ID);
                    if ~isempty(post_after) && node(node_ID).boundary_degree(post_after) < node(node_ID).no_cycle(post_after)
                        no_false1 = no_false1 + 1;
                    end
                end
                
                for k = 1 : length(cycle2)
                    node_ID = cycle2(k);
                    if isempty(node(node_ID).boundary_edge)
                        continue;
                    end
                    node_after_ID = cycle2(mod(k, length(cycle2)) + 1);                    

                    post_after = find(node(node_ID).boundary_edge == node_after_ID);
                    if ~isempty(post_after) && node(node_ID).boundary_degree(post_after) < node(node_ID).no_cycle(post_after)
                        no_false2 = no_false2 + 1;
                    end
                end
                
                if no_false1 > no_false2
                    if ~isAinB(i, delete_cycle)
                        delete_cycle = [delete_cycle, i];
                        break;
                    end
                else
                    if ~isAinB(j, delete_cycle)
                        delete_cycle = [delete_cycle, j];
                    end
                end
            end
        end
    end
    
    delete_cycle = sort(delete_cycle);
    if ~isempty(delete_cycle)
        for no_seq = 1: length(delete_cycle)
             no_seq_delete = delete_cycle(no_seq);
             cycle_temp = cycle_seq_delete{no_seq_delete - no_seq + 1};
             cycle_seq_delete(no_seq_delete - no_seq + 1) = [];
             for j = 1 : length(cycle_temp)
                node_ID = cycle_temp(j);
                if isempty(node(node_ID).boundary_edge)
                    continue;
                end

                index1 = j - 1;
                if index1 == 0
                    index1 = length(cycle_temp);
                end

                node_ID1 = cycle_temp(index1);
                post1 = find(node(node_ID).boundary_edge == node_ID1);

                if ~isempty(post1)
                    node(node_ID).no_cycle(post1) = node(node_ID).no_cycle(post1) - 1;
                end

                index2 = mod(j, length(cycle_temp)) + 1;
                node_ID2 = cycle_temp(index2);
                post2 = find(node(node_ID).boundary_edge == node_ID2);

                if ~isempty(post2)
                    node(node_ID).no_cycle(post2) = node(node_ID).no_cycle(post2) - 1;  
                end
            end
        end
    end
    
    no_cycle_left = size(cycle_seq_delete, 2);
    no_cycle = size(cycle, 2);
    
    if no_cycle_left == 0
        break;
    end

    for j = 1 : no_cycle_left
        cycle_add = cycle_seq_delete{j};
        cycle{no_cycle + j} = cycle_add;
    end

    % for each node, delete its neighbors or change the degree of its neighbors in the cycles found
    for i = 1 : size(cycle_seq_delete, 2)
        cycle_temp = cycle_seq_delete{i};

        for j = 1 : length(cycle_temp)
            node_ID = cycle_temp(j);
            if isempty(node(node_ID).boundary_edge)
                continue;
            end

            index1 = j - 1;
            if index1 == 0
                index1 = length(cycle_temp);
            end

            node_ID1 = cycle_temp(index1);
            post1 = find(node(node_ID).boundary_edge == node_ID1);

            if ~isempty(post1)
                if node(node_ID).boundary_degree(post1) == 1
                    node(node_ID).boundary_edge(post1) = [];
                    node(node_ID).boundary_degree(post1) = [];
                    node(node_ID).boundary_neighbor = setdiff(node(node_ID).boundary_neighbor, node_ID1);
                else
                    node(node_ID).boundary_degree(post1) = 1;
                end
                if node(node_ID).no_cycle(post1) ~= 0
                    node(node_ID).no_cycle(post1) = node(node_ID).no_cycle(post1) - 1;
                end
            end

            index2 = mod(j, length(cycle_temp)) + 1;
            node_ID2 = cycle_temp(index2);
            post2 = find(node(node_ID).boundary_edge == node_ID2);

            if ~isempty(post2)
                if node(node_ID).boundary_degree(post2) == 1
                    node(node_ID).boundary_edge(post2) = [];
                    node(node_ID).boundary_degree(post2) = [];
                    node(node_ID).boundary_neighbor = setdiff(node(node_ID).boundary_neighbor, node_ID2);
                else
                    node(node_ID).boundary_degree(post2) = 1;
                end
                if node(node_ID).no_cycle(post2) ~= 0
                    node(node_ID).no_cycle(post2) = node(node_ID).no_cycle(post2) - 1;
                end   
            end
        end
    end
    
    list_initiator = [];
    list_init_neighbor = [];

    % initialization
    for i = 1: length(node_left_seq)
        node_ID = node_left_seq(i);
        node(node_ID).init_flag = 0;
        node(node_ID).init = [];
        node(node_ID).father = [];
        node(node_ID).init_temp = [];
        node(node_ID).father_temp = [];
    end

    % choose initiators
    for i = 1 : length(boundary_node_set)  
        node_ID = boundary_node_set(i);
        no_boundary_edge = length(node(node_ID).boundary_edge);

        if no_boundary_edge == 2
            neigh_ID1 = node(node_ID).boundary_edge(1);
            neigh_ID2 = node(node_ID).boundary_edge(2);

            no_boundary_edge1 = length(node(neigh_ID1).boundary_edge);
            no_boundary_edge2 = length(node(neigh_ID2).boundary_edge);

            if no_boundary_edge1 ~= 2 && no_boundary_edge2 ~= 2
                node(node_ID).init_flag = 1;
                list_initiator = [list_initiator, node_ID];
                list_init_neighbor = [list_init_neighbor, neigh_ID1];
            elseif no_boundary_edge1 == 2 && no_boundary_edge2 == 2
                if node_ID < min(neigh_ID1, neigh_ID2)
                    node(node_ID).init_flag = 1;
                    list_initiator = [list_initiator, node_ID];
                    list_init_neighbor = [list_init_neighbor, neigh_ID1];
                end
            elseif no_boundary_edge1 == 2
                if no_boundary_edge2 == 1
                    node(node_ID).init_flag = 1;
                    list_initiator = [list_initiator, node_ID];
                    list_init_neighbor = [list_init_neighbor, neigh_ID2];
                elseif node_ID < neigh_ID1
                    node(node_ID).init_flag = 1;
                    list_initiator = [list_initiator, node_ID];
                    list_init_neighbor = [list_init_neighbor, neigh_ID1];
                end
            elseif no_boundary_edge2 == 2
                if no_boundary_edge1 == 1 || node_ID < neigh_ID2
                    node(node_ID).init_flag = 1;
                    list_initiator = [list_initiator, node_ID];
                    list_init_neighbor = [list_init_neighbor, neigh_ID1];
                end
            end
        end
    end

    if isempty(list_initiator)
        break;
    end
    
%     if flag_hex == 0
%         flag_special = 1;
%         flag_hex = 1;
%     else
%         flag_special = 0;
%     end
    
end

% for same cycles, only keep one
no_cycle = size(cycle, 2);
delete_cycle = [];
for i = 1 : no_cycle-1
    for j = i+1 : no_cycle
        cycle1 = sort(cycle{i});
        cycle2 = sort(cycle{j});

        cycle_inter = intersect(cycle1, cycle2);

        if length(cycle_inter) == length(cycle1) && length(cycle_inter) == length(cycle2)
            if ~isAinB(j, delete_cycle)
                delete_cycle = [delete_cycle, j];
            end
        end
    end
end

delete_cycle = sort(delete_cycle);
if ~isempty(delete_cycle)
    for no_seq = 1: length(delete_cycle)
         no_seq_delete = delete_cycle(no_seq);
         cycle(no_seq_delete - no_seq + 1) = [];
    end
end

% *************************************************************************
% *********************** minimization of cycles **************************
% *************************************************************************

% to find 2-hop nodes
for i = 1: length(node_coor)
    node(i).hop2 = [];
    if ~isempty(node(i).neighbors)
        for j = 1: length(node(i).neighbors)
            node_hop1 = node(i).neighbors(j);
            node(i).hop2 = union(node(i).hop2, node(node_hop1).neighbors);
        end
        node(i).hop2 = setdiff(node(i).hop2, node(i).neighbors);
        node(i).hop2 = setdiff(node(i).hop2, i);
    end
end

if ~isempty(cycle)  
    cycle_min = [];
    cycle_temp_set = cycle;
    i = 1;
    while(1)
        cycle_temp = cycle_temp_set{i};
        for n = 1: length(node_coor)
            node(n).visited = 0;
        end
        if length(cycle_temp) == 4
            flag = cycle_4_check(cycle_temp, node);
            if flag == 1
                no_min_cycle = size(cycle_min, 2);
                cycle_min{no_min_cycle + 1} = cycle_temp;
            end
        else
            % 1-hop check
            node_ID = cycle_temp(1);
            while(1)
                len = length(cycle_temp);
                if len >= 5
                    half_len = floor(len/2);
                    if mod(len, 2) == 0
                        half_len1 = floor(len/2) - 1;
                        half_len2 = floor(len/2);
                    else
                        half_len1 = floor(len/2);
                        half_len2 = floor(len/2);
                    end

                    sub_cycle1 = [node_ID];
                    sub_cycle2 = [node_ID];
                    j = find(cycle_temp == node_ID);
                    for ind = 1 : half_len1
                        node_ind = j - ind;
                        if node_ind == 0
                            node_ind = len;
                        else
                            node_ind = mod((j - ind), len);
                        end
                        sub_cycle1 = [sub_cycle1, cycle_temp(node_ind)];
                    end

                    for ind = 1 : half_len2
                        node_ind = mod((j+ind-1), len) + 1;
                        sub_cycle2 = [sub_cycle2, cycle_temp(node_ind)];
                    end

                    sub_sub_cycle1 = [];
                    sub_sub_cycle2 = [];
                    sub_sub_cycle3 = [];
                    flag1 = 0;
                    if half_len1 >= 2 % if the length of sub_cycle1 is larger or equal to 3
                        for pt1 = 1 : half_len1 - 1
                            node_temp_ID = sub_cycle1(half_len1 - pt1 + 2);
                            if ismember(node_temp_ID, node(node_ID).neighbors)
                                flag1 = 1;
                                break;
                            end
                        end
                        if flag1 == 1
                            sub_sub_cycle1 = sub_cycle1(1 : half_len1 - pt1 + 2);
                            sub_sub_cycle2_temp1 = [sub_cycle1(1), sub_cycle1(half_len1 - pt1 + 2 : half_len1 + 1)];
                        else
                            sub_sub_cycle2_temp1 = sub_cycle1;
                        end
                    else
                        sub_sub_cycle2_temp1 = sub_cycle1;
                    end

                    flag2 = 0;
                    if half_len2 >= 2 % if the length of sub_cycle1 is larger or equal to 3
                        for pt2 = 1 : half_len2 - 1
                            node_temp_ID = sub_cycle2(half_len2 - pt2 + 2);
                            if ismember(node_temp_ID, node(node_ID).neighbors)
                                flag2 = 1;
                                break;
                            end
                        end
                        if flag2 == 1
                            sub_sub_cycle3 = sub_cycle2(1 : half_len2 - pt2 + 2);
                            sub_sub_cycle2_temp2 = [sub_cycle2(1), sub_cycle2(half_len2 - pt2 + 2 : half_len2 + 1)];
                        else
                            sub_sub_cycle2_temp2 = sub_cycle2;
                        end
                    else
                        sub_sub_cycle2_temp2 = sub_cycle2;
                    end

                    sub_sub_cycle2 = [sub_sub_cycle2_temp2, fliplr(sub_sub_cycle2_temp1)];
                    sub_sub_cycle2(length(sub_sub_cycle2)) = [];

                    if length(sub_sub_cycle1) <= 3
                        sub_sub_cycle1 = [];
                    elseif length(sub_sub_cycle1) == 4
                        flag_hole = cycle_4_check(sub_sub_cycle1, node);
                        if flag_hole == 0
                            sub_sub_cycle1 = [];
                        else
                            no_min_cycle = size(cycle_min, 2);
                            cycle_min{no_min_cycle + 1} = sub_sub_cycle1;
                        end
                    else
                        no_cycle = size(cycle_temp_set, 2);
                        cycle_temp_set{no_cycle+1} = sub_sub_cycle1;
                    end

                    if length(sub_sub_cycle3) <= 3
                        sub_sub_cycle3 = [];
                    elseif length(sub_sub_cycle3) == 4
                        flag_hole = cycle_4_check(sub_sub_cycle3, node);
                        if flag_hole == 0
                            sub_sub_cycle3 = [];
                        else
                            no_min_cycle = size(cycle_min, 2);
                            cycle_min{no_min_cycle + 1} = sub_sub_cycle3;
                        end
                    else
                        no_cycle = size(cycle_temp_set, 2);
                        cycle_temp_set{no_cycle+1} = sub_sub_cycle3;
                    end

                    if length(sub_sub_cycle2) <= 3
                        sub_sub_cycle2 = [];
                    elseif length(sub_sub_cycle2) == 4
                        flag_hole = cycle_4_check(sub_sub_cycle2, node);
                        if flag_hole == 0
                            sub_sub_cycle2 = [];
                        end
                    end

                    if isempty(sub_sub_cycle2)
                        cycle_temp = [];
                        break;
                    end

                    % 2-hop check for sub_sub_cycle2
                    len2 = length(sub_sub_cycle2);
                    hop2_sub_cycle1 = [];
                    hop2_sub_cycle2 = [];
                    if len2 >= 6
                        for ind_temp = 1 : len2 - 5
                            if mod(ind_temp, 2) == 1
                                ind = floor(len2/2) + 1 - floor(ind_temp/2);
                            else
                                ind = floor(len2/2) + 1 + floor(ind_temp/2);
                            end

                            test_node_ID = sub_sub_cycle2(ind);
                            if ismember(test_node_ID, node(node_ID).hop2)
                                node_after_ID = sub_sub_cycle2(2);
                                node_before_ID = sub_sub_cycle2(len2);
                                if ismember(test_node_ID, node(node_after_ID).neighbors)
                                    hop2_sub_cycle1 = [node_ID, node_after_ID, sub_sub_cycle2(ind: len2)];
                                    hop2_sub_cycle2 = sub_sub_cycle2(2: ind);
                                elseif ismember(test_node_ID, node(node_before_ID).neighbors)
                                    hop2_sub_cycle1 = [sub_sub_cycle2(1: ind), node_before_ID];
                                    hop2_sub_cycle2 = sub_sub_cycle2(ind: len2);
                                else
                                    for no_neigh = 1 : length(node(node_ID).neighbors)
                                        neigh_ID = node(node_ID).neighbors(no_neigh);
                                        if ismember(test_node_ID, node(neigh_ID).neighbors)
                                            break;
                                        end
                                    end
                                    hop2_sub_cycle1 = [sub_sub_cycle2(1: ind), neigh_ID];
                                    hop2_sub_cycle2 = [node_ID, neigh_ID, sub_sub_cycle2(ind: len2)];
                                end
                                break;
                            end
                        end
                    end

                    if isempty(hop2_sub_cycle1) || isempty(hop2_sub_cycle2)
                        cycle_temp = sub_sub_cycle2;
                        node(node_ID).visited = 1;
                        node_ID = cycle_temp(2);
                        if node(node_ID).visited == 1
                            break;
                        end
                    else
                        if length(hop2_sub_cycle2) <= 3
                            hop2_sub_cycle2 = [];
                        elseif length(hop2_sub_cycle2) == 4
                            flag_hole = cycle_4_check(hop2_sub_cycle2, node);
                            if flag_hole == 0
                                hop2_sub_cycle2 = [];
                            else
                                no_min_cycle = size(cycle_min, 2);
                                cycle_min{no_min_cycle + 1} = hop2_sub_cycle2;
                            end
                        else
                            no_cycle = size(cycle_temp_set, 2);
                            cycle_temp_set{no_cycle+1} = hop2_sub_cycle2;
                        end

                        cycle_temp = hop2_sub_cycle1;
                        node(node_ID).visited = 1;
                        node_ID = cycle_temp(2);
                        if node(node_ID).visited == 1
                            break;
                        end
                    end
                else
                    break;
                end
            end
            
            if length(cycle_temp) >= 5
                sub_cycle_temp_temp1 = [];
                sub_cycle_temp_temp2 = [];
                for j = 1 : length(cycle_temp)
                    flag_shorter = 0;
                    node_ID = cycle_temp(j);
                    node_ID_after = cycle_temp(mod(j, length(cycle_temp))+1);
                    if j == 1
                        node_ID_before = cycle_temp(length(cycle_temp));
                    else
                        node_ID_before = cycle_temp(j-1);
                    end

                    common_neigh_set = intersect(node(node_ID).neighbors, node(node_ID_after).neighbors);
                    common_neigh_set = intersect(common_neigh_set, node(node_ID_before).neighbors);
                    if ~isempty(common_neigh_set)
                        for k = 1 : length(common_neigh_set)
                            node_test_ID = common_neigh_set(k);
                            cycle_temp_temp = cycle_temp;
                            cycle_temp_temp(j) = node_test_ID;
                            cycle_temp_temp = circshift(cycle_temp_temp, [1, -(j-1)]); % to make node_test_ID the first node of the cycle
                            % 1-hop check
                            for n = 3: length(cycle_temp_temp) - 1
                                node_temp_ID = cycle_temp_temp(n);
                                if ismember(node_temp_ID, node(node_test_ID).neighbors)
                                    sub_cycle_temp_temp1 = cycle_temp_temp(1: n);
                                    sub_cycle_temp_temp2 = [node_test_ID, cycle_temp_temp(n: length(cycle_temp_temp))];
                                    flag_shorter = 1;
                                    break;
                                end
                            end
                            
                            % 2-hop check
                            if flag_shorter == 1
                                break;
                            else
                                if length(cycle_temp_temp) == 5
                                    node_ID_2nd = cycle_temp_temp(2);
                                    node_ID_3rd = cycle_temp_temp(3);
                                    node_ID_4th = cycle_temp_temp(4);
                                    node_ID_5th = cycle_temp_temp(5);
                                    
                                    common_set1 = intersect(node(node_ID_3rd).neighbors, node(node_ID_4th).neighbors);
                                    common_set2 = intersect(node(node_test_ID).neighbors, node(node_ID_2nd).neighbors);
                                    common_set = intersect(common_set1, common_set2);
                                    if ~isempty(common_set)
                                        node_shorter_ID = common_set(1);
                                        sub_cycle_temp_temp1 = [node_test_ID, node_shorter_ID, node_ID_4th, node_ID_5th];
                                        flag_shorter = 1;
                                        break;
                                    else
                                        common_set2 = intersect(node(node_test_ID).neighbors, node(node_ID_5th).neighbors);
                                        common_set_another = intersect(common_set1, common_set2);
                                        if ~isempty(common_set_another)
                                            node_shorter_ID = common_set_another(1);
                                            sub_cycle_temp_temp1 = [cycle_temp_temp(1:3), node_shorter_ID];
                                            flag_shorter = 1;
                                            break;
                                        end
                                    end
                                else
                                    len2 = length(cycle_temp_temp);
                                    for ind_temp = 1 : len2 - 5
                                        if mod(ind_temp, 2) == 1
                                            ind = floor(len2/2) + 1 - floor(ind_temp/2);
                                        else
                                            ind = floor(len2/2) + 1 + floor(ind_temp/2);
                                        end

                                        node_shorter_ID = cycle_temp_temp(ind);
                                        if ismember(node_shorter_ID, node(node_test_ID).hop2)
                                            node_after_ID = cycle_temp_temp(2);
                                            node_before_ID = cycle_temp_temp(len2);
                                            if ismember(node_shorter_ID, node(node_after_ID).neighbors)
                                                sub_cycle_temp_temp1 = [node_test_ID, node_after_ID, cycle_temp_temp(ind: len2)];
                                                sub_cycle_temp_temp2 = cycle_temp_temp(2: ind);
                                            elseif ismember(node_shorter_ID, node(node_before_ID).neighbors)
                                                sub_cycle_temp_temp1 = [cycle_temp_temp(1: ind), node_before_ID];
                                                sub_cycle_temp_temp2 = cycle_temp_temp(ind: len2);
                                            else
                                                for no_neigh = 1 : length(node(node_test_ID).neighbors)
                                                    neigh_ID = node(node_test_ID).neighbors(no_neigh);
                                                    if ismember(node_shorter_ID, node(neigh_ID).neighbors)
                                                        break;
                                                    end
                                                end
                                                sub_cycle_temp_temp1 = [cycle_temp_temp(1: ind), neigh_ID];
                                                sub_cycle_temp_temp2 = [node_test_ID, neigh_ID, cycle_temp_temp(ind: len2)];
                                            end
                                            flag_shorter = 1;
                                            break;
                                        end
                                    end
                                    if flag_shorter == 1
                                        break;
                                    end
                                end
                            end
                        end
                    end
                    if flag_shorter == 1
                        break;
                    end
                end
            
                if ~isempty(sub_cycle_temp_temp1) || ~isempty(sub_cycle_temp_temp2)
                    cycle_temp = [];
                    if length(sub_cycle_temp_temp1) <= 3
                        sub_cycle_temp_temp1 = [];
                    else
                        no_cycle = size(cycle_temp_set, 2);
                        cycle_temp_set{no_cycle+1} = sub_cycle_temp_temp1;
                    end

                    if length(sub_cycle_temp_temp2) <= 3
                        sub_cycle_temp_temp2 = [];
                    else
                        no_cycle = size(cycle_temp_set, 2);
                        cycle_temp_set{no_cycle+1} = sub_cycle_temp_temp2;
                    end
                end
            end
                       
            if ~isempty(cycle_temp)
                no_min_cycle = size(cycle_min, 2);
                cycle_min{no_min_cycle + 1} = cycle_temp;
            end    
        end

        if i == size(cycle_temp_set, 2)
            break;
        else
            i = i + 1;
        end
    end
    
    % try to replace fence flag
    no_cycle = size(cycle_min, 2);
    for i = 1 : no_cycle
        cycle_temp = cycle_min{i};
        for j = 1 : length(cycle_temp)
            node_ID = cycle_temp(j);
            if node(node_ID).fence_flag == 1
                node_ID_after = cycle_temp(mod(j, length(cycle_temp))+1);
                if j == 1
                    node_ID_before = cycle_temp(length(cycle_temp));
                else
                    node_ID_before = cycle_temp(j-1);
                end

                common_neigh_set = intersect(node(node_ID).neighbors, node(node_ID_after).neighbors);
                common_neigh_set = intersect(common_neigh_set, node(node_ID_before).neighbors);
                if ~isempty(common_neigh_set)
                    node_ID_new = common_neigh_set(1);
                    cycle_temp(j) = node_ID_new;
                end
            end
        end
        cycle_min{i} = cycle_temp;
    end

    % for same cycles, only keep  one 
    delete_cycle = [];
    no_cycle_temp = size(cycle_min, 2);
    for i = 1 : no_cycle_temp-1
        for j = i+1 : no_cycle_temp
            cycle1 = sort(cycle_min{i});
            cycle2 = sort(cycle_min{j});

            cycle_inter = intersect(cycle1, cycle2);

            if length(cycle_inter) == length(cycle1) && length(cycle_inter) == length(cycle2)
                if ~isAinB(j, delete_cycle)
                    delete_cycle = [delete_cycle, j];
                end
            elseif length(cycle_inter)/length(cycle1) >= 3/4 && length(cycle_inter)/length(cycle2) >= 3/4
                if length(cycle1) > length(cycle2)
                    if ~isAinB(i, delete_cycle)
                        delete_cycle = [delete_cycle, i];
                        break;
                    end
                else
                    if ~isAinB(j, delete_cycle)
                        delete_cycle = [delete_cycle, j];
                    end
                end
            end
        end
    end

    delete_cycle = sort(delete_cycle);
    if ~isempty(delete_cycle)
        for no_seq = 1: length(delete_cycle)
             no_seq_delete = delete_cycle(no_seq);
             cycle_min(no_seq_delete - no_seq + 1) = [];
        end
    end
else
    cycle_min = [];
end
