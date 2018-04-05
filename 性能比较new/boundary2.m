%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to find all the boundary cycles of coverage holes based
% on location information. We use the method proposed in the paper "On 
% Discovering Sensing Coverage Holes in Large-Scale Sensor Networks".

function [cycle_seq_coor] = boundary2()
%clc;
%clear;
global node_coor node_x node_y

% set the fence_flag for all nodes to indicate whether they are fence node
% or not
for i = 1: length(node_coor) 
    if i > length(node_coor) - 24
        node(i).fence_flag = 1;
    else
        node(i).fence_flag = 0;
    end
end

% set the fence_flag for all nodes to indicate whether they are fence node
% or not
for i = 1: length(node_coor) 
    if i > length(node_coor) - 24
        node(i).fence_flag = 1;
    else
        node(i).fence_flag = 0;
    end
end

% to find neighbours in order to build 1-simplex
for i=1: length(node_coor) - 24
    node(i).simp1 = [];
    node(i).dist = [];
    node(i).angle = [];
    node(i).arc_start = [];
    node(i).arc_end = [];
    for j=1: length(node_coor)
        if (j==i) continue;
        else
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node(i).simp1 = [node(i).simp1, j];
                node(i).dist = [node(i).dist, distance];
                node(i).angle = [node(i).angle, 2*acos(distance/2)];
                
                delta_x = node_x(j) - node_x(i);
                delta_y = node_y(j) - node_y(i);
                center_angle = acos(delta_x/distance);
                if delta_y < 0
                    center_angle = 2*pi - center_angle;
                end
                
                start_angle = center_angle - acos(distance/2);
                end_angle = center_angle + acos(distance/2);
                
                if start_angle < 0 
                    start_angle = 2*pi + start_angle;
                end
                
                if end_angle > 2*pi
                    end_angle = end_angle - 2*pi;
                end
                
                node(i).arc_start = [node(i).arc_start, start_angle];
                node(i).arc_end = [node(i).arc_end, end_angle];
                if ~(node(j).fence_flag == 1)
                     % line([node_x(i),node_x(j)],[node_y(i),node_y(j)]);
                end
            end
        end
    end
end


% sort the arcs of each node
for i = 1: length(node_coor) - 24
    [arc_sort, index] = sort(node(i).arc_start);
    
    node(i).arc_start = arc_sort;
    node(i).arc_end =  node(i).arc_end(index);
    node(i).simp1 = node(i).simp1(index);
    node(i).dist = node(i).dist(index);
    node(i).angle = node(i).angle(index);
end

% delete the arcs that are subsumed by other arcs
for i = 1: length(node_coor) - 24
    neighbor_list = node(i).simp1;
    neighbor_delete = [];
    
    if ~isempty(neighbor_list)
        no_node = length(neighbor_list);
        
        for j = 1: no_node - 1
            for k = j+1 : no_node
                node_ID1 = neighbor_list(j);
                node_ID2 = neighbor_list(k);
                
                dist12 = dist2(node_coor(:,node_ID1),node_coor(:,node_ID2));
                dist01 = node(i).dist(j);
                dist02 = node(i).dist(k);
                
                if dist12 > 2 %the other two nodes are not neighbors
                    continue;
                else
                    delta_angle = acos((dist01^2 + dist02^2 - dist12^2)/(2*dist01*dist02));
                    
                    angle1 = node(i).angle(j);
                    angle2 = node(i).angle(k);
                    
                    if angle2 - angle1 >= 2*delta_angle
                        if ~isAinB(j, neighbor_delete)
                            neighbor_delete = [neighbor_delete, j];
                        end
                    elseif angle1 - angle2 >= 2*delta_angle
                        if ~isAinB(k, neighbor_delete)
                            neighbor_delete = [neighbor_delete, k];
                        end
                    end
                end
            end
        end
    end
    
    arc_start_list = node(i).arc_start;
    arc_end_list = node(i).arc_end;
    
    neighbor_delete = sort(neighbor_delete);
    for count = 1 : length(neighbor_delete)
        delete_ID = neighbor_delete(count);
        neighbor_list(delete_ID - count + 1) = [];
        arc_start_list(delete_ID - count + 1) = [];
        arc_end_list(delete_ID - count + 1) = [];
    end
    
    node(i).simp1_left = neighbor_list;
    node(i).arc_start_left = arc_start_list;
    node(i).arc_end_left = arc_end_list;    
end

for i = 1: length(node_coor) - 24
    node(i).arc_uncovered = [];
    
    if ~isempty(node(i).simp1_left)
        no_arc = length(node(i).arc_start_left);
        
        if no_arc == 1
            start_angle = node(i).arc_end_left(1);
            end_angle = node(i).arc_start_left(1);
            start_neighbor = node(i).simp1_left(1);
            end_neighbor = node(i).simp1_left(1);
            node(i).arc_uncovered = [node(i).arc_uncovered; start_angle, end_angle, start_neighbor, end_neighbor];
        else
            curEndangle = node(i).arc_end_left(1);
            
            for j = 1 : no_arc
                if j == no_arc
                    if node(i).arc_start_left(j) <= node(i).arc_end_left(j) || node(i).arc_end_left(j) < node(i).arc_start_left(1)
                        start_angle = node(i).arc_end_left(j);
                        end_angle = node(i).arc_start_left(1);
                        start_neighbor = node(i).simp1_left(j);
                        end_neighbor = node(i).simp1_left(1);
                        if end_angle - start_angle > 1e-12
                            node(i).arc_uncovered = [node(i).arc_uncovered; start_angle, end_angle, start_neighbor, end_neighbor];
                        elseif abs(2*pi + end_angle - start_angle) > 1e-12 
                            node(i).arc_uncovered = [node(i).arc_uncovered; start_angle, end_angle, start_neighbor, end_neighbor];
                        end
                    end
                else
                    if node(i).arc_start_left(j) <= node(i).arc_end_left(j)
                        if node(i).arc_start_left(j+1) > curEndangle + 1e-12
                            start_angle = curEndangle;
                            end_angle = node(i).arc_start_left(j+1);
                            start_neighbor = node(i).simp1_left(j);
                            end_neighbor = node(i).simp1_left(j+1);
                            if end_angle - start_angle > 1e-12
                                node(i).arc_uncovered = [node(i).arc_uncovered; start_angle, end_angle, start_neighbor, end_neighbor];
                            elseif abs(2*pi + end_angle - start_angle) > 1e-12 
                                node(i).arc_uncovered = [node(i).arc_uncovered; start_angle, end_angle, start_neighbor, end_neighbor];
                            end
                        end
                    end
                            
                    curEndangle = node(i).arc_end_left(j+1);
                end
            end
        end
    end
end

node_boundary_list = [];

for i = 1: length(node_coor) - 24
    if ~isempty(node(i).arc_uncovered)
        node_boundary_list = [node_boundary_list, i];
         plot(node_coor(1,i), node_coor(2,i), '*r'); hold on;  %%
    end
end

cycle_seq_coor = [];
        
while(node_boundary_list)
%     choose proper initiating node in order to avoid some extreme cases
    if length(node_boundary_list) < 3 
        break;
    end
    for i = 1 : length(node_boundary_list)
        flag_found = 0;
        flag_initiator = 0;
        initiator_ID = node_boundary_list(i);
        for j = 1 : size(node(initiator_ID).arc_uncovered, 1)
            start_neighbor_ID = node(initiator_ID).arc_uncovered(j, 3);
            end_neighbor_set = node(start_neighbor_ID).arc_uncovered(:, 4);
            index = end_neighbor_set == initiator_ID;
            
            node_temp = node(start_neighbor_ID).arc_uncovered(index, 3);
            
            if node_temp == initiator_ID
                flag_found =1;
                break;
            end
        end
        
        if flag_found
            continue;
        else
            for j = 1 : size(node(initiator_ID).arc_uncovered, 1)
                start_angle = node(initiator_ID).arc_uncovered(j, 1);
                end_angle = node(initiator_ID).arc_uncovered(j, 2);
                if start_angle <= end_angle
                    mid_angle = (node(initiator_ID).arc_uncovered(j, 1) + node(initiator_ID).arc_uncovered(j, 2))/2;
                else
                    mid_angle = (node(initiator_ID).arc_uncovered(j, 1) + node(initiator_ID).arc_uncovered(j, 2))/2;
                    if mid_angle < pi
                        mid_angle = mid_angle + pi;
                    else
                        mid_angle = mid_angle - pi;
                    end
                end
                
                if mid_angle >= pi
                    end_neighbor_ID = node(initiator_ID).arc_uncovered(j, 4);
                    start_neigh_set = node(end_neighbor_ID).arc_uncovered(:, 3);
                    
                    index = start_neigh_set == initiator_ID;
                    start_angle2 = node(end_neighbor_ID).arc_uncovered(index, 1);
                    end_angle2 = node(end_neighbor_ID).arc_uncovered(index, 2);
                    if start_angle2 <= end_angle2
                        mid_angle2 = (node(end_neighbor_ID).arc_uncovered(index, 1) + node(end_neighbor_ID).arc_uncovered(index, 2))/2;
                    else
                        mid_angle2 = (node(end_neighbor_ID).arc_uncovered(index, 1) + node(end_neighbor_ID).arc_uncovered(index, 2))/2;
                        if mid_angle2 < pi
                            mid_angle2 = mid_angle2 + pi;
                        else
                            mid_angle2 = mid_angle2 - pi;
                        end
                    end
                    
                    if mid_angle2 < pi
                        start_neighbor_ID = node(initiator_ID).arc_uncovered(j, 3);
                        beginning_ID = node(initiator_ID).arc_uncovered(j, 4);
                        flag_initiator = 1;
                        break;
                    end
                end
            end
            
            if flag_initiator
                break;
            else
                continue;
            end
        end
    end
    
    cycle_temp = [initiator_ID, start_neighbor_ID];
    line([node_x(initiator_ID),node_x(beginning_ID)],[node_y(initiator_ID),node_y(beginning_ID)], 'Color', 'r', 'LineWidth', 2); %%
     line([node_x(initiator_ID),node_x(start_neighbor_ID)],[node_y(initiator_ID),node_y(start_neighbor_ID)], 'Color', 'r', 'LineWidth', 2);  %%
    
    node(initiator_ID).arc_uncovered(j, :) = [];
    
    while(1)
        if ~isempty(node(start_neighbor_ID).arc_uncovered)
            end_neighbor_set = node(start_neighbor_ID).arc_uncovered(:, 4);
        
            index = find(end_neighbor_set == initiator_ID);

            initiator_ID = start_neighbor_ID;
            start_neighbor_ID = node(start_neighbor_ID).arc_uncovered(index, 3);
             line([node_x(initiator_ID),node_x(start_neighbor_ID)],[node_y(initiator_ID),node_y(start_neighbor_ID)], 'Color', 'r', 'LineWidth', 2);  %%

            node(initiator_ID).arc_uncovered(index, :) = [];

            if start_neighbor_ID == cycle_temp(1) && initiator_ID == beginning_ID
                end_neighbor_set = node(start_neighbor_ID).arc_uncovered(:, 4);
        
                index = find(end_neighbor_set == initiator_ID);
                node(start_neighbor_ID).arc_uncovered(index, :) = [];
                break;
            else
                cycle_temp = [cycle_temp, start_neighbor_ID];
            end
        else
            break;
        end
    end
    
    no_cycle = size(cycle_seq_coor, 2);
    cycle_seq_coor{no_cycle+1} = cycle_temp;  %cycleµÄ¼¯ºÏ
    
    node_boundary_list_temp = node_boundary_list;
    node_boundary_list = [];
    
    for i = 1: length(node_boundary_list_temp)
        node_boundary_left = node_boundary_list_temp(i);
        if ~isempty(node(node_boundary_left).arc_uncovered)
             node_boundary_list = [node_boundary_list, node_boundary_left];
        end
    end
    
end

% shorter the cycle containing two same nodes
no_cycle = size(cycle_seq_coor, 2);
for i = 1 : no_cycle
    cycle_temp = cycle_seq_coor{i};
    while(1)
        flag = 0;
        for j = 1 : length(cycle_temp)-1
            for k = j+1 : length(cycle_temp)
                if cycle_temp(j) == cycle_temp(k)
                    flag = 1;
                    break;
                end
            end
            if flag == 1
                break;
            end
        end
        if flag == 1
            if k-j <= length(cycle_temp)/2
                for no = 1 : k-j
                    cycle_seq_coor{i}(j+1) = [];
                end
            else
                cycle_seq_coor{i} = cycle_temp(j: k-1);
            end
            cycle_temp = cycle_seq_coor{i};
        else
            break;
        end
    end
end
