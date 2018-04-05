function [cycle_seq_coor] = find_cycle2()

global node_boundary_list node2 node_x node_y
% fiond cycle  node_boundary_list

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
        for j = 1 : size(node2(initiator_ID).arc_uncovered, 1)
            start_neighbor_ID = node2(initiator_ID).arc_uncovered(j, 3);
            end_neighbor_set = node2(start_neighbor_ID).arc_uncovered(:, 4);
            index = end_neighbor_set == initiator_ID;
            
            node_temp = node2(start_neighbor_ID).arc_uncovered(index, 3);
            
            if node_temp == initiator_ID
                flag_found =1;
                break;
            end
        end
        
        if flag_found
            continue;
        else
            for j = 1 : size(node2(initiator_ID).arc_uncovered, 1)
                start_angle = node2(initiator_ID).arc_uncovered(j, 1);
                end_angle = node2(initiator_ID).arc_uncovered(j, 2);
                if start_angle <= end_angle
                    mid_angle = (node2(initiator_ID).arc_uncovered(j, 1) + node2(initiator_ID).arc_uncovered(j, 2))/2;
                else
                    mid_angle = (node2(initiator_ID).arc_uncovered(j, 1) + node2(initiator_ID).arc_uncovered(j, 2))/2;
                    if mid_angle < pi
                        mid_angle = mid_angle + pi;
                    else
                        mid_angle = mid_angle - pi;
                    end
                end
                
                if mid_angle >= pi
                    end_neighbor_ID = node2(initiator_ID).arc_uncovered(j, 4);
                    start_neigh_set = node2(end_neighbor_ID).arc_uncovered(:, 3);
                    
                    index = start_neigh_set == initiator_ID;
                    start_angle2 = node2(end_neighbor_ID).arc_uncovered(index, 1);
                    end_angle2 = node2(end_neighbor_ID).arc_uncovered(index, 2);
                    if start_angle2 <= end_angle2
                        mid_angle2 = (node2(end_neighbor_ID).arc_uncovered(index, 1) + node2(end_neighbor_ID).arc_uncovered(index, 2))/2;
                    else
                        mid_angle2 = (node2(end_neighbor_ID).arc_uncovered(index, 1) + node2(end_neighbor_ID).arc_uncovered(index, 2))/2;
                        if mid_angle2 < pi
                            mid_angle2 = mid_angle2 + pi;
                        else
                            mid_angle2 = mid_angle2 - pi;
                        end
                    end
                    
                    if mid_angle2 < pi
                        start_neighbor_ID = node2(initiator_ID).arc_uncovered(j, 3);
                        beginning_ID = node2(initiator_ID).arc_uncovered(j, 4);
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
    
    node2(initiator_ID).arc_uncovered(j, :) = [];
    
    while(1)
        if ~isempty(node2(start_neighbor_ID).arc_uncovered)
            end_neighbor_set = node2(start_neighbor_ID).arc_uncovered(:, 4);
        
            index = find(end_neighbor_set == initiator_ID);

            initiator_ID = start_neighbor_ID;
            start_neighbor_ID = node2(start_neighbor_ID).arc_uncovered(index, 3);
             line([node_x(initiator_ID),node_x(start_neighbor_ID)],[node_y(initiator_ID),node_y(start_neighbor_ID)], 'Color', 'r', 'LineWidth', 2);  %%

            node2(initiator_ID).arc_uncovered(index, :) = [];

            if start_neighbor_ID == cycle_temp(1) && initiator_ID == beginning_ID
                end_neighbor_set = node2(start_neighbor_ID).arc_uncovered(:, 4);
        
                index = find(end_neighbor_set == initiator_ID);
                node2(start_neighbor_ID).arc_uncovered(index, :) = [];
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
        if ~isempty(node2(node_boundary_left).arc_uncovered)
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