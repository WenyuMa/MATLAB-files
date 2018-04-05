function [] = boundary2()
%clc;
%clear;
global node2 node_coor node_x node_y node_boundary_list k

% to find neighbours in order to build 1-simplex
for i=1: length(node_coor)   %- 24*k  %++
    node2(i).neighbour = [];
    node2(i).dist0 = [];
    node2(i).angle0 = [];
    node2(i).arc_start = [];
    node2(i).arc_end = [];
    
    if node2(i).status==1
      for j=1: length(node_coor)
        if (j==i) 
            continue;
        elseif (node2(j).status==1)
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node2(i).neighbour = [node2(i).neighbour, j];
                node2(i).dist0 = [node2(i).dist0, distance];
                node2(i).angle0 = [node2(i).angle0, 2*acos(distance/2)];
                
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
                
                node2(i).arc_start = [node2(i).arc_start, start_angle];
                node2(i).arc_end = [node2(i).arc_end, end_angle];
                %if ~(node2(j).fence_flag == 1)
                     % line([node_x(i),node_x(j)],[node_y(i),node_y(j)]);
                %end
            end
        end
      end 
    end     
end

% sort the arcs of each node
for i = 1: length(node_coor) - 24*k  %++
  if node2(i).status==1
    [arc_sort, index] = sort(node2(i).arc_start);
    
    node2(i).arc_start = arc_sort;
    node2(i).arc_end =  node2(i).arc_end(index);
    node2(i).simp1 = node2(i).neighbour(index);
    node2(i).dist = node2(i).dist0(index);
    node2(i).angle = node2(i).angle0(index);
  end
end

% delete the arcs that are subsumed by other arcs
for i = 1: length(node_coor) - 24*k   %++
  if node2(i).status==1
    neighbor_list = node2(i).simp1;
    neighbor_delete = [];
    
    if ~isempty(neighbor_list)
        no_node = length(neighbor_list);
        
        for j = 1: no_node - 1
            for m = j+1 : no_node
                node_ID1 = neighbor_list(j);
                node_ID2 = neighbor_list(m);
                
                dist12 = dist2(node_coor(:,node_ID1),node_coor(:,node_ID2));
                dist01 = node2(i).dist(j);
                dist02 = node2(i).dist(m);
                
                if dist12 > 2 %the other two nodes are not neighbors
                    continue;
                else
                    delta_angle = acos((dist01^2 + dist02^2 - dist12^2)/(2*dist01*dist02));
                    
                    angle1 = node2(i).angle(j);
                    angle2 = node2(i).angle(m);
                    
                    if angle2 - angle1 >= 2*delta_angle
                        if ~isAinB(j, neighbor_delete)
                            neighbor_delete = [neighbor_delete, j];
                        end
                    elseif angle1 - angle2 >= 2*delta_angle
                        if ~isAinB(m, neighbor_delete)
                            neighbor_delete = [neighbor_delete, m];
                        end
                    end
                end
            end
        end
    end
    
    arc_start_list = node2(i).arc_start;
    arc_end_list = node2(i).arc_end;
    
    neighbor_delete = sort(neighbor_delete);
    for count = 1 : length(neighbor_delete)
        delete_ID = neighbor_delete(count);
        neighbor_list(delete_ID - count + 1) = [];
        arc_start_list(delete_ID - count + 1) = [];
        arc_end_list(delete_ID - count + 1) = [];
    end
    
    node2(i).simp1_left = neighbor_list;
    node2(i).arc_start_left = arc_start_list;
    node2(i).arc_end_left = arc_end_list; 
  end 
end

for i = 1: length(node_coor) - 24*k   %++
  if node2(i).status==1
   
    node2(i).arc_uncovered = [];   
    if ~isempty(node2(i).simp1_left)
        no_arc = length(node2(i).arc_start_left);
        
        if no_arc == 1
            start_angle = node2(i).arc_end_left(1);
            end_angle = node2(i).arc_start_left(1);
            start_neighbor = node2(i).simp1_left(1);
            end_neighbor = node2(i).simp1_left(1);
            node2(i).arc_uncovered = [node2(i).arc_uncovered; start_angle, end_angle, start_neighbor, end_neighbor];
        else
            curEndangle = node2(i).arc_end_left(1);
            
            for j = 1 : no_arc
                if j == no_arc
                    if node2(i).arc_start_left(j) <= node2(i).arc_end_left(j) || node2(i).arc_end_left(j) < node2(i).arc_start_left(1)
                        start_angle = node2(i).arc_end_left(j);
                        end_angle = node2(i).arc_start_left(1);
                        start_neighbor = node2(i).simp1_left(j);
                        end_neighbor = node2(i).simp1_left(1);
                        if end_angle - start_angle > 1e-12
                            node2(i).arc_uncovered = [node2(i).arc_uncovered; start_angle, end_angle, start_neighbor, end_neighbor];
                        elseif abs(2*pi + end_angle - start_angle) > 1e-12 
                            node2(i).arc_uncovered = [node2(i).arc_uncovered; start_angle, end_angle, start_neighbor, end_neighbor];
                        end
                    end
                else
                    if node2(i).arc_start_left(j) <= node2(i).arc_end_left(j)
                        if node2(i).arc_start_left(j+1) > curEndangle + 1e-12
                            start_angle = curEndangle;
                            end_angle = node2(i).arc_start_left(j+1);
                            start_neighbor = node2(i).simp1_left(j);
                            end_neighbor = node2(i).simp1_left(j+1);
                            if end_angle - start_angle > 1e-12
                                node2(i).arc_uncovered = [node2(i).arc_uncovered; start_angle, end_angle, start_neighbor, end_neighbor];
                            elseif abs(2*pi + end_angle - start_angle) > 1e-12 
                                node2(i).arc_uncovered = [node2(i).arc_uncovered; start_angle, end_angle, start_neighbor, end_neighbor];
                            end
                        end
                    end
                            
                    curEndangle = node2(i).arc_end_left(j+1);
                end
            end
        end
    end
  end 
end

node_boundary_list = [];

for i = 1: length(node_coor) - 24*k    %++
  if node2(i).status==1
    if ~isempty(node2(i).arc_uncovered)
        node_boundary_list = [node_boundary_list, i];
        node2(i).bn=1;
        plot(node_coor(1,i), node_coor(2,i), '*r'); hold on;  %%
    end
  end
end

%--------------------------------------------------------------------------