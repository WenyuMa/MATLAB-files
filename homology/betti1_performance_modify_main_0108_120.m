% build points according to poisson point process
clc;
clear;

node_coor_set = [];
no_case = 0;

for count = 1 : 1000
    % lambda is the intensity 
    lambda=120;
    nmb=poissrnd(lambda);

    x = rand(1,nmb);
    y = rand(1,nmb);
    x = 10*x + 1;
    y = 10*y + 1;

    %coordinates of fence nodes
    fence_x = [0,0,0,0,0, 0, 0,2, 2,4, 4,6, 6,8, 8,10,10,12,12,12,12,12,12,12];
    fence_y = [0,2,4,6,8,10,12,0,12,0,12,0,12,0,12, 0,12, 0, 2, 4, 6, 8,10,12];

    inner_x = [1,1,1,1,1, 1,3, 3,5, 5,7, 7,9, 9,11,11,11,11,11,11];
    inner_y = [1,3,5,7,9,11,1,11,1,11,1,11,1,11, 1, 3, 5, 7, 9,11];

    %coordinates of all the nodes
    node_x_LBS = [x, inner_x, fence_x];
    node_y_LBS = [y, inner_y, fence_y];

    node_x_homology = [x, inner_x];
    node_y_homology = [y, inner_y];

    node_coor_LBS = [node_x_LBS; node_y_LBS];
    node_coor_homology = [node_x_homology; node_y_homology];
    
    dlmwrite('data120.txt', node_coor_LBS, 'delimiter', ' ');

    count
%     cycle_seq_coor{count} = betti1_coor_fun(node_coor);
    [cycle_seq_LBS] = betti1_coor_fun_080413(node_coor_LBS);
    [cycle_seq_homology, betti1_jplex] = betti1_homology_fun_1103_2(node_coor_homology);

    no_cycle_LBS = size(cycle_seq_LBS, 2);
    no_cycle_homology = size(cycle_seq_homology, 2);
    
    betti1_LBS(count) = no_cycle_LBS;
    betti1_Hom(count) = no_cycle_homology;
    betti1_jp(count) = betti1_jplex;
    
    if no_cycle_homology ~= betti1_jplex || no_cycle_LBS ~= betti1_jplex 
        no_case = size(node_coor_set, 2);
        node_coor_set{no_case + 1} = node_coor_LBS;
        continue;
    else
        % to find neighbours in order to build 1-simplex
        for i=1: length(node_coor_LBS)
            node(i).neighbors = [];
            for j=1: length(node_coor_LBS)
                if (j==i) continue;
                else
                    distance = dist2(node_coor_LBS(:,i),node_coor_LBS(:,j));
                    if distance <= 2 
                        node(i).neighbors = [node(i).neighbors, j];
                    end
                end
            end
        end

        cycle_seq_LBS_temp = cycle_seq_LBS;

        cycle_seq_homology_left = [];
        no_same_cycle = 0;
        for i = 1 : no_cycle_homology
            cycle_homology_temp = cycle_seq_homology{i};

            flag_same = 0;
            no_cycle_LBS_temp = size(cycle_seq_LBS_temp, 2);
            for j = 1 : no_cycle_LBS_temp
                cycle_LBS_temp = cycle_seq_LBS_temp{j};

                cycle_inter = intersect(cycle_homology_temp, cycle_LBS_temp);

                if length(cycle_inter)/length(cycle_homology_temp) >= 3/4  && length(cycle_inter)/length(cycle_LBS_temp) >= 3/4
                    flag_same = 1;
                    break;
                elseif length(cycle_inter) >= 1
                    node_ID_1st = cycle_inter(1);
                    post1 = find(cycle_LBS_temp == node_ID_1st);
                    post2 = find(cycle_homology_temp == node_ID_1st);
                    cycle_LBS_temp = circshift(cycle_LBS_temp, [1, -(post1-1)]);
                    cycle_homology_temp = circshift(cycle_homology_temp, [1, -(post2-1)]);

                    node_ID_left_LBS = cycle_LBS_temp(length(cycle_LBS_temp));
                    node_ID_right_LBS = cycle_LBS_temp(2);

                    node_ID_left_homology = cycle_homology_temp(length(cycle_homology_temp));
                    node_ID_right_homology = cycle_homology_temp(2);

                    if node_ID_left_LBS == node_ID_left_homology
                        cycle_LBS_temp = fliplr(circshift(cycle_LBS_temp, [1, -1]));
                        cycle_homology_temp = fliplr(circshift(cycle_homology_temp, [1, -1]));
                    elseif node_ID_left_LBS == node_ID_right_homology
                        cycle_LBS_temp = fliplr(circshift(cycle_LBS_temp, [1, -1]));
                    elseif node_ID_right_LBS == node_ID_left_homology
                        cycle_homology_temp = fliplr(circshift(cycle_homology_temp, [1, -1]));
                    elseif node_ID_right_LBS == node_ID_right_homology

                    elseif ismember(node_ID_left_LBS, node(node_ID_left_homology).neighbors)
                        cycle_LBS_temp = fliplr(circshift(cycle_LBS_temp, [1, -1]));
                        cycle_homology_temp = fliplr(circshift(cycle_homology_temp, [1, -1]));
                    elseif ismember(node_ID_left_LBS, node(node_ID_right_homology).neighbors)
                        cycle_LBS_temp = fliplr(circshift(cycle_LBS_temp, [1, -1]));
                    elseif ismember(node_ID_right_LBS, node(node_ID_left_homology).neighbors)
                        cycle_homology_temp = fliplr(circshift(cycle_homology_temp, [1, -1]));
                    elseif ismember(node_ID_right_LBS, node(node_ID_right_homology).neighbors)

                    else
                        continue;
                    end

                    pt_L = 2;
                    pt_H = 2;
                    while(1)
                        cycle_inter = intersect(cycle_homology_temp, cycle_LBS_temp);
                        if length(cycle_inter)/length(cycle_homology_temp) >= 3/4  && length(cycle_inter)/length(cycle_LBS_temp) >= 3/4
                            flag_same = 1;
                            break;
                        end

                        if pt_H > length(cycle_homology_temp) || pt_L > length(cycle_LBS_temp)
                            flag_same = 1;
                            break;
                        end

                        if cycle_homology_temp(pt_H) == cycle_LBS_temp(pt_L)
                            pt_L = pt_L + 1;
                            pt_H = pt_H + 1;
                        else
                            node_ID_L = cycle_LBS_temp(pt_L);
                            node_ID_H_mid = cycle_homology_temp(pt_H);
                            node_ID_H_after = cycle_homology_temp(mod(pt_H, length(cycle_homology_temp))+1);

                            if ismember(node_ID_H_mid, node(node_ID_L).neighbors) && ismember(node_ID_H_after, node(node_ID_L).neighbors)
                                cycle_homology_temp(pt_H) = cycle_LBS_temp(pt_L);
                                pt_L = pt_L + 1;
                                pt_H = pt_H + 1;
                            elseif ismember(node_ID_H_mid, node(node_ID_L).neighbors)
                                hom_seq_temp = [cycle_homology_temp(1 : pt_H-1), cycle_LBS_temp(pt_L), cycle_homology_temp(pt_H : length(cycle_homology_temp))];
                                cycle_homology_temp = hom_seq_temp;
                                pt_L = pt_L + 1;
                                pt_H = pt_H + 1;
                            else
                                node_ID_H_before = cycle_homology_temp(pt_H - 1);
                                inter1 = intersect(node(node_ID_L).neighbors, node(node_ID_H_before).neighbors);
                                inter2 = intersect(node(node_ID_H_mid).neighbors, node(node_ID_H_after).neighbors);
                                inter = intersect(inter1, inter2);

                                if ~isempty(inter)
                                    cycle_homology_temp(pt_H) = cycle_LBS_temp(pt_L);
                                    hom_seq_temp = [cycle_homology_temp(1 : pt_H), inter(1), cycle_homology_temp(pt_H+1 : length(cycle_homology_temp))];
                                    cycle_homology_temp = hom_seq_temp;
                                    pt_L = pt_L + 1;
                                    pt_H = pt_H + 1;
                                else
                                    break;
                                end
                            end
                            % check whether there is shorter path
                            node_ID_H_before = cycle_homology_temp(pt_H - 1);
                            if pt_H > 3 && pt_H <= length(cycle_homology_temp)
                                pt_temp = pt_H;
                                while(1)
                                    if pt_temp == length(cycle_homology_temp)
                                        node_ID_H_after = cycle_homology_temp(1);
                                        if ismember(node_ID_H_before, node(node_ID_H_after).neighbors)
                                            cycle_homology_temp = cycle_homology_temp(1: pt_H - 1);
                                        end
                                        break;
                                    else
                                        pt_temp = pt_temp + 1;
                                        node_ID_H_after = cycle_homology_temp(pt_temp);
                                        if ismember(node_ID_H_before, node(node_ID_H_after).neighbors)
                                            for no = 1 : pt_temp - pt_H
                                                cycle_homology_temp(pt_H) = [];
                                            end
                                            pt_temp = pt_H;
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if flag_same == 1
                        break;
                    else
                        cycle_homology_temp = cycle_seq_homology{i};
                    end
                end
            end
            if flag_same == 1
                no_same_cycle = no_same_cycle + 1;
                cycle_seq_LBS_temp(j) = [];
            else
                no_left = size(cycle_seq_homology_left, 2);
                cycle_seq_homology_left{no_left+1} = cycle_homology_temp;
            end
        end

        % this is for the case that two cycles bound the same hole but they have no
        % common node
        if ~isempty(cycle_seq_homology_left) && ~isempty(cycle_seq_LBS_temp)
            for i = 1 : size(cycle_seq_homology_left, 2);
                cycle_homology_temp = cycle_seq_homology_left{i};

                flag_same = 0;
                no_cycle_LBS_temp = size(cycle_seq_LBS_temp, 2);
                for j = 1 : no_cycle_LBS_temp
                    cycle_LBS_temp = cycle_seq_LBS_temp{j};

                    node_ID_1st_homology = cycle_homology_temp(1);
                    inter = intersect(cycle_LBS_temp, node(node_ID_1st_homology).neighbors);

                    if isempty(inter)
                        continue;
                    else
                        node_ID_1st_LBS = inter(1);
                        post = find(cycle_LBS_temp == node_ID_1st_LBS);
                        cycle_LBS_temp = circshift(cycle_LBS_temp, [1, -(post-1)]);

                        node_ID_left_LBS = cycle_LBS_temp(length(cycle_LBS_temp));
                        node_ID_right_LBS = cycle_LBS_temp(2);

                        node_ID_left_homology = cycle_homology_temp(length(cycle_homology_temp));
                        node_ID_right_homology = cycle_homology_temp(2);

                        if ismember(node_ID_1st_LBS, node(node_ID_left_homology).neighbors)
                            cycle_homology_temp = [node_ID_1st_LBS, cycle_homology_temp];
                        elseif ismember(node_ID_1st_LBS, node(node_ID_right_homology).neighbors)
                            cycle_homology_temp = [node_ID_1st_homology, node_ID_1st_LBS, cycle_homology_temp(2:length(cycle_homology_temp))];
                            cycle_homology_temp = circshift(cycle_homology_temp, [1, -1]);
                        end

                        node_ID_left_homology = cycle_homology_temp(length(cycle_homology_temp));
                        node_ID_right_homology = cycle_homology_temp(2);

                        if ismember(node_ID_left_LBS, node(node_ID_left_homology).neighbors)
                            cycle_LBS_temp = fliplr(circshift(cycle_LBS_temp, [1, -1]));
                            cycle_homology_temp = fliplr(circshift(cycle_homology_temp, [1, -1]));
                        elseif ismember(node_ID_left_LBS, node(node_ID_right_homology).neighbors)
                            cycle_LBS_temp = fliplr(circshift(cycle_LBS_temp, [1, -1]));
                        elseif ismember(node_ID_right_LBS, node(node_ID_left_homology).neighbors)
                            cycle_homology_temp = fliplr(circshift(cycle_homology_temp, [1, -1]));
                        elseif ismember(node_ID_right_LBS, node(node_ID_right_homology).neighbors)

                        else
                            continue;
                        end

                        pt_L = 2;
                        pt_H = 2;
                        while(1)
                            cycle_inter = intersect(cycle_homology_temp, cycle_LBS_temp);
                            if length(cycle_inter)/length(cycle_homology_temp) >= 3/4  && length(cycle_inter)/length(cycle_LBS_temp) >= 3/4
                                flag_same = 1;
                                break;
                            end

                            if pt_H > length(cycle_homology_temp) || pt_L > length(cycle_LBS_temp)
                                flag_same = 1;
                                break;
                            end

                            if cycle_homology_temp(pt_H) == cycle_LBS_temp(pt_L)
                                pt_L = pt_L + 1;
                                pt_H = pt_H + 1;
                            else
                                node_ID_L = cycle_LBS_temp(pt_L);
                                node_ID_H_mid = cycle_homology_temp(pt_H);
                                node_ID_H_after = cycle_homology_temp(mod(pt_H, length(cycle_homology_temp))+1);

                                if ismember(node_ID_H_mid, node(node_ID_L).neighbors) && ismember(node_ID_H_after, node(node_ID_L).neighbors)
                                    cycle_homology_temp(pt_H) = cycle_LBS_temp(pt_L);
                                    pt_L = pt_L + 1;
                                    pt_H = pt_H + 1;
                                elseif ismember(node_ID_H_mid, node(node_ID_L).neighbors)
                                    hom_seq_temp = [cycle_homology_temp(1 : pt_H-1), cycle_LBS_temp(pt_L), cycle_homology_temp(pt_H : length(cycle_homology_temp))];
                                    cycle_homology_temp = hom_seq_temp;
                                    pt_L = pt_L + 1;
                                    pt_H = pt_H + 1;
                                else
                                    node_ID_H_before = cycle_homology_temp(pt_H - 1);
                                    inter1 = intersect(node(node_ID_L).neighbors, node(node_ID_H_before).neighbors);
                                    inter2 = intersect(node(node_ID_H_mid).neighbors, node(node_ID_H_after).neighbors);
                                    inter = intersect(inter1, inter2);

                                    if ~isempty(inter)
                                        cycle_homology_temp(pt_H) = cycle_LBS_temp(pt_L);
                                        hom_seq_temp = [cycle_homology_temp(1 : pt_H), inter(1), cycle_homology_temp(pt_H+1 : length(cycle_homology_temp))];
                                        cycle_homology_temp = hom_seq_temp;
                                        pt_L = pt_L + 1;
                                        pt_H = pt_H + 1;
                                    else
                                        break;
                                    end
                                end
                                % check whether there is shorter path
                                node_ID_H_before = cycle_homology_temp(pt_H - 1);
                                if pt_H > 3 && pt_H <= length(cycle_homology_temp)
                                    pt_temp = pt_H;
                                    while(1)
                                        if pt_temp == length(cycle_homology_temp)
                                            node_ID_H_after = cycle_homology_temp(1);
                                            if ismember(node_ID_H_before, node(node_ID_H_after).neighbors)
                                                cycle_homology_temp = cycle_homology_temp(1: pt_H - 1);
                                            end
                                            break;
                                        else
                                            pt_temp = pt_temp + 1;
                                            node_ID_H_after = cycle_homology_temp(pt_temp);
                                            if ismember(node_ID_H_before, node(node_ID_H_after).neighbors)
                                                for no = 1 : pt_temp - pt_H
                                                    cycle_homology_temp(pt_H) = [];
                                                end
                                                pt_temp = pt_H;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if flag_same == 1
                            break;
                        else
                            cycle_homology_temp = cycle_seq_homology_left{i};
                        end
                    end
                end
                if flag_same == 1
                    no_same_cycle = no_same_cycle + 1;
                    cycle_seq_LBS_temp(j) = [];
                end
            end   
        end

        if no_same_cycle ~= betti1_jplex
            no_case = size(node_coor_set, 2);
            node_coor_set{no_case + 1} = node_coor_LBS;
        end
    end
end































