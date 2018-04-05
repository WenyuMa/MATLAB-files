% This is the function to check whether a cycle with length 4 can be
% triangulated or not.

function flag_hole = cycle_4_check(cycle, node)

node_ID1 = cycle(1);
node_ID2 = cycle(2);
node_ID3 = cycle(3);
node_ID4 = cycle(4);
if ismember(node_ID1, node(node_ID3).neighbors) || ismember(node_ID2, node(node_ID4).neighbors)
    flag_hole = 0;
else
    common_set1 = intersect(node(node_ID1).neighbors, node(node_ID2).neighbors);
    common_set2 = intersect(node(node_ID3).neighbors, node(node_ID4).neighbors);
    common_set = intersect(common_set1, common_set2);
    if isempty(common_set)
        flag_found = 0;
        for j = 1 : length(cycle)    
            common_neigh_set = intersect(node(cycle(j)).neighbors, node(cycle(mod(j, length(cycle))+1)).neighbors);
            if j == 1
                common_neigh_set = intersect(common_neigh_set, node(cycle(length(cycle))).neighbors);
            else
                common_neigh_set = intersect(common_neigh_set, node(cycle(j-1)).neighbors);
            end
            if ~isempty(common_neigh_set)
                for k = 1 : length(common_neigh_set)
                    node_temp_ID = common_neigh_set(k);
                    common_set = node(node_temp_ID).neighbors;
                    cycle_backup = cycle;
                    cycle_backup(j) = node_temp_ID;
                    cycle_temp = cycle;
                    cycle_temp(j) = [];
                    for n = 1 : length(cycle_temp)
                        common_set = intersect(common_set, node(cycle_temp(n)).neighbors);
                    end
                    if ~isempty(common_set)
                        flag_found = 1;
                        flag_hole = 0;
                        break;
                    else
                        for j2 = 1 : length(cycle_backup)
                            common_neigh_set2 = intersect(node(cycle_backup(j2)).neighbors, node(cycle_backup(mod(j2, length(cycle_backup))+1)).neighbors);
                            if j2 == 1
                                common_neigh_set2 = intersect(common_neigh_set2, node(cycle_backup(length(cycle_backup))).neighbors);
                            else
                                common_neigh_set2 = intersect(common_neigh_set2, node(cycle_backup(j2-1)).neighbors);
                            end
                            if ~isempty(common_neigh_set2)
                                for k2 = 1 : length(common_neigh_set2)
                                    node_temp_ID2 = common_neigh_set2(k2);
                                    common_set2 = node(node_temp_ID2).neighbors;
                                    cycle_temp2 = cycle_backup;
                                    cycle_temp2(j2) = [];
                                    for n2 = 1 : length(cycle_temp2)
                                        common_set2 = intersect(common_set2, node(cycle_temp2(n2)).neighbors);
                                    end
                                    if ~isempty(common_set2)
                                        flag_found = 1;
                                        flag_hole = 0;
                                        break;
                                    end
                                end
                            end
                            if flag_found == 1
                                break;
                            end
                        end
                        if flag_found == 1
                            break;
                        end
                    end
                end
            end
            if flag_found == 1
                break;
            end
        end
        if flag_found == 0
            flag_hole = 1;
        end
    else
        flag_hole = 0;
    end
end


                    