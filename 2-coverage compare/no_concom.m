% This is the function to compute the number of connected component in a 
% graph using adjacent matrix. 

function [no_cc] = no_concom(adj_matrix)

no_node = size(adj_matrix, 1);
flag = zeros(no_node, 1);
queue = [];
no_cc = 0;

for m = 1 : no_node
    if flag(m) == 0
        queue = [queue, m];
    end
    
    if ~isempty(queue)
        no_cc = no_cc + 1;

        while(~isempty(queue))
            vertex = queue(1);
            flag(vertex) = 1;
            queue(1) = [];

            for n = 1 : no_node
                if adj_matrix(vertex, n) == 1
                    if flag(n) == 0 && ~isAinB(n, queue)
                        queue = [queue, n];
                    end
                end
            end
        end
    end
end
