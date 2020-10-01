function new = builder(var, num, var_names, varargin)

[nPop, ~] = size(var);

parts = varargin;

pos_ind = 1:sum(num);
cumnum = cumsum(num);

next = 1;
part = 0;
while true
    
    if part == 0 || next > cumnum(part)
        
        part = part + 1;
    end
    
    current = var_names(next);
    con = current == var_names;
    
    vari = mat2cell(var(:,con), ones(nPop, 1), sum(con));
    
    [parts{part}.(current)] = vari{:};
    
    next = find(~con & pos_ind > next, 1);
    if isempty(next), break; end
end

%% Change to fully cell based
% Cannot do earlier due to struct access above
for i = numel(parts):-1:1

    new(:,i) = mat2cell(parts{i}, ones(nPop, 1));
end
end