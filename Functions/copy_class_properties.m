function obj2 = copy_class_properties(obj1, obj2)
%COPY_PROPERTIES copies common properties from one arbitrary class to
%another
%   Classes do not have to be the same type, if they contain common
%   properties then they will be copied
%
%   Inputs:
%   obj1 - the class to copy properties from
%   obj2 - the class to copy properties to

obj1_fn = fieldnames(obj1);
obj2_fn = fieldnames(obj2);

for i = 1:numel(obj1_fn)
    
    field = obj1_fn{i};
    
    if any(strcmp(field, obj2_fn))
        
        obj2.(field) = obj1.(field);
    end
end