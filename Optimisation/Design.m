classdef Design
    
    properties
        
        variables
        cost
        penalty
        penalties
    end
    methods
        function obj = design(varargin)
            %% Individual design class
            % Input is in-order class properties (variables, cost, etc)
            
            if nargin >= 1
                
                temp = varargin{1};
                
                [rows, ~] = size(temp);
                obj = obj.array(rows);
                fn = fieldnames(obj);
                
                for i = 1:nargin
                    for j = 1:rows
                        
                        obj(j).(fn{i}) = varargin{i}(j,:)';
                    end
                end
            end
        end
    end
    
    methods (Static)
        function obj = array(i, j)
            %% Create empty object array
            
            if nargin < 2, j = 1; end
            
            for i = i:-1:1
                for j = j:-1:1
                    
                    obj(i,j) = design();
                end
            end
        end
        function varargout = fromstruct(varargin)
            %% Turn struct/object of similar makeup to individual designs
            
            for arg = nargin:-1:1
                
                s = varargin{arg};
                
                fn_inp = fieldnames(s);
                [rows, ~] = size(s.(fn_inp{1}));
                
                obj = design.array(rows);
                
                fn_obj = fieldnames(obj);
                
                % Loop through input field names
                for i = 1:length(fn_inp)
                    % Loop through output (design object) field names
                    for j = 1:length(fn_obj)
                        
                        % Find matches
                        if strcmp(fn_inp{i}, fn_obj{j})
                            
                            % Struct matrices to individual (vector)
                            % designs
                            for k = 1:rows
                                
                                obj(k).(fn_obj{j}) = s.(fn_inp{i})(k,:);
                            end
                            break
                        end
                    end
                end
                
                varargout{arg} = obj;
            end
        end
    end
end