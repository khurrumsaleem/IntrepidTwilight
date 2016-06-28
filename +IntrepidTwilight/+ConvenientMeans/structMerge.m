function default = structMerge(default,new)
    
    if isstruct(default) && isstruct(new) && not(isempty(fieldnames(new)))
        default = merge(default,new);
    end
    
end

function default = merge(default,new)
    newFields = fieldnames(new);
    for k = 1:numel(newFields)
        newField = newFields{k};
        if isfield(default,newField)
            
            defaultFieldIsStruct = isstruct(default.(newField));
            newFieldIsStruct     = isstruct(new.(newField))    ;
            
            if not(defaultFieldIsStruct) && not(newFieldIsStruct),
                default.(newField) = new.(newField);

            elseif defaultFieldIsStruct && newFieldIsStruct
                default.(newField) = merge(default.(newField),new.(newField));
                
            end
        else
            
            default.(newField) = new.(newField);
            
        end
    end
end

