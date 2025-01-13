function firstpick(label::String)
    if isempty(label)
        return ""
    else
        return string(label[1])
    end
end