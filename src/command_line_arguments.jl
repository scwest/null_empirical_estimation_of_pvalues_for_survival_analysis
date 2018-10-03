"""
functions for parsing specific input arguments
"""
function upload_expression(input_filename, clinical_patient_order)
    expression_mat, element_order =
    open(input_filename) do infile
        patient_order = split(strip(readline(infile)), ",")[2:end]
        indx = sortperm(patient_order, by=i->clinical_patient_order[i])
        element_order = String[]
        expression_mat = Matrix(0, length(patient_order))
        for line in eachline(infile)
            line = split(strip(line), ",")
            append!(element_order, line[1])
            expressions = [parse(Float64, x) for x in line[2:end]]
            expressions = expressions[indx]
            expression_mat = [expression_mat; expressions']
        end
        (expression_mat, element_order)
    end
    return expression_mat, element_order
end

function upload_clinical(input_filename)
    days_to_event, event, clinical_patient_order =
    open(input_filename) do infile
        days_to_event = Int[]
        event = Int[]
        clinical_patient_order = String[]
        for line in eachline(infile)
            line = split(strip(line), ",")
            append!(clinical_patient_order, line[1])
            append!(days_to_event, parse(Float64, line[2]))
            append!(event, parse(Float64, line[3]))
        end
        combined = sort(collect(zip(days_to_event, event, clinical_patient_order)), by=x->x[1])
        days_to_event = getindex.(combined, 1)
        event = getindex.(combined, 2)
        clinical_patient_order = getindex.(combined, 3)
        (days_to_event, event, clinical_patient_order)
    end
    return days_to_event, event, clinical_patient_order
end


"""
defining variables from ARGS
"""
function get_input()
    clinical_patient_order, days_to_event, event = upload_clinical(ARGS[1])

    expression_mat, element_order =
    upload_expression(ARGS[2], clinical_patient_order)

    min_threshold = parse(Float64, ARGS[3])
    max_threshold = parse(Float64, ARGS[4])
    null_size = parse(Float64, ARGS[5])

    output_filename = ARGS[6]

    return clinical_patient_order, days_to_event, event, expression_mat,
    element_order, min_threshold, max_threshold, null_size, output_filename
end
