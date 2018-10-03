include("NEEPS.jl")
using NEEPS

clinical_patient_order, days_to_event, event, expression_mat,
element_order, min_threshold, max_threshold, null_size, output_filename =
get_input()

null_ps, lowest_pvals = parallel_null_and_curves(null_size, days_to_event,
event, min_threshold, max_threshold, expression_mat)

ordered_neep_adj_pvals = generate_neep_all(null_ps, lowest_pvals)

export_to_file(clinical_patient_order, ordered_neep_adj_pvals, output_filename)
