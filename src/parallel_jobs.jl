function parallel_null_and_curves(null_size, days_to_event, event, min_threshold,
    max_threshold, expression_mat)
    # moving existing variables to every process
    @everywhere null_size = null_size
    @everywhere days_to_event = days_to_event
    @everywhere event = event
    @everywhere min_threshold = min_threshold
    @everywhere max_threshold = max_threshold
    @everywhere expression_mat = expression_mat

    @everywhere include("survival_log_rank_pvals.jl")
    include("transformation.jl")

    # create a channel with the jobs for the null distribution
    @sync begin
        null_ps = zeros(null_size)
        for i in 1:null_size
            @spawn null_ps[i] = null_run(days_to_event, event, min_threshold, max_threshold)
        end

        lowest_pvals = lowest_pvals = zeros(size(expression_mat)[1])
        for i in 1:size(expression_mat)[1]
            @spawn lowest_pvals[i] = lowest_logrank_p(days_to_event, event, expression_mat[i], min_threshold, max_threshold)
        end
    end

    return null_ps, lowest_pvals
end
