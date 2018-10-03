using Distributions


"""
Simple Functions
"""
# calculate the stepwise e value for generating log-rank statistic
get_k(m1, m2, n1, n2) = m1 - n1 * (m1+m2) / (n1+n2)

# calculate variance from sum and sum of squares, and total
get_var(s::Float64, ss::Float64, total::Int16) = (ss - (s^2 / total)) / (total - 1)

# return ordered list of indices to change from group = 0 to group = 1
# also initializes group array
function get_ordered_indices(expression, min_threshold, max_threshold)
    group = zeros(length(expression))
    lowest_0_index = convert(Int16, floor(length(expression)*min_threshold))
    largest_threshold_index = convert(Int16, floor(expression)*max_threshold)
    p = sortperm(expression)[largest_0_index:largest_threshold_index]
    group[p] = 1
    p, group
end

"""
Survival Functions
"""
function lowest_logrank_p(days_to_event, event, expression,
    min_threshold, max_threshold)
    # days_to_event is sorted
    # indices of event correspond to indices of days_to_event
    # indices of expression correspond to indices of days_to_event
    # min/max threshold are for splitting KM curves

    # prepare for iterations
    p, group = get_ordered_indices(expression, min_threshold, max_threshold)
    # reverse p since pop is faster than shift
    p = flipdim(p, 1)
    lowest_pval = 1
    while !isempty(p)
        # move to the next threshold
        # find the index of the next lowest expression value
        group[pop!(p)] = 0  # set that spot to equal 0 (be in KM curve 'low')
        test_statistic = get_test_statistic(days_to_event, event, group)
        pval = ccdf(Chisq(1), test_statistic)
        lowest_pval = ifelse(pval < lowest_pval, pval, lowest_pval)

    end
end

function get_test_statistic(days_to_event, event, group)
    total = 0
    n1s = sum(group)
    n = [length(group)-n1s, n1s]
    m = [0,0]
    s = 0
    ss = 0
    prev_days = 0

    for i in 1:length(days_to_event)
        if e[i] == 0
            n[group[i]+1] -= 1
        else
            if days_to_event[i] == prev_days
                m[group[i]+1] += 1
                n[group[i]+1] -= 1
            else
                total += 1
                k = get_k(m[1], m[2], n[1], n[2])
                s += k
                s += k^2
                m[group[i]+1] = 1
                m[abs(group[i]-1)+1] = 0
                n[group[i]+1] -= 1
                prev_days = days_to_event[i]
            end
        end
    end
    return (s/total)^2 / get_var(s, ss, total)
end

function null_run(days_to_event, event, min_threshold,
    max_threshold)
    expression = rand(length(days_to_event))
    llp = lowest_logrank_p(days_to_event, event, expression, min_threshold,
    max_threshold)
    return llp
end
