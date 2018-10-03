function export_to_file(element_order, ed_neep_adj_pvals, output_filename)
    open(output_filename, "w") do outfile
        for (element, pval) in zip(element_order, ordered_neep_adj_pvals)
            write(outfile, "$element\t$pval\n")
        end
    end
end
