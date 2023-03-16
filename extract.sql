-- length: end - start
-- midpoint: start + (end-start) / 2
-- distance from 5' of sense strand: mod.start - gene.start
-- distance from 3' of sense strand: gene.end - mod.end
-- distance from 5' cutoff: mod.start - (gene.start - (gene.end - gene.start)
-- distance from 3' cutoff: (gene.end + (gene.end - gene.start) - mod.end 
-- inverted for antisense strand
with percentage_windows as (
    select Count(*) as sites,
        h.type as mod,
        FLOOR(
            1000 * (
                CASE
                    WHEN a.strand = '+' then (h.start - a.start)
                    else (a.end - h.
                end
            )
        end
) * 1.0 / (a.end - a.start
)
) as percentage
from annotations a,
    mods h
where h.chromosome = a.chromosome
    and h.start >= a.start - (a.end - a.start
)
and h.
end <= a.
end + (a.end - a.start
)
group by percentage,
    h.type
)
select sites,
    --	SUM(sites) OVER (ORDER BY percentage ROWS BETWEEN CURRENT ROW AND 4 FOLLOWING) AS sites_sliding_sum,
    Round(percentage / 10, 1) as percentile,
    mod
from percentage_windows