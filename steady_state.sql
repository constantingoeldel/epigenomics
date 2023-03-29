-- length: end - start
-- midpoint: start + (end-start) / 2
-- distance from 5' of sense strand: mod.start - gene.start
-- distance from 3' of sense strand: gene.end - mod.end
-- distance from 5' cutoff: mod.start - (gene.start - (gene.end - gene.start)
-- distance from 3' cutoff: (gene.end + (gene.end - gene.start) - mod.end 
-- inverted for antisense strand
select Count(*) as sites,
    AVG(h.meth_lvl) as avg_meth_lvl,
    AVG(
        h.count_methylated * 1.0 / (h.count_total + 0.0001)
    ) as avg_meth_sites,
    AVG(posteriormax) as avg_posteriormax,
    h.generation,
    h.line
from methylome h
    join sites_annotation s on h.id = s.site
    join annotations a on s.annotation = a.id
where a.annotation = 'gbM'
    and h.sample = 'MA3'
group by h.generation,
    h.line