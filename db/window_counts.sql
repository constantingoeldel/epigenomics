with percentage_windows as (
select  
	Count(*) as sites, 
	FLOOR((CASE WHEN c.strand then (c.location -a.start) else (a.end - c.location) end * 100) / (a.end - a.start))  as percentile ,
	c.generation,
	c.layer 
from cg_within_gbM_genes c, annotation a
where c.strand = a.strand and c.chromosome = a.chromosome and c.location >= a.start and c.location <= a.end -- and c.status = 'I'
group by  c.generation, c.layer,percentile
limit 100
)

select 
	*, 
	SUM(sites) OVER (ORDER BY percentile ROWS BETWEEN CURRENT ROW AND 4 FOLLOWING) AS sites_sliding_sum
from percentage_windows

