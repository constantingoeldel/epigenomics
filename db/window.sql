select  c.*
from methylome c, annotations a
where c.strand = a.strand and c.chromosome = a.chromosome and c.location >= a.start and c.location <= a.end and 
	FLOOR((CASE WHEN c.strand = '+'::strandness then (c.location -a.start) else (a.end - c.location) end * 100) / (a.end - a.start))  = $1 and generation = $2 and layer =  $3
limit 100
