select * from methylome m, annotations a, sites_annotation t
where m.id = t.site and a.id = t.annotation
