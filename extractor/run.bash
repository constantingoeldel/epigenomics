# # WT
# cargo run --release -- \
# -m ../../methylome/within_gbM_genes/ \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -o ../../windows \
# -e ../../methylome/edgelist.txt \
# -n ../../methylome/nodelist.txt \
# -s 1 -w 5 

# # CMT3
# cargo run --release -- \
# -m /mnt/nas/zhilin/others/constantin-sergio/CMT3/total_original_methylome \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -o ../../windows \
# -e /home/constantin/methylome/cmt3_edgelist.txt \
# -n /home/constantin/methylome/cmt3_nodelist.txt \
# --name cmt3 \
# -s 1 -w 5 

# # suv 4/5/6
# cargo run --release -- \
# -m /mnt/nas/zhilin/others/constantin-sergio/SUV456/total_original_methylome \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -o ../../windows \
# -e /mnt/nas/zhilin/others/constantin-sergio/SUV456/SUV456_edgelist.txt \
# -n /mnt/nas/zhilin/others/constantin-sergio/SUV456/SUV456_nodelist.txt \
# --name suv \
# -s 1 -w 5 

# ros 
cargo run --release -- \
-m /mnt/nas/zhilin/others/constantin-sergio/biostress-data \
-g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
-o ../../windows \
-e /mnt/nas/zhilin/others/constantin-sergio/biostress-data/edgelist_ros_mock.tsv \
-n /mnt/nas/zhilin/others/constantin-sergio/biostress-data/nodelist_ros_mock.tsv \
--name ros \
--alphabeta \
-s 1 -w 5 

# # nrpe 
# cargo run --release -- \
# -m /mnt/nas/zhilin/others/constantin-sergio/biostress-data \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -o ../../windows \
# -e /mnt/nas/zhilin/others/constantin-sergio/biostress-data/edgelist_nrpe_mock.tsv \
# -n /mnt/nas/zhilin/others/constantin-sergio/biostress-data/nodelist_nrpe_mock.tsv \
# --name nrpe \
# -s 1 -w 5 


# # col 
# cargo run --release -- \
# -m /mnt/nas/zhilin/others/constantin-sergio/biostress-data \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -o ../../windows \
# -e /mnt/nas/zhilin/others/constantin-sergio/biostress-data/edgelist_col_mock.tsv \
# -n /mnt/nas/zhilin/others/constantin-sergio/biostress-data/nodelist_col_mock.tsv \
# --name col \
# -s 1 -w 5 

