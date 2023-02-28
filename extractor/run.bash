# WT
# cargo run --release -- \
# -m ../../methylome/within_gbM_genes/ \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -o /mnt/extStorage/workingDir/constantin_not_owned_by_postgres/windows/wt \
# -e ../../methylome/edgelist.txt \
# -n ../../methylome/nodelist.txt \
# --alphabeta \
# --name wildtype \
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
#cargo run --release -- \
#-m /mnt/nas/zhilin/others/constantin-sergio/biostress-data \
#-g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
#-o ../../windows \
#-e /mnt/nas/zhilin/others/constantin-sergio/biostress-data/edgelist_ros_mock.tsv \
#-n /mnt/nas/zhilin/others/constantin-sergio/biostress-data/nodelist_ros_mock.tsv \
#--name ros \
#--alphabeta \
#-s 1 -w 5 

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
 cargo run --release -- \
 -m /mnt/nas/zhilin/others/constantin-sergio/biostress-data \
 -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
 -o /mnt/extStorage/workingDir/constantin_not_owned_by_postgres/windows/windows_col \
 -e /mnt/nas/zhilin/others/constantin-sergio/biostress-data/edgelist_col_mock.tsv \
 -n ../../methylome/nodelist_col_mock.tsv \
 --name col \
 -s 1 -w 5 \
--alphabeta

# mods
# cargo run --release -- \
# -m ../../modifications/bed \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -s 1 -w 1 \
# -o ../../windows \
# --force
# --cutoff-gene-length


# # chromatin states
# cargo run --release -- \
# -m ../../chr_states \
# -g ../../methylome/gbM_gene_anotation_extract_Arabidopsis.bed \
# -s 1 -w 1 \
# -o ../../windows \

# # chromatin states within red CS && gbM genes
# cargo run  -- \
# -m ../../chr_states \
# -g ../../methylome/redCS_SPMRs_in_gbM_genes.txt \
# -s 1 -w 1 \
# -o ../../windows \

# # chromatin states within red CS
# cargo run --release -- \
# -m ../../chr_states \
# -g ../../methylome/redCS-SPMRs.bed \
# -s 1 -w 1 \
# -o ../../windows \