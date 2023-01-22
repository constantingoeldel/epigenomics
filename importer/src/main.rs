use db::*;
use dotenv::dotenv;
#[tokio::main]
async fn main() {
    dotenv().ok();

    let db = db::connect()
        .await
        .expect("Connection could not be established");
    // import_genes(
    //     String::from("UM"),
    //     String::from("/mnt/nas/zhilin/others/UM_gene_anotation_extract_Arabidopsis.bed"),
    //     db,
    // )
    // .await
    // .expect("Import failed");

    let files = [
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_6.txt.change.bed",
            2,
            1,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_51.txt.change.bed",
            6,
            1,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_96.txt.change.bed",
            11,
            1,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_7.txt.change.bed",
            2,
            3,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_52.txt.change.bed",
            6,
            3,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_97.txt.change.bed",
            11,
            3,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_8.txt.change.bed",
            2,
            4,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_53.txt.change.bed",
            6,
            4,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_98.txt.change.bed",
            11,
            4,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_9.txt.change.bed",
            2,
            5,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_54.txt.change.bed",
            6,
            5,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_99.txt.change.bed",
            11,
            5,
            "ros",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_11.txt.change.bed",
            2,
            1,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_56.txt.change.bed",
            6,
            1,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_101.txt.change.bed",
            11,
            1,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_12.txt.change.bed",
            2,
            3,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_57.txt.change.bed",
            6,
            3,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_102.txt.change.bed",
            11,
            3,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_13.txt.change.bed",
            2,
            4,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_58.txt.change.bed",
            6,
            4,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_103.txt.change.bed",
            11,
            4,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_14.txt.change.bed",
            2,
            5,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_59.txt.change.bed",
            6,
            5,
            "nrpe",
        ),
        (
            "/mnt/nas/zhilin/MA-lines/CG/methylome_CG_bed/methylome_104.txt.change.bed",
            11,
            5,
            "nrpe",
        ),
    ];
    for (file, generation, layer, sample) in files {
        import_sites(sample, generation, layer, file, &db)
            .await
            .expect("Import failed")
    }
}
