use db::*;
use dotenv::dotenv;
use lib::*;
#[tokio::main]
async fn main() {
    dotenv().ok();

    let db = db::connect()
        .await
        .expect("Connection could not be established");
    import_genes(
        String::from("UM"),
        String::from("/mnt/nas/zhilin/others/UM_gene_anotation_extract_Arabidopsis.bed"),
        db,
    )
    .await
    .expect("Import failed");

    let dir = std::fs::read_dir("/mnt/nas/zhilin/others/constantin-sergio/biostress-data/")
        .expect("Directory must exist");
    // for file in dir {
    //     import_sites(
    //         sample,
    //         generation,
    //         layer,
    //         String::from(
    //             file.expect("File must be present")
    //                 .path()
    //                 .to_str()
    //                 .expect("File must be convertible to String"),
    //         ),
    //         db,
    //     )
    //     .await
    //     .expect("Import failed")
    // }
}
