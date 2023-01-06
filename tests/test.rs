use extractor::arguments::Args;
use extractor::extract;
use serial_test::serial;

// #[test]
// #[serial]
// fn run_relative() {
//     let args = Args {
//         methylome: "/home/constantin/methylome/within_gbM_genes".to_string(),
//         genome: "/home/constantin/methylome/gbM_gene_anotation_extract_Arabidopsis.bed".to_string(),
//         window_size: 5,
//         window_step: 1,

//         output_dir: "/home/constantin/windows".to_string(),
//         absolute: false,
//         cutoff: 2048,
//     };
//     match extract(args) {
//         Ok(_) => println!("Done!"),
//         Err(e) => panic!("Error: {}", e),
//     }
// }
// #[test]
// #[serial]
// fn run_absolute() {
//     let args = Args {
//         methylome: "/home/constantin/methylome/within_gbM_genes".to_string(),
//         genome: "/home/constantin/methylome/gbM_gene_anotation_extract_Arabidopsis.bed".to_string(),
//         window_size: 512,
//         window_step: 256,

//         output_dir: "/home/constantin/windows".to_string(),
//         absolute: true,
//         cutoff: 2048,
//     };
//     match extract(args) {
//         Ok(_) => println!("Done!"),
//         Err(e) => panic!("Error: {}", e),
//     }
// }
