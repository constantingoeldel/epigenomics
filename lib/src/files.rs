use std::{ffi::OsStr, path::Path};

use crate::*;

pub fn open_file(path: &PathBuf) -> Result<File> {
    let file = File::open(path).or(Err(Error::File(path.to_owned())))?;
    Ok(file)
}

pub fn lines_from_file(path: &PathBuf) -> Result<io::Lines<io::BufReader<File>>> {
    let file = File::open(path).or(Err(Error::File(path.to_owned())))?;
    Ok(io::BufReader::new(file).lines())
}

pub fn load_methylome(methylome: &PathBuf) -> Result<Vec<PathBuf>> {
    let methylome_dir = fs::read_dir(methylome).or(Err(Error::File(methylome.to_owned())))?;
    let methylome_files = methylome_dir
        .map(|f| f.as_ref().unwrap().path())
        .filter(|path| !path.extension().contains(&"tsv") && !path.extension().contains(&"fn")) // Filter out tsv and fn files, which are often nodelist/edgelist files.
        .collect();
    Ok(methylome_files)
}
/// Highly sacriligous
pub fn file_name(path: &Path) -> String {
    let name = path.components().nth_back(0).unwrap();
    name.as_os_str().to_string_lossy().into()
}
