use crate::*;

pub fn open_file(path: &PathBuf, filename: &OsString) -> Result<File> {
    let file = File::open(path).map_err(|_| {
        Error::File(
            String::from("methylome file"),
            String::from(filename.to_str().unwrap()),
        )
    })?;
    Ok(file)
}

pub fn lines_from_file(filename: &str) -> Result<io::Lines<io::BufReader<File>>> {
    let file = File::open(filename)
        .map_err(|_| Error::File(String::from("Annotation file"), String::from(filename)))?;
    Ok(io::BufReader::new(file).lines())
}

pub fn load_methylome(methylome: &str) -> Result<Vec<(PathBuf, OsString)>> {
    let methylome_dir = fs::read_dir(methylome)
        .map_err(|_| Error::File(String::from("Methylome directory"), String::from(methylome)))?;
    let methylome_files = methylome_dir
        .map(|f| (f.as_ref().unwrap().path(), f.unwrap().file_name()))
        .collect();
    Ok(methylome_files)
}
