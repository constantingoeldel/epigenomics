use crate::*;

pub fn set_up_output_dir(args: Args, max_gene_length: u32) -> Result<()> {
    let output_dir = PathBuf::from(&args.output_dir).canonicalize()?;
    fs::read_dir(&output_dir).map_err(|_| {
        Error::File(
            String::from("Output directory"),
            String::from(&args.output_dir),
        )
    })?; // Throw error if base output dir does not exist
         // Replace existing content of output dir
    fs::remove_dir_all(&args.output_dir).unwrap();
    fs::create_dir(&args.output_dir).unwrap();

    let edgelist = fs::read_to_string(args.edges)?;
    let nodes = fs::read_to_string(args.nodes)?;

    let sides = vec![
        ("upstream", args.cutoff),
        ("gene", max_gene_length),
        ("downstream", args.cutoff),
    ];

    for side in sides {
        let max = if args.absolute { side.1 } else { 100 };
        let side = side.0;

        for window in (0..=max).step_by(args.window_step as usize) {
            let mut nodelist = String::new();
            let lines = nodes.split('\n');
            for line in lines {
                if line.starts_with('/') {
                    let old_file = line.split('\t').next().unwrap();
                    let filename = old_file.split('/').last().unwrap();
                    let file = format!(
                        "{}/{}/{}/{}",
                        &output_dir.to_string_lossy(),
                        side,
                        window,
                        filename
                    );
                    nodelist += &line.replace(old_file, &file);
                } else {
                    nodelist += line;
                }
                nodelist += "\n";
            }

            let path = format!("{}/{}/{}", &output_dir.to_string_lossy(), side, window);
            let window_dir = fs::read_dir(&path);

            match window_dir {
                Ok(_) => continue,
                Err(_) => {
                    fs::create_dir_all(&path).unwrap();
                    fs::write(path.to_owned() + "/nodelist.fn", nodelist)
                        .expect("Nodelist not writable at ");
                    fs::write(path.to_owned() + "/edgelist.fn", &edgelist).expect("msg");
                }
            };
        }
    }
    Ok(())
}
