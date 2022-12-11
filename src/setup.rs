use std::fs;

use crate::*;

pub fn set_up_output_dir(max_gene_length: i32, args: Args) -> Result<()> {
    fs::read_dir(&args.output_dir).map_err(|_| {
        Error::File(
            String::from("Output directory"),
            String::from(&args.output_dir),
        )
    })?; // Throw error if base output dir does not exist

    if args.force {
        fs::remove_dir_all(&args.output_dir).unwrap();
        fs::create_dir(&args.output_dir).unwrap();
    }
    let edgelist = String::from(
        "from to
0_0 1_2
1_2 2_2
2_2 3_2
3_2 4_2
4_2 5_2
5_2 6_2
6_2 7_2
7_2 8_2
8_2 9_2
9_2 10_2
10_2 11_2
0_0 1_8
1_8 2_8
2_8 3_8
3_8 4_8
4_8 5_8
5_8 6_8
6_8 7_8
7_8 8_8
8_8 9_8
9_8 10_8
10_8 11_8",
    );

    let sides = vec![
        ("upstream", args.cutoff),
        ("gene", max_gene_length),
        ("downstream", args.cutoff),
    ];

    for side in sides {
        let max = if args.absolute { side.1 } else { 100 };
        let side = side.0;

        for window in (0..=max).step_by(args.window_step as usize) {
            let nodelist = format!(
                "filename,node,gen,meth
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G0_All.txt,0_0,0,Y
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G1_L2_All.txt,1_2,1,Y
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G1_L8_All.txt,1_8,1,Y
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G2_L2_All.txt,2_2,2,Y
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G2_L8_All.txt,2_8,2,Y
-,3_2,3,N
-,3_8,3,N
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G4_L2_All.txt,4_2,4,Y
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G4_L8_All.txt,4_8,4,Y
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G5_L2_All.txt,5_2,5,Y
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G5_L8_All.txt,5_8,5,Y
-,6_2,6,N
-,6_8,6,N
-,7_2,7,N
-,7_8,7,N
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G8_L2_All.txt,8_2,8,Y
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G8_L8_All.txt,8_8,8,Y
-,9_2,9,N
-,9_8,9,N
-,10_2,10,N
-,10_8,10,N
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G11_L2_All.txt,11_2,11,Y
/mnt/extStorage/constantin/windows/{side}/{window}/methylome_Col0_G11_L8_All.txt,11_8,11,Y
"
            );

            let path = format!("{}/{}/{}", args.output_dir, side, window);
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
