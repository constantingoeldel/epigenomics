use std::fs;

use crate::*;

pub fn set_up_output_dir(max_gene_length: i32, args: Args) -> Result<()> {
    fs::read_dir(&args.output_dir).map_err(|_| {
        Error::File(
            String::from("Output directory"),
            String::from(&args.output_dir),
        )
    })?; // Throw error if base output dir does not exist
         // Replace existing content of output dir
    fs::remove_dir_all(&args.output_dir).unwrap();
    fs::create_dir(&args.output_dir).unwrap();

    let edgelist = String::from(
        "from	to	gendiff	group
G0	G1_A_MOCK	1	M1
G1_A_MOCK	G2_A_MOCK	1	M1
G1_A_MOCK	G2_A_MOCK_11	1	M1
G2_A_MOCK	G5_A_MOCK	1	M1
G5_A_MOCK	G6_A_MOCK	1	M1
G5_A_MOCK	G6_A_MOCK_56	1	M1
G6_A_MOCK	G10_A_MOCK	1	M1
G10_A_MOCK	G11_A_MOCK	1	M1
G10_A_MOCK	G11_A_MOCK_101	1	M1
G0	G1_C_MOCK	1	M2
G1_C_MOCK	G2_C_MOCK	1	M2
G1_C_MOCK	G2_C_MOCK_12	1	M2
G2_C_MOCK	G5_C_MOCK	1	M2
G5_C_MOCK	G6_C_MOCK	1	M2
G5_C_MOCK	G6_C_MOCK_57	1	M2
G6_C_MOCK	G10_C_MOCK	1	M2
G10_C_MOCK	G11_C_MOCK	1	M2
G10_C_MOCK	G11_C_MOCK_102	1	M2
G0	G1_D_MOCK	1	M3
G1_D_MOCK	G2_D_MOCK	1	M3
G1_D_MOCK	G2_D_MOCK_13	1	M3
G2_D_MOCK	G5_D_MOCK	1	M3
G5_D_MOCK	G6_D_MOCK	1	M3
G5_D_MOCK	G6_D_MOCK_58	1	M3
G6_D_MOCK	G10_D_MOCK	1	M3
G10_D_MOCK	G11_D_MOCK	1	M3
G10_D_MOCK	G11_D_MOCK_103	1	M3
G0	G1_E_MOCK	1	M4
G1_E_MOCK	G2_E_MOCK	1	M4
G1_E_MOCK	G2_E_MOCK_14	1	M4
G2_E_MOCK	G5_E_MOCK	1	M4
G5_E_MOCK	G6_E_MOCK	1	M4
G5_E_MOCK	G6_E_MOCK_59	1	M4
G6_E_MOCK	G10_E_MOCK	1	M4
G10_E_MOCK	G11_E_MOCK	1	M4
G10_E_MOCK	G11_E_MOCK_104	1	M4
",
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
                "filename	node	gen	meth
-	G0	0	N
-	G1_A_MOCK	1	N
-	G2_A_MOCK	2	N
/home/constantin/windows/{side}/{window}/methylome_11_All.txt	G2_A_MOCK_11	2	Y
-	G5_A_MOCK	5	N
-	G6_A_MOCK	6	N
/home/constantin/windows/{side}/{window}/methylome_56_All.txt	G6_A_MOCK_56	6	Y
-	G10_A_MOCK	10	N
-	G11_A_MOCK	11	N
/home/constantin/windows/{side}/{window}/methylome_101_All.txt	G11_A_MOCK_101	11	Y
-	G1_C_MOCK	1	N
-	G2_C_MOCK	2	N
/home/constantin/windows/{side}/{window}/methylome_12_All.txt	G2_C_MOCK_12	2	Y
-	G5_C_MOCK	5	N
-	G6_C_MOCK	6	N
/home/constantin/windows/{side}/{window}/methylome_57_All.txt	G6_C_MOCK_57	6	Y
-	G10_C_MOCK	10	N
-	G11_C_MOCK	11	N
/home/constantin/windows/{side}/{window}/methylome_102_All.txt	G11_C_MOCK_102	11	Y
-	G1_D_MOCK	1	N
-	G2_D_MOCK	2	N
/home/constantin/windows/{side}/{window}/methylome_13_All.txt	G2_D_MOCK_13	2	Y
-	G5_D_MOCK	5	N
-	G6_D_MOCK	6	N
/home/constantin/windows/{side}/{window}/methylome_58_All.txt	G6_D_MOCK_58	6	Y
-	G10_D_MOCK	10	N
-	G11_D_MOCK	11	N
/home/constantin/windows/{side}/{window}/methylome_103_All.txt	G11_D_MOCK_103	11	Y
-	G1_E_MOCK	1	N
-	G2_E_MOCK	2	N
/home/constantin/windows/{side}/{window}/methylome_14_All.txt	G2_E_MOCK_14	2	Y
-	G5_E_MOCK	5	N
-	G6_E_MOCK	6	N
/home/constantin/windows/{side}/{window}/methylome_59_All.txt	G6_E_MOCK_59	6	Y
-	G10_E_MOCK	10	N
-	G11_E_MOCK	11	N
/home/constantin/windows/{side}/{window}/methylome_104_All.txt	G11_E_MOCK_104	11	Y

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
