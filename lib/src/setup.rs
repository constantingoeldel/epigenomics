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
G1_A_MOCK	G2_A_MOCK_6	1	M1
G2_A_MOCK	G5_A_MOCK	1	M1
G5_A_MOCK	G6_A_MOCK	1	M1
G5_A_MOCK	G6_A_MOCK_51	1	M1
G6_A_MOCK	G10_A_MOCK	1	M1
G10_A_MOCK	G11_A_MOCK	1	M1
G10_A_MOCK	G11_A_MOCK_96	1	M1
G0	G1_C_MOCK	1	M2
G1_C_MOCK	G2_C_MOCK	1	M2
G1_C_MOCK	G2_C_MOCK_7	1	M2
G2_C_MOCK	G5_C_MOCK	1	M2
G5_C_MOCK	G6_C_MOCK	1	M2
G5_C_MOCK	G6_C_MOCK_52	1	M2
G6_C_MOCK	G10_C_MOCK	1	M2
G10_C_MOCK	G11_C_MOCK	1	M2
G10_C_MOCK	G11_C_MOCK_97	1	M2
G0	G1_D_MOCK	1	M3
G1_D_MOCK	G2_D_MOCK	1	M3
G1_D_MOCK	G2_D_MOCK_8	1	M3
G2_D_MOCK	G5_D_MOCK	1	M3
G5_D_MOCK	G6_D_MOCK	1	M3
G5_D_MOCK	G6_D_MOCK_53	1	M3
G6_D_MOCK	G10_D_MOCK	1	M3
G10_D_MOCK	G11_D_MOCK	1	M3
G10_D_MOCK	G11_D_MOCK_98	1	M3
G0	G1_E_MOCK	1	M4
G1_E_MOCK	G2_E_MOCK	1	M4
G1_E_MOCK	G2_E_MOCK_9	1	M4
G2_E_MOCK	G5_E_MOCK	1	M4
G5_E_MOCK	G6_E_MOCK	1	M4
G5_E_MOCK	G6_E_MOCK_54	1	M4
G6_E_MOCK	G10_E_MOCK	1	M4
G10_E_MOCK	G11_E_MOCK	1	M4
G10_E_MOCK	G11_E_MOCK_99	1	M4
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
/home/constantin/windows/{side}/{window}/methylome_6_All.txt	G2_A_MOCK_6	2	Y
-	G5_A_MOCK	5	N
-	G6_A_MOCK	6	N
/home/constantin/windows/{side}/{window}/methylome_61_All.txt	G6_A_MOCK_51	6	Y
-	G10_A_MOCK	10	N
-	G11_A_MOCK	11	N
/home/constantin/windows/{side}/{window}/methylome_96_All.txt	G11_A_MOCK_96	11	Y
-	G1_C_MOCK	1	N
-	G2_C_MOCK	2	N
/home/constantin/windows/{side}/{window}/methylome_7_All.txt	G2_C_MOCK_7	2	Y
-	G5_C_MOCK	5	N
-	G6_C_MOCK	6	N
/home/constantin/windows/{side}/{window}/methylome_52_All.txt	G6_C_MOCK_52	6	Y
-	G10_C_MOCK	10	N
-	G11_C_MOCK	11	N
/home/constantin/windows/{side}/{window}/methylome_97_All.txt	G11_C_MOCK_97	11	Y
-	G1_D_MOCK	1	N
-	G2_D_MOCK	2	N
/home/constantin/windows/{side}/{window}/methylome_8_All.txt	G2_D_MOCK_8	2	Y
-	G5_D_MOCK	5	N
-	G6_D_MOCK	6	N
/home/constantin/windows/{side}/{window}/methylome_53_All.txt	G6_D_MOCK_53	6	Y
-	G10_D_MOCK	10	N
-	G11_D_MOCK	11	N
/home/constantin/windows/{side}/{window}/methylome_98_All.txt	G11_D_MOCK_98	11	Y
-	G1_E_MOCK	1	N
-	G2_E_MOCK	2	N
/home/constantin/windows/{side}/{window}/methylome_9_All.txt	G2_E_MOCK_9	2	Y
-	G5_E_MOCK	5	N
-	G6_E_MOCK	6	N
/home/constantin/windows/{side}/{window}/methylome_54_All.txt	G6_E_MOCK_5	6	Y
-	G10_E_MOCK	10	N
-	G11_E_MOCK	11	N
/home/constantin/windows/{side}/{window}/methylome_99_All.txt	G11_E_MOCK_104	11	Y

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
