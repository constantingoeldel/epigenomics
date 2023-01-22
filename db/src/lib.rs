use lib::{
    arguments::Args,
    files::lines_from_file,
    methylation_site::MethylationSite,
    structs::{Gene, Result, Strand},
    windows::Windows,
};
use sqlx::{postgres::PgPoolOptions, Pool, Postgres};
use std::{collections::HashMap, env};

pub async fn connect() -> Result<Pool<Postgres>> {
    let pool = PgPoolOptions::new()
        .max_connections(50)
        .connect(
            &env::var("POSTGRES_CONNECTION_STRING")
                .expect("Environment variable POSTGRES_CONNECTION_STRING must be set"),
        )
        .await?;
    Ok(pool)
}

pub async fn windows(args: Args, db: Pool<Postgres>) -> Result<HashMap<String, Windows>> {
    let mut results: HashMap<String, Windows> = HashMap::new();
    for generation in [0, 1, 2, 4, 5, 8, 11] {
        for layer in [0, 2, 8] {
            if (generation != 0 && layer == 0) || (generation == 0 && layer != 0) {
                continue;
            }
            for percentile in 0..100 {
                let filename = format!("{}", &args.output_dir)
                    + &format!("/{percentile}/methylome_Col0_G{generation}_L{layer}_All.txt");

                let mut windows = Windows::new(100, &args);
                // let result: Vec<MethylationSite> = sqlx::query_file_as!(
                //     MethylationSite,
                //     "./window.sql",
                //     percentile,
                //     generation,
                //     layer
                // )
                // .fetch_all(&db)
                // .await?;
                // results.insert(filename, windows);
            }
        }
    }
    Ok(results)
}

pub async fn import_genes(note: String, filename: String, db: Pool<Postgres>) -> Result<()> {
    let lines = lines_from_file(&filename)?;
    let mut rows_affected = 0;
    for line in lines {
        if line.is_err() {
            continue;
        }
        let annotation = Gene::from_annotation_file_line(
            &line.expect("Line should not be an error after the check"),
            false,
        );
        if let Some(a) = annotation {
            rows_affected += sqlx::query!(
                r#"INSERT INTO annotations (chromosome, start, "end", name, strand, annotation) VALUES($1,$2,$3,$4,$5, $6)"#,
                a.chromosome as i32,
                a.start,
                a.end,
                a.name,
                a.strand as Strand,
                note
            )
            .execute(&db)
            .await?.rows_affected();
        }
    }
    println!("{rows_affected} annotations were anned to the database!");
    Ok(())
}

pub async fn import_sites(
    sample: &str,
    generation: i32,
    layer: i32,
    filename: &str,
    db: &Pool<Postgres>,
) -> Result<()> {
    let lines = lines_from_file(&filename)?;
    let mut rows_affected = 0;
    for line in lines {
        if line.is_err() {
            continue;
        }
        let site = MethylationSite::from_methylome_file_line(
            &line
                .as_ref()
                .expect("Line should not be an error after the check"),
            false,
        );
        match site {
            Err(e) => println!("{e} when parsing: {}", line.unwrap()),
            Ok(s) => {
                rows_affected += sqlx::query!(
                    r#"INSERT INTO methylome (chromosome, location, context, count_methylated, strand, count_total, posteriormax, status, meth_lvl,generation, layer, sample) VALUES($1,$2,$3,$4,$5,$6, $7,$8,$9, $10, $11, $12)"#,
                    s.chromosome as i32,
                    s.location,
                    s.context,
                    s.count_methylated as i32,
                    s.strand as Strand,
                    s.count_total as i32,
                    s.posteriormax,
                    s.status as i8,
                    s.meth_lvl,
                    generation,
                    layer,
                    sample
                )
                .execute(db)
                .await?.rows_affected();
            }
        }
    }
    println!("{rows_affected} annotations were anned to the database!");
    Ok(())
}
