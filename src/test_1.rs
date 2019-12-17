const USAGE: &str = "
Usage: test_1 [--dt N] [--plots]
       test_1 --help

An example template of simulation

Options:
    --dt N          Time step of the simulation [default: 5e-5]
    --tf N          Runtime of the simulation [default: 0.1]
    --plots         Show the plots
    -h, --help      Show this message.
";

// -------------------------------------------------

// std imports
use std::f64::consts::PI;
use std::fs;
use std::fs::OpenOptions;
use std::io::Write;

// external crate imports
use docopt::Docopt;
use gnuplot::*;
use multiphysics::prelude::*;

// -------------------------
// local imports
// -------------------------
// dem granular
pub use dem_rust::basic_equations::*;
pub use dem_rust::dem_3d::*;
pub use dem_rust::wall::*;
pub use dem_rust::prelude::*;

// rigid body imports
pub use dem_rust::rb::rb_2d::Rigidbody2D;

// external crate imports
// use gnuplot::*;
// use simple_shapes::prelude::*;

// -------------------------------------------------


#[derive(Deserialize, Debug)]
pub struct Args {
    flag_dt: f64,
    flag_plots: bool,
}


pub fn main(args: &[String]) {
    // --------------------------------------
    // GET THE COMMAND LINE VARIABLES
    // --------------------------------------
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.argv(args).deserialize())
        .unwrap_or_else(|e| e.exit());
    // println!("{:?}", args);
    // --------------------------------------
    // GET THE COMMAND LINE VARIABLES ENDS
    // --------------------------------------

    // --------------------------------------
    // Create particle arrays
    // --------------------------------------
    // Uncomment here to create particle arrays
    // --------------------------------------
    // Create particle arrays
    // --------------------------------------

    // ---------------------------------------
    // setup corresponding nbs nnps
    // ----------------------------------------
    // Uncomment here to create neighbours list data structure
    // ---------------------------------------
    // setup corresponding nbs nnps
    // ----------------------------------------

    // ----------------------------------------
    // solver data
    // ----------------------------------------
    let tf = args.tf;
    let dt = args.dt;
    let total_steps = (tf / dt) as u64;
    let mut t = 0.;
    let mut step_no = 0;
    let total_output_file = 1000;
    let pfreq = total_steps / total_output_file;
    // ----------------------------------------
    // solver data
    // ----------------------------------------

    // ----------------------------------------
    // output directory
    // ----------------------------------------
    let project_root = env!("CARGO_MANIFEST_DIR");
    let dir_name = project_root.to_owned() + "/test_1_output";
    let _p = fs::create_dir(&dir_name);
    // ----------------------------------------
    // output directory
    // ----------------------------------------

    // ----------------------------------------
    // some constants
    // ----------------------------------------
    let stage1 = 1;
    let stage2 = 2;
    // coefficient of friction
    let mu = 0.5;
    // ----------------------------------------
    // some constants
    // ----------------------------------------

    // create a progress bar
    let pb = setup_progress_bar(total_steps);

    while t < tf {
        // ----------------------
        // Rk2 initialize
        // ----------------------

        // ----------------------
        // stage 1 equations
        // ----------------------

        // ----------------------
        // stage 1 stepper
        // ----------------------

        // ----------------------
        // stage 2 equations
        // ----------------------

        // ----------------------
        // stage 2 stepper
        // ----------------------

        // ----------------------
        // output the vtk file
        // ----------------------
        if step_no % pfreq == 0 && step_no != 0{
            // write_to_vtk!(sand, format!("{}/sand_{}.vtk", &dir_name, step_no));
        }

        step_no += 1;
        t += dt;

        // progress bar increment
        pb.inc(1);
    }

    pb.finish_with_message("Simulation succesfully completed");


    // ---------------------------------------
    // Plotting
    // ---------------------------------------
    // UNCOMMENT AND USE THE PLOTTING FACILITY
    // let (incident_angle_experiment_kharaz, rebound_angle_experiment_kharaz) =
    //     read_xy_pairs(&format!(
    //         "{}/data/chung_test_4_incident_angle_vs_rebound_angle_experiment_kharaz.txt",
    //         &project_root
    //     ));

    // let mut fg = Figure::new();
    // fg.axes2d()
    //     .set_x_label("Incident angle (degree)", &[])
    //     .set_y_label("Rebound angle (degree)", &[])
    //     // .set_x_range(Fix(0.), Fix(90.))
    //     // .set_y_range(Fix(-800.), Fix(0.))
    //     .lines(
    //         &incident_angle_experiment_kharaz,
    //         &rebound_angle_experiment_kharaz,
    //         &[Caption("Kharaz experiment"), Color("black")],
    //     )
    //     .lines(
    //         &incident_angle_paper_simulated,
    //         &rebound_angle_paper_simulated,
    //         &[Caption("Paper simulated"), Color("blue")],
    //     )
    //     .points(
    //         &incident_angle,
    //         &rebound_angle_al_alloy,
    //         &[Caption("Al alloy"), Color("black")],
    //     )
    //     .points(
    //         &incident_angle,
    //         &rebound_angle_al_oxide,
    //         &[Caption("Al oxide"), Color("blue")],
    //     );

    // let px = 1000;
    // fg.save_to_png(
    //     &format!(
    //         "{}/chung_test_4_incident_angle_vs_rebound_angle.png",
    //         &dir_name
    //     ),
    //     px,
    //     px,
    // )
    // .unwrap();
    // if args.flag_plots {
    //     fg.show().unwrap();
    // }

}