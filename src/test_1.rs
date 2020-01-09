const USAGE: &str = "
Usage: test_1 [--dt N] [--plots]
       test_1 --help

A template for DEM spherical particles modelling.


This benchmark is taken from

author: XXXX
paper: XXXX
link: XXXX

Description:
-----------

Describe the test in a few lines.


Method:
--------

What specific DEM model is used to model the current simulation.


Input:
------------

========Any inputs==============


Analysis:
-----------

How is this benchmark validated? What plots are used? What conclusions are been
made?



New features of dem_rust:
-------------------------

Any new features implemented in this current example, which can be later used in
another examples?


Options:
    --dt N          Time step of the simulation [default: 1e-4]
    --tf N          Runtime of the simulation [default: 1.]
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

// // rigid body imports
// pub use dem_rust::rb::rb_2d::Rigidbody2D;

// for reading data from file (comma separated)
use crate::read_xy_pairs;

// external crate imports
// use gnuplot::*;
// use simple_shapes::prelude::*;

// -------------------------------------------------


#[derive(Deserialize, Debug)]
pub struct Args {
    flag_tf: f64,
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
    // CREATE PARTICLE ARRAYS
    // --------------------------------------
    let x = vec![0.];
    let y = vec![50. * 1e-2];
    let z = vec![0.; x.len()];
    let radius = vec![10. * 1e-2; x.len()];
    let mut sand = DEM3DParticleArray::from_xyz_rad(&x, &y, &z, &radius);
    let rho = 2600.;
    sand.rho = vec![rho; x.len()];
    sand.m = vec![rho * PI * radius[0] * radius[0]; x.len()];
    sand.m_inv = vec![1. / sand.m[0]; x.len()];

    let inertia = 4. * (2. * radius[0]) * (2. * radius[0]) / 10.;
    sand.mi = vec![inertia; x.len()];
    sand.mi_inv = vec![1. / sand.mi[0]; x.len()];

    let stiffness = 5. * 1e6;
    sand.kn = vec![stiffness; x.len()];
    sand.kt = vec![stiffness; x.len()];
    // set some dummy Young's modulus for linear DEM case; change these values if
    // you are updating this example for nonlinear DEM
    let yng = 1.;
    let nu = 0.2;
    let shr = yng / (2. * (1. + nu));
    sand.young_mod = vec![yng; x.len()];
    sand.poisson = vec![nu; x.len()];
    sand.shear_mod = vec![shr; x.len()];

    // this is nice step for debugging
    sand.validate_particle_array();

    // -------------------------
    // create an infinite wall
    // -------------------------
    let x = vec![0.];
    let y = vec![0.];
    let z = vec![0.; x.len()];
    let wall_points = vec![[Vector3::new(-1., 0., 0.), Vector3::new(1., 0., 0.)]];

    let nx = vec![0.];
    let ny = vec![1.];
    let nz = vec![0.; x.len()];
    let mut wall =
        WallDEMParticleArray::from_xyz_wall_points_normals(&x, &y, &z, &wall_points, &nx, &ny, &nz);
    wall.kn = vec![stiffness; x.len()];
    wall.kt = vec![stiffness; x.len()];
    // --------------------------------------
    // CREATE PARTICLE ARRAYS ENDS
    // --------------------------------------

    // ---------------------------------------
    // SETUP CORRESPONDING NBS NNPS
    // ----------------------------------------
    let max_coordinate = 50. * 1e-2;
    let max_size = 2. * radius[0];
    let mut nbs2d_sand =
        NBS2D::from_maximum_and_no_of_particles(max_coordinate, max_size, sand.x.len());
    // ---------------------------------------
    // SETUP CORRESPONDING NBS NNPS ENDS
    // ----------------------------------------

    // ----------------------------------------
    // SOLVER DATA
    // ----------------------------------------
    let tf = args.flag_tf;
    let dt = args.flag_dt;
    let total_steps = (tf / dt) as u64;
    let mut t = 0.;
    let mut step_no = 0;
    let total_output_file = 1000;
    let pfreq = if total_steps < total_output_file {
        1
    } else {
        total_steps / total_output_file
    };
    // ----------------------------------------
    // SOLVER DATA ENDS
    // ----------------------------------------

    // ----------------------------------------
    // OUTPUT DIRECTORY
    // ----------------------------------------
    let project_root = env!("CARGO_MANIFEST_DIR");
    let dir_name = project_root.to_owned() + "/test_1_output";
    let _p = fs::create_dir(&dir_name);
    // ----------------------------------------
    // OUTPUT DIRECTORY ENDS
    // ----------------------------------------

    // ----------------------------------------
    // SOME CONSTANTS
    // ----------------------------------------
    let stage1 = 1;
    let stage2 = 2;
    // coefficient of friction
    let mu = 0.5;
    let kn = 5. * 1e5;
    let kt = 5. * 1e5;

    let en = 0.9;
    // ----------------------------------------
    // SOME CONSTANTS ENDS
    // ----------------------------------------

    // create a progress bar
    let pb = setup_progress_bar(total_steps);

    // -------------------------------
    // FOR PARAVIEW VISUALIZATION
    // -------------------------------
    // define variables for automatic visualization

    // particle array files
    let mut sand_files = "".to_string();
    write_to_vtk!(sand, format!("{}/sand_{}.vtk", &dir_name, step_no));
    sand_files.push_str(&format!("'{}/sand_{}.vtk', ", &dir_name, step_no));

    let mut wall_files = "".to_string();
    write_wall_to_vtk!(wall, format!("{}/wall_{}.vtk", &dir_name, step_no));
    wall_files.push_str(&format!("'{}/wall_{}.vtk', ", &dir_name, step_no));

    // neighbours files
    let mut nnps_files = "".to_string();
    write_nnps_2d_to_vtk!(nbs2d_sand, format!("{}/nnps.vtk", &dir_name));
    nnps_files.push_str(&format!("'{}/nnps.vtk'", &dir_name));
    // -------------------------------
    // FOR PARAVIEW VISUALIZATION ENDS
    // -------------------------------

    // -------------------------------
    // POST PROCESSING VARIABLES
    // -------------------------------
    // uncomment and initialize the variables
    // -------------------------------
    // POST PROCESSING VARIABLES ENDS
    // -------------------------------

    while t < tf {
        // ----------------------
        // RK2 INITIALIZE
        // ----------------------
        rk2_initialize_copy_contacts_3d!((sand));
        rk2_initialize_dem!((sand));

        // ----------------------
        // STAGE 1 EQUATIONS
        // ----------------------
        nbs2d_sand.register_particles_to_nnps(&sand.x, &sand.y, &sand.z);
        make_forces_torques_zero!((sand));
        body_force!((sand), 0., -9.81, 0.);
        // dem_3d_force_linear_pp!(sand, (sand), (nbs2d_sand), (kn), (kt), (en), (mu), stage1, dt);
        dem_3d_force_linear_pw!(sand, (wall), (kn), (kt), (en), (mu), stage1, dt);

        // ----------------------
        // STAGE 1 STEPPER
        // ----------------------
        rk2_stage1_dem!((sand), dt);

        // ----------------------
        // STAGE 2 EQUATIONS
        // ----------------------
        nbs2d_sand.register_particles_to_nnps(&sand.x, &sand.y, &sand.z);
        make_forces_torques_zero!((sand));
        body_force!((sand), 0., -9.81, 0.);
        // dem_3d_force_linear_pp!(sand, (sand), (nbs2d_sand), (kn), (kt), (en), (mu), stage2, dt);
        dem_3d_force_linear_pw!(sand, (wall), (kn), (kt), (en), (mu), stage2, dt);

        // ----------------------
        // STAGE 2 STEPPER
        // ----------------------
        rk2_stage2_dem!((sand), dt);

        update_contacts_pp!(sand, (sand));
        // ----------------------
        // OUTPUT THE VTK FILE AND WRITE ALL FILES DATA
        // ----------------------
        if step_no % pfreq == 0 && step_no != 0 {
            // ----------------------
            // WRITE DATA TO VTK FILE
            // ----------------------
            write_to_vtk!(sand, format!("{}/sand_{}.vtk", &dir_name, step_no));
            write_wall_to_vtk!(wall, format!("{}/wall_{}.vtk", &dir_name, step_no));

            // ----------------------
            // FOR PARAVIEW AUTOMATION
            // ----------------------
            sand_files.push_str(&format!("'{}/sand_{}.vtk', ", &dir_name, step_no));
            wall_files.push_str(&format!("'{}/wall_{}.vtk', ", &dir_name, step_no));
        }

        step_no += 1;
        t += dt;

        // progress bar increment
        pb.inc(1);
    }

    pb.finish_with_message("Simulation succesfully completed");

    // ---------------------------------------
    // write an easy paraview visualization file
    // ---------------------------------------
    // truncate the extra part of the string
    sand_files.truncate(sand_files.len() - 2);

    write_vis_file(
        format!("{}/vis_paraview.py", &dir_name),
        vec!["sand", "wall"],
        vec![sand_files, wall_files],
        vec![true, false],
        vec!["nnps"],
        vec![nnps_files],
    );
    // ---------------------------------------
    // write an easy paraview visualization file ends
    // ---------------------------------------

    // ---------------------------------------
    // PLOTTING
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
    // ---------------------------------------
    // PLOTTING ENDS
    // ---------------------------------------

}
