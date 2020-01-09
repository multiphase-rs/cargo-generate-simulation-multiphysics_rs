const USAGE: &str = "
Usage: test_2 [--dt N] [--plots]
       test_2 --help

A template for rigid bodies interacting with wall and among themselves

Options:
    --dt N          Time step of the simulation [default: 5e-5]
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
// use gnuplot::*;
use multiphysics::prelude::*;
use simple_shapes::{benchmarks::create_zhang_geometry};

// -------------------------
// local imports
// -------------------------
// dem granular
pub use dem_rust::basic_equations::*;
pub use dem_rust::dem_3d::*;
pub use dem_rust::prelude::*;
pub use dem_rust::wall::*;

// rigid body imports
pub use dem_rust::rb::rb_2d::Rigidbody2D;
pub use dem_rust::rb::rb_utils::get_body_id_limits_from_body_id;

// -------------------------------------------------

fn get_center_of_mass(com: &[Vector3<f64>]) -> Vector3<f64>{
    let mut avg_com = Vector3::new(0., 0., 0.);
    for i in 0..com.len(){
        avg_com += com[i];
    }
    avg_com = avg_com / 33.;
    return avg_com;
}

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
    // CREATE PARTICLE ARRAYS
    // --------------------------------------
    // let spacing = 0.0005;
    let spacing = 0.001;

    // create cylinders
    let (xc, yc, bid, xt, yt) = create_zhang_geometry(spacing, 2, true);
    let (mut x, mut y) = (xc, yc);
    // x.iter_mut().for_each(|t| *t += spacing);
    // y.iter_mut().for_each(|t| *t += spacing);
    let body_id_limits = get_body_id_limits_from_body_id(&bid);
    let dem_id = bid.clone();

    let z = vec![0.; x.len()];
    let radius = vec![spacing / 2.; x.len()];

    let mut bodies = Rigidbody2D::from_xyz_rad(&x, &y, &z, &radius);
    // println!("I am here");
    let m = 2.7 * 1e3 * PI * (spacing / 2.).powf(2.);
    bodies.m = vec![m; x.len()];
    bodies.body_id = bid;
    bodies.body_id_limits = body_id_limits;
    // println!("{:?}", bodies.body_id_limits);
    bodies.dem_id = dem_id;
    // println!("{:?}", bodies.dem_id);

    // set material properties of the particle
    let yng = 1e10;
    let nu = 0.3;
    let shr = yng / (2. * (1. + nu));
    bodies.young_mod = vec![yng; x.len()];
    bodies.poisson = vec![nu; x.len()];
    bodies.shear_mod = vec![shr; x.len()];

    setup_2d_rigid_body!(bodies);

    // set some velocity to the center of mass of the body
    // bodies.velocity_com[0][0] = 2.;
    // bodies.velocity_com[0][1] = 0.;
    // bodies.angular_velocity[0][2] = 0.;

    // bodies.velocity_com[1][0] = -2.;
    // bodies.velocity_com[1][1] = 0.;
    // bodies.angular_velocity[1][2] = 0.;

    // let stiffness_kn = 5. * 1e6;
    // let nu = 0.3;
    // let stiffness_kt = 2. * (1. - nu) / (2. - nu) * stiffness_kn;
    // body_1.kn = vec![stiffness_kn, stiffness_kn];
    // body_1.kt = vec![stiffness_kt, stiffness_kt];

    // this is nice step for debugging
    bodies.validate_particle_array();

    // create a wall
    let x = vec![0., 12. * 1e-2, 26. * 1e-2, 6. * 1e-2];
    let y = vec![10. * 1e-2, 0., 10. * 1e-2, 10. * 1e-2];
    let z = vec![0., 0., 0., 0.];
    let wall_points = vec![
        [Vector3::new(0., 0., 0.), Vector3::new(0., 26. * 1e-2, 0.)],
        [Vector3::new(0., 0., 0.), Vector3::new(26. * 1e-2, 0., 0.)],
        [Vector3::new(26. * 1e-2, 0., 0.), Vector3::new(26. * 1e-2, 26. * 1e-2, 0.)],
        [Vector3::new(6. * 1e-2, 0., 0.), Vector3::new(6. * 1e-2, 26. * 1e-2, 0.)],
    ];

    let nx = vec![1., 0., -1., -1.];
    let ny = vec![0., 1., 0., 0.];
    let nz = vec![0., 0., 0., 0.];
    let mut wall =
        WallDEMParticleArray::from_xyz_wall_points_normals(&x, &y, &z, &wall_points, &nx, &ny, &nz);
    wall.young_mod = vec![yng; x.len()];
    wall.poisson = vec![nu; x.len()];
    wall.shear_mod = vec![shr; x.len()];
    wall.dem_id = vec![2; x.len()];
    // --------------------------------------
    // CREATE PARTICLE ARRAYS ENDS
    // --------------------------------------


    // ---------------------------------------
    // SETUP CORRESPONDING NBS NNPS
    // ----------------------------------------
    let max_coordinate = 27. * 1e-2;
    let max_size = 2. * radius[0];
    let mut nbs2d_bodies =
        NBS2D::from_maximum_and_no_of_particles(max_coordinate, max_size, bodies.x.len());
    // ---------------------------------------
    // SETUP CORRESPONDING NBS NNPS ENDS
    // ----------------------------------------


    // ----------------------------------------
    // SOLVER DATA
    // ----------------------------------------
    let release_gate_time = 0.05;
    let tf = release_gate_time + 0.4;
    let dt = args.flag_dt;
    let total_steps = (tf / dt) as u64;
    let mut t = 0.;
    let mut step_no = 0;
    let total_output_file = 100;
    let pfreq = if total_steps < total_output_file {
        1
    } else {
        total_steps / total_output_file
    };
    // let pfreq = 100;
    // ----------------------------------------
    // SOLVER DATA ENDS
    // ----------------------------------------

    // ----------------------------------------
    // OUTPUT DIRECTORY
    // ----------------------------------------
    let project_root = env!("CARGO_MANIFEST_DIR");
    let dir_name = project_root.to_owned() + "/test_2_output";
    let _p = fs::create_dir(&dir_name);
    // ----------------------------------------
    // OUTPUT DIRECTORY ENDS
    // ----------------------------------------


    // ----------------------------------------
    // SOME CONSTANTS
    // ----------------------------------------
    // some constants
    let stage1 = 1;
    let stage2 = 2;
    // coefficient of friction
    let mu = 0.5;
    // ----------------------------------------
    // SOME CONSTANTS ENDS
    // ----------------------------------------

    // create a progress bar
    let pb = setup_progress_bar(total_steps);

    // -------------------------------
    // FOR PARAVIEW VISUALIZATION
    // -------------------------------
    // for paraview visualization
    let mut bodies_files = "".to_string();
    write_to_vtk!(bodies, format!("{}/bodies_{}.vtk", &dir_name, step_no));
    bodies_files.push_str(&format!("'{}/bodies_{}.vtk', ", &dir_name, step_no));

    let mut wall_files = "".to_string();
    write_wall_to_vtk!(wall, format!("{}/wall_{}.vtk", &dir_name, step_no));
    wall_files.push_str(&format!("'{}/wall_{}.vtk', ", &dir_name, step_no));

    let mut nnps_files = "".to_string();
    write_nnps_2d_to_vtk!(nbs2d_bodies, format!("{}/nnps.vtk", &dir_name));
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
        rk2_initialize_copy_contacts_3d!((bodies));
        rk2_rb2d_initialize!((bodies));


        // ----------------------
        // STAGE 1 EQUATIONS
        // ----------------------
        nbs2d_bodies.register_particles_to_nnps(&bodies.x, &bodies.y, &bodies.z);
        make_forces_torques_zero!((bodies));
        body_force!((bodies), 0., -9.81, 0.);
        rb_dem_interaction_force_hertz_pw!(bodies, (wall), (0.5), (0.4), stage1, dt);
        rb_dem_interaction_force_hertz_pp!(
            bodies,
            (bodies),
            (nbs2d_bodies),
            (0.5),
            (0.4),
            stage1,
            dt
        );
        sumUpExternalForces!((bodies));


        // ----------------------
        // STAGE 1 STEPPER
        // ----------------------
        rk2_rb2d_stage_1!((bodies), dt);


        // ----------------------
        // STAGE 2 EQUATIONS
        // ----------------------
        nbs2d_bodies.register_particles_to_nnps(&bodies.x, &bodies.y, &bodies.z);
        make_forces_torques_zero!((bodies));
        body_force!((bodies), 0., -9.81, 0.);
        rb_dem_interaction_force_hertz_pw!(bodies, (wall), (0.5), (0.4), stage2, dt);
        rb_dem_interaction_force_hertz_pp!(
            bodies,
            (bodies),
            (nbs2d_bodies),
            (0.5),
            (0.4),
            stage2,
            dt
        );
        sumUpExternalForces!((bodies));

        // ----------------------
        // STAGE 2 STEPPER
        // ----------------------
        rk2_rb2d_stage_2!((bodies), dt);

        update_contacts_pp!(bodies, (bodies));
        update_contacts_pw!(bodies, (wall));

        // ----------------------
        // OUTPUT THE VTK FILE AND WRITE ALL FILES DATA
        // ----------------------
        if step_no % pfreq == 0 && step_no != 0 {
            // ----------------------
            // WRITE DATA TO VTK FILE
            // ----------------------
            write_to_vtk!(bodies, format!("{}/bodies_{}.vtk", &dir_name, step_no));
            write_wall_to_vtk!(wall, format!("{}/wall_{}.vtk", &dir_name, step_no));

            // ----------------------
            // FOR PARAVIEW AUTOMATION
            // ----------------------
            wall_files.push_str(&format!("'{}/wall_{}.vtk', ", &dir_name, step_no));
            bodies_files.push_str(&format!("'{}/bodies_{}.vtk', ", &dir_name, step_no));
        }

        step_no += 1;
        t += dt;

        // progress bar increment
        pb.inc(1);

        // move a wall after some time
        if t < release_gate_time + dt / 2. && t > release_gate_time - dt / 2. {
            wall.x[3] = 27. * 1e-2;
            wall.wall_points[3][0] = Vector3::new(27. * 1e-2, 0., 0.);
            wall.wall_points[3][1] = Vector3::new(27. * 1e-2, 26. * 1e-2, 0.);

            // print the center of mass of the system
            let com = get_center_of_mass(&bodies.com);
            println!("------------------------------");
            println!("------------------------------");
            println!("------------------------------");
            println!("------------------------------");
            println!("center of mass is {:?}", com / (26. * 1e-2));
            println!("------------------------------");
            println!("------------------------------");
            println!("------------------------------");
            println!("------------------------------");
        }
    }

    pb.finish_with_message("Simulation succesfully completed");

    // ---------------------------------------
    // write an easy paraview visualization file
    // ---------------------------------------
    // truncate the extra part of the string
    bodies_files.truncate(bodies_files.len() - 2);
    wall_files.truncate(wall_files.len() - 2);

    write_vis_file(
        format!("{}/vis_paraview.py", &dir_name),
        vec!["bodies", "wall"],
        vec![bodies_files, wall_files],
        vec![true, false],
        vec!["nnps"],
        vec![nnps_files],
    );
}
