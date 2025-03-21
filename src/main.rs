use std::{f32, vec};

use bevy::{prelude::*, window::PresentMode};
use bevy_asset_loader::prelude::*;
use bevy_fps_counter::FpsCounterPlugin;
use bevy_inspector_egui::quick::WorldInspectorPlugin;
use particle::{Position, SmoothingLength, Velocity, WaterDepth};

mod particle;

const SCALE: u32 = 30;
const PARTICLES_PER_PIXEL: u32 = 3;
const DAM_COLOR: Color = Color::srgba(1.0, 0.0, 0.0, 1.0);
const MAX_PARTICLES_IN_BUCKET: usize = 32;
const TOTAL_VOLUME: f32 = 12700000.0;
const GRAVITY: f32 = 9.8;
const CFL: f32 = 0.1;
const ROUGHNESS_COEFF: f32 = 0.0162;

#[derive(AssetCollection, Resource)]
struct MapAssets {
    #[asset(path = "dam.png")]
    dam: Handle<Image>,
    #[asset(path = "elevation.png")]
    elevation: Handle<Image>,
}

#[derive(Clone, Eq, PartialEq, Debug, Hash, Default, States)]
enum AppState {
    #[default]
    AssetLoading,
    Initializing,
    Simulation,
    Paused,
}

#[derive(Resource, Default)]
struct Bounds {
    render: Rect,
    simulation: Rect,
}

#[derive(Resource, Default)]
struct PixelSize(f32);

#[derive(Resource)]
struct MaxSmoothingLength(f32);

#[derive(Resource, Default)]
struct Bucket {
    data: Vec<[Entity; MAX_PARTICLES_IN_BUCKET]>,
    size: UVec2,
}

#[derive(Resource)]
struct Volume(f32);

#[derive(Resource)]
struct TimeStep(f32);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Dam Break Simulation".into(),
                name: Some("twod_sph_swe".into()),
                present_mode: PresentMode::Immediate,
                ..default()
            }),
            ..default()
        }))
        .add_plugins(WorldInspectorPlugin::new())
        .add_plugins(FpsCounterPlugin)
        .init_state::<AppState>()
        .insert_resource(TimeStep(1.0 / 60.0))
        .add_loading_state(
            LoadingState::new(AppState::AssetLoading)
                .continue_to_state(AppState::Initializing)
                .load_collection::<MapAssets>(),
        )
        .add_systems(
            OnEnter(AppState::Initializing),
            (setup, initialize_particles).chain(),
        )
        .add_systems(
            Update,
            (
                generate_bucket,
                update_water_depth,
                update_smoothing_length,
                update_time_step,
                update_velocity_and_position,
                check_particle_values,
            )
                .chain()
                .run_if(in_state(AppState::Simulation)),
        )
        .run();
}

fn setup(mut commands: Commands, assets: Res<Assets<Image>>, map_assets: Res<MapAssets>) {
    commands.spawn(Camera2d);

    let scale = 3.5;
    commands.insert_resource(PixelSize(scale));

    commands.spawn((
        Sprite::from_image(map_assets.dam.clone_weak()),
        Transform {
            scale: Vec3::splat(scale),
            ..Default::default()
        },
        Name::new("Map"),
    ));

    let image_dimentions = assets.get(&map_assets.dam).unwrap().size().as_vec2();
    let scaled_image_dimentions = image_dimentions * scale;

    let bounds = Bounds {
        render: Rect::from_center_size(Vec2::ZERO, scaled_image_dimentions),
        simulation: Rect::from_center_size(Vec2::ZERO, image_dimentions * SCALE as f32),
    };

    commands.insert_resource(bounds);
}

fn initialize_particles(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    images: Res<Assets<Image>>,
    bounds: Res<Bounds>,
    pixel_size: Res<PixelSize>,
    map_assets: Res<MapAssets>,
    mut next_state: ResMut<NextState<AppState>>,
) {
    commands.insert_resource(Bucket::default());

    let shape = meshes.add(Circle::new(0.5));
    let color = Color::srgb(0.0, 0.0, 1.0);

    let initial_smoothing_length = SCALE as f32 / PARTICLES_PER_PIXEL as f32 * 2.5;
    commands.insert_resource(MaxSmoothingLength(initial_smoothing_length));

    // Fix the panic by safely getting the image
    let image = match images.get(&map_assets.dam) {
        Some(img) => img,
        None => {
            error!("Failed to get dam image - asset not loaded");
            next_state.set(AppState::AssetLoading);
            return;
        }
    };

    // Rest of the function remains the same
    // let scale = pixel_size.0 / PARTICLES_PER_PIXEL as f32 * 0.5;

    let [width, height] = image.size().to_array();

    let mut number_particles = 0;
    for x in 0..width {
        for y in 0..height {
            if !is_dam_pixel(image, x, y) {
                continue;
            }

            let start_i = if x == 0 || PARTICLES_PER_PIXEL == 1 || !is_dam_pixel(image, x - 1, y) {
                0
            } else {
                1
            };

            let start_j = if y == 0 || PARTICLES_PER_PIXEL == 1 || !is_dam_pixel(image, x, y - 1) {
                0
            } else {
                1
            };

            for i in start_i..PARTICLES_PER_PIXEL {
                for j in start_j..PARTICLES_PER_PIXEL {
                    let mut position = Vec3::new(
                        x as f32 * pixel_size.0
                            + i as f32 * pixel_size.0 / PARTICLES_PER_PIXEL as f32,
                        y as f32 * pixel_size.0
                            + j as f32 * pixel_size.0 / PARTICLES_PER_PIXEL as f32,
                        0.0,
                    );

                    let render_bounds_size = bounds.render.max - bounds.render.min;

                    position.y = render_bounds_size.y - position.y;

                    let normalized_position = position.truncate() / render_bounds_size;

                    let real_world_position = bounds.simulation.min
                        + normalized_position * (bounds.simulation.max - bounds.simulation.min);

                    let map_position = bounds.render.min + normalized_position * render_bounds_size;

                    // println!(
                    //     "Spawning particle at image position: ({}, {}), sub-pixel: ({}, {}), position: {:?}, world position: {:?}",
                    //     x, y, i, j, map_position, real_world_position
                    // );

                    commands.spawn((
                        Mesh2d(shape.clone()),
                        MeshMaterial2d(materials.add(color)),
                        Transform {
                            translation: map_position.extend(1.0),
                            scale: Vec3::splat(1.0),
                            ..default()
                        },
                        Position(real_world_position),
                        Velocity::default(),
                        WaterDepth::default(),
                        SmoothingLength::new(initial_smoothing_length),
                    ));

                    number_particles += 1;
                }
            }
        }
    }

    commands.insert_resource(Volume(TOTAL_VOLUME / number_particles as f32));
    next_state.set(AppState::Simulation);
}

fn generate_bucket(
    q_particles: Query<(Entity, &Position)>,
    max_smoothing_length: Res<MaxSmoothingLength>,
    mut bucket: ResMut<Bucket>,
    bounds: Res<Bounds>,
) {
    let sim_bounds_size = bounds.simulation.max - bounds.simulation.min;
    bucket.size = (sim_bounds_size / max_smoothing_length.0).ceil().as_uvec2();
    bucket.data = vec![
        [Entity::from_raw(u32::MAX); MAX_PARTICLES_IN_BUCKET];
        (bucket.size.x * bucket.size.y) as usize
    ];

    for (i, p_i) in q_particles.iter() {
        let bucket_pos = ((p_i.0 - bounds.simulation.min)
            / (bounds.simulation.max - bounds.simulation.min)
            * (bucket.size.as_vec2() - Vec2::ONE))
            .as_uvec2();

        let index = (bucket_pos.x + bucket_pos.y * bucket.size.x) as usize;

        for j in 0..MAX_PARTICLES_IN_BUCKET {
            if bucket.data[index][j].index() == u32::MAX {
                bucket.data[index][j] = i;
                break;
            }
        }
    }
}

fn update_water_depth(
    mut query_set: ParamSet<(
        Query<(Entity, &Position, &SmoothingLength)>,
        Query<(Entity, &mut WaterDepth)>,
    )>,
    bucket: Res<Bucket>,
    bounds: Res<Bounds>,
    volume: Res<Volume>,
) {
    // First stage: Calculate water depths
    let mut water_depths: Vec<(Entity, f32)> = Vec::new();

    let q_positions = query_set.p0();

    for (i, p_i, l_i) in q_positions.iter() {
        let mut new_water_depth = 0.0;

        let bucket_pos = ((p_i.0 - bounds.simulation.min)
            / (bounds.simulation.max - bounds.simulation.min)
            * (bucket.size.as_vec2() - Vec2::ONE))
            .as_uvec2();

        for x_offset in -1..=1 {
            for y_offset in -1..=1 {
                let x_pos = bucket_pos.x as i32 + x_offset;
                let y_pos = bucket_pos.y as i32 + y_offset;

                if x_pos < 0
                    || x_pos >= bucket.size.x as i32
                    || y_pos < 0
                    || y_pos >= bucket.size.y as i32
                {
                    continue;
                }

                let bucket_pos_j = UVec2::new(x_pos as u32, y_pos as u32);
                let index_j = (bucket_pos_j.x + bucket_pos_j.y * bucket.size.x) as usize;

                for j2 in 0..MAX_PARTICLES_IN_BUCKET {
                    let j = bucket.data[index_j][j2];
                    if j.index() == u32::MAX {
                        continue;
                    }

                    if let Ok((_, p_j, _)) = q_positions.get(j) {
                        let r_ij = p_i.0.distance(p_j.0);
                        new_water_depth += volume.0 * w(r_ij, l_i.current);
                    }
                }
            }
        }

        water_depths.push((i, new_water_depth));
    }

    // Second stage: Update water depths
    let mut q_water_depths = query_set.p1();

    for (entity, new_water_depth) in water_depths.iter() {
        if let Ok((_, mut h_i)) = q_water_depths.get_mut(*entity) {
            if *new_water_depth > 0.01 {
                if h_i.first == -1.0 {
                    h_i.first = *new_water_depth;
                }
                h_i.current = *new_water_depth;

                h_i.is_dry = false;
            } else {
                h_i.is_dry = true;
            }
        }
    }
}

fn update_smoothing_length(
    mut q_particles: Query<(&WaterDepth, &mut SmoothingLength)>,
    mut max_smoothing_length: ResMut<MaxSmoothingLength>,
) {
    let mut new_max_smoothing_length = 0.0;

    for (h_i, mut l_i) in &mut q_particles {
        if h_i.is_dry {
            continue;
        }

        // Avoid division by zero
        let ratio = h_i.current / h_i.first;
        if ratio > 0.0 {
            // Avoid negative values inside sqrt
            l_i.current = l_i.first * f32::sqrt(ratio);

            // Update max smoothing length
            new_max_smoothing_length = f32::max(new_max_smoothing_length, l_i.current);
        }
    }

    max_smoothing_length.0 = new_max_smoothing_length;
}

fn update_time_step(
    q_particles: Query<(&Velocity, &WaterDepth, &SmoothingLength)>,
    mut time_step: ResMut<TimeStep>,
) {
    let mut min_time_step = f32::INFINITY;

    for (u_i, h_i, l_i) in q_particles.iter() {
        if h_i.is_dry || u_i.0.x <= f32::EPSILON || u_i.0.y <= f32::EPSILON {
            continue;
        }

        // Gravity wave speed of propagation
        let c_i = f32::sqrt(GRAVITY * h_i.current);

        min_time_step = f32::min(min_time_step, l_i.current / (c_i + u_i.0.length()))
    }

    if min_time_step > f32::EPSILON && min_time_step.is_finite() {
        time_step.0 = CFL * min_time_step;
    }
}

fn update_velocity_and_position(
    mut query_set: ParamSet<(
        Query<(Entity, &Position, &Velocity, &WaterDepth, &SmoothingLength)>,
        Query<(Entity, &mut Transform, &mut Position, &mut Velocity)>,
    )>,
    images: Res<Assets<Image>>,
    bounds: Res<Bounds>,
    map_assets: Res<MapAssets>,
    bucket: Res<Bucket>,
    volume: Res<Volume>,
    time_step: Res<TimeStep>,
) {
    let elevation_image_handle = &map_assets.elevation;
    let elevation = images.get(elevation_image_handle).unwrap();

    let mut accelerations: Vec<(Entity, Vec2)> = Vec::new();

    let q_immutable = query_set.p0();

    // Then process each particle
    for (i, p_i, u_i, h_i, l_i) in q_immutable.iter() {
        if h_i.is_dry {
            continue;
        }

        let mut sum = Vec2::ZERO;

        let bucket_pos = ((p_i.0 - bounds.simulation.min)
            / (bounds.simulation.max - bounds.simulation.min)
            * (bucket.size.as_vec2() - Vec2::ONE))
            .as_uvec2();

        let z_i = sample_elevation(elevation, p_i.0, bounds.simulation);

        for x_offset in -1..=1 {
            for y_offset in -1..=1 {
                let x_pos = bucket_pos.x as i32 + x_offset;
                let y_pos = bucket_pos.y as i32 + y_offset;

                if x_pos < 0
                    || x_pos >= bucket.size.x as i32
                    || y_pos < 0
                    || y_pos >= bucket.size.y as i32
                {
                    continue;
                }

                let bucket_pos_j = UVec2::new(x_pos as u32, y_pos as u32);
                let index_j = (bucket_pos_j.x + bucket_pos_j.y * bucket.size.x) as usize;

                for j2 in 0..MAX_PARTICLES_IN_BUCKET {
                    let j = bucket.data[index_j][j2];
                    if j.index() == u32::MAX {
                        break;
                    }

                    if let Ok((_, p_j, _, _, _)) = q_immutable.get(j) {
                        let z_j = sample_elevation(elevation, p_j.0, bounds.simulation);

                        // Calculate distance vector
                        let r_vec = p_j.0 - p_i.0;
                        let r_ij = r_vec.length();

                        let weight = w(r_ij, l_i.current);

                        // Apply elevation gradient to the sum

                        let bkp = sum;

                        let elevation_diff = z_j - z_i;

                        if elevation_diff.abs() < f32::EPSILON {
                            continue;
                        }

                        sum +=
                            (volume.0 / h_i.current) * elevation_diff * weight * r_vec.normalize();

                        if sum.x.is_nan()
                            || sum.y.is_nan()
                            || sum.x.is_infinite()
                            || sum.y.is_infinite()
                        {
                            warn!(
                                "Sum is NaN/Inf: before={:?}, v_j = {}, h_i = {}, z_i = {}, z_j = {}, w = {}, r_ij = {:?}",
                                bkp, volume.0, h_i.current, z_i, z_j, weight, r_vec
                            );

                            return;
                        }
                    }
                }
            }
        }

        let friction =
            -ROUGHNESS_COEFF.powi(2) * u_i.0 * u_i.0.length() / h_i.current.powf(4.0 / 3.0);

        accelerations.push((i, GRAVITY * (sum + friction)));

        if friction.x.is_nan()
            || friction.y.is_nan()
            || friction.x.is_infinite()
            || friction.y.is_infinite()
            || friction.x.abs() > 9999.0
            || friction.y.abs() > 9999.0
        {
            warn!(
                "Friction is {:?}: u_i = {:?}, h_i = {}",
                friction, u_i.0, h_i.current
            );
        }
    }

    let mut q_mutable = query_set.p1();

    for (entity, acceleration) in accelerations.iter() {
        if let Ok((_, mut transform, mut p_i, mut u_i)) = q_mutable.get_mut(*entity) {
            u_i.0 += acceleration * time_step.0;

            if u_i.0.x.abs().is_nan()
                || u_i.0.y.abs().is_nan()
                || u_i.0.x.abs().is_infinite()
                || u_i.0.y.abs().is_infinite()
                || u_i.0.x.abs() > 9999.0
                || u_i.0.y.abs() > 9999.0
            {
                warn!(
                    "Velocity is {:?}: acceleration = {:?}, timestep = {}",
                    u_i.0, acceleration, time_step.0
                );

                return;
            }

            let new_position = p_i.0 + u_i.0 * time_step.0;

            // Clamp position to simulation bounds
            let clamped_position = Vec2::new(
                new_position
                    .x
                    .clamp(bounds.simulation.min.x, bounds.simulation.max.x),
                new_position
                    .y
                    .clamp(bounds.simulation.min.y, bounds.simulation.max.y),
            );

            // Update positions
            p_i.0 = clamped_position;

            // Update the render position by mapping from simulation to render space
            let normalized_position =
                (p_i.0 - bounds.simulation.min) / (bounds.simulation.max - bounds.simulation.min);
            let render_position =
                bounds.render.min + normalized_position * (bounds.render.max - bounds.render.min);
            transform.translation = render_position.extend(1.0);
        }
    }
}

fn sample_elevation(image: &Image, position: Vec2, simulation_bounds: Rect) -> f32 {
    let uv = (position - simulation_bounds.min) / (simulation_bounds.max - simulation_bounds.min);

    let size = image.size().as_vec2();

    bilinear_sample(image, uv) * f32::max(size.x, size.y) * SCALE as f32
}

fn bilinear_sample(image: &Image, uv: Vec2) -> f32 {
    let size = image.size().as_vec2();

    // Convert UV coordinates to pixel coordinates
    let pixel_pos = uv * (size - Vec2::ONE);

    // Get the four closest pixels
    let x0 = pixel_pos.x.floor() as u32;
    let y0 = pixel_pos.y.floor() as u32;
    let x1 = (x0 + 1).min(size.x as u32 - 1);
    let y1 = (y0 + 1).min(size.y as u32 - 1);

    // Calculate the fractional part for interpolation weights
    let fx = pixel_pos.x - x0 as f32;
    let fy = pixel_pos.y - y0 as f32;

    // Get the four neighboring pixels
    let c00 = image
        .get_color_at(x0, y0)
        .unwrap_or(Color::BLACK)
        .to_srgba()
        .red;
    let c10 = image
        .get_color_at(x1, y0)
        .unwrap_or(Color::BLACK)
        .to_srgba()
        .red;
    let c01 = image
        .get_color_at(x0, y1)
        .unwrap_or(Color::BLACK)
        .to_srgba()
        .red;
    let c11 = image
        .get_color_at(x1, y1)
        .unwrap_or(Color::BLACK)
        .to_srgba()
        .red;

    // Bilinear interpolation
    let top = c00.lerp(c10, fx);
    let bottom = c01.lerp(c11, fx);
    top.lerp(bottom, fy)
}

fn w(distance: f32, smoothing_length: f32) -> f32 {
    let q = distance / smoothing_length;

    if q > 2.0 || q < 0.0 {
        return 0.0;
    }

    let alpha = 15.0 / (7.0 * f32::consts::PI * smoothing_length * smoothing_length);

    if q < 1.0 {
        alpha * (2.0 / 3.0 - q * q + (q * q * q) * 0.5)
    } else {
        let value = 2.0 - q;
        alpha * value * value * value / 6.0
    }
}

fn is_dam_pixel(image: &Image, x: u32, y: u32) -> bool {
    let pixel = image.get_color_at(x, y).unwrap();
    pixel == DAM_COLOR
}

fn check_particle_values(
    q_particles: Query<(Entity, &Position, &Velocity, &WaterDepth, &SmoothingLength)>,
    mut next_state: ResMut<NextState<AppState>>,
    time_step: Res<TimeStep>,
) {
    for (entity, position, velocity, water_depth, smoothing_length) in &q_particles {
        // Check for NaN or Infinity in position
        if position.0.x.is_nan()
            || position.0.x.is_infinite()
            || position.0.y.is_nan()
            || position.0.y.is_infinite()
        {
            warn!(
                "Entity {:?} has invalid position value: {:?}\n\tvelocity: {:?}\n\twater depth: {:?}\n\tis dry: {}\n\tsmoothing length: {:?}\n\ttime step: {}",
                entity.index(),
                position.0,
                velocity.0,
                water_depth.current,
                water_depth.is_dry,
                smoothing_length.current,
                time_step.0
            );

            next_state.set(AppState::Paused);
            break;
        }

        // Check for NaN or Infinity in velocity
        if velocity.0.x.is_nan()
            || velocity.0.x.is_infinite()
            || velocity.0.y.is_nan()
            || velocity.0.y.is_infinite()
        {
            warn!(
                "Entity {:?} has invalid velocity value: {:?}\n\tposition: {:?}\n\twater depth: {:?}\n\tis dry: {}\n\tsmoothing length: {:?}\n\ttime step: {}",
                entity.index(),
                velocity.0,
                position.0,
                water_depth.current,
                water_depth.is_dry,
                smoothing_length.current,
                time_step.0
            );

            next_state.set(AppState::Paused);
            break;
        }
    }
}
