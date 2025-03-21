use bevy::prelude::*;

#[derive(Component, Default, Clone)]
pub struct Position(pub Vec2);

#[derive(Component, Default)]
pub struct Velocity(pub Vec2);

#[derive(Component)]
pub struct WaterDepth {
    pub current: f32,
    pub first: f32,
    pub is_dry: bool,
}

impl Default for WaterDepth {
    fn default() -> Self {
        Self {
            current: 0.0,
            first: -1.0,
            is_dry: true,
        }
    }
}

#[derive(Component)]
pub struct SmoothingLength {
    pub current: f32,
    pub first: f32,
}

impl SmoothingLength {
    pub fn new(initial: f32) -> Self {
        Self {
            current: initial,
            first: initial,
        }
    }
}
