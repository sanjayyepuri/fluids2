use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern {
    pub fn alert(s: &str);
}

#[wasm_bindgen]
pub fn add(x: f32, y: f32) -> f32 {
    x + y
}

// the first step to create a pic simulation is to get a datastructure to represent the particles
// we then need to ensure this datastructure is friendly to be read directly from memory in javascript