mod barnes_hut;

fn main() {
    let mut sim = barnes_hut::BarnesHut::new(1024, 0.5, 4, 4, 10.0);

    sim.init_cube(1024);
    println!("{:?}", sim.get_particle(0));

    for _ in 1..256 {
        sim.step(1.0/60.0);
        println!("{:?}", sim.get_particle(128));
    }

    println!("Hello, world");
}
