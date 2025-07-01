// main.rs now uses items from lib.rs

use anyhow::Result;
use log::error;
use orion_kmer::{cli, commands}; // Use items from the library part of the crate

fn main() -> Result<()> {
    let matches = cli::parse_cli();

    if let Err(e) = commands::dispatch_command(matches.command, matches.threads, matches.verbose) {
        error!("Error: {}", e);
        std::process::exit(1);
    }

    Ok(())
}
