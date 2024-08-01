use clap::{Parser, Subcommand};

use crate::variant_key::{decode_variant_key, parse_variant};

#[cfg(test)]
mod test;

mod variant_key;

#[derive(Parser)]
#[command(name = "vkrs", version = "0.1.0")]
#[command(bin_name = "vkrs")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Converts a variant string into a variantkey
    #[command(arg_required_else_help = true,
        name = "vk",
        about = "Returns the variantkey for the input variant",
    )]
    ToKey {
        /// Variant in the form chr-pos-ref-alt e.g. MT-12345-A-T, 1-23345-C-CT
        variant_string: String,
        #[arg(short = 'x', long = "hex")]
        convert_hex: bool,
    },

    /// Converts a variantkey into a variant string
    #[command(arg_required_else_help = true,
        name = "kv",
        about = "Returns the variant for the input key",
    )]
    ToVariant {
        /// Variant key
        variant_key: String,
        #[arg(short = 'x', long = "hex")]
        convert_hex: bool,
    },
}

fn main() {
    let args = Cli::parse();
    match args.command {
        Commands::ToKey { variant_string, convert_hex } => {
            run_to_key_command(&variant_string, convert_hex);
        }
        Commands::ToVariant { variant_key, convert_hex } => {
            run_to_variant_command(&variant_key, convert_hex)
        }
    }
}

fn run_to_key_command(variant_string: &str, convert_hex: bool) {
    let variant = match parse_variant(&variant_string) {
        Some(variant) => variant,
        None => return eprintln!("Illegal variant string! {}", &variant_string),
    };
    let variant_key = &variant.to_variant_key();
    if convert_hex {
        return println!("{:x}", variant_key);
    }
    println!("{}", variant_key);
}

fn run_to_variant_command(variant_key_string: &str, convert_hex: bool) {
    let variant_key = if convert_hex { u64::from_str_radix(variant_key_string, 16).unwrap() } else { u64::from_str_radix(variant_key_string, 10).unwrap() };
    let variant = decode_variant_key(&variant_key);
    println!("{}", variant);
}

