use std::cmp::min;
use std::collections::HashMap;
use std::iter::zip;
use anyhow::{anyhow, Result};
use log::debug;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PointND {
    dimensions: usize,
    coordinates: Vec<i128>,
}
impl PointND {
    pub fn new(coordinates: Vec<i128>) -> Self {
        Self {
            dimensions: coordinates.len(),
            coordinates,
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
struct Fraction {
    num: i128,
    denom: i128,
}

impl Fraction {
    /// Get a normalised Fraction.
    pub fn new_normalised(numerator: i128, denominator: i128) -> Self {
        let gcf = gcd(numerator, denominator);
        let normalised_numerator = numerator / gcf;
        let normalised_denominator = denominator / gcf;
        if normalised_denominator < 0 {
            Self {
                num: -normalised_numerator,
                denom: -normalised_denominator,
            }
        } else {
            Self {
                num: normalised_numerator,
                denom: normalised_denominator,
            }
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
struct LineND {
    dimensions: usize,
    point: Vec<Fraction>,
    direction: Vec<Fraction>,
}

impl LineND {

    /// Get a normalised 3-dimensional line object from two points.
    ///
    /// Normalise the line by getting a normalised direction and normalised point.
    pub fn new_normalised(point_1: &PointND, point_2: &PointND) -> Self {
        if point_1.dimensions != point_2.dimensions {
            panic!("Points provided to line constructor have different dimensions.")
        }
        let mut point_numerators: Vec<i128> = point_1.coordinates.clone();
        let mut point_denominators: Vec<i128> = vec![1; point_1.dimensions];
        let mut direction_numerators: Vec<i128> = zip(
            &point_1.coordinates, &point_2.coordinates
        ).map(|(p1_value, p2_value)| p2_value - p1_value).collect();
        let mut direction_denominators: Vec<i128> = vec![1; point_1.dimensions];
        let mut normalisation_coordinate_index = 0;
        while direction_numerators[normalisation_coordinate_index] == 0 {
            normalisation_coordinate_index += 1;
        }
        let direction = LineND::normalise_direction(
            &mut direction_numerators,
            &mut direction_denominators,
            normalisation_coordinate_index
        );
        let point = LineND::normalise_point(
            &mut point_numerators,
            &mut point_denominators,
            &direction_numerators,
            &direction_denominators,
            normalisation_coordinate_index,
        );
        Self {
            dimensions: point_1.dimensions,
            point,
            direction,
        }
    }

    /// The direction is normalised by making the first nonzero coordinate for
    /// the direction vector equal to one and representing the remaining
    /// coordinates in reduced fractions of integers.
    fn normalise_direction(
        direction_numerators: &mut Vec<i128>,
        direction_denominators: &mut Vec<i128>,
        normalisation_coordinate_index: usize
    ) -> Vec<Fraction> {
        let dimensions = direction_numerators.len();
        let normalisation_multiplier = direction_numerators[normalisation_coordinate_index];
        for index in normalisation_coordinate_index..dimensions {
            direction_denominators[index] *= normalisation_multiplier;
            if direction_denominators[index] < 0 {
                direction_numerators[index] *= -1;
                direction_denominators[index] *= -1;
            }
        }
        for index in 0..dimensions {
            if direction_numerators[index] == 0 {
                direction_denominators[index] = 1;
                continue;
            }
            let div = gcd(direction_numerators[index], direction_denominators[index]);
            direction_numerators[index] /= div;
            direction_denominators[index] /= div;
        }
        zip(direction_numerators, direction_denominators).map(|(num, denom)| Fraction::new_normalised(*num, *denom)).collect()
    }

    /// The point is normalised by shifting it along the line such that the coordinate
    /// corresponding to the first nonzero coordinate of the direction vector is 0.
    /// The remaining coordinates of the point are represented as reduced fractions
    /// of integers.
    fn normalise_point(
        point_numerators: &mut Vec<i128>,
        point_denominators: &mut Vec<i128>,
        direction_numerators: &Vec<i128>,
        direction_denominators: &Vec<i128>,
        normalisation_coordinate_index: usize,
    ) -> Vec<Fraction> {
        let dimensions = point_numerators.len();
        let normalisation_multiplier = point_numerators[normalisation_coordinate_index];
        for index in 0..dimensions {
            point_numerators[index] *= direction_denominators[index];
            point_denominators[index] *= direction_denominators[index];
            point_numerators[index] -= normalisation_multiplier * direction_numerators[index];
        }
        zip(point_numerators, point_denominators).map(|(num, denom)| Fraction::new_normalised(*num, *denom)).collect()
    }
}

/// Gets the greatest common divisor of two integers.
///
/// The result returned will always be a nonnegative integer. If
/// one of the arguments is 0, then the result will be the absolute
/// value of the other argument.
pub fn gcd(a: i128, b: i128) -> i128 {
    if a < 0 {
        return gcd(-a, b);
    }
    if b == 0 {
        return a;
    }
    gcd(b, a % b)
}


/// Get the largest number of points in the sequence intersected by a single line.
///
/// Only lines going through points with indices in [`start_index`, `end_index`) are
/// considered and only points after `start_index` are considered.
/// Only lines whose point indices are separated by at most `window_size` are
/// considered.
pub fn count_collinear_points(
    point_sequence: &[PointND],
    start_index: usize,
    end_index: usize,
    window_size: usize,
) -> i32 {
    if window_size == 0 {
        return 1;
    }
    let mut count_max = 0;
    for i in start_index..end_index {
        debug!(
            "Progress: considering lines through {}, current max collinear is: {}",
            i, count_max
        );
        let mut line_counts = HashMap::new();
        let window_end = min(point_sequence.len(), i + window_size + 1);
        for j in i + 1..window_end {
            let line = LineND::new_normalised(&point_sequence[i], &point_sequence[j]);
            let count = line_counts.entry(line).or_insert(1);
            *count += 1;
            if *count > count_max {
                count_max = *count;
            }
        }
    }
    count_max
}

fn build_point_sequence(char_string: &str, dimensions: usize) -> Result<Vec<PointND>> {
    let mut points: Vec<PointND> = vec![];
    let mut curr_position = PointND {
        dimensions,
        coordinates: vec![0; dimensions],
    };
    points.push(curr_position.clone());
    for ch in char_string.chars() {
        let digit = ch.to_digit(10).ok_or(anyhow!("Invalid character {}", ch))? as usize;
        if digit >= dimensions {
            return Err(anyhow!("Digit {} too large for {} dimensions", digit, dimensions))
        }
        let mut unit_vector = vec![0; dimensions];
        unit_vector[digit] = 1;
        for index in 0..dimensions {
            curr_position.coordinates[index] += unit_vector[index];
        }
        points.push(curr_position.clone());
    }
    Ok(points)
}

fn main() {
    env_logger::init();
    let point_sequence = build_point_sequence("0101", 4).unwrap();
    println!("Sequence length: {}", point_sequence.len());
    println!(
        "Maximum number of collinear points: {}",
        count_collinear_points(
            &point_sequence, 0, point_sequence.len(), point_sequence.len()
        )
    );
}
