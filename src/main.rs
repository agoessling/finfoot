use uom::fmt::DisplayStyle::Abbreviation;
use uom::si::acceleration::meter_per_second_squared;
use uom::si::length::meter;
use uom::si::velocity::meter_per_second;
use uom::si::{
    f64::{Acceleration, Length, Velocity},
    ISQ,
};
use uom::si::{Dimension, Quantity, Units};
use uom::Conversion;

use typenum::operator_aliases::Diff;
use typenum::{Integer, P1};

use uom::num_traits::Num;

trait DivideByTime {
    type Type;
}

impl<D, U, V> DivideByTime for Quantity<D, U, V>
where
    D: Dimension + ?Sized,
    U: Units<V> + ?Sized,
    V: Num + Conversion<V>,
    D::T: std::ops::Sub<P1>,
    Diff<D::T, P1>: Integer,
{
    type Type = Quantity<ISQ<D::L, D::M, Diff<D::T, P1>, D::I, D::Th, D::N, D::J, D::Kind>, U, V>;
}

type TimeDerivative<T> = <T as DivideByTime>::Type;

fn main() {
    let position = Length::new::<meter>(1.0);
    let velocity = Velocity::new::<meter_per_second>(1.0);
    let test = TimeDerivative::<Velocity>::new::<meter_per_second_squared>(2.0);
    println!(
        "{}, {}, {}",
        Length::format_args(meter, Abbreviation).with(position),
        Velocity::format_args(meter_per_second, Abbreviation).with(velocity),
        Acceleration::format_args(meter_per_second_squared, Abbreviation).with(test),
    );
}
