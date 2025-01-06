use anyhow::{anyhow, Result};
use colorgrad::Gradient;
use lammps_util_rust::range_f64;
use nalgebra::{point, DMatrix, Point2};
use plotters::{prelude::*, style::BLACK};

pub struct Domain {
    lo: Point2<f64>,
    hi: Point2<f64>,
}

impl Domain {
    pub fn new(p1: Point2<f64>, p2: Point2<f64>) -> Self {
        Domain {
            lo: point![p1.x.min(p2.x), p1.y.min(p2.y)],
            hi: point![p1.x.max(p2.x), p1.y.max(p2.y)],
        }
    }

    pub fn lo(&self) -> Point2<f64> {
        self.lo
    }

    pub fn hi(&self) -> Point2<f64> {
        self.hi
    }

    pub fn width(&self) -> f64 {
        self.hi.x - self.lo.x
    }

    pub fn height(&self) -> f64 {
        self.hi.y - self.lo.y
    }

    pub fn area(&self) -> f64 {
        self.width() * self.height()
    }
}

fn filled_style<C: Into<RGBAColor>>(color: C) -> ShapeStyle {
    ShapeStyle {
        color: color.into(),
        filled: true,
        stroke_width: 0,
    }
}

pub struct Colorbar<T: Gradient> {
    min: f64,
    max: f64,
    gradient: T,
}

impl<T: Gradient> Colorbar<T> {
    pub fn new(min: f64, max: f64, gradient: T) -> Self {
        Self {
            min: min.min(max),
            max: min.max(max),
            gradient,
        }
    }

    pub fn color(&self, value: f64) -> RGBColor {
        let value = self.min.max(value).min(self.max);
        let scaled = (value - self.min) / (self.max - self.min);
        let rgba = self.gradient.at(scaled as f32).to_rgba8();
        RGBColor(rgba[0], rgba[1], rgba[2])
    }

    pub fn draw<DB: DrawingBackend>(&self, mut chart_builder: ChartBuilder<DB>) {
        let &Self { min, max, .. } = self;
        let step = (max - min) / 256.0;
        let mut chart_context = chart_builder
            .margin_top(10)
            .x_label_area_size(25)
            .y_label_area_size(40)
            .build_cartesian_2d(0.0..1.0, min..max)
            .unwrap();
        chart_context
            .configure_mesh()
            .set_all_tick_mark_size(5)
            .disable_x_axis()
            .disable_x_mesh()
            .disable_y_mesh()
            .axis_style(BLACK)
            .label_style("sans-serif".into_font().color(&BLACK))
            .draw()
            .unwrap();
        let plotting_area = chart_context.plotting_area();
        range_f64(min, max, 256).iter().for_each(|&value| {
            let color = self.color(value);
            let rectangle = Rectangle::new(
                [(0.0, value - step / 2.0), (1.0, value + step / 2.0)],
                filled_style(color),
            );
            plotting_area.draw(&rectangle).unwrap();
        });
    }
}

pub fn heatmap<DB: DrawingBackend, T: Gradient>(
    data: &DMatrix<f64>,
    domain: &Domain,
    colorbar: &Colorbar<T>,
    mut chart_builder: ChartBuilder<DB>,
) -> Result<()> {
    if !domain.area().is_normal() {
        return Err(anyhow!("Invalid domain"));
    }

    let mut chart_context = chart_builder
        .margin_top(10)
        .x_label_area_size(25)
        .y_label_area_size(40)
        .build_cartesian_2d(domain.lo().x..domain.hi().x, domain.lo().y..domain.hi().y)
        .unwrap();

    chart_context
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .set_all_tick_mark_size(5)
        .axis_style(BLACK)
        .label_style("sans-serif".into_font().color(&BLACK))
        .draw()
        .unwrap();

    let plotting_area = chart_context.plotting_area();

    let x_count = data.shape().1;
    let x_values = (0..x_count).collect::<Vec<_>>();
    let x_step = domain.width() / x_count as f64;
    let y_count = data.shape().0;
    let y_values = (0..y_count).collect::<Vec<_>>();
    let y_step = domain.height() / y_count as f64;

    x_values.iter().for_each(|&x| {
        let r_x = domain.lo().x + x as f64 * x_step;
        y_values.iter().for_each(|&y| {
            let r_y = domain.lo().y + y as f64 * y_step;
            let rectangle = Rectangle::new(
                [(r_x, r_y), (r_x + x_step, r_y + y_step)],
                filled_style(colorbar.color(data[(x, y)])),
            );
            plotting_area.draw(&rectangle).unwrap();
        });
    });

    Ok(())
}
