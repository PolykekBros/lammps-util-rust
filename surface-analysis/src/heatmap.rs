use anyhow::Result;
use colorgrad::Gradient;
use geomutil_util::Point2;
use lammps_util_rust::range_f32;
use plotters::{prelude::*, style::BLACK};

pub struct Domain {
    lo: Point2,
    hi: Point2,
}

impl Domain {
    pub fn new(a: Point2, b: Point2) -> Self {
        Domain {
            lo: Point2::from([a.x.min(b.x), a.y.min(b.y)]),
            hi: Point2::from([a.x.max(b.x), a.y.max(b.y)]),
        }
    }

    pub fn lo(&self) -> Point2 {
        self.lo
    }

    pub fn hi(&self) -> Point2 {
        self.hi
    }

    pub fn width(&self) -> f32 {
        self.hi.x - self.lo.x
    }

    pub fn height(&self) -> f32 {
        self.hi.y - self.lo.y
    }

    pub fn area(&self) -> f32 {
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
    min: f32,
    max: f32,
    gradient: T,
}

impl<T: Gradient> Colorbar<T> {
    pub fn new(min: f32, max: f32, gradient: T) -> Self {
        Self {
            min: min.min(max),
            max: min.max(max),
            gradient,
        }
    }

    pub fn color(&self, value: f32) -> RGBColor {
        let value = self.min.max(value).min(self.max);
        let scaled = (value - self.min) / (self.max - self.min);
        let rgba = self.gradient.at(scaled as f32).to_rgba8();
        RGBColor(rgba[0], rgba[1], rgba[2])
    }

    pub fn draw<DB: DrawingBackend>(&self, mut chart_builder: ChartBuilder<DB>) {
        let &Self { min, max, .. } = self;
        let step = (max - min) / 255.0;
        let mut chart_context = chart_builder
            .margin_top(10)
            .margin_right(15)
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
            .label_style(("sans-serif", 16).into_font().color(&BLACK))
            .draw()
            .unwrap();
        let plotting_area = chart_context.plotting_area();
        range_f32(min, max, 256).iter().for_each(|&value| {
            let color = self.color(value);
            let rectangle =
                Rectangle::new([(0.0, value), (1.0, value + step)], filled_style(color));
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
    assert!(domain.area().is_normal());

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
        .label_style(("sans-serif", 16).into_font().color(&BLACK))
        .draw()
        .unwrap();

    let plotting_area = chart_context.plotting_area();

    let x_count = data.shape().0;
    let x_values = (0..x_count).collect::<Vec<_>>();
    let x_step = domain.width() / x_count as f64;
    let y_count = data.shape().1;
    let y_values = (0..y_count).collect::<Vec<_>>();
    let y_step = domain.height() / y_count as f64;

    x_values.iter().for_each(|&x| {
        let r_x = domain.lo().x + x as f32 * x_step;
        y_values.iter().for_each(|&y| {
            let r_y = domain.lo().y + y as f32 * y_step;
            let rectangle = Rectangle::new(
                [(r_x, r_y), (r_x + x_step, r_y + y_step)],
                filled_style(colorbar.color(data[(x, y)])),
            );
            plotting_area.draw(&rectangle).unwrap();
        });
    });

    Ok(())
}
