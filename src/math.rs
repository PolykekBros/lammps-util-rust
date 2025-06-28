pub fn range_f32(begin: f32, end: f32, count: usize) -> Vec<f32> {
    assert!(begin <= end);
    let step = (end - begin) / 1.max(count - 1) as f32;
    (0..count).map(|n| begin + n as f32 * step).collect()
}

pub fn range_f64(begin: f64, end: f64, count: usize) -> Vec<f64> {
    assert!(begin <= end);
    let step = (end - begin) / 1.max(count - 1) as f64;
    (0..count).map(|n| begin + n as f64 * step).collect()
}

pub trait IteratorAvg: Iterator {
    fn avg(self) -> Option<f64>
    where
        Self::Item: Into<f64> + Copy;

    fn avg_with_std(self) -> Option<(f64, f64)>
    where
        Self::Item: Into<f64> + Copy;

    fn std(self) -> Option<f64>
    where
        Self::Item: Into<f64> + Copy;
}

impl<I> IteratorAvg for I
where
    I: Iterator,
    I::Item: Into<f64> + Copy,
{
    fn avg(self) -> Option<f64> {
        self.map(|v| (v.into(), 1))
            .reduce(|(sum, cnt), (next, next_cnt)| (sum + next, cnt + next_cnt))
            .map(|(sum, cnt)| sum / cnt as f64)
    }

    fn avg_with_std(self) -> Option<(f64, f64)> {
        let values = self.collect::<Vec<_>>();
        let avg = values.iter().copied().avg()?;
        values
            .iter()
            .map(|&v| (v.into() - avg).powi(2))
            .avg()
            .map(|x| (avg, x.sqrt()))
    }

    fn std(self) -> Option<f64> {
        self.avg_with_std().map(|(_, x)| x)
    }
}
