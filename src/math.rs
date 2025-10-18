use std::iter;

macro_rules! impl_range {
    ($($a:ident)*) => ($(
        pub fn $a(begin: $a, end: $a, count: usize) -> impl Iterator<Item = $a> {
            let step = (end - begin) / 1.max(1.max(count) - 1)  as $a;
            (0..count).map(move |i| (i as $a).mul_add(step, begin))
        }
    )*)
}

pub mod range {
    impl_range! { f32 f64 }
}

pub trait IteratorAvg<T>: Iterator<Item = T> {
    fn avg(self) -> Option<T>;
    fn avg_with_std(self) -> Option<(T, T)>;
    fn std(self) -> Option<T>;
}

macro_rules! impl_avg {
    ($($a:ident)*) => ($(
        impl<I> IteratorAvg<$a> for I
        where
            I: Iterator<Item = $a>,
        {
            fn avg(self) -> Option<$a> {
                iter::zip(self, 1usize..)
                    .reduce(|(sum, _), (next, cnt)| (sum + next, cnt))
                    .map(|(sum, cnt)| sum / cnt as $a)
            }

            fn avg_with_std(self) -> Option<($a, $a)> {
                let values = self.collect::<Vec<_>>();
                let avg = values.iter().copied().avg()?;
                values
                    .into_iter()
                    .map(|x| (x - avg).powi(2))
                    .avg()
                    .map(|std| (avg, std.sqrt()))
            }

            fn std(self) -> Option<$a> {
                self.avg_with_std().map(|(_, std)| std)
            }
        }
    )*)
}

impl_avg! { f32 f64 }

#[cfg(test)]
mod tests {
    use super::*;
    use assert_float_eq::assert_f32_near;

    #[test]
    fn test_iterator_avg_f32() {
        let data: Vec<f32> = vec![];
        let avg = data.iter().copied().avg();
        let std = data.iter().copied().std();
        assert!(avg.is_none());
        assert!(std.is_none());
        let data: Vec<f32> = vec![5.0];
        let avg = data.iter().copied().avg();
        let std = data.iter().copied().std();
        assert!(avg.is_some());
        assert_f32_near!(avg.unwrap(), 5.0, 6);
        assert!(std.is_some());
        assert_f32_near!(std.unwrap(), 0.0, 6);
        let data: Vec<f32> = vec![1.0, 2.0, 3.0];
        let avg = data.iter().copied().avg();
        let std = data.iter().copied().std();
        assert!(avg.is_some());
        assert_f32_near!(avg.unwrap(), 2.0, 6);
        assert!(std.is_some());
        assert_f32_near!(std.unwrap(), 0.81649658, 6);
    }

    #[test]
    fn test_range_f32() {
        let range = range::f32(0.0, 10.0, 5).collect::<Vec<_>>();
        let values = [0.0, 2.5, 5.0, 7.5, 10.0];
        assert_eq!(range.len(), values.len());
        iter::zip(range, values).for_each(|(got, want)| {
            assert_f32_near!(got, want, 6);
        });
        let range = range::f32(4.0, 5.0, 1).collect::<Vec<_>>();
        let values = [4.0];
        assert_eq!(range.len(), values.len());
        iter::zip(range, values).for_each(|(got, want)| {
            assert_f32_near!(got, want, 6);
        });
        assert_eq!(range::f32(4.0, 5.0, 0).count(), 0);
    }
}
