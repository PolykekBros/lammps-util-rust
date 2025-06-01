pub fn range_f64(begin: f64, end: f64, count: usize) -> Vec<f64> {
    assert!(begin <= end);
    let step = (end - begin) / 1.max(count - 1) as f64;
    (0..count).map(|n| begin + n as f64 * step).collect()
}

pub fn get_avg(data: &[f64]) -> Option<f64> {
    if data.is_empty() {
        return None;
    }
    Some(data.iter().sum::<f64>() / data.len() as f64)
}

pub fn get_std(data: &[f64]) -> Option<f64> {
    let n = data.len();
    match n {
        0 => None,
        1 => Some(0.0),
        2.. => {
            let avg = get_avg(data)?;
            Some(data.iter().copied().map(|x| (x - avg).powi(2)).sum::<f64>() / (n - 1) as f64)
        }
    }
}

pub fn get_avg_with_std(data: &[f64]) -> Option<(f64, f64)> {
    get_avg(data).zip(get_std(data))
}
