//! Dense linear solver: Gaussian elimination with partial pivoting.

/// Solves `a * x = b`. Returns `None` if the matrix is singular.
pub fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        // Find the row with the largest value in this column.
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1e-14 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        let (upper, lower) = a.split_at_mut(col + 1);
        let pivot_row = &upper[col];
        for (offset, row) in lower.iter_mut().enumerate() {
            let factor = row[col] / pivot_row[col];
            if factor == 0.0 {
                continue;
            }
            for (value, pivot_value) in row[col..].iter_mut().zip(&pivot_row[col..]) {
                *value -= factor * pivot_value;
            }
            b[col + 1 + offset] -= factor * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let sum: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - sum) / a[row][row];
    }
    Some(x)
}

/// LU factorization with partial pivoting. Use it to solve the same
/// matrix with many right-hand sides.
pub struct Lu {
    /// L (below the diagonal, unit diagonal) and U (diagonal and above).
    factors: Vec<Vec<f64>>,
    /// Row order after the row exchanges.
    order: Vec<usize>,
}

impl Lu {
    /// Returns `None` if the matrix is singular.
    pub fn new(mut a: Vec<Vec<f64>>) -> Option<Self> {
        let n = a.len();
        let mut order: Vec<usize> = (0..n).collect();
        for col in 0..n {
            let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
            if a[pivot][col].abs() < 1e-14 {
                return None;
            }
            a.swap(col, pivot);
            order.swap(col, pivot);
            let (upper, lower) = a.split_at_mut(col + 1);
            let pivot_row = &upper[col];
            for row in lower.iter_mut() {
                let factor = row[col] / pivot_row[col];
                row[col] = factor;
                if factor == 0.0 {
                    continue;
                }
                for (value, pivot_value) in row[col + 1..].iter_mut().zip(&pivot_row[col + 1..]) {
                    *value -= factor * pivot_value;
                }
            }
        }
        Some(Self { factors: a, order })
    }

    /// Solves `a * x = b`.
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.order.len();
        let mut x: Vec<f64> = self.order.iter().map(|&i| b[i]).collect();
        for row in 1..n {
            let sum: f64 = (0..row).map(|k| self.factors[row][k] * x[k]).sum();
            x[row] -= sum;
        }
        for row in (0..n).rev() {
            let sum: f64 = (row + 1..n).map(|k| self.factors[row][k] * x[k]).sum();
            x[row] = (x[row] - sum) / self.factors[row][row];
        }
        x
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn solves_small_system() {
        let a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let x = super::solve(a, vec![3.0, 5.0]).unwrap();
        assert!((x[0] - 0.8).abs() < 1e-12 && (x[1] - 1.4).abs() < 1e-12);
    }

    #[test]
    fn needs_row_exchange() {
        // The first pivot is zero, thus the solver must exchange rows.
        let a = vec![
            vec![0.0, 2.0, 1.0],
            vec![1.0, 1.0, 1.0],
            vec![2.0, 0.0, 3.0],
        ];
        let x = super::solve(a, vec![7.0, 6.0, 11.0]).unwrap();
        for (value, expected) in x.iter().zip([1.0, 2.0, 3.0]) {
            assert!((value - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn lu_gives_the_same_result_as_direct_solution() {
        let a = vec![
            vec![0.0, 2.0, 1.0, -1.0],
            vec![1.0, 1.0, 1.0, 3.0],
            vec![2.0, 0.0, 3.0, 0.5],
            vec![-1.0, 4.0, 0.0, 2.0],
        ];
        let lu = super::Lu::new(a.clone()).unwrap();
        for b in [vec![7.0, 6.0, 11.0, 1.0], vec![0.0, 1.0, 0.0, -2.0]] {
            let direct = super::solve(a.clone(), b.clone()).unwrap();
            for (x, y) in lu.solve(&b).iter().zip(&direct) {
                assert!((x - y).abs() < 1e-12);
            }
        }
        assert!(super::Lu::new(vec![vec![1.0, 2.0], vec![2.0, 4.0]]).is_none());
    }

    #[test]
    fn singular_matrix_gives_none() {
        let a = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert!(super::solve(a, vec![1.0, 2.0]).is_none());
    }
}
