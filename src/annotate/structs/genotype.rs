#[derive(Debug, Clone)]
pub struct Genotype {
    pub alleles: Vec<usize>,
    pub phased: bool,
}

impl Genotype {
    pub fn parse(gt: &str) -> Option<Self> {
        if gt == "." {
            return None;
        }

        let phased = gt.contains('|');
        let sep = if phased { '|' } else { '/' };

        let alleles = gt
            .split(sep)
            .map(|a| a.parse::<usize>().ok())
            .collect::<Option<Vec<_>>>()?;

        Some(Self { alleles, phased })
    }

    pub fn ploidy(&self) -> usize {
        self.alleles.len()
    }

    pub fn validate_against_map(&self, allele_map: &[Option<usize>]) -> bool {
        self.alleles.iter().all(|&a| {
            if a == 0 {
                true
            } else {
                allele_map.get(a - 1).and_then(|x| *x).is_some()
            }
        })
    }

    pub fn remap(&self, allele_map: &[Option<usize>]) -> Option<Self> {
        if !self.validate_against_map(allele_map) {
            return None;
        }

        let mut out = Vec::with_capacity(self.alleles.len());
        for &a in &self.alleles {
            if a == 0 {
                out.push(0);
            } else {
                out.push(allele_map[a - 1]? + 1);
            }
        }

        Some(Self {
            alleles: out,
            phased: self.phased,
        })
    }

    pub fn to_string(&self) -> String {
        let sep = if self.phased { "|" } else { "/" };
        self.alleles
            .iter()
            .map(|a| a.to_string())
            .collect::<Vec<_>>()
            .join(sep)
    }

    /// G-index (bcftools-compatible)
    pub fn g_index(&self) -> usize {
        match self.ploidy() {
            1 => self.alleles[0],
            2 => {
                let a = self.alleles[0];
                let b = self.alleles[1];
                if a <= b {
                    b * (b + 1) / 2 + a
                } else {
                    a * (a + 1) / 2 + b
                }
            }
            _ => 0, // fallback for higher ploidy
        }
    }
}
