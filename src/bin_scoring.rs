pub struct BinScorer {
    pub contamination_weight: f64,
    pub completion_weight: f64,


}

impl BinScorer {

    pub fn initialise_bin_scorer(contamination_weight: f64, completion_weight: f64) -> BinScorer {
        BinScorer {contamination_weight: contamination_weight, completion_weight: completion_weight}
    }

    pub fn score_bin(&self, completion_value: f64, contamination_value: f64) -> f64 {

        let bin_score = ((&self.completion_weight * completion_value) - (&self.contamination_weight * contamination_value)) / 100.0;
        bin_score

    }
}