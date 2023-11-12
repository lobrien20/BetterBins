use std::cmp::Ordering;
use std::collections::HashMap;

use log::debug;
use log::info;
use log::warn;
use ndarray::ShapeBuilder;
use non_dominated_sort::DominanceOrd;
use non_dominated_sort::non_dominated_sort;
use rand::Rng;
use rand::seq::SliceRandom;
use crate::bin_info_storing::BinSet;

pub fn get_parent_picking_approach_from_string(parent_picking_string: String, parents_required: usize) -> Box<dyn ParentPicking> {

    // ----NOTE---- here the parents required mean that it will generate x number of parents, if required number of permutations is greater than parents required, then elite parents generated.


  //  if parent_picking_string == "Elite" {
   //     Box::new(EliteParentPicker {parents_required: elite_parents})
  //  }  else if parent_picking_string == "ParentoElite" {


    Box::new(nsga2 {parents_required: parents_required, tournament_size: 2})

}


pub trait ParentPicking {
    fn pick_parents(&self, all_possible_parents: &mut Vec<BinSet>) -> Vec<BinSet>;
}

struct EliteParentPicker {
    parents_required: usize,
}
impl ParentPicking for EliteParentPicker {

    fn pick_parents(&self, all_possible_parents: &mut Vec<BinSet>) -> Vec<BinSet> {

        debug!("Running elite parent picking!");
        all_possible_parents.sort_by(|a, b| b.overall_score.unwrap().partial_cmp(&a.overall_score.unwrap()).unwrap());
        let mut elite_parents = Vec::new();
        for parent in all_possible_parents {
            elite_parents.push(parent.clone());
            if elite_parents.len() == self.parents_required {
                break
            }
        }

        elite_parents
    }
}

pub struct ParentoFrontEliteParents {
    pub parents_required: usize
}
impl DominanceOrd for ParentoFrontEliteParents {
    type T = BinSet;

    fn dominance_ord(&self, a: &Self::T, b: &Self::T) -> Ordering {
        best_bins_dominance_ord(a, b)
    }
}


pub fn best_bins_dominance_ord(a: &BinSet, b: &BinSet) -> Ordering {
    let a_comp = a.overall_completeness;
    let b_comp = b.overall_completeness;
    let a_cont = a.overall_contamination;
    let b_cont = b.overall_contamination;
    let a_bins = a.best_bins.len();
    let b_bins = b.best_bins.len();

    if (a_comp > b_comp && a_cont <= b_cont && a_bins <= b_bins) || (a_comp >= b_comp && a_cont < b_cont && a_bins <= b_bins) || (a_comp >= b_comp && a_cont <= b_cont && a_bins < b_bins) {
        Ordering::Greater
    } else if (a_comp < b_comp && a_cont >= b_cont && a_bins >= b_bins) || (a_comp <= b_comp && a_cont > b_cont && a_bins >= b_bins) || (a_comp <= b_comp && a_cont >= b_cont && a_bins > b_bins) {
        Ordering::Less
    } else {
        Ordering::Equal
    }
}


impl ParentPicking for ParentoFrontEliteParents {

    fn pick_parents(&self, all_possible_parents: &mut Vec<BinSet>) -> Vec<BinSet> {
        let result = non_dominated_sort(&all_possible_parents, self);
        let mut elite_parents = Vec::new();
        while elite_parents.len() != self.parents_required {
            let parento_front: Vec<BinSet> = result.iter().map(|x| x.0.clone()).collect();
            for bin_set in parento_front {
                if elite_parents.len() == self.parents_required {
                    break
                } else {
                    elite_parents.push(bin_set);
                }

            }
        }
        elite_parents

    }
}




pub trait OffSpringCreator {
    fn create_offspring_from_parents(&self, vector_of_permutations: &Vec<Vec<usize>>, disallowed_permutations: &mut Vec<Vec<usize>>, required_permutations: usize) -> Vec<Vec<usize>>;
}


pub struct Mutator {
    pub num_of_mutations: usize
}
impl Mutator {

    pub fn initialise(perc_mutation: f64, number_of_potential_bins: usize) -> Mutator {
        let mut num_of_mutations: usize = (number_of_potential_bins as f64 * (perc_mutation / 100.0) ).round() as usize;
        info!("Mutation percentage is {}, amount is {}", perc_mutation, num_of_mutations);
        if num_of_mutations == 0 {
            info!("Warning: mutation percentage so low that number of mutations (rounded) would equal 0, so made number of mutations 1.");
            num_of_mutations = 1;
        }
        Mutator { num_of_mutations: num_of_mutations }
    }
    fn mutate_permutation(&self, permutation: &Vec<usize>, disallowed_permutations: &Vec<Vec<usize>> ) -> Vec<usize>{
        let mut new_permutation = permutation.clone();


        // Swap the values at idx1 and idx2
        for i in 0..self.num_of_mutations as i16 {

            let mut rng = rand::thread_rng();
            let idx1 = rng.gen_range(0..permutation.len());
            let idx2 = rng.gen_range(0..permutation.len());
            new_permutation.swap(idx1, idx2);

        }

        while disallowed_permutations.contains(&&new_permutation) {
            new_permutation = self.mutate_permutation(permutation, disallowed_permutations);
        }

        new_permutation
    }
}

impl OffSpringCreator for Mutator {
    fn create_offspring_from_parents(&self, vector_of_permutations: &Vec<Vec<usize>>, disallowed_permutations: &mut Vec<Vec<usize>>, required_permutations: usize ) -> Vec<Vec<usize>> {
        let mut offspring = Vec::new();
        let mut rng = rand::thread_rng();


        loop {

            let permutation_to_test = vector_of_permutations[rng.gen_range(0..(vector_of_permutations.len() -1))].clone();
            let mutated_child = self.mutate_permutation(&permutation_to_test, &offspring);

            offspring.push(mutated_child);

            if offspring.len() == required_permutations {
                break
            }

        }
        offspring
    }



}

pub struct OrderedCrossover {
    crossover_amount: usize
}

impl OrderedCrossover {

    pub fn initialise(perc_crossover: f64, number_of_potential_bins: usize) -> OrderedCrossover {

        let crossover_amount: usize = (number_of_potential_bins as f64 * (perc_crossover / 100.0) ).round() as usize;
        if crossover_amount == number_of_potential_bins {
            panic!("Crossover percentage is too high, (percentage picked is resulting in size of crossover being equal to the total number of bins");
        }
        if crossover_amount == 0 {
            panic!("Crossover percentage is too low ((percentage picked is resulting in size of crossover being 0)");
        }
        info!("Crossover percentage is {}, amount is {}", perc_crossover, crossover_amount);

        OrderedCrossover {crossover_amount: crossover_amount}


    }
    pub fn crossover_from_list_of_parents(&self, mut list_of_parents: Vec<Vec<usize>> ) -> Vec<Vec<usize>>{
        let mut rng = rand::thread_rng();
        list_of_parents.shuffle(&mut rng);
        let mut parent_iter = list_of_parents.iter();
        let mut all_crossed_children = Vec::new();
        let mut parent_1 = parent_iter.next().unwrap().clone();
        let mut parent_2 = parent_iter.next().unwrap().clone();
        let mut loop_count = 0;
        loop {
            loop_count += 1;

            let (crossed_child_1, crossed_child_2) = self.crossover_parents(parent_1.clone(), parent_2.clone());

            if (!all_crossed_children.contains(&crossed_child_1) && !all_crossed_children.contains(&crossed_child_2)) || loop_count == 100000 {
                if loop_count == 100000 {
                    warn!("Large amount of looping has occured for crossover of parents.");
                    loop_count = 0;
                }
                all_crossed_children.push(crossed_child_1);
                all_crossed_children.push(crossed_child_2);
                if all_crossed_children.len() == list_of_parents.len() {
                    break;
                }
                parent_1 = parent_iter.next().unwrap().clone();
                parent_2 = parent_iter.next().unwrap().clone();


            }

        };
        all_crossed_children
    }
    fn crossover_parents(&self, parent_1: Vec<usize>, parent_2: Vec<usize>) -> (Vec<usize>, Vec<usize>) {

        let (crossover_start, crossover_end) = self.generate_crossover_region(parent_1.len());
        let crossed_child_1 = self.create_crossed_child(&parent_1, (crossover_start, crossover_end), &parent_2);
        let crossed_child_2 = self.create_crossed_child(&parent_2, (crossover_start, crossover_end), &parent_1);
        (crossed_child_1, crossed_child_2)

    }
    fn generate_crossover_region(&self, length_of_permutation: usize) -> (usize, usize) {
        let mut rng = rand::thread_rng();
        let crossover_end = rng.gen_range(self.crossover_amount..length_of_permutation);
        let crossover_start = crossover_end - self.crossover_amount;
        (crossover_start, crossover_end)

    }
    fn create_crossed_child(&self, parent_to_be_crossed: &Vec<usize>, (crossover_start, crossover_end): (usize, usize), crossing_parent: &Vec<usize>) -> Vec<usize>{
        let other_parent_crossing_region = &crossing_parent[crossover_start..crossover_end];
        let mut parent_remaining_unpicked = parent_to_be_crossed.iter().filter(|x| !other_parent_crossing_region.contains(x));
        let mut child = Vec::new();
        let mut other_parent_cross_iter = other_parent_crossing_region.iter();
        for i in 0..parent_to_be_crossed.len() {
            if i >= crossover_start && i < crossover_end {
                child.push(other_parent_cross_iter.next().unwrap().clone());
            } else {
                child.push(parent_remaining_unpicked.next().unwrap().clone());
            }
        }
        child
    }

}

pub struct ChildGenerator {
    pub crossover: Option<OrderedCrossover>,
    pub mutator_struct: Option<Mutator>
}
impl ChildGenerator {
    pub fn initialise(perc_mutation: f64, number_of_potential_bins: usize, perc_crossover: f64) -> ChildGenerator {
        let mut crossover_struct = None;
        if perc_crossover > 0.0 {
            info!("Running genetic algorithm with ordered crossover, percentage: {}", perc_crossover);
            crossover_struct = Some(OrderedCrossover::initialise(perc_crossover, number_of_potential_bins));
        }
        let mut mutator_struct = None;
        if perc_mutation > 0.0 {
            info!("Running genetic algorithm with mutator, percentage: {}", perc_mutation);
            mutator_struct = Some(Mutator::initialise(perc_mutation, number_of_potential_bins));

        }

        ChildGenerator {crossover: crossover_struct, mutator_struct: mutator_struct}
    }

    pub fn generate_children(&mut self, list_of_parents: &Vec<Vec<usize>>) -> Vec<Vec<usize>> {
        let mut offspring = list_of_parents.clone();
        match &self.crossover {
            Some(crossover_struct) => {
                info!("Crossing over...");
                debug!("Offspring equals {}", offspring.len());
                offspring = crossover_struct.crossover_from_list_of_parents(offspring);
                info!("Crossover complete!");
                debug!("Offspring equals {}", offspring.len());
            },
            None => ()
        }
        match &self.mutator_struct {
            Some(mutator_struct) => {
                info!("Mutating...");
                debug!("Offspring equals {}", offspring.len());
                mutator_struct.create_offspring_from_parents(&offspring, &mut Vec::new(), offspring.len());
                info!("Mutating complete!");
                debug!("Offspring equals {}", offspring.len());
            }
            None => ()
        }
        offspring
    }

}


pub fn generate_random_permutation(vec_of_bin_ids: &Vec<usize>, used_permutations: &Vec<Vec<usize>>) -> Vec<usize> {
    let mut rng = rand::thread_rng();
    let mut new_permutation = vec_of_bin_ids.clone();
    new_permutation.shuffle(&mut rng);
    while used_permutations.contains(&vec_of_bin_ids) {
        new_permutation.shuffle(&mut rng);
    }

    new_permutation
}


pub fn mutate_permutation(permutation: &Vec<usize>, number_of_mutations: usize) -> Vec<usize>{
    let mut new_permutation = permutation.clone();
    let mut rng = rand::thread_rng();
    let idx1 = rng.gen_range(0..permutation.len());
    let idx2 = rng.gen_range(0..permutation.len());

    // Swap the values at idx1 and idx2
    new_permutation.swap(idx1, idx2);

    new_permutation

}

pub fn create_new_permutations_from_best_bin_permut(vector_of_all_bins_starting_with_best_bins_in_order: Vec<usize>, original_best_bins: Vec<usize>, required_permutations: usize, number_of_mutations: usize) -> Vec<Vec<usize>> {
    // method done for when have the best bins first with the remains added on at the end. Mutates the order but specifically requires that the best bin region has been targetd n times.
    let num_of_original_best_bins = &original_best_bins.len();
    let mut new_permutations = Vec::new();
    loop {
        let new_permutation = mutate_permutation(&vector_of_all_bins_starting_with_best_bins_in_order, number_of_mutations);
        let new_permut_original_best_bin_region = &new_permutation[0..*num_of_original_best_bins];
        println!("NEW PERMUT IS {:?} region, and {:?}", &new_permut_original_best_bin_region, new_permutation);
        println!("TESTING AGAINST ORIGINAL BEST BIN REG: {:?}", original_best_bins);
        if new_permut_original_best_bin_region != original_best_bins { // checks if the bin region is different to the original region where the best bins where decided
            new_permutations.push(new_permutation);
            if new_permutations.len() == required_permutations {
                break
            }
        }
    }
    new_permutations

}




struct nsga2 {
    tournament_size: usize,
    parents_required: usize
}

impl nsga2  {

    fn get_bin_set_parento_window_and_crowd_dists(&self, all_possible_parents: &mut Vec<BinSet>) ->  Vec<(BinSet, f64, usize)> {

            let mut result = non_dominated_sort(&all_possible_parents, self);
            let mut all_bin_set_results_with_crowd_dists = Vec::new();
            let mut front_count = 0;

            loop {
                if result.is_empty() {
                    break
                }

                let parento_front: Vec<BinSet> = result.iter().map(|x| x.0.clone()).collect();
                let mut front_crowd_dist_result = self.generate_crowd_distances_for_binsets(&parento_front, front_count);
                all_bin_set_results_with_crowd_dists.extend(front_crowd_dist_result);
                front_count += 1;
                result = result.next_front();
            }
            all_bin_set_results_with_crowd_dists

    }

    fn generate_tournament_options<'b>(&self, all_bin_set_results_with_crowd_dists: &'b Vec<(BinSet, f64, usize)>) -> Vec<&'b (BinSet, f64, usize)> {
        let mut tournament_bin_sets = Vec::new();
        let mut rng = rand::thread_rng();

        while tournament_bin_sets.len() < self.tournament_size { // caution, need to add something to stop infinite loop
            let bin_set_indice_for_tourn = rng.gen_range(0..all_bin_set_results_with_crowd_dists.len());
            let bin_set = &all_bin_set_results_with_crowd_dists[bin_set_indice_for_tourn];
            if !tournament_bin_sets.contains(&bin_set) {
                tournament_bin_sets.push(bin_set);
            }
        }
        tournament_bin_sets
    }

    fn run_a_tournament(&self, all_bin_set_results_with_crowd_dists: &Vec<(BinSet, f64, usize)>) -> (BinSet, f64, usize) {

        let tournament_bin_sets = self.generate_tournament_options(all_bin_set_results_with_crowd_dists);
        let mut rng = rand::thread_rng();

        let mut best_window = 0;
        let mut best_window_bin_sets = Vec::new();

        for bin_set in tournament_bin_sets {

            if bin_set.2 >= best_window {

                best_window = bin_set.2;
                best_window_bin_sets = vec![bin_set];

            } else if bin_set.2 == best_window {

                best_window_bin_sets.push(bin_set);

            } else {

                continue;

            }
        }
        let mut best_bin_set = best_window_bin_sets.remove(0);
        for bin_set in best_window_bin_sets {

            if bin_set.1 > best_bin_set.1 {
                best_bin_set = bin_set;
            }

            if bin_set.1 == best_bin_set.1 {

                if rng.gen_bool(0.5) == true {
                    best_bin_set = bin_set;

                }
            }
        }
        best_bin_set.clone()


    }



    fn generate_crowd_distances_for_binsets(&self, vector_of_bin_sets: &Vec<BinSet>, parento_front_num: usize) -> Vec<(BinSet, f64, usize)> {
        // note, for the crowding distance calculation it doesn't matter that higher completeness is better and lower contamination is better. It's just about getting the range
        // just in case of confusion for later

        if vector_of_bin_sets.len() == 1 {

            return vec![(vector_of_bin_sets[0].clone(), 9999.9, parento_front_num)]

        } else if vector_of_bin_sets.len() == 2 {
            return vec![(vector_of_bin_sets[0].clone(), 9999.9, parento_front_num), (vector_of_bin_sets[1].clone(), 9999.9, parento_front_num)];
        } else if vector_of_bin_sets.len() == 3 {
            return vec![(vector_of_bin_sets[0].clone(), 9999.9, parento_front_num), (vector_of_bin_sets[1].clone(), 0.0, parento_front_num), (vector_of_bin_sets[2].clone(), 9999.9, parento_front_num)];
        }


        println!("Length of this binset window is... {}", vector_of_bin_sets.len());
        let mut best_bins_sorted_by_comp = vector_of_bin_sets.clone();
        best_bins_sorted_by_comp.sort_by(|a, b| a.overall_completeness.partial_cmp(&b.overall_completeness).unwrap());
        let mut ranked_comp_values: Vec<(f64, Vec<usize>)> = best_bins_sorted_by_comp.iter().map(|x| (x.overall_completeness, x.bin_set_order.clone().unwrap())).collect();
        let comp_crowd_dist_val_hashmap = self.get_crowding_distance_values_for_objective(ranked_comp_values);

        let mut best_bins_sorted_by_cont = vector_of_bin_sets.clone();
        best_bins_sorted_by_cont.sort_by(|a, b| a.overall_contamination.partial_cmp(&b.overall_contamination).unwrap());
        let mut ranked_cont_values: Vec<(f64, Vec<usize>)> = best_bins_sorted_by_cont.iter().map(|x| (x.overall_contamination, x.bin_set_order.clone().unwrap())).collect();
        let cont_crowd_dist_val_hashmap = self.get_crowding_distance_values_for_objective(ranked_cont_values);

        let mut best_bins_sorted_by_num_of_bins = vector_of_bin_sets.clone();
        let mut ranked_num_of_bin_values: Vec<(f64, Vec<usize>)> = best_bins_sorted_by_num_of_bins.iter().map(|x| (x.best_bins.len() as f64, x.bin_set_order.clone().unwrap())).collect();
        let bin_num_val_hashmap = self.get_crowding_distance_values_for_objective(ranked_num_of_bin_values);


        let mut permutation_and_crowd_dists = Vec::new();

       for bin_set in vector_of_bin_sets {

            let order = bin_set.bin_set_order.clone().unwrap();
            let sum_crowd_dist = comp_crowd_dist_val_hashmap.get(&order).unwrap() + cont_crowd_dist_val_hashmap.get(&order).unwrap() + bin_num_val_hashmap.get(&order).unwrap();
            permutation_and_crowd_dists.push((bin_set.clone(), sum_crowd_dist, parento_front_num));

        }

        permutation_and_crowd_dists



    }


    fn get_crowding_distance_values_for_objective(&self, all_values: Vec<(f64, Vec<usize>)>) -> HashMap<Vec<usize>, f64>{
        let mut all_values = self.scale_objective_values(all_values);
        let lowest_value = all_values.remove(0);
        let biggest_value = all_values.pop().unwrap();
        let val_difference = biggest_value.0 - lowest_value.0;
        let mut vals_hashmap = HashMap::new();
        vals_hashmap.insert(lowest_value.1.clone(), 9999.9);
        vals_hashmap.insert(biggest_value.1.clone(), 9999.9);

        println!("Lowest val is {}, Biggest val is {}, val difference is: {}", lowest_value.0, biggest_value.0, val_difference);
        let mut incremented_distance = 0.0;
        for (index, val) in all_values.iter().enumerate() {
            println!("index is {}", index);
            let mut crowding_dist_val = 0.0;
            if index == 0 {

                crowding_dist_val = all_values[(index + 1)].0 - lowest_value.0;

            } else if index == (all_values.len() - 1) {

                crowding_dist_val = biggest_value.0 - all_values[(index - 1)].0;

            } else {

                crowding_dist_val = all_values[(index + 1)].0 - all_values[(index - 1)].0;

            }

            crowding_dist_val = crowding_dist_val / val_difference;
            incremented_distance = incremented_distance + crowding_dist_val;
            vals_hashmap.insert(val.1.clone(), incremented_distance);



        }

        vals_hashmap
    }
    fn scale_objective_values(&self, all_values: Vec<(f64, Vec<usize>)>) -> Vec<(f64, Vec<usize>)>{
        let min_value = all_values[0].0;
        let max_value = all_values[all_values.len() - 1].0;
        let mut scaled_values = Vec::new();
        for value_set in all_values.into_iter() {
            let scaled_value = (value_set.0 - min_value) / (max_value - min_value);
            scaled_values.push((scaled_value, value_set.1));
        }

        scaled_values
    }


}

impl ParentPicking for nsga2 {

    fn pick_parents(&self, all_possible_parents: &mut Vec<BinSet>) -> Vec<BinSet> {
        info!("nsgaII - getting bin set parento windows and crowd distances");
        debug!("All possible parents equal: {}", all_possible_parents.len());
        let parents_crowd_dist_and_window = self.get_bin_set_parento_window_and_crowd_dists(all_possible_parents);
        debug!("All possible parents after bin set parento window and crowd dists equal: {}. This should be the same as the all possible parents equal line above.", parents_crowd_dist_and_window.len());
        info!("nsgaII - bin set parento windows and crowd distanced obtained. Running tournaments...");
        let mut created_parents = Vec::new();
        debug!("nsgaII - parents required = {}", self.parents_required);
        while created_parents.len() != self.parents_required {

            created_parents.push(self.run_a_tournament(&parents_crowd_dist_and_window).0);
            debug!("nsgaII - one parent produced. Created parents equals: {}", created_parents.len());
        }

        created_parents

    }

}
impl DominanceOrd for nsga2 {
    type T = BinSet;

    fn dominance_ord(&self, a: &Self::T, b: &Self::T) -> Ordering {
        best_bins_dominance_ord(a, b)
    }
}


mod tests {
    use super::*;
    use std::{fs, env, path::Path};

    fn create_some_test_permutations() -> Vec<Vec<usize>> {

        let permut_1 = vec![0, 1, 2, 3, 4, 6, 5];
        let permut_2 = vec![2, 4, 3, 1, 5, 0, 6];
        let permut_3 = vec![6, 4, 1, 5, 3, 0, 2];
        let permut_4 = vec![1, 4, 2, 5, 0, 3, 6];

        vec![permut_1, permut_2, permut_3, permut_4]
    }
    fn create_some_test_best_bins() -> Vec<BinSet> {
        // creates four test best bins.
        // overall score has a difference of 10 between each one (for weighted picking)
        // overall completeness has the same, but the third best bin has a higher contamination than the rest.
        let test_permutations = create_some_test_permutations();
        let mut test_best_bins = Vec::new();
        let mut score_to_add = 0.0;

        for permut in test_permutations {

            let mut best_bins_t= BinSet::constructnew();

            best_bins_t.overall_completeness = 50.0;
            best_bins_t.overall_contamination = 0.0;
            best_bins_t.overall_score = Some(score_to_add);
            best_bins_t.bin_set_order = Some(permut);
            test_best_bins.push(best_bins_t);
            score_to_add += 10.0;


        }
        test_best_bins[3].overall_contamination = 10.0;
        test_best_bins[2].overall_completeness = 90.0;
        test_best_bins
    }

    fn create_nsga_ii_struct() -> nsga2 {
        nsga2 { tournament_size: 2, parents_required: 10 }
    }
    #[test]
    fn test_get_crowding_distance_values_for_objective() {

        let nsga_obj = create_nsga_ii_struct();
        let test_permuts = create_some_test_permutations();
        let mut all_values = vec![(10.0, test_permuts[0].clone()), (55.0, test_permuts[1].clone()), (60.0, test_permuts[2].clone()), (90.0, test_permuts[3].clone())];

        let test_hash_result = nsga_obj.get_crowding_distance_values_for_objective(all_values);

        let mut expected_hashmap = HashMap::new();
        expected_hashmap.insert(test_permuts[0].clone(), 9999.9);
        expected_hashmap.insert(test_permuts[1].clone(), 0.625);
        expected_hashmap.insert(test_permuts[2].clone(), 1.0625);
        expected_hashmap.insert(test_permuts[3].clone(), 9999.9);
        assert_eq!(expected_hashmap, test_hash_result);


    }
    #[test]
    fn test_generate_bin_set_parento_window_and_crowd_dists_for_binsets() {
        let mut test_best_bins = create_some_test_best_bins();
        let nsga_obj = create_nsga_ii_struct();
        let test_result = nsga_obj.get_bin_set_parento_window_and_crowd_dists(&mut test_best_bins);
        assert_eq!(test_result.len(), 4);
        let expected_windows = 3;
        let mut unique_windows = Vec::new();
        for test_res in test_result {
            println!("{}, {}, comp is: {} contam is {}", test_res.1, test_res.2, test_res.0.overall_completeness, test_res.0.overall_contamination);
            if !unique_windows.contains(&test_res.2) {
                unique_windows.push(test_res.2);
            }
        }

        assert_eq!(expected_windows, unique_windows.len());
    }
    #[test]
    fn test_run_a_tournament() {
        let mut test_best_bins = create_some_test_best_bins();
        let nsga_obj = create_nsga_ii_struct();
        let test_result = nsga_obj.get_bin_set_parento_window_and_crowd_dists(&mut test_best_bins);
        println!("{:?}", test_result);
        for i in 0..50 {
            let tournament_result = nsga_obj.run_a_tournament(&test_result);
            if tournament_result == test_result[0] {
                println!("{:?}", test_result[0]);
                panic!("Somehow the worst item has won the tournament, should never happen.")
            }
        }
        let additional_test = vec![test_result[0].clone(), test_result[2].clone()];
        let tournament_result_2 = nsga_obj.run_a_tournament(&additional_test);
        assert_eq!(test_result[2], tournament_result_2);

    }
}
