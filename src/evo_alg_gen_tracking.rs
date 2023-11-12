use std::{collections::HashMap, sync::Arc};

use log::{debug, info};
use rand::Rng;

use crate::{evo_alg_utils::{ParentPicking, ParentoFrontEliteParents, OffSpringCreator, Mutator, ChildGenerator}, contigs::Contig, bin_info_storing::BinSet, evo_alg_utils::{generate_random_permutation, create_new_permutations_from_best_bin_permut}};



pub struct EvoAlgRunner {

    pub permutations_per_generation: usize,
    pub max_generations: usize,
    pub evo_alg_state: EvoAlgState,
    pub child_generator: ChildGenerator,
    pub bin_to_contig_set_hashmap: HashMap<usize, Vec<Arc<Contig>>>,
    pub parent_picking_approach: Box<dyn ParentPicking>


}

impl EvoAlgRunner {

    pub fn initialise(initial_potential_bins: Vec<Vec<Arc<Contig>>>, parent_picking_approach: Box<dyn ParentPicking>, permutations_per_generation: usize, max_generations: usize, child_generator: ChildGenerator) -> EvoAlgRunner { // default options
        let bin_to_contig_set_hashmap: HashMap<usize, Vec<Arc<Contig>>> = initial_potential_bins.into_iter().enumerate().map(|(i, contig_set)| (i, contig_set)).collect();
        let evo_alg_state = EvoAlgState::initialise();
        EvoAlgRunner {

            permutations_per_generation: permutations_per_generation,
            max_generations: max_generations,
            evo_alg_state: evo_alg_state,
            child_generator: child_generator,
            bin_to_contig_set_hashmap: bin_to_contig_set_hashmap,
            parent_picking_approach: parent_picking_approach,


        }
    }




    pub fn create_initial_permutations(&mut self) -> Vec<Vec<usize>>{
        let all_bin_ids: Vec<usize> = self.bin_to_contig_set_hashmap.keys().cloned().collect();
        let mut initial_permutations = Vec::new();
        while initial_permutations.len() < self.permutations_per_generation {
        initial_permutations.push(generate_random_permutation(&all_bin_ids, &initial_permutations));
        }

        initial_permutations
    }

    pub fn create_new_bin_set_permutations_from_original(&mut self, initial_mutation_val: usize, initial_best_bins: &mut BinSet) -> Vec<Vec<usize>> {
        // takes the bin set found in the classic best bins initially. Which will only contain a fragment of the bins, since it creates a filtered non redundant bin set.
        // It then creates a full list of all the bins with the best bins first (so that is the order first tested) and then mutates, but the mutations this time target the best bin region
        // since mutating the indices of bins not picked in the first place would do nothing.
        let best_bin_contig_sets: Vec<&Vec<Arc<Contig>>> = initial_best_bins.best_bins.iter().map(|x| &x.0).collect();
        let initial_best_bin_nums: Vec<usize> = self.bin_to_contig_set_hashmap.iter()
        .filter_map(|(x, y)| if best_bin_contig_sets.contains(&y) { Some(*x) } else { None })
        .collect();
        debug!("The following nums for the initial best bin set: {:?}", &initial_best_bin_nums);
        let unused_bin_nums_for_permutation: Vec<usize> = self.bin_to_contig_set_hashmap.iter()
        .filter_map(|(x, y)| if !best_bin_contig_sets.contains(&y) { Some(*x) } else { None })
        .collect();

        let mut initial_best_bin_nums_remaining_unused_bins_added = initial_best_bin_nums.clone();
        initial_best_bin_nums_remaining_unused_bins_added.extend(unused_bin_nums_for_permutation);
        initial_best_bins.bin_set_order = Some(initial_best_bin_nums_remaining_unused_bins_added.clone());

        let initial_permutations_to_test_by_number = create_new_permutations_from_best_bin_permut(initial_best_bin_nums_remaining_unused_bins_added, initial_best_bin_nums, self.permutations_per_generation, initial_mutation_val);

        self.evo_alg_state.add_permutation(initial_best_bins.clone());
        initial_permutations_to_test_by_number

    }

    pub fn generate_contig_set_order_from_permutation(&self, permutation_order: &Vec<usize>) -> Vec<Vec<Arc<Contig>>> {

        permutation_order.iter().map(|x| self.bin_to_contig_set_hashmap.get(x).unwrap().clone()).collect()

    }


    pub fn create_new_generation(&mut self, last_generation: Vec<BinSet>) -> Vec<Vec<usize>>

    {
        info!("Picking parents...");
        let parents = &self.parent_picking_approach.pick_parents(&mut last_generation.clone());
        let parent_orders: Vec<Vec<usize>> = parents.iter().map(|x| x.bin_set_order.clone().unwrap()).collect();

        debug!("creating new generation, current generation is: {}", &self.evo_alg_state.generation);
        let new_orders = self.child_generator.generate_children(&parent_orders);

        debug!("New permutation of {} created.", &new_orders.len());
        new_orders



    }


}

pub struct EvoAlgState {

    pub generation: usize,
    pub used_permutations: Vec<BinSet>,
    pub best_permutation: Option<BinSet>,

}

impl EvoAlgState {


    pub fn initialise() -> EvoAlgState{

        EvoAlgState {
            generation: 0,
            used_permutations: Vec::new(),
            best_permutation: None
        }
    }


    fn initialise_with_a_best_permutation(permutation: BinSet) -> EvoAlgState {

        EvoAlgState {
            generation: 0,
            used_permutations: vec![permutation.clone()],
            best_permutation: Some(permutation)
        }


    }

    pub fn next_generation(&mut self, completed_bin_sets: &Vec<BinSet> ) {
        for bin_set in completed_bin_sets {
            self.add_permutation(bin_set.clone());
        }
        self.generation += 1;

    }

    pub fn add_permutation(&mut self, permutation: BinSet) {

        self.used_permutations.push(permutation);

    }


    pub fn give_used_permutations_by_vector_of_ids(&self) -> Vec<Vec<usize>> {
        println!("{:?}", &self.used_permutations.len());
        let bin_set_order: Vec<Vec<usize>> = self.used_permutations.iter().map(|x| x.bin_set_order.clone().unwrap()).collect();
        bin_set_order

    }
}

/*
struct EvoAlgBuilder {

    permutations_per_generation: usize,
    mutation_val: usize,
    max_generations: usize,
    evo_alg_state: EvoAlgState,
    bin_to_contig_set_hashmap: HashMap<usize, Vec<Arc<Contig>>>,
    parent_picking_approach: Box<dyn ParentPicking>

}

impl EvoAlgBuilder {
    pub fn new(&self,  initial_potential_bins: Vec<Vec<Arc<Contig>>>) -> EvoAlgBuilder {
        let bin_to_contig_set_hashmap: HashMap<usize, Vec<Arc<Contig>>> = initial_potential_bins.into_iter().enumerate().map(|(i, contig_set)| (i, contig_set)).collect();

        EvoAlgBuilder {

            permutations_per_generation: 100,
            mutation_val: 5,
            max_generations: 100,
            evo_alg_state: EvoAlgState::initialise(),
            bin_to_contig_set_hashmap: bin_to_contig_set_hashmap,
            parent_picking_approach: Box::new(ParentoFrontEliteParents {parents_required: 20})

        }


    }

    pub fn build(self) -> EvoAlgRunner {
        EvoAlgRunner {

            permutations_per_generation: self.permutations_per_generation,
            mutation_val: self.mutation_val,
            max_generations: self.max_generations,
            evo_alg_state: self.evo_alg_state,
            bin_to_contig_set_hashmap:self.bin_to_contig_set_hashmap,
            parent_picking_approach: self.parent_picking_approach

        }

    }
}
*/
