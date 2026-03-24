#!/usr/bin/env python
import os
import sys
import argparse

# Add the current directory to sys.path to allow imports from ht_screening
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from ht_screening.ht_workflow import HighThroughputWorkflow

def main():
    parser = argparse.ArgumentParser(description="High Throughput Workflow for FFBuilder")
    parser.add_argument("-c", "--candidates", default="candidates.txt", help="Path to candidates SMILES file")
    parser.add_argument("-t", "--training", default="training_data.csv", help="Path to training data CSV (smiles, property)")
    parser.add_argument("-n", "--nsuggest", type=int, default=5, help="Number of molecules to suggest")
    parser.add_argument("-r", "--repeats", type=int, default=1, help="Number of repeat units for polymers")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.candidates):
        print(f"Candidates file {args.candidates} not found. Creating a sample one.")
        with open(args.candidates, "w") as f:
            f.write("CCCC\nCCCCCO\nc1ccccc1\nC1CCCCC1\nCC(=O)O\n")
            
    workflow = HighThroughputWorkflow(args.candidates, args.training)
    workflow.run_iteration(n_suggest=args.nsuggest, repeats=args.repeats)

if __name__ == "__main__":
    main()
