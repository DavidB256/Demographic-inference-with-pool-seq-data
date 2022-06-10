import moments

def control_model(params, ns):
    # params = [nu, T]
    return moments.Demographics1D.growth(params, ns)

def two_pop_split_model():
    pass

# This script should iterate through pipeline_instructions
instructions = "/scratch/djb3ve/data/pipeline_instructions.txt"

with open(instructions, "r") as f:
