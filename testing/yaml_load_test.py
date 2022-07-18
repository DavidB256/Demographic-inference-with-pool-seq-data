import yaml

with open("config.yaml", "r") as f:
    yd = yaml.safe_load(f)

print(yd)
