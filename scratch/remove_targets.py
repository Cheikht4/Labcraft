import yaml
with open("validation/DENV2.yaml", "r") as f:
    data = yaml.safe_load(f)

if "targets" in data:
    del data["targets"]

with open("validation/DENV2.yaml", "w") as f:
    yaml.dump(data, f, default_flow_style=False, sort_keys=False)
