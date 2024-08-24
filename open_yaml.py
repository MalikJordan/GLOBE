import yaml
import os

file = os.getcwd() + '/model_description.yaml'
with open(file, 'r') as f:
    data = yaml.full_load(f)

print()