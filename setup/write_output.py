import json
import numpy as np
import os

# Open json file containing tracer information
file_path = os.getcwd() + "/model_description.json"
with open(file_path) as read_model:
    model_info = json.load(read_model)
    parameters = model_info["parameters"]

class OutputData():

    def __init__(self,parameters):
        self.count = 0
        self.day = 0
        self.month = 0
        self.daily_ave = np.zeros(())
        self.monthly_ave = np.zeros(())